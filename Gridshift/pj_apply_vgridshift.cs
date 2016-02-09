//*****************************************************************************
//
// Project:	PROJ.4
// Purpose:	Apply vertical datum shifts based on grid shift files, normally
//			geoid grids mapping WGS84 to NAVD88 or something similar.
// Author:	Frank Warmerdam, warmerdam@pobox.com
//
//*****************************************************************************
// Copyright (c) 2010, Frank Warmerdam <warmerdam@pobox.com>
// Copyright (c) 2009-2015 by the Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//****************************************************************************

using System;
using System.Text;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Gridshift
{
	public static partial class Grid
	{
		//**********************************************************************
		//						pj_apply_vgridshift()
		//
		//		This implmentation takes uses the gridlist from a coordinate
		//		system definition. If the gridlist has not yet been
		//		populated in the coordinate system definition we set it up
		//		now.
		//**********************************************************************
		static int debug_count_apply_vgridshift=0;
		public static int pj_apply_vgridshift(PJ defn, string listname, ref PJ_GRIDINFO[] gridlist_p, ref int gridlist_count_p, bool inverse, int point_count, int point_offset, double[] x, double[] y, double[] z)
		{
			if(gridlist_p==null)
			{
				gridlist_p=pj_gridlist_from_nadgrids(Proj.pj_get_ctx(defn), Proj.pj_param_s(defn.ctx, defn.parameters, listname), out gridlist_count_p);
				if(gridlist_p==null||gridlist_count_p==0) return defn.ctx.last_errno;
			}

			if(gridlist_p==null||gridlist_count_p==0)
			{
				Proj.pj_ctx_set_errno(defn.ctx, -38);
				return -38;
			}

			PJ_GRIDINFO[] tables=gridlist_p;
			defn.ctx.last_errno=0;

			for(int i=0; i<point_count; i++)
			{
				int io=i*point_offset;
				double value=Libc.HUGE_VAL;

				LP input;
				input.phi=y[io];
				input.lam=x[io];

				// keep trying till we find a table that works
				for(int itable=0; itable<gridlist_count_p; itable++)
				{
					PJ_GRIDINFO gi=tables[itable];
					CTABLE ct=gi.ct;
					double grid_x, grid_y;
					int grid_ix, grid_iy;

					// skip tables that don't match our point at all.
					if(ct.ll.phi>input.phi||ct.ll.lam>input.lam||ct.ll.phi+(ct.lim.phi-1)*ct.del.phi<input.phi||ct.ll.lam+(ct.lim.lam-1)*ct.del.lam<input.lam)
						continue;

					// If we have child nodes, check to see if any of them apply.
					while(gi.child!=null)
					{
						PJ_GRIDINFO child=gi.child;

						for(; child!=null; child=child.next)
						{
							CTABLE ct1=child.ct;

							if(ct1.ll.phi>input.phi||ct1.ll.lam>input.lam||ct1.ll.phi+(ct1.lim.phi-1)*ct1.del.phi<input.phi||ct1.ll.lam+(ct1.lim.lam-1)*ct1.del.lam<input.lam)
								continue;

							break;
						}

						// If we didn't find a child then nothing more to do,
						if(child==null) break;

						// otherwise use the child, first checking it's children.
						gi=child;
						ct=child.ct;
					}

					// load the grid shift info if we don't have it.
					if(ct.cvs==null&&!pj_gridinfo_load(Proj.pj_get_ctx(defn), gi))
					{
						Proj.pj_ctx_set_errno(defn.ctx, -38);
						return -38;
					}

					// Interpolation a location within the grid
					grid_x=(input.lam-ct.ll.lam)/ct.del.lam;
					grid_y=(input.phi-ct.ll.phi)/ct.del.phi;
					grid_ix=(int)Math.Floor(grid_x);
					grid_iy=(int)Math.Floor(grid_y);
					grid_x-=grid_ix;
					grid_y-=grid_iy;

					LP cvs1=ct.cvs[(grid_ix+grid_iy*ct.lim.lam)/2];
					LP cvs2=ct.cvs[(grid_ix+(grid_iy+1)*ct.lim.lam)/2];
					value=cvs1.lam*(1.0-grid_x)*(1.0-grid_y)+cvs1.phi*grid_x*(1.0-grid_y)+
							cvs2.lam*(1.0-grid_x)*grid_y+cvs2.phi*grid_x*grid_y;

					// nodata?
					// GTX official nodata value is -88.88880f, but some grids also
					// use other big values for nodata (e.g naptrans2008.gtx has
					// nodata values like -2147479936), so test them too
					if (value > 1000 || value < -1000 || Math.Abs(value+88.8888)<0.0001) value=Libc.HUGE_VAL;
					else
					{
						if(inverse) z[io]-=value;
						else z[io]+=value;
					}

					if(value!=Libc.HUGE_VAL)
					{
						if(debug_count_apply_vgridshift++<20) Proj.pj_log(defn.ctx, PJ_LOG.DEBUG_MINOR, "pj_apply_gridshift(): used {0}", ct.id);
						break;
					}
				}

				if(value==Libc.HUGE_VAL)
				{
					Proj.pj_log(defn.ctx, PJ_LOG.DEBUG_MAJOR, "pj_apply_vgridshift(): failed to find a grid shift table for\n\t\t\tlocation ({0}dW,{1}dN)", x[io]*Proj.RAD_TO_DEG, y[io]*Proj.RAD_TO_DEG);
					StringBuilder gridlist=new StringBuilder();
					for(int itable=0; itable<gridlist_count_p; itable++)
					{
						PJ_GRIDINFO gi=tables[itable];
						if(itable==0) gridlist.AppendFormat("   tried: {0}", gi.gridname);
						else gridlist.AppendFormat(",{0}", gi.gridname);
					}
					Proj.pj_log(defn.ctx, PJ_LOG.DEBUG_MAJOR, "{0}", gridlist.ToString());

					Proj.pj_ctx_set_errno(defn.ctx, (int)PJD_ERR.GRID_AREA);
					return (int)PJD_ERR.GRID_AREA;
				}
			}

			return 0;
		}
	}
}
