//*****************************************************************************
//
// Project:	PROJ.4
// Purpose:	Code to manage the list of currently loaded (cached) PJ_GRIDINFOs
//			See pj_gridinfo.cs for details of loading individual grids.
// Author:	Frank Warmerdam, warmerdam@pobox.com
//
//*****************************************************************************
// Copyright (c) 2000, Frank Warmerdam <warmerdam@pobox.com>
// Copyright (c) 2009-2011 by the Authors
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
//*****************************************************************************

using System;
using System.Collections.Generic;
using System.Text;

namespace Free.Ports.Proj4.Gridshift
{
	public static partial class Grid
	{
		static PJ_GRIDINFO grid_list=null;

		//**********************************************************************
		//						pj_deallocate_grids()
		//
		//		Deallocate all loaded grids.
		//**********************************************************************
		public static void pj_deallocate_grids()
		{
			while(grid_list!=null)
			{
				PJ_GRIDINFO item=grid_list;
				grid_list=grid_list.next;
				item.next=null;
				pj_gridinfo_free(Proj.pj_get_default_ctx(), item);
				item=null;
			}
		}

		//**********************************************************************
		//						pj_gridlist_merge_grid()
		//
		//		Find/load the named gridfile and merge it into the
		//		last_nadgrids_list.
		//**********************************************************************
		static bool pj_gridlist_merge_gridfile(projCtx ctx, string gridname, ref PJ_GRIDINFO[] p_gridlist, ref int p_gridcount, ref int p_gridmax)
		{
			bool got_match=false;
			PJ_GRIDINFO this_grid, tail=null;

			// --------------------------------------------------------------------
			//		Try to find in the existing list of loaded grids. Add all
			//		matching grids as with NTv2 we can get many grids from one
			//		file (one shared gridname).
			// --------------------------------------------------------------------
			for(this_grid=grid_list; this_grid!=null; this_grid=this_grid.next)
			{
				if(this_grid.gridname==gridname)
				{
					got_match=true;

					// dont add to the list if it is invalid.
					if(this_grid.ct==null) return false;

					// do we need to grow the list?
					if(p_gridcount>=p_gridmax-2)
					{
						PJ_GRIDINFO[] new_list;
						int new_max=p_gridmax+20;

						new_list=new PJ_GRIDINFO[new_max];
						if(p_gridlist!=null)
						{
							Array.Copy(p_gridlist, new_list, p_gridmax);
							p_gridlist=null;
						}

						p_gridlist=new_list;
						p_gridmax=new_max;
					}

					// add to the list
					p_gridlist[p_gridcount++]=this_grid;
					p_gridlist[p_gridcount]=null;
				}

				tail=this_grid;
			}

			if(got_match) return true;

			// --------------------------------------------------------------------
			//		Try to load the named grid.
			// --------------------------------------------------------------------
			this_grid=pj_gridinfo_init(ctx, gridname);

			// we should get at least a stub grid with a missing "ct" member
			if(this_grid==null) return false;

			if(tail!=null) tail.next=this_grid;
			else grid_list=this_grid;

			// --------------------------------------------------------------------
			//		Recurse to add the grid now that it is loaded.
			// --------------------------------------------------------------------
			return pj_gridlist_merge_gridfile(ctx, gridname, ref p_gridlist, ref p_gridcount, ref p_gridmax);
		}

		//**********************************************************************
		//						pj_gridlist_from_nadgrids()
		//
		//		This functions loads the list of grids corresponding to a
		//		particular nadgrids string into a list, and returns it. The
		//		list is kept around till a request is made with a different
		//		string in order to cut down on the string parsing cost, and
		//		the cost of building the list of tables each time.
		//**********************************************************************
		static object staticLockDummy=new object();
		public static PJ_GRIDINFO[] pj_gridlist_from_nadgrids(projCtx ctx, string nadgrids, out int grid_count)
		{
			PJ_GRIDINFO[] gridlist=null;
			int grid_max=0;

			Proj.pj_errno=0;
			grid_count=0;

			lock(staticLockDummy)
			{
				// --------------------------------------------------------------------
				//		Loop processing names out of nadgrids one at a time.
				// --------------------------------------------------------------------
				string s=nadgrids;

				for(int i=0; i<s.Length; )
				{
					bool required=true;
					string name;

					if(s[i]=='@')
					{
						required=false;
						i++;
					}

					int end_char=i;
					while(end_char<s.Length&&s[end_char]!=',') end_char++;
					name=s.Substring(i, end_char-i);

					i=end_char;
					if(i<s.Length&&s[i]==',') i++;

					if(!pj_gridlist_merge_gridfile(ctx, name, ref gridlist, ref grid_count, ref grid_max)&&required)
					{
						Proj.pj_ctx_set_errno(ctx, -38);
						return null;
					}

					Proj.pj_errno=0;
				}

				return gridlist;
			}
		}
	}
}
