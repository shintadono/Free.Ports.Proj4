//*****************************************************************************
//
// Project:	PROJ.4
// Purpose:	Perform overall coordinate system to coordinate system 
//			transformations (pj_transform() function) including reprojection
//			and datum shifting.
// Author:	Frank Warmerdam, warmerdam@pobox.com
//
//*****************************************************************************
// Copyright (c) 2000, Frank Warmerdam
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
using Free.Ports.Proj4.Geocentric;
using Free.Ports.Proj4.Gridshift;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		// This table is intended to indicate for any given error code in
		// the range 0 to -49, whether that error will occur for all locations (ie.
		// it is a problem with the coordinate system as a whole) in which case the
		// value would be 0, or if the problem is with the point being transformed
		// in which case the value is 1.
		//
		// At some point we might want to move this array in with the error message
		// list or something, but while experimenting with it this should be fine.
		static readonly int[] transient_error={
		//	0  1  2  3  4  5  6  7  8  9
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0,	// 0 to 9
			0, 0, 0, 0, 1, 1, 0, 1, 1, 1,	// 10 to 19
			1, 0, 0, 0, 0, 0, 0, 1, 0, 0,	// 20 to 29
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0,	// 30 to 39
			0, 0, 0, 0, 0, 0, 0, 0, 1, 0};	// 40 to 49

		//**********************************************************************
		//							pj_transform()
		//
		//		Currently this function doesn't recognise if two projections
		//		are identical (to short circuit reprojection) because it is
		//		difficult to compare PJ structures (since there are some
		//		projection specific components).
		//**********************************************************************
		public static int pj_transform(PJ srcdefn, PJ dstdefn, ref double x, ref double y, ref double z)
		{
			double[] xa=new double[1];
			double[] ya=new double[1];
			double[] za=new double[1];

			xa[0]=x;
			ya[0]=y;
			za[0]=z;

			int ret=pj_transform(srcdefn, dstdefn, 1, 0, xa, ya, za);

			x=xa[0];
			y=ya[0];
			z=za[0];

			return ret;
		}

		public static int pj_transform(PJ srcdefn, PJ dstdefn, double[] x, double[] y, double[] z)
		{
			if(x.Length==y.Length)
			{
				if(x.Length==z.Length)
				{
					int ret=pj_transform(srcdefn, dstdefn, x.Length, 0, x, y, z);
					return ret;
				}
				else
				{
					throw new Exception("Coordinate Arrays had not the same length!");
				}
			}
			else
			{
				throw new Exception("Coordinate Arrays had not the same length!");
			}
		}

		public static int pj_transform(PJ srcdefn, PJ dstdefn, int point_count, int point_offset, double[] x, double[] y, double[] z)
		{
			srcdefn.ctx.last_errno=0;
			dstdefn.ctx.last_errno=0;

			if(point_offset==0) point_offset=1;

			// --------------------------------------------------------------------
			//		Transform unusual input coordinate axis orientation to
			//		standard form if needed.
			// --------------------------------------------------------------------
			if(srcdefn.axis!="enu")
			{
				int err=pj_adjust_axis(srcdefn.ctx, srcdefn.axis, false, point_count, point_offset, x, y, z);
				if(err!=0) return err;
			}

			// --------------------------------------------------------------------
			//		Transform Z to meters if it isn't already.
			// --------------------------------------------------------------------
			if(srcdefn.vto_meter!=1.0&&z!=null)
			{
				for(int i=0; i<point_count; i++) z[point_offset*i]*=srcdefn.vto_meter;
			}

			// --------------------------------------------------------------------
			//		Transform geocentric source coordinates to lat/long.
			// --------------------------------------------------------------------
			if(srcdefn.is_geocent)
			{
				if(z==null)
				{
					pj_ctx_set_errno(pj_get_ctx(srcdefn), (int)PJD_ERR.GEOCENTRIC);
					return (int)PJD_ERR.GEOCENTRIC;
				}

				if(srcdefn.to_meter!=1.0)
				{
					for(int i=0; i<point_count; i++)
					{
						if(x[point_offset*i]!=Libc.HUGE_VAL)
						{
							x[point_offset*i]*=srcdefn.to_meter;
							y[point_offset*i]*=srcdefn.to_meter;
						}
					}
				}

				int err=pj_geocentric_to_geodetic(srcdefn.a_orig, srcdefn.es_orig, point_count, point_offset, x, y, z);
				if(err!=0) return err;
			}
			// --------------------------------------------------------------------
			//		Transform source points to lat/long, if they aren't
			//		already.
			// --------------------------------------------------------------------
			else if(!srcdefn.is_latlong)
			{
				if(srcdefn.inv==null)
				{
					pj_ctx_set_errno(pj_get_ctx(srcdefn), -17);
#if DEBUG
					Console.Error.WriteLine("pj_transform(): source projection not invertable");
#endif

					return -17;
				}

				for(int i=0; i<point_count; i++)
				{
					XY projected_loc;
					LP geodetic_loc;

					projected_loc.x=x[point_offset*i];
					projected_loc.y=y[point_offset*i];

					if(projected_loc.x==Libc.HUGE_VAL) continue;

					geodetic_loc=pj_inv(projected_loc, srcdefn);
					if(srcdefn.ctx.last_errno!=0)
					{
						if((srcdefn.ctx.last_errno!=(int)ERRORNUMBER.EDOM&&srcdefn.ctx.last_errno!=(int)ERRORNUMBER.ERANGE)&&
							(srcdefn.ctx.last_errno>0||srcdefn.ctx.last_errno<-49||point_count==1||transient_error[-srcdefn.ctx.last_errno]==0))
							return srcdefn.ctx.last_errno;

						geodetic_loc.lam=Libc.HUGE_VAL;
						geodetic_loc.phi=Libc.HUGE_VAL;
					}

					x[point_offset*i]=geodetic_loc.lam;
					y[point_offset*i]=geodetic_loc.phi;
				}
			}

			// --------------------------------------------------------------------
			//		But if they are already lat long, adjust for the prime
			//		meridian if there is one in effect.
			// --------------------------------------------------------------------
			if(srcdefn.from_greenwich!=0.0)
			{
				for(int i=0; i<point_count; i++)
				{
					if(x[point_offset*i]!=Libc.HUGE_VAL) x[point_offset*i]+=srcdefn.from_greenwich;
				}
			}

			// --------------------------------------------------------------------
			//		Do we need to translate from geoid to ellipsoidal vertical
			//		datum?
			// --------------------------------------------------------------------
			if(srcdefn.has_geoid_vgrids&&z!=null)
			{
				if(Grid.pj_apply_vgridshift(srcdefn, "geoidgrids", ref srcdefn.vgridlist_geoid, ref srcdefn.vgridlist_geoid_count, false, point_count, point_offset, x, y, z)!=0)
					return pj_ctx_get_errno(srcdefn.ctx);
			}

			// --------------------------------------------------------------------
			//		Convert datums if needed, and possible.
			// --------------------------------------------------------------------
			if(pj_datum_transform(srcdefn, dstdefn, point_count, point_offset, x, y, z)!=0)
			{
				if(srcdefn.ctx.last_errno!=0) return srcdefn.ctx.last_errno;
				else return dstdefn.ctx.last_errno;
			}

			// --------------------------------------------------------------------
			//		Do we need to translate from geoid to ellipsoidal vertical
			//		datum?
			// --------------------------------------------------------------------
			if(dstdefn.has_geoid_vgrids&&z!=null)
			{
				if(Grid.pj_apply_vgridshift(dstdefn, "geoidgrids", ref dstdefn.vgridlist_geoid, ref dstdefn.vgridlist_geoid_count, true, point_count, point_offset, x, y, z)!=0)
					return dstdefn.ctx.last_errno;
			}

			// --------------------------------------------------------------------
			//		But if they are staying lat long, adjust for the prime
			//		meridian if there is one in effect.
			// --------------------------------------------------------------------
			if(dstdefn.from_greenwich!=0.0)
			{
				for(int i=0; i<point_count; i++)
				{
					if(x[point_offset*i]!=Libc.HUGE_VAL)
						x[point_offset*i]-=dstdefn.from_greenwich;
				}
			}

			// --------------------------------------------------------------------
			//		Transform destination latlong to geocentric if required.
			// --------------------------------------------------------------------
			if(dstdefn.is_geocent)
			{
				if(z==null)
				{
					pj_ctx_set_errno(dstdefn.ctx, (int)PJD_ERR.GEOCENTRIC);
					return (int)PJD_ERR.GEOCENTRIC;
				}

				pj_geodetic_to_geocentric(dstdefn.a_orig, dstdefn.es_orig, point_count, point_offset, x, y, z);

				if(dstdefn.fr_meter!=1.0)
				{
					for(int i=0; i<point_count; i++)
					{
						if(x[point_offset*i]!=Libc.HUGE_VAL)
						{
							x[point_offset*i]*=dstdefn.fr_meter;
							y[point_offset*i]*=dstdefn.fr_meter;
						}
					}
				}
			}
			// --------------------------------------------------------------------
			//		Transform destination points to projection coordinates, if
			//		desired.
			// --------------------------------------------------------------------
			else if(!dstdefn.is_latlong)
			{
				for(int i=0; i<point_count; i++)
				{
					XY projected_loc;
					LP geodetic_loc;

					geodetic_loc.lam=x[point_offset*i];
					geodetic_loc.phi=y[point_offset*i];

					if(geodetic_loc.lam==Libc.HUGE_VAL) continue;

					projected_loc=pj_fwd(geodetic_loc, dstdefn);
					if(dstdefn.ctx.last_errno!=0)
					{
						if((dstdefn.ctx.last_errno!=(int)ERRORNUMBER.EDOM&&dstdefn.ctx.last_errno!=(int)ERRORNUMBER.ERANGE)&&
							(dstdefn.ctx.last_errno>0||dstdefn.ctx.last_errno<-49||point_count==1||transient_error[-dstdefn.ctx.last_errno]==0))
							return dstdefn.ctx.last_errno;

						projected_loc.x=Libc.HUGE_VAL;
						projected_loc.y=Libc.HUGE_VAL;

					}

					x[point_offset*i]=projected_loc.x;
					y[point_offset*i]=projected_loc.y;
				}
			}

			// --------------------------------------------------------------------
			//		If a wrapping center other than 0 is provided, rewrap around
			//		the suggested center (for latlong coordinate systems only).
			// --------------------------------------------------------------------
			else if(dstdefn.is_latlong&&dstdefn.is_long_wrap_set)
			{
				for(int i=0; i<point_count; i++)
				{
					if(x[point_offset*i]==Libc.HUGE_VAL) continue;

					while(x[point_offset*i]<dstdefn.long_wrap_center-PI) x[point_offset*i]+=TWOPI;
					while(x[point_offset*i]>dstdefn.long_wrap_center+PI) x[point_offset*i]-=TWOPI;
				}
			}

			// --------------------------------------------------------------------
			//		Transform Z from meters if needed.
			// --------------------------------------------------------------------
			if(dstdefn.vto_meter!=1.0&&z!=null)
			{
				for(int i=0; i<point_count; i++) z[point_offset*i]*=dstdefn.vfr_meter;
			}

			// --------------------------------------------------------------------
			//		Transform normalized axes into unusual output coordinate axis
			//		orientation if needed.
			// --------------------------------------------------------------------
			if(dstdefn.axis!="enu")
			{

				int err=pj_adjust_axis(dstdefn.ctx, dstdefn.axis, true, point_count, point_offset, x, y, z);
				if(err!=0) return err;
			}

			return 0;
		}

		//**********************************************************************
		//						pj_geodetic_to_geocentric()
		//**********************************************************************
		public static int pj_geodetic_to_geocentric(double a, double es, int point_count, int point_offset, double[] x, double[] y, double[] z)
		{
			GeocentricInfo gi=new GeocentricInfo();
			int ret_errno=0;

			double b;
			if(es==0.0) b=a;
			else b=a*Math.Sqrt(1-es);

			if(gi.pj_Set_Geocentric_Parameters(a, b)!=0) return (int)PJD_ERR.GEOCENTRIC;

			for(int i=0; i<point_count; i++)
			{
				int io=i*point_offset;

				if(x[io]==Libc.HUGE_VAL) continue;

				if(gi.pj_Convert_Geodetic_To_Geocentric(y[io], x[io], z[io], out x[io], out y[io], out z[io])!=0)
				{
					ret_errno=-14;
					x[io]=y[io]=Libc.HUGE_VAL;
					// but keep processing points!
				}
			}

			return ret_errno;
		}

		//**********************************************************************
		//						pj_geodetic_to_geocentric()
		//**********************************************************************
		public static int pj_geocentric_to_geodetic(double a, double es, int point_count, int point_offset, double[] x, double[] y, double[] z)
		{
			GeocentricInfo gi=new GeocentricInfo();

			double b;
			if(es==0.0) b=a;
			else b=a*Math.Sqrt(1-es);

			if(gi.pj_Set_Geocentric_Parameters(a, b)!=0) return (int)PJD_ERR.GEOCENTRIC;

			for(int i=0; i<point_count; i++)
			{
				int io=i*point_offset;

				if(x[io]==Libc.HUGE_VAL) continue;

				gi.pj_Convert_Geocentric_To_Geodetic(x[io], y[io], z[io], out y[io], out x[io], out z[io]);
			}

			return 0;
		}

		//**********************************************************************
		//							pj_compare_datums()
		//
		//		Returns true if the two datums are identical, otherwise false
		//**********************************************************************
		public static bool pj_compare_datums(PJ srcdefn, PJ dstdefn)
		{
			if(srcdefn.datum_type!=dstdefn.datum_type) return false;

			// the tolerence for es is to ensure that GRS80 and WGS84 are considered identical
			if(srcdefn.a_orig!=dstdefn.a_orig||Math.Abs(srcdefn.es_orig-dstdefn.es_orig)>0.000000000050) return false;

			if(srcdefn.datum_type==PJD._3PARAM)
				return (srcdefn.datum_params[0]==dstdefn.datum_params[0]&&
						srcdefn.datum_params[1]==dstdefn.datum_params[1]&&
						srcdefn.datum_params[2]==dstdefn.datum_params[2]);

			if(srcdefn.datum_type==PJD._7PARAM)
				return (srcdefn.datum_params[0]==dstdefn.datum_params[0]&&
						srcdefn.datum_params[1]==dstdefn.datum_params[1]&&
						srcdefn.datum_params[2]==dstdefn.datum_params[2]&&
						srcdefn.datum_params[3]==dstdefn.datum_params[3]&&
						srcdefn.datum_params[4]==dstdefn.datum_params[4]&&
						srcdefn.datum_params[5]==dstdefn.datum_params[5]&&
						srcdefn.datum_params[6]==dstdefn.datum_params[6]);

			if(srcdefn.datum_type==PJD.GRIDSHIFT)
				return pj_param_s(srcdefn.ctx, srcdefn.parameters, "nadgrids")==pj_param_s(dstdefn.ctx, dstdefn.parameters, "nadgrids");

			return true;
		}

		//**********************************************************************
		//						pj_geocentic_to_wgs84()
		//**********************************************************************
		static int pj_geocentric_to_wgs84(PJ defn, int point_count, int point_offset, double[] x, double[] y, double[] z)
		{
			if(defn.datum_type==PJD._3PARAM)
			{
				for(int i=0; i<point_count; i++)
				{
					int io=i*point_offset;

					if(x[io]==Libc.HUGE_VAL) continue;

					x[io]=x[io]+defn.datum_params[0];
					y[io]=y[io]+defn.datum_params[1];
					z[io]=z[io]+defn.datum_params[2];
				}
			}
			else if(defn.datum_type==PJD._7PARAM)
			{
				for(int i=0; i<point_count; i++)
				{
					int io=i*point_offset;
					double x_out, y_out, z_out;

					if(x[io]==Libc.HUGE_VAL) continue;

					x_out=defn.datum_params[6]*(x[io]-defn.datum_params[5]*y[io]+defn.datum_params[4]*z[io])+defn.datum_params[0];
					y_out=defn.datum_params[6]*(defn.datum_params[5]*x[io]+y[io]-defn.datum_params[3]*z[io])+defn.datum_params[1];
					z_out=defn.datum_params[6]*(-defn.datum_params[4]*x[io]+defn.datum_params[3]*y[io]+z[io])+defn.datum_params[2];

					x[io]=x_out;
					y[io]=y_out;
					z[io]=z_out;
				}
			}

			return 0;
		}

		//**********************************************************************
		//						pj_geocentic_from_wgs84()
		//**********************************************************************
		static int pj_geocentric_from_wgs84(PJ defn, int point_count, int point_offset, double[] x, double[] y, double[] z)
		{
			if(defn.datum_type==PJD._3PARAM)
			{
				for(int i=0; i<point_count; i++)
				{
					long io=i*point_offset;

					if(x[io]==Libc.HUGE_VAL) continue;

					x[io]=x[io]-defn.datum_params[0];
					y[io]=y[io]-defn.datum_params[1];
					z[io]=z[io]-defn.datum_params[2];
				}
			}
			else if(defn.datum_type==PJD._7PARAM)
			{
				for(int i=0; i<point_count; i++)
				{
					long io=i*point_offset;
					double x_tmp, y_tmp, z_tmp;

					if(x[io]==Libc.HUGE_VAL) continue;

					x_tmp=(x[io]-defn.datum_params[0])/defn.datum_params[6];
					y_tmp=(y[io]-defn.datum_params[1])/defn.datum_params[6];
					z_tmp=(z[io]-defn.datum_params[2])/defn.datum_params[6];

					x[io]=x_tmp+defn.datum_params[5]*y_tmp-defn.datum_params[4]*z_tmp;
					y[io]=-defn.datum_params[5]*x_tmp+y_tmp+defn.datum_params[3]*z_tmp;
					z[io]=defn.datum_params[4]*x_tmp-defn.datum_params[3]*y_tmp+z_tmp;
				}
			}

			return 0;
		}

		//**********************************************************************
		//							pj_datum_transform()
		//
		//		The input should be long/lat/z coordinates in radians in the
		//		source datum, and the output should be long/lat/z
		//		coordinates in radians in the destination datum.
		//**********************************************************************
		public static int pj_datum_transform(PJ srcdefn, PJ dstdefn, int point_count, int point_offset, double[] x, double[] y, double[] z)
		{
			const double SRS_WGS84_SEMIMAJOR=6378137.0;
			const double SRS_WGS84_ESQUARED=0.0066943799901413165;

			// --------------------------------------------------------------------
			//		We cannot do any meaningful datum transformation if either
			//		the source or destination are of an unknown datum type
			//		(ie. only a +ellps declaration, no +datum). This is new
			//		behavior for PROJ 4.6.0.
			// --------------------------------------------------------------------
			if(srcdefn.datum_type==PJD.UNKNOWN||dstdefn.datum_type==PJD.UNKNOWN) return 0;

			// --------------------------------------------------------------------
			//		Short cut if the datums are identical.
			// --------------------------------------------------------------------
			if(pj_compare_datums(srcdefn, dstdefn)) return 0;

			double src_a=srcdefn.a_orig;
			double src_es=srcdefn.es_orig;

			double dst_a=dstdefn.a_orig;
			double dst_es=dstdefn.es_orig;

			// --------------------------------------------------------------------
			//		Create a temporary Z array if one is not provided.
			// --------------------------------------------------------------------
			if(z==null) z=new double[point_count*point_offset];

			// --------------------------------------------------------------------
			//		If this datum requires grid shifts, then apply it to geodetic
			//		coordinates.
			// --------------------------------------------------------------------
			if(srcdefn.datum_type==PJD.GRIDSHIFT)
			{
				Grid.pj_apply_gridshift_2(srcdefn, false, point_count, point_offset, x, y, z);
				if(srcdefn.ctx.last_errno!=0&&(srcdefn.ctx.last_errno>0||transient_error[-srcdefn.ctx.last_errno]==0)) return srcdefn.ctx.last_errno;

				src_a=SRS_WGS84_SEMIMAJOR;
				src_es=SRS_WGS84_ESQUARED;
			}

			if(dstdefn.datum_type==PJD.GRIDSHIFT)
			{
				dst_a=SRS_WGS84_SEMIMAJOR;
				dst_es=SRS_WGS84_ESQUARED;
			}

			// ====================================================================
			//		Do we need to go through geocentric coordinates?
			// ====================================================================
			if(src_es!=dst_es||src_a!=dst_a||srcdefn.datum_type==PJD._3PARAM||srcdefn.datum_type==PJD._7PARAM||
				dstdefn.datum_type==PJD._3PARAM||dstdefn.datum_type==PJD._7PARAM)
			{
				// --------------------------------------------------------------------
				//		Convert to geocentric coordinates.
				// --------------------------------------------------------------------
				srcdefn.ctx.last_errno=pj_geodetic_to_geocentric(src_a, src_es, point_count, point_offset, x, y, z);
				if(srcdefn.ctx.last_errno!=0&&(srcdefn.ctx.last_errno>0||transient_error[-srcdefn.ctx.last_errno]==0)) return srcdefn.ctx.last_errno;

				// --------------------------------------------------------------------
				//		Convert between datums.
				// --------------------------------------------------------------------
				if(srcdefn.datum_type==PJD._3PARAM||srcdefn.datum_type==PJD._7PARAM)
				{
					pj_geocentric_to_wgs84(srcdefn, point_count, point_offset, x, y, z);
					if(srcdefn.ctx.last_errno!=0&&(srcdefn.ctx.last_errno>0||transient_error[-srcdefn.ctx.last_errno]==0)) return srcdefn.ctx.last_errno;
				}

				if(dstdefn.datum_type==PJD._3PARAM||dstdefn.datum_type==PJD._7PARAM)
				{
					pj_geocentric_from_wgs84(dstdefn, point_count, point_offset, x, y, z);
					if(dstdefn.ctx.last_errno!=0&&(dstdefn.ctx.last_errno>0||transient_error[-dstdefn.ctx.last_errno]==0)) return dstdefn.ctx.last_errno;
				}

				// --------------------------------------------------------------------
				//		Convert back to geodetic coordinates.
				// --------------------------------------------------------------------
				dstdefn.ctx.last_errno=pj_geocentric_to_geodetic(dst_a, dst_es, point_count, point_offset, x, y, z);
				if(dstdefn.ctx.last_errno!=0&&(dstdefn.ctx.last_errno>0||transient_error[-dstdefn.ctx.last_errno]==0)) return dstdefn.ctx.last_errno;
			}

			// --------------------------------------------------------------------
			//		Apply grid shift to destination if required.
			// --------------------------------------------------------------------
			if(dstdefn.datum_type==PJD.GRIDSHIFT)
			{
				Grid.pj_apply_gridshift_2(dstdefn, true, point_count, point_offset, x, y, z);
				if(dstdefn.ctx.last_errno!=0&&(dstdefn.ctx.last_errno>0||transient_error[-dstdefn.ctx.last_errno]==0)) return dstdefn.ctx.last_errno;
			}

			return 0;
		}

		//***********************************************************************
		//*							pj_adjust_axis()							*
		//*																		*
		//*		Normalize or de-normalized the x/y/z axes. The normal form		*
		//*		is "enu" (easting, northing, up).								*
		//***********************************************************************
		public static int pj_adjust_axis(projCtx ctx, string axis, bool denormalize_flag, long point_count, int point_offset, double[] x, double[] y, double[] z)
		{
			double x_in, y_in, z_in=0.0;
			int i, i_axis;

			if(!denormalize_flag)
			{
				for(i=0; i<point_count; i++)
				{
					x_in=x[point_offset*i];
					y_in=y[point_offset*i];
					if(z!=null) z_in=z[point_offset*i];

					for(i_axis=0; i_axis<3; i_axis++)
					{
						double value;

						if(i_axis==0) value=x_in;
						else if(i_axis==1) value=y_in;
						else value=z_in;

						switch(axis[i_axis])
						{
							case 'e': x[point_offset*i]=value; break;
							case 'w': x[point_offset*i]=-value; break;
							case 'n': y[point_offset*i]=value; break;
							case 's': y[point_offset*i]=-value; break;
							case 'u': if(z!=null) z[point_offset*i]=value; break;
							case 'd': if(z!=null) z[point_offset*i]=-value; break;
							default: pj_ctx_set_errno(ctx, (int)PJD_ERR.AXIS); return (int)PJD_ERR.AXIS;
						}
					} // i_axis
				} // i (point)
			}
			else // denormalize
			{
				for(i=0; i<point_count; i++)
				{
					x_in=x[point_offset*i];
					y_in=y[point_offset*i];
					if(z!=null) z_in=z[point_offset*i];

					for(i_axis=0; i_axis<3; i_axis++)
					{
						double[] target;

						if(i_axis==2&&z==null) continue;

						if(i_axis==0) target=x;
						else if(i_axis==1) target=y;
						else target=z;

						switch(axis[i_axis])
						{
							case 'e': target[point_offset*i]=x_in; break;
							case 'w': target[point_offset*i]=-x_in; break;
							case 'n': target[point_offset*i]=y_in; break;
							case 's': target[point_offset*i]=-y_in; break;
							case 'u': target[point_offset*i]=z_in; break;
							case 'd': target[point_offset*i]=-z_in; break;
							default: pj_ctx_set_errno(ctx, (int)PJD_ERR.AXIS); return (int)PJD_ERR.AXIS;
						}
					} // i_axis
				} // i (point)
			}

			return 0;
		}
	}
}
