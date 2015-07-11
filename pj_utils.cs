//*****************************************************************************
//
// Project:	PROJ.4
// Purpose:	Some utility functions we don't want to bother putting in
//			their own source files.
// Author:	Frank Warmerdam, warmerdam@pobox.com
//
//*****************************************************************************
// Copyright (c) 2001, Frank Warmerdam
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

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		//**********************************************************************
		//							pj_is_latlong()
		//
		//		Returns TRUE if this coordinate system object is geographic.
		//**********************************************************************
		public static bool pj_is_latlong(PJ pj)
		{
			return pj==null||pj.is_latlong;
		}

		//**********************************************************************
		//							pj_is_geocent()
		//
		//		Returns TRUE if this coordinate system object is geocentric.
		//**********************************************************************
		public static bool pj_is_geocent(PJ pj)
		{
			return pj!=null&&pj.is_geocent;
		}

		//**********************************************************************
		//							pj_latlong_from_proj()
		//
		//		Return a PJ* definition defining the lat/long coordinate
		//		system on which a projection is based. If the coordinate
		//		system passed in is latlong, a clone of the same will be
		//		returned.
		//**********************************************************************
		public static PJ pj_latlong_from_proj(PJ pj_in)
		{
			string defn;
			bool got_datum=false;

			pj_errno=0;
			defn="+proj=latlong";

			if(pj_param_t(pj_in.ctx, pj_in.parameters, "datum"))
			{
				got_datum=true;
				defn+=" +datum="+pj_param_s(pj_in.ctx, pj_in.parameters, "datum");
			}
			else if(pj_param_t(pj_in.ctx, pj_in.parameters, "ellps"))
			{
				defn+=" +ellps="+pj_param_s(pj_in.ctx, pj_in.parameters, "ellps");
			}
			else if(pj_param_t(pj_in.ctx, pj_in.parameters, "a"))
			{
				defn+=" +a="+pj_param_s(pj_in.ctx, pj_in.parameters, "a");

				if(pj_param_t(pj_in.ctx, pj_in.parameters, "b")) defn+=" +b="+pj_param_s(pj_in.ctx, pj_in.parameters, "b");
				else if(pj_param_t(pj_in.ctx, pj_in.parameters, "es")) defn+=" +es="+pj_param_s(pj_in.ctx, pj_in.parameters, "es");
				else if(pj_param_t(pj_in.ctx, pj_in.parameters, "f")) defn+=" +f="+pj_param_s(pj_in.ctx, pj_in.parameters, "f");
				else defn+=" +es="+pj_in.es.ToString("G16", nc);
			}
			else
			{
				pj_ctx_set_errno(pj_in.ctx, -13);
				return null;
			}

			if(!got_datum)
			{
				if(pj_param_t(pj_in.ctx, pj_in.parameters, "towgs84")) defn+=" +towgs84="+pj_param_s(pj_in.ctx, pj_in.parameters, "towgs84");
				if(pj_param_t(pj_in.ctx, pj_in.parameters, "nadgrids")) defn+=" +nadgrids="+pj_param_s(pj_in.ctx, pj_in.parameters, "nadgrids");
			}

			// copy over some other information related to ellipsoid
			if(pj_param_t(pj_in.ctx, pj_in.parameters, "R")) defn+=" +R="+pj_param_s(pj_in.ctx, pj_in.parameters, "R");

			if(pj_param_t(pj_in.ctx, pj_in.parameters, "R_A")) defn+=" +R_A";
			if(pj_param_t(pj_in.ctx, pj_in.parameters, "R_V")) defn+=" +R_V";
			if(pj_param_t(pj_in.ctx, pj_in.parameters, "R_a")) defn+=" +R_a";
			if(pj_param_t(pj_in.ctx, pj_in.parameters, "R_lat_a")) defn+=" +R_lat_a="+pj_param_s(pj_in.ctx, pj_in.parameters, "R_lat_a");
			if(pj_param_t(pj_in.ctx, pj_in.parameters, "R_lat_g")) defn+=" +R_lat_g="+pj_param_s(pj_in.ctx, pj_in.parameters, "R_lat_g");

			// copy over prime meridian
			if(pj_param_t(pj_in.ctx, pj_in.parameters, "pm")) defn+=" +pm="+pj_param_s(pj_in.ctx, pj_in.parameters, "pm");

			return pj_init_plus_ctx(pj_in.ctx, defn);
		}

		//**********************************************************************
		//						pj_get_spheroid_defn()
		//
		//		Fetch the internal definition of the spheroid. Note that
		//		you can compute "b" from eccentricity_squared as:
		//
		//		b=a*sqrt(1-es)
		//**********************************************************************
		public static void pj_get_spheroid_defn(PJ defn, out double major_axis, out double eccentricity_squared)
		{
			major_axis=defn.a;
			eccentricity_squared=defn.es;
		}
	}
}
