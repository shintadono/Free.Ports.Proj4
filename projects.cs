//*****************************************************************************
//
// Project:	PROJ.4
// Purpose:	Primary (private) include file for PROJ.4 library.
// Author:	Gerald Evenden
//
//*****************************************************************************
// Copyright (c) 2000, Frank Warmerdam
// Copyright (c) 2008-2011 by the Authors
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
using System.Globalization;
using System.Text;
using Free.Ports.Proj4.Gridshift;

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		// some useful constants
		internal const double HALFPI = 1.5707963267948966;
		internal const double FORTPI = 0.78539816339744833;
		internal const double PI = 3.14159265358979323846;
		internal const double TWOPI = 6.2831853071795864769;

		// environment parameter name
		const string PROJ_LIB = "PROJ_LIB";

		// directory delimiter for DOS support
		const char DIR_CHAR = '\\';
		//const char DIR_CHAR='/';
	}

	// proj thread context
	public class projCtx
	{
		public delegate void Logger(object appdata, PJ_LOG level, string msg);

		public int last_errno;
		public PJ_LOG debug_level;
		public Logger logger;
		public object app_data;
	}

	// datum_type values
	public enum PJD
	{
		UNKNOWN = 0,
		_3PARAM = 1,
		_7PARAM = 2,
		GRIDSHIFT = 3,
		WGS84 = 4,		// WGS84 (or anything considered equivelent)
	}

	enum PJD_ERR
	{
		// library errors
		GEOCENTRIC = -45,
		AXIS = -47,
		GRID_AREA = -48,
		CATALOG = -49,
	}

	public struct projUV
	{
		public double u, v;

		public static implicit operator XY(projUV a)
		{
			XY ret; ret.x = a.u; ret.y = a.v;
			return ret;
		}

		public static implicit operator LP(projUV a)
		{
			LP ret; ret.lam = a.u; ret.phi = a.v;
			return ret;
		}
	}

	public struct projUVW
	{
		public double u, v, w;

		public static implicit operator XYZ(projUVW a)
		{
			XYZ ret; ret.x = a.u; ret.y = a.v; ret.z = a.w;
			return ret;
		}

		public static implicit operator LPZ(projUVW a)
		{
			LPZ ret; ret.lam = a.u; ret.phi = a.v; ret.z = a.w;
			return ret;
		}
	}

	public struct XY
	{
		public double x, y;

		public static implicit operator projUV(XY a)
		{
			projUV ret; ret.u = a.x; ret.v = a.y;
			return ret;
		}
	}

	public struct XYZ
	{
		public double x, y, z;

		public static implicit operator projUVW(XYZ a)
		{
			projUVW ret; ret.u = a.x; ret.v = a.y; ret.w = a.z;
			return ret;
		}
	}

	public struct LP
	{
		/// <summary>
		/// Longitude
		/// </summary>
		public double lam;

		/// <summary>
		/// Latitude
		/// </summary>
		public double phi;

		public static implicit operator projUV(LP a)
		{
			projUV ret; ret.u = a.lam; ret.v = a.phi;
			return ret;
		}
	}

	public struct LPZ
	{
		/// <summary>
		/// Longitude
		/// </summary>
		public double lam;

		/// <summary>
		/// Latitude
		/// </summary>
		public double phi;

		public double z;

		public static implicit operator projUVW(LPZ a)
		{
			projUVW ret; ret.u = a.lam; ret.v = a.phi; ret.w = a.z;
			return ret;
		}
	}

	// base projection data structure
	public abstract class PJ
	{
		internal static CultureInfo nc = new CultureInfo("");

		public projCtx ctx;
		public Func<LP, XY> fwd;
		public Func<XY, LP> inv;
		public Func<LPZ, XYZ> fwd3d;
		public Func<XYZ, LPZ> inv3d;

		public abstract string Name { get; }
		public abstract string DescriptionName { get; }
		public abstract string DescriptionType { get; }
		public abstract string DescriptionParameters { get; }
		public abstract bool Invertible { get; }

		protected virtual string Proj4ParameterString { get { return ""; } }
		public virtual string ToProj4String()
		{
			StringBuilder ret = new StringBuilder();
			ret.Append("+proj=").Append(Name);

			ret.AppendFormat(nc, " +lon_0={0}", lam0 * Proj.RAD_TO_DEG);
			ret.AppendFormat(nc, " +lat_0={0}", phi0 * Proj.RAD_TO_DEG);
			ret.AppendFormat(nc, " +x_0={0}", x0);
			ret.AppendFormat(nc, " +y_0={0}", y0);
			if (k0 != 1.0) ret.AppendFormat(nc, " +k_0={0}", k0);

			ret.Append(Proj4ParameterString);

			ret.AppendFormat(nc, " +a={0}", a_orig);
			ret.AppendFormat(nc, " +es={0}", es_orig);

			const double SEC_TO_RAD = 4.84813681109535993589914102357e-6;

			switch (datum_type)
			{
				default:
				case PJD.UNKNOWN: break;
				case PJD.WGS84:
					ret.Append(" +datum=WGS84");
					break;
				case PJD._7PARAM:
					ret.AppendFormat(nc, " +towgs84={0},{1},{2},{3},{4},{5},{6}",
						datum_params[0], datum_params[1], datum_params[2],
						datum_params[3] / SEC_TO_RAD, datum_params[4] / SEC_TO_RAD, datum_params[5] / SEC_TO_RAD,
						(datum_params[6] - 1) * 1000000);
					break;
				case PJD._3PARAM:
					ret.AppendFormat(nc, " +towgs84={0},{1},{2}", datum_params[0], datum_params[1], datum_params[2]);
					break;
				case PJD.GRIDSHIFT:
					foreach (string p in parameters)
					{
						if (p.Length > 9 && p.StartsWith("nadgrids="))
						{
							ret.Append(" +nadgrids=").Append(p.Substring(9));
							break;
						}
					}
					break;
			}

			if (to_meter != 1.0) ret.AppendFormat(nc, " +to_meter={0}", to_meter);
			if (vto_meter != 1.0) ret.AppendFormat(nc, " +vto_meter={0}", vto_meter);
			if (from_greenwich != 0.0) ret.AppendFormat(nc, " +pm={0}", from_greenwich);

			if (is_long_wrap_set) ret.AppendFormat(nc, " +lon_wrap={0}", long_wrap_center * Proj.RAD_TO_DEG);
			if (axis != "enu") ret.Append(" +axis=").Append(axis);

			if (geoc) ret.Append(" +geoc");
			if (over) ret.Append(" +over");

			return ret.ToString();
		}

		internal List<string> parameters;	// parameter list
		public bool over;				// over-range flag
		public bool geoc;				// geocentric latitude flag
		public bool is_latlong;			// proj=latlong ... not really a projection at all
		public bool is_geocent;			// proj=geocent ... not really a projection at all
		public double
			a,					// major axis or radius if es==0
			a_orig,				// major axis before any +proj related adjustment
			es,					// e^2
			es_orig,			// es before any +proj related adjustment
			e,					// eccentricity
			ra,					// 1/A
			one_es,				// 1 - e^2
			rone_es,			// 1/one_es
			lam0, phi0,			// central longitude, latitude
			x0, y0,				// easting and northing
			k0,					// general scaling factor
			to_meter, fr_meter;	// cartesian scaling

		public PJD datum_type;			// PJD.UNKNOWN/_3PARAM/_7PARAM/GRIDSHIFT/WGS84
		public double[] datum_params = new double[7];
		public PJ_GRIDINFO[] gridlist;
		public int gridlist_count;
		public bool has_geoid_vgrids;
		public PJ_GRIDINFO[] vgridlist_geoid;
		public int vgridlist_geoid_count;
		public double vto_meter, vfr_meter;
		public double from_greenwich;	// prime meridian offset (in radians)
		public double long_wrap_center;	// 0.0 for -180 to 180, actually in radians
		public bool is_long_wrap_set;
		public string axis;

		// New Datum Shift Grid Catalogs
		public string catalog_name;
		public PJ_GridCatalog catalog;

		public double datum_date;

		public PJ_GRIDINFO last_before_grid;
		public PJ_Region last_before_region;
		public double last_before_date;

		public PJ_GRIDINFO last_after_grid;
		public PJ_Region last_after_region;
		public double last_after_date;

		public abstract PJ Init();
	}
}
