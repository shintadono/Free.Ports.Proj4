//*****************************************************************************
//
// Project:  SCH Coordinate system
// Purpose:  Implementation of SCH Coordinate system
// References :
//      1. Hensley. Scott. SCH Coordinates and various transformations. June 15, 2000.
//      2. Buckley, Sean Monroe. Radar interferometry measurement of land subsidence. 2000..
//         PhD Thesis. UT Austin. (Appendix)
//      3. Hensley, Scott, Elaine Chapin, and T. Michel. "Improved processing of AIRSAR
//         data based on the GeoSAR processor." Airsar earth science and applications
//         workshop, March. 2002. (http://airsar.jpl.nasa.gov/documents/workshop2002/papers/T3.pdf)
//
// Author:   Piyush Agram (piyush.agram@jpl.nasa.gov)
// Copyright (c) 2015 California Institute of Technology.
// Government sponsorship acknowledged.
// Copyright (c) 2015 by the Authors
//
// NOTE:  The SCH coordinate system is a sensor aligned coordinate system
// developed at JPL for radar mapping missions. Details pertaining to the
// coordinate system have been release in the public domain (see references above).
// This code is an independent implementation of the SCH coordinate system
// that conforms to the PROJ.4 conventions and uses the details presented in these
// publicly released documents. All credit for the development of the coordinate
// system and its use should be directed towards the original developers at JPL.
//*****************************************************************************
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//***************************************************************************/

using System;
using System.Text;
using Free.Ports.Proj4.Geocentric;

namespace Free.Ports.Proj4
{
	class PJ_sch : PJ
	{
		protected double plat; // Peg Latitude
		protected double plon; // Peg Longitude
		protected double phdg; // Peg heading
		protected double h0;   // Average altitude
		protected double[] transMat = new double[9];
		protected double[] xyzoff = new double[3];
		protected double rcurv;
		protected GeocentricInfo sph=new GeocentricInfo();
		protected GeocentricInfo elp_0=new GeocentricInfo();

		public override string Name { get { return "sch"; } }
		public override string DescriptionName { get { return "Spherical Cross-track Height"; } }
		public override string DescriptionType { get { return "tMisc"; } }
		public override string DescriptionParameters { get { return "plat_0= plon_0= phdg_0= [h_0=]"; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret = new StringBuilder();
				ret.AppendFormat(nc, " +plat_0={0}", plat);
				ret.AppendFormat(nc, " +plon_0={0}", plon);
				ret.AppendFormat(nc, " +phdg_0={0}", phdg);
				if(h0!=0.0) ret.AppendFormat(nc, " +h_0={0}", h0);
				return ret.ToString();
			}
		}

		LPZ inverse3d(XYZ xyz)
		{
			LPZ lpz;
			lpz.lam = lpz.phi = lpz.z = 0;

			// Local lat, lon using radius
			double pxyz0 = xyz.y * a / rcurv;
			double pxyz1 = xyz.x * a / rcurv;
			double pxyz2 = xyz.z;

			double temp0, temp1, temp2;
			if (sph.pj_Convert_Geodetic_To_Geocentric(pxyz0, pxyz1, pxyz2, out temp0, out temp1, out temp2) != 0) { Proj.pj_ctx_set_errno(ctx, -20); return lpz; }

			// Apply rotation
			pxyz0 = transMat[0] * temp0 + transMat[1] * temp1 + transMat[2] * temp2;
			pxyz1 = transMat[3] * temp0 + transMat[4] * temp1 + transMat[5] * temp2;
			pxyz2 = transMat[6] * temp0 + transMat[7] * temp1 + transMat[8] * temp2;

			// Apply offset
			pxyz0 += xyzoff[0];
			pxyz1 += xyzoff[1];
			pxyz2 += xyzoff[2];

			// Convert geocentric coordinates to lat lon
			elp_0.pj_Convert_Geocentric_To_Geodetic(pxyz0, pxyz1, pxyz2, out temp0, out temp1, out temp2);

			lpz.lam = temp1;
			lpz.phi = temp0;
			lpz.z = temp2;

			//printf("INVERSE: \n");
			//printf("XYZ: %f %f %f \n", xyz.x, xyz.y, xyz.z);
			//printf("LPZ: %f %f %f \n", lpz.lam, lpz.phi, lpz.z);
			return lpz;
		}

		XYZ forward3d(LPZ lpz)
		{
			XYZ xyz;
			xyz.x = xyz.y = xyz.z = 0;

			double temp0, temp1, temp2;

			// Convert lat lon to geocentric coordinates
			if (elp_0.pj_Convert_Geodetic_To_Geocentric(lpz.phi, lpz.lam, lpz.z, out temp0, out temp1, out temp2) != 0) { Proj.pj_ctx_set_errno(ctx, -20); return xyz; }

			// Adjust for offset
			temp0 -= xyzoff[0];
			temp1 -= xyzoff[1];
			temp2 -= xyzoff[2];

			// Apply rotation
			double pxyz0 = transMat[0] * temp0 + transMat[3] * temp1 + transMat[6] * temp2;
			double pxyz1 = transMat[1] * temp0 + transMat[4] * temp1 + transMat[7] * temp2;
			double pxyz2 = transMat[2] * temp0 + transMat[5] * temp1 + transMat[8] * temp2;

			// Convert to local lat,lon
			sph.pj_Convert_Geocentric_To_Geodetic(pxyz0, pxyz1, pxyz2, out temp0, out temp1, out temp2);

			// Scale by radius
			xyz.x = temp1 * rcurv / a;
			xyz.y = temp0 * rcurv / a;
			xyz.z = temp2;

			//printf("FORWARD: \n");
			//printf("LPZ: %f %f %f \n", lpz.lam, lpz.phi, lpz.z);
			//printf("XYZ: %f %f %f \n", xyz.x, xyz.y, xyz.z);
			return xyz;
		}

		// general initialization
		PJ setup()
		{
			double temp = a * Math.Sqrt(1.0 - es);

			// Setup original geocentric system
			if (elp_0.pj_Set_Geocentric_Parameters(a, temp) != 0) { Proj.pj_ctx_set_errno(ctx, -37); return null; }

			double clt = Math.Cos(plat);
			double slt = Math.Sin(plat);
			double clo = Math.Cos(plon);
			double slo = Math.Sin(plon);

			// Estimate the radius of curvature for given peg
			temp = Math.Sqrt(1.0 - (es) * slt * slt);
			double reast = a / temp;
			double rnorth = a * (1.0 - (es)) / Math.Pow(temp, 3);

			double chdg = Math.Cos(phdg);
			double shdg = Math.Sin(phdg);

			rcurv = h0 + (reast * rnorth) / (reast * chdg * chdg + rnorth * shdg * shdg);

			//printf("North Radius: %f \n", rnorth);
			//printf("East Radius: %f \n", reast);
			//printf("Effective Radius: %f \n", rcurv);

			// Set up local sphere at the given peg point
			if (sph.pj_Set_Geocentric_Parameters(rcurv, rcurv) != 0) { Proj.pj_ctx_set_errno(ctx, -37); return null; }

			// Set up the transformation matrices
			transMat[0] = clt * clo;
			transMat[1] = -shdg * slo - slt * clo * chdg;
			transMat[2] = slo * chdg - slt * clo * shdg;
			transMat[3] = clt * slo;
			transMat[4] = clo * shdg - slt * slo * chdg;
			transMat[5] = -clo * chdg - slt * slo * shdg;
			transMat[6] = slt;
			transMat[7] = clt * chdg;
			transMat[8] = clt * shdg;

			double pxyz0, pxyz1, pxyz2;

			if (elp_0.pj_Convert_Geodetic_To_Geocentric(plat, plon, h0, out pxyz0, out pxyz1, out pxyz2) != 0)
			{
				Proj.pj_ctx_set_errno(ctx, -14); return null;
			}

			xyzoff[0] = pxyz0 - (rcurv) * clt * clo;
			xyzoff[1] = pxyz1 - (rcurv) * clt * slo;
			xyzoff[2] = pxyz2 - (rcurv) * slt;

			//printf("Offset: %f %f %f \n", xyzoff[0], xyzoff[1], xyzoff[2]);

			fwd3d = forward3d;
			inv3d = inverse3d;

			return this;
		}

		public override PJ Init()
		{
			h0 = 0.0;

			// Check if peg latitude was defined
			if (Proj.pj_param_t(ctx, parameters, "plat_0")) plat = Proj.pj_param_r(ctx, parameters, "plat_0");
			else { Proj.pj_ctx_set_errno(ctx, -37); return null; }

			// Check if peg longitude was defined
			if (Proj.pj_param_t(ctx, parameters, "plon_0")) plon = Proj.pj_param_r(ctx, parameters, "plon_0");
			else { Proj.pj_ctx_set_errno(ctx, -37); return null; }

			// Check if peg latitude is defined
			if (Proj.pj_param_t(ctx, parameters, "phdg_0")) phdg = Proj.pj_param_r(ctx, parameters, "phdg_0");
			else { Proj.pj_ctx_set_errno(ctx, -37); return null; }

			// Check if average height was defined
			// If so read it in
			if (Proj.pj_param_t(ctx, parameters, "h_0")) h0 = Proj.pj_param_d(ctx, parameters, "h_0");

			// Completed reading in the projection parameters
			//printf("PSA: Lat = %f Lon = %f Hdg = %f \n", plat, plon, phdg);
			return setup();
		}
	}
}
