using System;
using System.Text;

namespace Free.Ports.Proj4.Projections
{
	class PJ_ocea : PJ
	{
		protected double rok, rtk, sinphi, cosphi, singam, cosgam;

		public override string Name { get { return "ocea"; } }
		public override string DescriptionName { get { return "Oblique Cylindrical Equal Area"; } }
		public override string DescriptionType { get { return "Cyl, Sph"; } }
		public override string DescriptionParameters { get { return "lonc= alpha= or lat_1= lat_2= lon_1= lon_2="; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret = new StringBuilder();
				ret.AppendFormat(nc, " +lonc={0}", Math.Atan2(singam, cosgam) * Proj.RAD_TO_DEG - 90.0);
				ret.AppendFormat(nc, " +alpha={0}", Math.Asin(sinphi) * Proj.RAD_TO_DEG);
				return ret.ToString();
			}
		}

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x = xy.y = 0;

			xy.y = Math.Sin(lp.lam);
			//xy.x=Math.Atan2((Math.Tan(lp.phi)*cosphi+sinphi*xy.y), Math.Cos(lp.lam));

			double t = Math.Cos(lp.lam);
			xy.x = Math.Atan((Math.Tan(lp.phi) * cosphi + sinphi * xy.y) / t);
			if (t < 0.0) xy.x += Proj.PI;
			xy.x *= rtk;
			xy.y = rok * (sinphi * Math.Sin(lp.phi) - cosphi * Math.Cos(lp.phi) * xy.y);

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam = lp.phi = 0;

			xy.y /= rok;
			xy.x /= rtk;
			double t = Math.Sqrt(1.0 - xy.y * xy.y);
			double s = Math.Sin(xy.x);
			lp.phi = Math.Asin(xy.y * sinphi + t * cosphi * s);
			lp.lam = Math.Atan2(t * sinphi * s - xy.y * cosphi, t * Math.Cos(xy.x));

			return lp;
		}

		public override PJ Init()
		{
			rok = a / k0;
			rtk = a * k0;
			// If the keyword "alpha" is found in the sentence then use 1point+1azimuth
			if (Proj.pj_param_t(ctx, parameters, "alpha"))
			{
				double phi_0 = 0.0;
				// Define Pole of oblique transformation from 1 point & 1 azimuth
				double alpha = Proj.pj_param_r(ctx, parameters, "alpha");
				double lonz = Proj.pj_param_r(ctx, parameters, "lonc");
				// Equation 9-8 page 80 (http://pubs.usgs.gov/pp/1395/report.pdf)
				singam = Math.Atan(-Math.Cos(alpha) / (-Math.Sin(phi_0) * Math.Sin(alpha))) + lonz;
				// Equation 9-7 page 80 (http://pubs.usgs.gov/pp/1395/report.pdf)
				sinphi = Math.Asin(Math.Cos(phi_0) * Math.Sin(alpha));
			}
			else // If the keyword "alpha" is NOT found in the sentence then use 2points
			{
				// Define Pole of oblique transformation from 2 points
				double phi_1 = Proj.pj_param_r(ctx, parameters, "lat_1");
				double phi_2 = Proj.pj_param_r(ctx, parameters, "lat_2");
				double lam_1 = Proj.pj_param_r(ctx, parameters, "lon_1");
				double lam_2 = Proj.pj_param_r(ctx, parameters, "lon_2");
				// Equation 9-1 page 80 (http://pubs.usgs.gov/pp/1395/report.pdf)
				singam = Math.Atan2(Math.Cos(phi_1) * Math.Sin(phi_2) * Math.Cos(lam_1) - Math.Sin(phi_1) * Math.Cos(phi_2) * Math.Cos(lam_2),
					Math.Sin(phi_1) * Math.Cos(phi_2) * Math.Sin(lam_2) - Math.Cos(phi_1) * Math.Sin(phi_2) * Math.Sin(lam_1));
				// Equation 9-2 page 80 (http://pubs.usgs.gov/pp/1395/report.pdf)
				sinphi = Math.Atan(-Math.Cos(singam - lam_1) / Math.Tan(phi_1));
			}

			lam0 = singam + Proj.HALFPI;
			cosphi = Math.Cos(sinphi);
			sinphi = Math.Sin(sinphi);
			cosgam = Math.Cos(singam);
			singam = Math.Sin(singam);
			inv = s_inverse;
			fwd = s_forward;
			es = 0.0;

			return this;
		}
	}
}
