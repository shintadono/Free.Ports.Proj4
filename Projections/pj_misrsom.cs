//*****************************************************************************
// This implements Space Oblique Mercator (SOM) projection, used by the
// Multi-angle Imaging SpectroRadiometer (MISR) products, from the NASA EOS Terra
// platform.
//
// The code is identical to that of Landsat SOM (PJ_lsat.c) with the following
// parameter changes:
//
//   inclination angle = 98.30382 degrees
//   period of revolution = 98.88 minutes
//   ascending longitude = 129.3056 degrees - (360 / 233) * path_number
//
// and the following code change:
//
//   P->rlm = PI * (1. / 248. + .5161290322580645);
//
// changed to:
//
//   P->rlm = 0
//
//*****************************************************************************
// based upon Snyder and Linck, USGS-NMD

using System;
using System.Text;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_misrsom : PJ
	{
		protected double a2, a4, b, c1, c3;
		protected double q, t, u, w, p22, sa, ca, xj, rlm, rlm2;

		const double TOL = 1e-7;
		const double PI_HALFPI = 4.71238898038468985766;
		const double TWOPI_HALFPI = 7.85398163397448309610;

		public override string Name { get { return "misrsom"; } }
		public override string DescriptionName { get { return "Space oblique for MISR"; } }
		public override string DescriptionType { get { return "Cyl, Sph&Ell"; } }
		public override string DescriptionParameters { get { return "path="; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				int path = (int)Math.Round((Proj.DEG_TO_RAD * 129.3056 - lam0) * 233.0 / Proj.TWOPI);

				StringBuilder ret = new StringBuilder();
				ret.AppendFormat(nc, " +path={0}", path);
				return ret.ToString();
			}
		}

		static void seraz0(double lam, double mult, PJ_misrsom P)
		{
			lam *= Proj.DEG_TO_RAD;
			double sd = Math.Sin(lam);
			double sdsq = sd * sd;
			double s = P.p22 * P.sa * Math.Cos(lam) * Math.Sqrt((1.0 + P.t * sdsq) / ((1.0 + P.w * sdsq) * (1.0 + P.q * sdsq)));
			double d__1 = 1.0 + P.q * sdsq;
			double h = Math.Sqrt((1.0 + P.q * sdsq) / (1.0 + P.w * sdsq)) * ((1.0 + P.w * sdsq) / (d__1 * d__1) - P.p22 * P.ca);
			double sq = Math.Sqrt(P.xj * P.xj + s * s);
			double fc = mult * (h * P.xj - s * s) / sq;
			P.b += fc;
			P.a2 += fc * Math.Cos(lam + lam);
			P.a4 += fc * Math.Cos(lam * 4.0);
			fc = mult * s * (h + P.xj) / sq;
			P.c1 += fc * Math.Cos(lam);
			P.c3 += fc * Math.Cos(lam * 3.0);
		}

		// ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x = xy.y = 0;

			if (lp.phi > Proj.HALFPI) lp.phi = Proj.HALFPI;
			else if (lp.phi < -Proj.HALFPI) lp.phi = -Proj.HALFPI;

			double lampp = lp.phi >= 0.0 ? Proj.HALFPI : PI_HALFPI;
			double tanphi = Math.Tan(lp.phi);

			int l;
			double lamt = 0, lamdp = 0;
			for (int nn = 0; ;)
			{
				double sav = lampp;

				double lamtp = lp.lam + p22 * lampp;
				double cl = Math.Cos(lamtp);
				if (Math.Abs(cl) < TOL) lamtp -= TOL;
				double fac = lampp - Math.Sin(lampp) * (cl < 0.0 ? -Proj.HALFPI : Proj.HALFPI);

				for (l = 50; l != 0; --l)
				{
					lamt = lp.lam + p22 * sav;
					double c = Math.Cos(lamt);
					if (Math.Abs(c) < TOL) lamt -= TOL;
					double xlam = (one_es * tanphi * sa + Math.Sin(lamt) * ca) / c;
					lamdp = Math.Atan(xlam) + fac;
					if (Math.Abs(Math.Abs(sav) - Math.Abs(lamdp)) < TOL) break;

					sav = lamdp;
				}

				if (l == 0 || ++nn >= 3 || (lamdp > rlm && lamdp < rlm2)) break;

				if (lamdp <= rlm) lampp = TWOPI_HALFPI;
				else if (lamdp >= rlm2) lampp = Proj.HALFPI;
			}

			if (l != 0)
			{
				double sp = Math.Sin(lp.phi);
				double phidp = Proj.aasin(ctx, (one_es * ca * sp - sa * Math.Cos(lp.phi) * Math.Sin(lamt)) / Math.Sqrt(1.0 - es * sp * sp));
				double tanph = Math.Log(Math.Tan(Proj.FORTPI + 0.5 * phidp));
				double sd = Math.Sin(lamdp);
				double sdsq = sd * sd;
				double s = p22 * sa * Math.Cos(lamdp) * Math.Sqrt((1.0 + t * sdsq) / ((1.0 + w * sdsq) * (1.0 + q * sdsq)));
				double d = Math.Sqrt(xj * xj + s * s);
				xy.x = b * lamdp + a2 * Math.Sin(2.0 * lamdp) + a4 * Math.Sin(lamdp * 4.0) - tanph * s / d;
				xy.y = c1 * sd + c3 * Math.Sin(lamdp * 3.0) + tanph * xj / d;
			}
			else xy.x = xy.y = Libc.HUGE_VAL;

			return xy;
		}

		// ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam = lp.phi = 0;

			double lamdp = xy.x / b;

			double s, sav;
			int nn = 50;
			do
			{
				sav = lamdp;
				double sd = Math.Sin(lamdp);
				double sdsq = sd * sd;
				s = p22 * sa * Math.Cos(lamdp) * Math.Sqrt((1.0 + t * sdsq) / ((1.0 + w * sdsq) * (1.0 + q * sdsq)));
				lamdp = xy.x + xy.y * s / xj - a2 * Math.Sin(2.0 * lamdp) - a4 * Math.Sin(lamdp * 4.0) - s / xj * (c1 * Math.Sin(lamdp) + c3 * Math.Sin(lamdp * 3.0));
				lamdp /= b;
			} while (Math.Abs(lamdp - sav) >= TOL && --nn != 0);

			double sl = Math.Sin(lamdp);
			double fac = Math.Exp(Math.Sqrt(1.0 + s * s / xj / xj) * (xy.y - c1 * sl - c3 * Math.Sin(lamdp * 3.0)));
			double phidp = 2.0 * (Math.Atan(fac) - Proj.FORTPI);
			double dd = sl * sl;
			if (Math.Abs(Math.Cos(lamdp)) < TOL) lamdp -= TOL;

			double spp = Math.Sin(phidp);
			double sppsq = spp * spp;
			double lamt = Math.Atan(((1.0 - sppsq * rone_es) * Math.Tan(lamdp) * ca - spp * sa * Math.Sqrt((1.0 + q * dd) * (1.0 - sppsq) - sppsq * u) / Math.Cos(lamdp)) / (1.0 - sppsq * (1.0 + u)));
			sl = lamt >= 0.0 ? 1.0 : -1.0;
			double scl = Math.Cos(lamdp) >= 0.0 ? 1.0 : -1;
			lamt -= Proj.HALFPI * (1.0 - scl) * sl;
			lp.lam = lamt - p22 * lamdp;

			if (Math.Abs(sa) < TOL) lp.phi = Proj.aasin(ctx, spp / Math.Sqrt(one_es * one_es + es * sppsq));
			else lp.phi = Math.Atan((Math.Tan(lamdp) * Math.Cos(lamt) - ca * Math.Sin(lamt)) / (one_es * sa));

			return lp;
		}

		public override PJ Init()
		{
			int path = Proj.pj_param_i(ctx, parameters, "path");
			if (path <= 0 || path > 233) { Proj.pj_ctx_set_errno(ctx, -29); return null; }

			lam0 = Proj.DEG_TO_RAD * 129.3056 - Proj.TWOPI / 233.0 * path;
			double alf = 98.30382 * Proj.DEG_TO_RAD;
			p22 = 98.88 / 1440.0;

			sa = Math.Sin(alf);
			ca = Math.Cos(alf);
			if (Math.Abs(ca) < 1e-9) ca = 1e-9;

			double esc = es * ca * ca;
			double ess = es * sa * sa;
			w = (1.0 - esc) * rone_es;
			w = w * w - 1.0;
			q = ess * rone_es;
			t = ess * (2.0 - es) * rone_es * rone_es;
			u = esc * rone_es;
			xj = one_es * one_es * one_es;
			rlm = 0;
			rlm2 = rlm + Proj.TWOPI;
			a2 = a4 = b = c1 = c3 = 0.0;

			seraz0(0.0, 1.0, this);
			for (double lam = 9.0; lam <= 81.0001; lam += 18.0) seraz0(lam, 4.0, this);
			for (double lam = 18; lam <= 72.0001; lam += 18.0) seraz0(lam, 2.0, this);

			seraz0(90.0, 1.0, this);
			a2 /= 30.0;
			a4 /= 60.0;
			b /= 30.0;
			c1 /= 15.0;
			c3 /= 45.0;
			inv = e_inverse;
			fwd = e_forward;

			return this;
		}
	}
}
