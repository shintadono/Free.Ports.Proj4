using System;
using System.Text;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_hammer : PJ
	{
		protected double w, m, rm;

		const double EPS10 = 1.0e-10;

		public override string Name { get { return "hammer"; } }
		public override string DescriptionName { get { return "Hammer & Eckert-Greifendorff"; } }
		public override string DescriptionType { get { return "Misc Sph"; } }
		public override string DescriptionParameters { get { return "W= M="; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret = new StringBuilder();
				if (w != 0.5) ret.AppendFormat(nc, " +W={0}", w);
				if (rm != 1.0) ret.AppendFormat(nc, " +M={0}", 1 / rm);
				return ret.ToString();
			}
		}

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x = xy.y = 0;

			double cosphi = Math.Cos(lp.phi);
			lp.lam *= w;
			double d = Math.Sqrt(2.0 / (1.0 + cosphi * Math.Cos(lp.lam)));
			xy.x = m * d * cosphi * Math.Sin(lp.lam);
			xy.y = rm * d * Math.Sin(lp.phi);

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam = lp.phi = 0;

			double z = Math.Sqrt(1.0 - 0.25 * w * w * xy.x * xy.x - 0.25 * xy.y * xy.y);
			if (Math.Abs(2.0 * z * z - 1.0) < EPS10)
			{
				lp.lam = Libc.HUGE_VAL;
				lp.phi = Libc.HUGE_VAL;
				Proj.pj_ctx_set_errno(ctx, -14);
			}
			else
			{
				lp.lam = Proj.aatan2(w * xy.x * z, 2.0 * z * z - 1) / w;
				lp.phi = Proj.aasin(ctx, z * xy.y);
			}

			return lp;
		}

		public override PJ Init()
		{
			if (Proj.pj_param_t(ctx, parameters, "W"))
			{
				w = Math.Abs(Proj.pj_param_d(ctx, parameters, "W"));
				if (w <= 0.0) { Proj.pj_ctx_set_errno(ctx, -27); return null; }
			}
			else w = 0.5;

			if (Proj.pj_param_t(ctx, parameters, "M"))
			{
				m = Math.Abs(Proj.pj_param_d(ctx, parameters, "M"));
				if (m <= 0.0) { Proj.pj_ctx_set_errno(ctx, -27); return null; }
			}
			else m = 1.0;

			rm = 1.0 / m;
			m /= w;
			es = 0.0;
			fwd = s_forward;
			inv = s_inverse;

			return this;
		}
	}
}
