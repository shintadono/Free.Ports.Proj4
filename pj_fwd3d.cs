using System;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		// general forward projection

		// forward projection entry
		public static XYZ pj_fwd3d(LPZ lpz, PJ P)
		{
			XYZ xyz;

			// check for forward and latitude or longitude overange
			double t = Math.Abs(lpz.phi) - HALFPI;
			if (t > EPS12 || Math.Abs(lpz.lam) > 10.0)
			{
				xyz.x = xyz.y = xyz.z = Libc.HUGE_VAL;
				pj_ctx_set_errno(P.ctx, -14);
				return xyz;
			}

			// proceed with projection
			P.ctx.last_errno = 0;
			pj_errno = 0;
			Libc.errno = 0;

			if (Math.Abs(t) <= EPS12) lpz.phi = lpz.phi < 0.0 ? -HALFPI : HALFPI;
			else if (P.geoc) lpz.phi = Math.Atan(P.rone_es * Math.Tan(lpz.phi)); // Maybe redundant and never used.
			lpz.lam -= P.lam0; // compute del lp.lam
			if (!P.over) lpz.lam = adjlon(lpz.lam); // adjust del longitude

			// Check for NULL pointer
			if (P.fwd3d != null)
			{
				xyz = P.fwd3d(lpz); // project
				if (P.ctx.last_errno != 0) xyz.x = xyz.y = xyz.z = Libc.HUGE_VAL;
				// adjust for major axis and easting/northings
				else
				{
					xyz.x = P.fr_meter * (P.a * xyz.x + P.x0);
					xyz.y = P.fr_meter * (P.a * xyz.y + P.y0);
					// z is not scaled since this handled by vto_meter outside
				}
			}
			else
			{
				xyz.x = xyz.y = xyz.z = Libc.HUGE_VAL;
			}

			return xyz;
		}
	}
}
