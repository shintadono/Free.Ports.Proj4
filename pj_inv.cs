using System;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		// general inverse projection

		// inverse projection entry
		public static LP pj_inv(XY xy, PJ P)
		{
			LP lp;

			// can't do as much preliminary checking as with forward
			if (xy.x == Libc.HUGE_VAL || xy.y == Libc.HUGE_VAL)
			{
				lp.lam = lp.phi = Libc.HUGE_VAL;
				pj_ctx_set_errno(P.ctx, -15);

				return lp;
			}

			Libc.errno = pj_errno = 0;
			P.ctx.last_errno = 0;
			xy.x = (xy.x * P.to_meter - P.x0) * P.ra; // descale and de-offset
			xy.y = (xy.y * P.to_meter - P.y0) * P.ra;

			// Check for NULL pointer
			if (P.inv != null)
			{
				lp = P.inv(xy); // inverse project
				if (P.ctx.last_errno != 0)
				{
					lp.lam = lp.phi = Libc.HUGE_VAL;
				}
				else
				{
					lp.lam += P.lam0; // reduce from del lp.lam
					if (!P.over) lp.lam = adjlon(lp.lam); // adjust longitude to CM
					if (P.geoc && Math.Abs(Math.Abs(lp.phi) - HALFPI) > EPS12)
						lp.phi = Math.Atan(P.one_es * Math.Tan(lp.phi));
				}
			}
			else
			{
				lp.lam = lp.phi = Libc.HUGE_VAL;
			}

			return lp;
		}
	}
}
