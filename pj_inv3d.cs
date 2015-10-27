using System;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		// general inverse projection

		// inverse projection entry
		public static LPZ pj_inv3d(XYZ xyz, PJ P)
		{
			LPZ lpz;

			// can't do as much preliminary checking as with forward
			if (xyz.x == Libc.HUGE_VAL || xyz.y == Libc.HUGE_VAL || xyz.z == Libc.HUGE_VAL)
			{
				lpz.lam = lpz.phi = lpz.z = Libc.HUGE_VAL;
				pj_ctx_set_errno(P.ctx, -15);
				return lpz;
			}

			Libc.errno = pj_errno = 0;
			P.ctx.last_errno = 0;

			xyz.x = (xyz.x * P.to_meter - P.x0) * P.ra; // descale and de-offset
			xyz.y = (xyz.y * P.to_meter - P.y0) * P.ra;
			// z is not scaled since that is handled by vto_meter before we get here

			// Check for NULL pointer
			if (P.inv3d != null)
			{
				lpz = P.inv3d(xyz); // inverse project
				if (P.ctx.last_errno != 0) lpz.lam = lpz.phi = lpz.z = Libc.HUGE_VAL;
				else
				{
					lpz.lam += P.lam0; // reduce from del lp.lam
					if (!P.over) lpz.lam = adjlon(lpz.lam); // adjust longitude to CM

					// This maybe redundant and never used
					if (P.geoc && Math.Abs(Math.Abs(lpz.phi) - HALFPI) > EPS12) lpz.phi = Math.Atan(P.one_es * Math.Tan(lpz.phi));
				}
			}
			else
			{
				lpz.lam = lpz.phi = lpz.z = Libc.HUGE_VAL;
			}

			return lpz;
		}
	}
}
