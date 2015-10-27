using System;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		// general forward projection

		// forward projection entry
		public static XY pj_fwd(LP lp, PJ P)
		{
			XY xy;

			// check for forward and latitude or longitude overange
			double t=Math.Abs(lp.phi)-HALFPI;
			if(t>EPS12||Math.Abs(lp.lam)>10.0)
			{
				xy.x=xy.y=Libc.HUGE_VAL;
				pj_ctx_set_errno(P.ctx, -14);
				return xy;
			}

			// proceed with projection
			Libc.errno=pj_errno=0;
			P.ctx.last_errno=0;
			if(Math.Abs(t)<=EPS12) lp.phi=lp.phi<0.0?-HALFPI:HALFPI;
			else if(P.geoc) lp.phi=Math.Atan(P.rone_es*Math.Tan(lp.phi));

			lp.lam-=P.lam0; // compute del lp.lam
			if(!P.over) lp.lam=adjlon(lp.lam); // adjust del longitude

			// Check for NULL pointer
			if (P.fwd != null)
			{
				xy = P.fwd(lp); // project
				if (P.ctx.last_errno != 0)
				{
					xy.x = xy.y = Libc.HUGE_VAL;
				}
				else // adjust for major axis and easting/northings
				{
					xy.x = P.fr_meter * (P.a * xy.x + P.x0);
					xy.y = P.fr_meter * (P.a * xy.y + P.y0);
				}
			}
			else
			{
				xy.x = xy.y = Libc.HUGE_VAL;
			}

			if(pj_errno==0) pj_errno=Libc.errno;

			return xy;
		}
	}
}
