using System;
using System.Text;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_tmerc : PJ
	{
		protected double esp, ml0;
		protected double[] en;

		public override string Name { get { return "tmerc"; } }
		public override string DescriptionName { get { return "Transverse Mercator"; } }
		public override string DescriptionType { get { return "Cyl, Sph&Ell"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double FC1=1.0;
		const double FC2=0.5;
		const double FC3=0.16666666666666666666;
		const double FC4=0.08333333333333333333;
		const double FC5=0.05;
		const double FC6=0.03333333333333333333;
		const double FC7=0.02380952380952380952;
		const double FC8=0.01785714285714285714;

		// ellipse
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			// Fail if our longitude is more than 90 degrees from the
			// central meridian since the results are essentially garbage.
			// Is error -20 really an appropriate return value?
			//
			// http://trac.osgeo.org/proj/ticket/5
			if(lp.lam<-Proj.HALFPI||lp.lam>Proj.HALFPI)
			{
				xy.x=Libc.HUGE_VAL;
				xy.y=Libc.HUGE_VAL;
				Proj.pj_ctx_set_errno(ctx, -14);

				return xy;
			}

			double sinphi=Math.Sin(lp.phi);
			double cosphi=Math.Cos(lp.phi);
			double t=Math.Abs(cosphi)>1e-10?sinphi/cosphi:0.0;
			t*=t;
			double al=cosphi*lp.lam;
			double als=al*al;
			al/=Math.Sqrt(1.0-es*sinphi*sinphi);
			double n=esp*cosphi*cosphi;

			xy.x=k0*al*(FC1+FC3*als*(1.0-t+n+FC5*als*(5.0+t*(t-18.0)+n*(14.0-58.0*t)+FC7*als*(61.0+t*(t*(179.0-t)-479.0)))));

			xy.y=k0*(Proj.pj_mlfn(lp.phi, sinphi, cosphi, en)-ml0+sinphi*al*lp.lam*FC2*(1.0+FC4*als*(5.0-t+n*(9.0+4.0*n)+
				FC6*als*(61.0+t*(t-58.0)+n*(270.0-330*t)+FC8*als*(1385.0+t*(t*(543.0-t)-3111.0))))));

			return xy;
		}

		// sphere
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			// Fail if our longitude is more than 90 degrees from the
			// central meridian since the results are essentially garbage.
			// Is error -20 really an appropriate return value?
			//
			// http://trac.osgeo.org/proj/ticket/5
			if(lp.lam<-Proj.HALFPI||lp.lam>Proj.HALFPI)
			{
				xy.x=Libc.HUGE_VAL;
				xy.y=Libc.HUGE_VAL;
				Proj.pj_ctx_set_errno(ctx, -14);

				return xy;
			}

			double cosphi=Math.Cos(lp.phi);

			double b=cosphi*Math.Sin(lp.lam);
			if(Math.Abs(Math.Abs(b)-1.0)<=1.0e-10) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }

			xy.x=ml0*Math.Log((1.0+b)/(1.0-b));
			if((b=Math.Abs(xy.y=cosphi*Math.Cos(lp.lam)/Math.Sqrt(1.0-b*b)))>=1.0)
			{
				if((b-1.0)>1.0e-10) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }
				else xy.y=0.0;
			}
			else xy.y=Math.Acos(xy.y);

			if(lp.phi<0.0) xy.y=-xy.y;
			xy.y=esp*(xy.y-phi0);

			return xy;
		}

		// ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=Proj.pj_inv_mlfn(ctx, ml0+xy.y/k0, es, en);
			if(Math.Abs(lp.phi)>=Proj.HALFPI)
			{
				lp.phi=xy.y<0.0?-Proj.HALFPI:Proj.HALFPI;
				lp.lam=0.0;
			}
			else
			{
				double sinphi=Math.Sin(lp.phi);
				double cosphi=Math.Cos(lp.phi);
				double t=Math.Abs(cosphi)>1e-10?sinphi/cosphi:0.0;
				double n=esp*cosphi*cosphi;
				double con=1.0-es*sinphi*sinphi;
				double d=xy.x*Math.Sqrt(con)/k0;
				con*=t;
				t*=t;
				double ds=d*d;

				lp.phi-=(con*ds/(1.0-es))*FC2*(1.0-ds*FC4*(5.0+t*(3.0-9.0*n)+n*(1.0-4*n)-ds*FC6*(61.0+t*(90.0-252.0*n+45.0*t)+46.0*n
					-ds*FC8*(1385.0+t*(3633.0+t*(4095.0+1574.0*t))))));

				lp.lam=d*(FC1-ds*FC3*(1.0+2.0*t+n-ds*FC5*(5.0+t*(28.0+24.0*t+8.0*n)+6.0*n-ds*FC7*(61.0+t*(662.0+t*(1320.0+720.0*t))))))/cosphi;
			}

			return lp;
		}

		// sphere
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double h=Math.Exp(xy.x/esp);
			double g=0.5*(h-1.0/h);
			h=Math.Cos(phi0+xy.y/esp);
			lp.phi=Math.Asin(Math.Sqrt((1.0-h*h)/(1.0+g*g)));
			if(xy.y<0.0) lp.phi=-lp.phi;
			lp.lam=(g!=0||h!=0)?Math.Atan2(g, h):0.0;

			return lp;
		}

		// general initialization
		protected PJ_tmerc setup()
		{
			if(es!=0)
			{
				en=Proj.pj_enfn(es);
				if(en==null) return null;

				ml0=Proj.pj_mlfn(phi0, Math.Sin(phi0), Math.Cos(phi0), en);
				esp=es/(1.0-es);
				inv=e_inverse;
				fwd=e_forward;
			}
			else
			{
				esp=k0;
				ml0=0.5*esp;
				inv=s_inverse;
				fwd=s_forward;
			}

			return this;
		}

		public override PJ Init()
		{
			return setup();
		}
	}
}
