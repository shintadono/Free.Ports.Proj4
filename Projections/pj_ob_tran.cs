using System;
using System.Text;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_ob_tran : PJ
	{
		protected PJ link;
		protected double lamp, cphip, sphip;

		public override string Name { get { return "ob_tran"; } }
		public override string DescriptionName { get { return "General Oblique Transformation"; } }
		public override string DescriptionType { get { return "Misc Sph"; } }
		public override string DescriptionParameters { get { return "o_proj= plus parameters for projection o_lat_p= o_lon_p= (new pole) or o_alpha= o_lon_c= o_lat_c= or o_lon_1= o_lat_1= o_lon_2= o_lat_2="; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +o_proj={0}", link.Name);
				ret.AppendFormat(nc, " +o_lon_p={0}", lamp*Proj.RAD_TO_DEG);
				ret.AppendFormat(nc, " +o_lat_p={0}", Math.Atan2(sphip, cphip)*Proj.RAD_TO_DEG);
				return ret.ToString();
			}
		}

		const double TOL=1.0e-10;

		// spheroid
		XY o_forward(LP lp)
		{
			double coslam=Math.Cos(lp.lam);
			double sinphi=Math.Sin(lp.phi);
			double cosphi=Math.Cos(lp.phi);
			lp.lam=Proj.adjlon(Proj.aatan2(cosphi*Math.Sin(lp.lam), sphip*cosphi*coslam+cphip*sinphi)+lamp);
			lp.phi=Proj.aasin(ctx, sphip*sinphi-cphip*cosphi*coslam);

			return link.fwd(lp);
		}

		// spheroid
		XY t_forward(LP lp)
		{
			double cosphi=Math.Cos(lp.phi);
			double coslam=Math.Cos(lp.lam);
			lp.lam=Proj.adjlon(Proj.aatan2(cosphi*Math.Sin(lp.lam), Math.Sin(lp.phi))+lamp);
			lp.phi=Proj.aasin(ctx, -cosphi*coslam);

			return link.fwd(lp);
		}

		// spheroid
		LP o_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp=link.inv(xy);
			if(lp.lam!=Libc.HUGE_VAL)
			{
				lp.lam-=lamp;
				double coslam=Math.Cos(lp.lam);
				double sinphi=Math.Sin(lp.phi);
				double cosphi=Math.Cos(lp.phi);
				lp.phi=Proj.aasin(ctx, sphip*sinphi+cphip*cosphi*coslam);
				lp.lam=Proj.aatan2(cosphi*Math.Sin(lp.lam), sphip*cosphi*coslam-cphip*sinphi);
			}

			return lp;
		}

		// spheroid
		LP t_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp=link.inv(xy);
			if(lp.lam!=Libc.HUGE_VAL)
			{
				double cosphi=Math.Cos(lp.phi);
				double t=lp.lam-lamp;
				lp.lam=Proj.aatan2(cosphi*Math.Sin(t), -Math.Sin(lp.phi));
				lp.phi=Proj.aasin(ctx, cosphi*Math.Cos(t));
			}

			return lp;
		}

		public override PJ Init()
		{
			// get name of projection to be translated
			string name=Proj.pj_param_s(ctx, parameters, "o_proj");
			if(name==null||name=="") { Proj.pj_ctx_set_errno(ctx, -26); return null; }

			link=Proj.GetPJ(name);
			if(link==null) { Proj.pj_ctx_set_errno(ctx, -5); return null; }

			// copy existing header into new
			es=0.0; // force to spherical
			link.parameters=parameters;
			link.ctx=ctx;
			link.over=over;
			link.geoc=geoc;
			link.a=a;
			link.es=es;
			link.ra=ra;
			link.lam0=lam0;
			link.phi0=phi0;
			link.x0=x0;
			link.y0=y0;
			link.k0=k0;

			// force spherical earth
			link.one_es=link.rone_es=1.0;
			link.es=link.e=0.0;
			link=link.Init();
			if(link==null) return null;

			double phip=0;

			if(Proj.pj_param_t(ctx, parameters, "o_alpha"))
			{
				double lamc=Proj.pj_param_r(ctx, parameters, "o_lon_c");
				double phic=Proj.pj_param_r(ctx, parameters, "o_lat_c");
				double alpha=Proj.pj_param_r(ctx, parameters, "o_alpha");
				
				//if(Math.Abs(phic)<=TOL10||Math.Abs(Math.Abs(phic)-HALFPI)<=TOL10||Math.Abs(Math.Abs(alpha)-HALFPI)<=TOL10)

				if(Math.Abs(Math.Abs(phic)-Proj.HALFPI)<=TOL) { Proj.pj_ctx_set_errno(ctx, -32); return null; }

				lamp=lamc+Proj.aatan2(-Math.Cos(alpha), -Math.Sin(alpha)*Math.Sin(phic));
				phip=Proj.aasin(ctx, Math.Cos(phic)*Math.Sin(alpha));
			}
			else if(Proj.pj_param_t(ctx, parameters, "o_lat_p"))
			{ // specified new pole
				lamp=Proj.pj_param_r(ctx, parameters, "o_lon_p");
				phip=Proj.pj_param_r(ctx, parameters, "o_lat_p");
			}
			else
			{ // specified new "equator" points
				double lam1=Proj.pj_param_r(ctx, parameters, "o_lon_1");
				double phi1=Proj.pj_param_r(ctx, parameters, "o_lat_1");
				double lam2=Proj.pj_param_r(ctx, parameters, "o_lon_2");
				double phi2=Proj.pj_param_r(ctx, parameters, "o_lat_2");

				double con=Math.Abs(phi1);
				if(Math.Abs(phi1-phi2)<=TOL||con<=TOL||
					Math.Abs(con-Proj.HALFPI)<=TOL||Math.Abs(Math.Abs(phi2)-Proj.HALFPI)<=TOL) { Proj.pj_ctx_set_errno(ctx, -33); return null; }
				lamp=Math.Atan2(Math.Cos(phi1)*Math.Sin(phi2)*Math.Cos(lam1)-Math.Sin(phi1)*Math.Cos(phi2)*Math.Cos(lam2),
					Math.Sin(phi1)*Math.Cos(phi2)*Math.Sin(lam2)-Math.Cos(phi1)*Math.Sin(phi2)*Math.Sin(lam1));
				phip=Math.Atan(-Math.Cos(lamp-lam1)/Math.Tan(phi1));
			}

			if(Math.Abs(phip)>TOL)
			{ // oblique
				cphip=Math.Cos(phip);
				sphip=Math.Sin(phip);
				fwd=o_forward;
				inv=link.inv!=null?o_inverse:(Func<XY, LP>)null;
			}
			else
			{ // transverse
				fwd=t_forward;
				inv=link.inv!=null?t_inverse:(Func<XY, LP>)null;
			}

			return this;
		}
	}
}
