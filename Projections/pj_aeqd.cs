//*****************************************************************************
//
// Project:	PROJ.4
// Purpose:	Implementation of the aeqd (Azimuthal Equidistant) projection.
// Author:	Gerald Evenden
//
//*****************************************************************************
// Copyright (c) 1995, Gerald Evenden
// Copyright (c) 2008-2015 by the Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//*****************************************************************************

using System;
using Free.Ports.Proj4.LibcStuff;
using Free.Ports.Proj4.Geodesic;

namespace Free.Ports.Proj4.Projections
{
	class PJ_aeqd : PJ
	{
		protected double sinph0;
		protected double cosph0;
		protected double[] en;
		protected double M1;
		protected double N1;
		protected double Mp;
		protected double He;
		protected double G;
		protected aeqd_mode mode;
		protected geod_geodesic g;

		const double EPS10=1.0e-10;
		const double TOL=1.0e-14;

		public override string Name { get { return "aeqd"; } }
		public override string DescriptionName { get { return "Azimuthal Equidistant"; } }
		public override string DescriptionType { get { return "Azi, Sph&Ell"; } }
		public override string DescriptionParameters { get { return "lat_0= guam"; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				if(inv==e_guam_inv) return "guam";
				return "";
			}
		}

		protected enum aeqd_mode
		{
			N_POLE=0,
			S_POLE=1,
			EQUIT=2,
			OBLIQ=3
		}

		// Guam elliptical
		XY e_guam_fwd(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double cosphi=Math.Cos(lp.phi);
			double sinphi=Math.Sin(lp.phi);
			double t=1.0/Math.Sqrt(1.0-es*sinphi*sinphi);
			xy.x=lp.lam*cosphi*t;
			xy.y=Proj.pj_mlfn(lp.phi, sinphi, cosphi, en)-M1+0.5*lp.lam*lp.lam*cosphi*sinphi*t;

			return xy;
		}

		// elliptical
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double coslam=Math.Cos(lp.lam);
			double cosphi=Math.Cos(lp.phi);
			double sinphi=Math.Sin(lp.phi);

			switch(mode)
			{
				case aeqd_mode.N_POLE: coslam=-coslam; goto case aeqd_mode.S_POLE;
				case aeqd_mode.S_POLE:
					double rho=Math.Abs(Mp-Proj.pj_mlfn(lp.phi, sinphi, cosphi, en));
					xy.x=rho*Math.Sin(lp.lam);
					xy.y=rho*coslam;
					break;
				case aeqd_mode.EQUIT:
				case aeqd_mode.OBLIQ:
					if(Math.Abs(lp.lam)<EPS10&&Math.Abs(lp.phi-phi0)<EPS10)
					{
						xy.x=xy.y=0.0;
						break;
					}

					double phi1=phi0*Proj.DEG_TO_RAD;
					double lam1=lam0*Proj.DEG_TO_RAD;
					double phi2=lp.phi*Proj.DEG_TO_RAD;
					double lam2=(lp.lam+lam0)*Proj.DEG_TO_RAD;

					double azi1, azi2, s12;
					g.geod_inverse(phi1, lam1, phi2, lam2, out s12, out azi1, out azi2);
					
					azi1/=Proj.DEG_TO_RAD;
					xy.x=s12*Math.Sin(azi1)/a;
					xy.y=s12*Math.Cos(azi1)/a;
					break;
			}

			return xy;
		}

		// spherical
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double sinphi=Math.Sin(lp.phi);
			double cosphi=Math.Cos(lp.phi);
			double coslam=Math.Cos(lp.lam);

			switch(mode)
			{
				case aeqd_mode.EQUIT:
				case aeqd_mode.OBLIQ:
					if(mode==aeqd_mode.EQUIT) xy.y=cosphi*coslam;
					else xy.y=sinph0*sinphi+cosph0*cosphi*coslam;
					if(Math.Abs(Math.Abs(xy.y)-1.0)<TOL)
					{
						if(xy.y<0.0)
						{
							Proj.pj_ctx_set_errno(ctx, -20);
							return xy;
						}
						else xy.x=xy.y=0.0;
					}
					else
					{
						xy.y=Math.Acos(xy.y);
						xy.y/=Math.Sin(xy.y);
						xy.x=xy.y*cosphi*Math.Sin(lp.lam);
						xy.y*=(mode==aeqd_mode.EQUIT)?sinphi:cosph0*sinphi-sinph0*cosphi*coslam;
					}
					break;
				case aeqd_mode.N_POLE:
					lp.phi=-lp.phi;
					coslam=-coslam;
					goto case aeqd_mode.S_POLE;
				case aeqd_mode.S_POLE:
					if(Math.Abs(lp.phi-Proj.HALFPI)<EPS10)
					{
						Proj.pj_ctx_set_errno(ctx, -20);
						return xy;
					}
					xy.y=Proj.HALFPI+lp.phi;
					xy.x=xy.y*Math.Sin(lp.lam);
					xy.y*=coslam;
					break;
			}
			return xy;
		}

		// Guam elliptical
		LP e_guam_inv(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double x2, t=0;

			x2=0.5*xy.x*xy.x;
			lp.phi=phi0;
			for(int i=0; i<3; i++)
			{
				t=e*Math.Sin(lp.phi);
				t=Math.Sqrt(1.0-t*t);
				lp.phi=Proj.pj_inv_mlfn(ctx, M1+xy.y-x2*Math.Tan(lp.phi)*t, es, en);
			}

			lp.lam=xy.x*t/Math.Cos(lp.phi);
			return lp;
		}

		// elliptical
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double c=Libc.hypot(xy.x, xy.y);
			if(c<EPS10)
			{
				lp.phi=phi0;
				lp.lam=0.0;
				return lp;
			}

			if(mode==aeqd_mode.OBLIQ||mode==aeqd_mode.EQUIT)
			{
				double x2=xy.x*a;
				double y2=xy.y*a;
				double lat1=phi0*Proj.DEG_TO_RAD;
				double lon1=lam0*Proj.DEG_TO_RAD;
				double azi1=Math.Atan2(x2, y2)*Proj.DEG_TO_RAD;
				double s12=Math.Sqrt(x2*x2+y2*y2);

				double lat2, lon2, azi2;
				g.geod_direct(lat1, lon1, azi1, s12, out lat2, out lon2, out azi2);
				
				lp.phi=lat2/Proj.DEG_TO_RAD;
				lp.lam=lon2/Proj.DEG_TO_RAD;
				lp.lam-=lam0;
			}
			else
			{ // Polar
				lp.phi=Proj.pj_inv_mlfn(ctx, mode==aeqd_mode.N_POLE?Mp-c:Mp+c, es, en);
				lp.lam=Math.Atan2(xy.x, mode==aeqd_mode.N_POLE?-xy.y:xy.y);
			}
			return lp;
		}

		// spherical
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double c_rh=Libc.hypot(xy.x, xy.y);
			if(c_rh>Proj.PI)
			{
				if(c_rh-EPS10>Proj.PI)
				{
					Proj.pj_ctx_set_errno(ctx, -20);
					return lp;
				}
				c_rh=Proj.PI;
			}
			else if(c_rh<EPS10)
			{
				lp.phi=phi0;
				lp.lam=0.0;
				return lp;
			}

			if(mode==aeqd_mode.OBLIQ||mode==aeqd_mode.EQUIT)
			{
				double sinc=Math.Sin(c_rh);
				double cosc=Math.Cos(c_rh);
				if(mode==aeqd_mode.EQUIT)
				{
					lp.phi=Proj.aasin(ctx, xy.y*sinc/c_rh);
					xy.x*=sinc;
					xy.y=cosc*c_rh;
				}
				else
				{
					lp.phi=Proj.aasin(ctx, cosc*sinph0+xy.y*sinc*cosph0/c_rh);
					xy.y=(cosc-sinph0*Math.Sin(lp.phi))*c_rh;
					xy.x*=sinc*cosph0;
				}
				lp.lam=Math.Atan2(xy.x, xy.y);
			}
			else if(mode==aeqd_mode.N_POLE)
			{
				lp.phi=Proj.HALFPI-c_rh;
				lp.lam=Math.Atan2(xy.x, -xy.y);
			}
			else
			{
				lp.phi=c_rh-Proj.HALFPI;
				lp.lam=Math.Atan2(xy.x, xy.y);
			}

			return lp;
		}

		public override PJ Init()
		{
			g=new geod_geodesic(a, es/(1+Math.Sqrt(one_es)));

			phi0=Proj.pj_param_r(ctx, parameters, "lat_0");
			if(Math.Abs(Math.Abs(phi0)-Proj.HALFPI)<EPS10)
			{
				mode=phi0<0.0?aeqd_mode.S_POLE:aeqd_mode.N_POLE;
				sinph0=phi0<0.0?-1.0:1.0;
				cosph0=0.0;
			}
			else if(Math.Abs(phi0)<EPS10)
			{
				mode=aeqd_mode.EQUIT;
				sinph0=0.0;
				cosph0=1.0;
			}
			else
			{
				mode=aeqd_mode.OBLIQ;
				sinph0=Math.Sin(phi0);
				cosph0=Math.Cos(phi0);
			}

			if(es==0)
			{
				inv=s_inverse;
				fwd=s_forward;
			}
			else
			{
				en=Proj.pj_enfn(es);
				if(en==null) return null;
				if(Proj.pj_param_b(ctx, parameters, "guam"))
				{
					M1=Proj.pj_mlfn(phi0, sinph0, cosph0, en);
					inv=e_guam_inv;
					fwd=e_guam_fwd;
				}
				else
				{
					switch(mode)
					{
						case aeqd_mode.N_POLE:
							Mp=Proj.pj_mlfn(Proj.HALFPI, 1.0, 0.0, en);
							break;
						case aeqd_mode.S_POLE:
							Mp=Proj.pj_mlfn(-Proj.HALFPI, -1.0, 0.0, en);
							break;
						case aeqd_mode.EQUIT:
						case aeqd_mode.OBLIQ:
							inv=e_inverse;
							fwd=e_forward;
							N1=1.0/Math.Sqrt(1.0-es*sinph0*sinph0);
							He=e/Math.Sqrt(one_es);
							G=sinph0*He;
							He*=cosph0;
							break;
					}
					inv=e_inverse;
					fwd=e_forward;
				}
			}

			return this;
		}
	}
}
