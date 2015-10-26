// libproj -- library of cartographic projections
//
// Copyright (c) 2008   Gerald I. Evenden
// Copyright (c) 2009-2015 by the Authors
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// The code in this file is largly based upon procedures:
//
// Written by: Knud Poder and Karsten Engsager
//
// Base upon math from: R.Koenig and K.H. Weise, "Mathematische
// Grundlagen der hoeheren Geodaesie und Kartographie,
// Springer-Verlag, Berlin/Goettingen" Heidelberg, 1951.
//
// Modified and used here by permission of Reference Networks
// Division, Kort og Matrikelstyrelsen (KMS), Copenhagen, Denmark

using System;
using System.Text;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_etmerc : PJ
	{
		protected double Qn;					// Merid. quad., scaled to the projection
		protected double Zb;					// Radius vector in polar coord. systems
		protected double[] cgb=new double[6];	// Constants for Gauss => Geo lat
		protected double[] cbg=new double[6];	// Constants for Geo lat => Gauss
		protected double[] utg=new double[6];	// Constants for transv. merc. => geo
		protected double[] gtu=new double[6];	// Constants for geo > transv. merc.

		public override string Name { get { return "etmerc"; } }
		public override string DescriptionName { get { return "Extended Transverse Mercator"; } }
		public override string DescriptionType { get { return "Cyl, Ell"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double M_PI_4=Math.PI/4;

		const int PROJ_ETMERC_ORDER=6;

		// Compute log(1+x) accurately
		static double log1py(double x)
		{
			double y=1+x, z=y-1;

			// Here's the explanation for this magic: y = 1 + z, exactly, and z
			// approx x, thus log(y)/z (which is nearly constant near z = 0) returns
			// a good approximation to the true log(1 + x)/x. The multiplication x *
			// (log(y)/z) introduces little additional error.
			return z==0?x:x*Math.Log(y)/z;
		}

		// Compute asinh(x) accurately
		static double asinhy(double x)
		{
			double y=Math.Abs(x); // Enforce odd parity

			y=log1py(y*(1+y/(Libc.hypot(1.0, y)+1)));

			return x<0?-y:y;
		}

		static double gatg(double[] p1, int len_p1, double B)
		{
			double h=0, h2=0;
			double cos_2B=2*Math.Cos(2*B);
			int p=len_p1;
			double h1=p1[--p];
			while(p>0)
			{
				h=-h2+cos_2B*h1+p1[--p];
				h2=h1;
				h1=h;
			}

			return B+h*Math.Sin(2*B);
		}

		static double clenS(double[] a, int size, double arg_r, double arg_i, out double R, out double I)
		{
			// arguments
			int p=size;
			double sin_arg_r=Math.Sin(arg_r);
			double cos_arg_r=Math.Cos(arg_r);
			double sinh_arg_i=Math.Sinh(arg_i);
			double cosh_arg_i=Math.Cosh(arg_i);
			double r=2*cos_arg_r*cosh_arg_i;
			double i=-2*sin_arg_r*sinh_arg_i;

			// summation loop
			double hr=a[--p], hr1=0, hi=0, hi1=0;
			double hr2, hi2;
			while(p>0)
			{
				hr2=hr1;
				hi2=hi1;
				hr1=hr;
				hi1=hi;
				hr=-hr2+r*hr1-i*hi1+a[--p];
				hi=-hi2+i*hr1+r*hi1;
			}
			r=sin_arg_r*cosh_arg_i;
			i=cos_arg_r*sinh_arg_i;
			R=r*hr-i*hi;
			I=r*hi+i*hr;
			return R;
		}

		static double clens(double[] a, int size, double arg_r)
		{
			int p=size;
			double cos_arg_r=Math.Cos(arg_r);
			double r=2*cos_arg_r;

			// summation loop
			double hr=a[--p], hr1=0, hr2;
			while(p>0)
			{
				hr2=hr1;
				hr1=hr;
				hr=-hr2+r*hr1+a[--p];
			}
			return Math.Sin(arg_r)*hr;
		}

		// ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double Cn=lp.phi, Ce=lp.lam;

			// ell. LAT, LNG => Gaussian LAT, LNG
			Cn=gatg(cbg, PROJ_ETMERC_ORDER, Cn);

			// Gaussian LAT, LNG => compl. sph. LAT
			double sin_Cn=Math.Sin(Cn);
			double cos_Cn=Math.Cos(Cn);
			double cos_Ce=Math.Cos(Ce);
			Cn=Math.Atan2(sin_Cn, cos_Ce*cos_Cn);
			Ce=Math.Atan2(Math.Sin(Ce)*cos_Cn, Libc.hypot(sin_Cn, cos_Cn*cos_Ce));

			// compl. sph. N, E -> ell. norm. N, E
			Ce=asinhy(Math.Tan(Ce)); // Replaces: Ce=Math.Log(Math.Tan(M_PI_4+Ce*0.5));
			double dCn, dCe;
			Cn+=clenS(gtu, PROJ_ETMERC_ORDER, 2*Cn, 2*Ce, out dCn, out dCe);
			Ce+=dCe;
			if(Math.Abs(Ce)<=2.623395162778)
			{ // 150 degrees
				xy.y=Qn*Cn+Zb;	// Northing
				xy.x=Qn*Ce;		// Easting
			}
			else xy.x=xy.y=Libc.HUGE_VAL;

			return xy;
		}

		// ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double Cn=xy.y, Ce=xy.x;

			// normalize N, E
			Cn=(Cn-Zb)/Qn;
			Ce=Ce/Qn;
			if(Math.Abs(Ce)<=2.623395162778)
			{ // 150 degrees
				// norm. N, E -> compl. sph. LAT, LNG
				double dCn, dCe;
				Cn+=clenS(utg, PROJ_ETMERC_ORDER, 2*Cn, 2*Ce, out dCn, out dCe);
				Ce+=dCe;
				Ce=Math.Atan(Math.Sinh(Ce)); // Replaces: Ce=2*(Math.Atan(Math.Exp(Ce))-M_PI_4);

				// compl. sph. LAT -> Gaussian LAT, LNG
				double cos_Cn=Math.Cos(Cn);
				double sin_Ce=Math.Sin(Ce);
				double cos_Ce=Math.Cos(Ce);
				Ce=Math.Atan2(sin_Ce, cos_Ce*cos_Cn);
				Cn=Math.Atan2(Math.Sin(Cn)*cos_Ce, Libc.hypot(sin_Ce, cos_Ce*cos_Cn));
	
				// Gaussian LAT, LNG -> ell. LAT, LNG
				lp.phi=gatg(cgb, PROJ_ETMERC_ORDER, Cn);
				lp.lam=Ce;
			}
			else lp.phi=lp.lam=Libc.HUGE_VAL;

			return lp;
		}

		// general initialization
		protected PJ_etmerc setup()
		{
			if(es<=0) { Proj.pj_ctx_set_errno(ctx, -34); return null; }

			double f=es/(1+Math.Sqrt(1-es)); // Replaces: f=1-Math.Sqrt(1-es);

			// third flattening
			double n=f/(2-f);
			double np=n;

			// COEF. OF TRIG SERIES GEO <-> GAUSS
			// cgb := Gaussian -> Geodetic, KW p190 - 191 (61) - (62)
			// cbg := Geodetic -> Gaussian, KW p186 - 187 (51) - (52)
			// PROJ_ETMERC_ORDER = 6th degree : Engsager and Poder: ICC2007
			cgb[0]=n*(2+n*(-2/3.0+n*(-2+n*(116/45.0+n*(26/45.0+n*(-2854/675.0))))));
			cbg[0]=n*(-2+n*(2/3.0+n*(4/3.0+n*(-82/45.0+n*(32/45.0+n*(4642/4725.0))))));
			np*=n;
			cgb[1]=np*(7/3.0+n*(-8/5.0+n*(-227/45.0+n*(2704/315.0+n*(2323/945.0)))));
			cbg[1]=np*(5/3.0+n*(-16/15.0+n*(-13/9.0+n*(904/315.0+n*(-1522/945.0)))));
			np*=n;

			// n^5 coeff corrected from 1262/105 -> -1262/105
			cgb[2]=np*(56/15.0+n*(-136/35.0+n*(-1262/105.0+n*(73814/2835.0))));
			cbg[2]=np*(-26/15.0+n*(34/21.0+n*(8/5.0+n*(-12686/2835.0))));
			np*=n;

			// n^5 coeff corrected from 322/35 -> 332/35
			cgb[3]=np*(4279/630.0+n*(-332/35.0+n*(-399572/14175.0)));
			cbg[3]=np*(1237/630.0+n*(-12/5.0+n*(-24832/14175.0)));
			np*=n;
			cgb[4]=np*(4174/315.0+n*(-144838/6237.0));
			cbg[4]=np*(-734/315.0+n*(109598/31185.0));
			np*=n;
			cgb[5]=np*(601676/22275.0);
			cbg[5]=np*(444337/155925.0);

			// Constants of the projections
			// Transverse Mercator (UTM, ITM, etc)
			np=n*n;

			// Norm. mer. quad, K&W p.50 (96), p.19 (38b), p.5 (2)
			Qn=k0/(1+n)*(1+np*(1/4.0+np*(1/64.0+np/256.0)));

			// coef of trig series
			// utg := ell. N, E -> sph. N, E,  KW p194 (65)
			// gtu := sph. N, E -> ell. N, E,  KW p196 (69)
			utg[0]=n*(-0.5+n*(2/3.0+n*(-37/96.0+n*(1/360.0+n*(81/512.0+n*(-96199/604800.0))))));
			gtu[0]=n*(0.5+n*(-2/3.0+n*(5/16.0+n*(41/180.0+n*(-127/288.0+n*(7891/37800.0))))));
			utg[1]=np*(-1/48.0+n*(-1/15.0+n*(437/1440.0+n*(-46/105.0+n*(1118711/3870720.0)))));
			gtu[1]=np*(13/48.0+n*(-3/5.0+n*(557/1440.0+n*(281/630.0+n*(-1983433/1935360.0)))));
			np*=n;
			utg[2]=np*(-17/480.0+n*(37/840.0+n*(209/4480.0+n*(-5569/90720.0))));
			gtu[2]=np*(61/240.0+n*(-103/140.0+n*(15061/26880.0+n*(167603/181440.0))));
			np*=n;
			utg[3]=np*(-4397/161280.0+n*(11/504.0+n*(830251/7257600.0)));
			gtu[3]=np*(49561/161280.0+n*(-179/168.0+n*(6601661/7257600.0)));
			np*=n;
			utg[4]=np*(-4583/161280.0+n*(108847/3991680.0));
			gtu[4]=np*(34729/80640.0+n*(-3418889/1995840.0));
			np*=n;
			utg[5]=np*(-20648693/638668800.0);
			gtu[5]=np*(212378941/319334400.0);

			// Gaussian latitude value of the origin latitude
			double Z=gatg(cbg, PROJ_ETMERC_ORDER, phi0);

			// Origin northing minus true northing at the origin latitude
			// i.e. true northing = N - P->Zb
			Zb=-Qn*(Z+clens(gtu, PROJ_ETMERC_ORDER, 2*Z));

			inv=e_inverse;
			fwd=e_forward;

			return this;
		}

		public override PJ Init()
		{
			return setup();
		}
	}

	class PJ_utm : PJ_etmerc
	{
		public override string Name { get { return "utm"; } }
		public override string DescriptionName { get { return "Universal Transverse Mercator (UTM)"; } }
		public override string DescriptionParameters { get { return "zone= south"; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret = new StringBuilder();

				int zone = (int)Math.Floor((Proj.adjlon(lam0) + Proj.PI) * 30.0 / Proj.PI);
				if (zone < 0) zone = 0;
				else if (zone >= 60) zone = 59;

				ret.AppendFormat(" +zone={0}", zone + 1);

				if (y0 == 10000000.0) ret.Append(" +south");
				return ret.ToString();
			}
		}

		// utm uses etmerc for the underlying projection
		public override PJ Init()
		{
			if (es == 0) { Proj.pj_ctx_set_errno(ctx, -34); return null; }

			y0 = Proj.pj_param_b(ctx, parameters, "south") ? 10000000.0 : 0.0;
			x0 = 500000.0;

			int zone;

			if (Proj.pj_param_t(ctx, parameters, "zone")) // zone input ?
			{
				zone = Proj.pj_param_i(ctx, parameters, "zone");
				if (zone > 0 && zone <= 60) zone--;
				else { Proj.pj_ctx_set_errno(ctx, -35); return null; }
			}
			else // nearest central meridian input
			{
				zone = (int)Math.Floor((Proj.adjlon(lam0) + Proj.PI) * 30.0 / Proj.PI);
				if (zone < 0) zone = 0;
				else if (zone >= 60) zone = 59;
			}

			lam0 = (zone + 0.5) * Proj.PI / 30.0 - Proj.PI;
			k0 = 0.9996;
			phi0 = 0.0;

			return setup();
		}
	}
}
