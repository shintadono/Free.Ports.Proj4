using System;

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		// Convert radian argument to DMS ascii format

		// RES is fractional second figures
		// RES60 = 60 * RES
		// CONV = 180 * 3600 * RES / PI (radians to RES seconds)
		static double RES=1000.0, RES60=60000.0, CONV=206264806.24709635516;

		static string format="{0}d{1}'{2:0.###}\"{3}";
		static bool dolong=false;

		public static void set_rtodms(int fract, bool con_w)
		{
			if(fract>=0&&fract<9)
			{
				RES=1.0;

				// following not very elegant, but used infrequently
				for(int i=0; i<fract; ++i) RES*=10.0;

				RES60=RES*60.0;
				CONV=180.0*3600.0*RES/PI;

				if(fract==0)
				{
					if(!con_w) format=string.Format("{{0}}d{{1}}'{{2}}\"{{3}}");
					else format=string.Format("{{0}}d{{1:00}}'{{2:00}}\"{{3}}");
				}
				else
				{
					string rauten="".PadRight(fract, '#');

					if(!con_w) format=string.Format("{{0}}d{{1}}'{{2:0.{0}}}\"{{3}}", rauten);
					else format=string.Format("{{0}}d{{1:00}}'{{2:00.{0}}}\"{{3}}", rauten);
				}

				dolong=con_w;
			}
		}

		public static string rtodms(double r, char pos, char neg)
		{
			string sign="";
			int deg, min;
			double sec;

			string ss="";

			if(r<0)
			{
				r=-r;
				if(pos!='\0'&&neg!='\0') sign+=neg;
				else ss+='-';
			}
			else
			{
				if(pos!='\0') sign+=pos;
			}

			r=Math.Floor(r*CONV+0.5);
			sec=(r/RES)%60.0;
			r=Math.Floor(r/RES60);
			min=(int)(r%60.0);
			r=Math.Floor(r/60.0);
			deg=(int)r;

			if(dolong) ss+=string.Format(nc, format, deg, min, sec, sign);
			else if(sec!=0.0) ss+=string.Format(nc, format, deg, min, sec, sign);
			else if(min!=0) ss+=string.Format(nc, "{0}d{1}'{2}", deg, min, sign);
			else ss+=string.Format(nc, "{0}d{1}", deg, sign);

			return ss;
		}
	}
}