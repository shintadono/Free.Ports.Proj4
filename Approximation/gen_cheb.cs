#if USE_NOT_LIB_STUFF
using System;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Approximation
{
	public static partial class Approximation
	{
		// generates 'T' option output
		delegate double InputFunc(string in_string, out string rs);

		public static void gen_cheb(bool inverse, ProjFunc proj, string s, PJ P, int iargc, string[] iargv)
		{
			InputFunc input = Proj.dmstor;
			if (inverse) input = Libc.strtod;

			int errin = 0;
			projUV low, upp;
			low.u = low.v = upp.u = upp.v = 0;
			if (s.Length > 0) low.u = input(s, out s); else { low.u = 0; errin++; }
			if (s.Length > 1 && s[0] == ',') upp.u = input(s.Substring(1), out s); else { upp.u = 0; errin++; }
			if (s.Length > 1 && s[0] == ',') low.v = input(s.Substring(1), out s); else { low.v = 0; errin++; }
			if (s.Length > 1 && s[0] == ',') upp.v = input(s.Substring(1), out s); else { upp.v = 0; errin++; }

			if (errin != 0) emess(16, "null or absent -T parameters");

			int NU = 15, NV = 15, res = -1;

			if (s.Length > 1 && s[0] == ',')
			{
				s = s.Substring(1);
				if (s[0] != ',') res = Libc.strtol(s, out s, 10);
			}
			if (s.Length > 1 && s[0] == ',')
			{
				s = s.Substring(1);
				if (s[0] != ',') NU = Libc.strtol(s, out s, 10);
			}
			if (s.Length > 1 && s[0] == ',')
			{
				s = s.Substring(1);
				if (s[0] != ',') NV = Libc.strtol(s, out s, 10);
			}

			bool pwr = s.Length > 0 && s.StartsWith(",P");
			Console.WriteLine("#proj_{0}\n#\trun-line:", pwr ? "Power" : "Chebyshev");

			// proj execution audit trail
			if (iargc > 0)
			{
				int n = 0;

				for (int i = 0; i < iargc; i++)
				{
					string arg = iargv[i];
					if (arg[0] != '+')
					{
						if (n == 0)
						{
							Console.Write("#");
							n++;
						}
						Console.Write(" " + arg);
						n += arg.Length + 1;
						if (n > 50)
						{
							Console.WriteLine(); n = 0;
						}
					}
				}
				if (n != 0) Console.WriteLine();
			}

			Console.WriteLine("# projection parameters");
			Console.WriteLine("#" + P.DescriptionName);
			Console.WriteLine("#" + P.DescriptionParameters);
			Console.WriteLine("#" + P.DescriptionType);
			Console.WriteLine("#" + P.ToProj4String());

			if (low.u == upp.u || low.v >= upp.v) emess(16, "approx. argument range error");
			if (low.u > upp.u) low.u -= Proj.TWOPI;
			if (NU < 2 || NV < 2) emess(16, "approx. work dimensions ({0} {1}) too small", NU, NV);

			projUV resid = new projUV();
			Tseries F = mk_cheby(low, upp, Math.Pow(10.0, res) * 0.5, ref resid, proj, NU, NV, pwr);
			if (F == null) emess(16, "generation of approx failed\nreason: {0}\n", Proj.pj_strerrno(Libc.errno));

			Console.WriteLine("{0},{1:G12},{2:G12},{3:G12},{4:G12},{5:G12}", inverse ? 'I' : 'F', P.lam0 * Proj.RAD_TO_DEG,
				low.u * (inverse ? 1.0 : Proj.RAD_TO_DEG), upp.u * (inverse ? 1.0 : Proj.RAD_TO_DEG),
				low.v * (inverse ? 1.0 : Proj.RAD_TO_DEG), upp.v * (inverse ? 1.0 : Proj.RAD_TO_DEG));

			string fmt;
			if (pwr) fmt = "G15";
			else if (res <= 0) fmt = string.Format("F{0}", -res + 1);
			else fmt = "F0";

			p_series(F, Console.Out, fmt);
			Console.WriteLine("# |u,v| sums {0} {1}\n#end_proj_{2}", resid.u, resid.v, pwr ? "Power" : "Chebyshev");
		}
	}
}
#endif
