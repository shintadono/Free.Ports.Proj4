using System;
using System.Collections.Generic;

namespace Free.Ports.Proj4.Geodesic
{
	public partial class Geod
	{
		static string emess(string fmt, params object[] args)
		{
			return Proj.pj_get_release() + "\r\n" + string.Format(fmt, args);
		}

		public void Set(string[] argv)
		{
			List<string> start = new List<string>();

			// put arguments into internal linked list
			if (argv.Length <= 0) throw new Exception(emess("no arguments in initialization list"));
			for (int i = 0; i < argv.Length; ++i) start.Add(Proj.pj_mkparam(argv[i]));

			// set elliptical parameters
			double es;
			if (Proj.pj_ell_set(Proj.pj_get_default_ctx(), start, out geod_a, out es)) throw new Exception(emess("ellipse setup failure"));

			// set units
			string name = Proj.pj_param_s(null, start, "units");
			if (name != null && name != "")
			{
				to_meter = Proj.GetUnitFactor(name);
				if (double.IsNaN(to_meter)) throw new Exception(emess("{0} unknown unit conversion id", name));
				fr_meter = 1.0 / to_meter;
			}
			else to_meter = fr_meter = 1.0;

			geod_f = es / (1 + Math.Sqrt(1 - es));
			Ini();

			// check if line or arc mode
			if (Proj.pj_param_t(null, start, "lat_1"))
			{
				phi1 = Proj.pj_param_r(null, start, "lat_1");
				lam1 = Proj.pj_param_r(null, start, "lon_1");
				if (Proj.pj_param_t(null, start, "lat_2"))
				{
					phi2 = Proj.pj_param_r(null, start, "lat_2");
					lam2 = Proj.pj_param_r(null, start, "lon_2");
					Inv();
					Pre();
				}
				else
				{
					geod_S = Proj.pj_param_d(null, start, "S");
					if (geod_S != 0.0)
					{
						al12 = Proj.pj_param_r(null, start, "A");
						Pre();
						For();
					}
					else throw new Exception(emess("incomplete geodesic/arc info"));
				}

				n_alpha = Proj.pj_param_i(null, start, "n_A");
				if (n_alpha > 0)
				{
					del_alpha = Proj.pj_param_r(null, start, "del_A");
					if (del_alpha == 0) throw new Exception(emess("del azimuth == 0"));
				}
				else
				{
					double del_S = Math.Abs(Proj.pj_param_d(null, start, "del_S"));
					if (del_S != 0)
					{
						n_S = (int)(geod_S / del_S + 0.5);
					}
					else
					{
						n_S = Proj.pj_param_i(null, start, "n_S");
						if (n_S <= 0) throw new Exception(emess("no interval divisor selected"));
					}
				}
			}

			start.Clear();
		}
	}
}
