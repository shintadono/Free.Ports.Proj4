using System.Collections.Generic;
using System.Globalization;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		//***********************************************************************
		//								pj_parsedate()
		//
		//		Parse a date into a floating point year value. Acceptable
		//		values are "yyyy.fraction" and "yyyy-mm-dd". Anything else
		//		returns 0.0.
		//***********************************************************************
		public static double pj_parsedate(projCtx ctx, string date_string)
		{
			try
			{
				if (date_string.Length == 10 && date_string[4] == '-' && date_string[7] == '-')
				{
					int year = int.Parse(date_string.Substring(0, 4));
					int month = int.Parse(date_string.Substring(5, 2));
					int day = int.Parse(date_string.Substring(8, 2));

					// simplified calculation so we don't need to know all about months
					return year + ((month - 1) * 31 + (day - 1)) / 372.0;
				}

				return double.Parse(date_string, nc);
			}
			catch
			{
				return 0;
			}
		}
	}
}
