using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		static readonly double[] vm = new double[]
		{
			DEG_TO_RAD,
			0.0002908882086657216,
			0.0000048481368110953599
		};

		public static double dmstor(string nptr, int stuff) // for: ...=dmstor(s, 0);
		{
			string dummy="";
			return dmstor_ctx(pj_get_default_ctx(), nptr, out dummy);
		}

		public static double dmstor_ctx(projCtx ctx, string nptr, int stuff) // for: ...=dmstor(s, 0);
		{
			string dummy="";
			return dmstor_ctx(ctx, nptr, out dummy);
		}

		public static double dmstor(string nptr)
		{
			string dummy="";
			return dmstor_ctx(pj_get_default_ctx(), nptr, out dummy);
		}

		public static double dmstor_ctx(projCtx ctx, string nptr)
		{
			string dummy="";
			return dmstor_ctx(ctx, nptr, out dummy);
		}

		public static double dmstor(string nptr, out string newptr)
		{
			return dmstor_ctx(pj_get_default_ctx(), nptr, out newptr);
		}

		public static double dmstor_ctx(projCtx ctx, string nptr, out string newptr)
		{
			const string sym="NnEeSsWw";

			newptr=nptr;
			if(nptr==null)
			{
				pj_ctx_set_errno(ctx, -16);
				return Libc.HUGE_VAL;
			}

			string s=nptr.TrimStart();
			if(s.Length==0)
			{
				newptr="";
				pj_ctx_set_errno(ctx, -16);
				return Libc.HUGE_VAL;
			}

			newptr=s;

			string work="";
			foreach(char c in s)
			{
				if((c>='0'&&c<='9')||(c>='a'&&c<='z')||(c>='A'&&c<='Z')||c=='.'||c=='°'||c=='\''||c=='"'||c=='+'||c=='-')
				{
					work+=c;
				}
				else break;
			}

			s=work;
			char sign=work[0];
			if(sign=='+'||sign=='-') s=s.Substring(1);
			else sign='+';

			if((s.Length==0)||(!((s[0]>='0'&&s[0]<='9')||s[0]=='.')))
			{
				pj_ctx_set_errno(ctx, -16);
				return Libc.HUGE_VAL;
			}

			int n=0, p;
			double v=0, tv;

			for(int nl=0; nl<3; nl=n+1)
			{
				if(s.Length==0) break;

				char c=s[0];
				if(!((c>='0'&&c<='9')||c=='.')) break;

				tv=proj_strtod(s, out s);
				if(tv==Libc.HUGE_VAL) return tv;

				if(s.Length==0)
				{
					v+=tv*vm[nl];
					n=4;
					continue;
				}

				switch(s[0])
				{
					case '°':
					case 'D':
					case 'd':
						n=0; break;
					case '\'':
						n=1; break;
					case '"':
						n=2; break;
					case 'r':
					case 'R':
						if(nl!=0)
						{
							pj_ctx_set_errno(ctx, -16);
							return Libc.HUGE_VAL;
						}
						s=s.Substring(1);
						v=tv;
						n=4;
						continue;
					default:
						v+=tv*vm[nl];
						n=4;
						continue;
				}
				if(n<nl)
				{
					pj_ctx_set_errno(ctx, -16);
					return Libc.HUGE_VAL;
				}

				v+=tv*vm[n];
				s=s.Substring(1);
			}

			// postfix sign
			if(s.Length>0)
			{
				p=sym.IndexOf(s[0]);
				if(p>-1)
				{
					sign=p>=4?(sign=='-'?'+':'-'):sign;
					s=s.Substring(1);
				}
			}

			if(sign=='-') v=-v;

			// return point of next char after valid string
			if(newptr.Length<=(work.Length-s.Length)) newptr="";
			else newptr=newptr.Substring(work.Length-s.Length);
			return v;
		}

		static double proj_strtod(string nptr, out string newptr)
		{
			// Scan for characters which cause problems with VC++ strtod()
			int ind=nptr.IndexOfAny("dD".ToCharArray());
			if(ind!=-1)
			{
				string cp=nptr.Replace('d', '\0').Replace('D', '\0');
				double result=Libc.strtod(cp, out newptr);
				newptr=nptr.Substring(nptr.Length-newptr.Length);
				return result;
			}

			// no offending characters, just handle normally
			return Libc.strtod(nptr, out newptr);
		}
	}
}
