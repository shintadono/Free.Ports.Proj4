namespace Free.Ports.Proj4.ComplexPoly
{
	public struct COMPLEX
	{
		public double r, i;

		public COMPLEX(double r, double i)
		{
			this.r=r;
			this.i=i;
		}

		// evaluate complex polynomial

		// note: coefficients are always from C_1 to C_n
		// i.e. C_0 == (0., 0)
		// n should always be >= 1 though no checks are made
		public COMPLEX pj_zpoly1(COMPLEX[] C, int n)
		{
			int C_ind=n;
			COMPLEX a=C[C_ind];

			double t;
			while((n--)>0)
			{
				t=a.r;
				a.r=C[--C_ind].r+r*t-i*a.i;
				a.i=C[C_ind].i+r*a.i+i*t;
			}

			t=a.r;
			a.r=r*t-i*a.i;
			a.i=r*a.i+i*t;

			return a;
		}

		// evaluate complex polynomial and derivative
		public COMPLEX pj_zpolyd1(COMPLEX[] C, int n, out COMPLEX der)
		{
			int C_ind=n;
			double t;
			bool first = true;

			COMPLEX a = C[C_ind];
			COMPLEX b = a;
			while ((n--) > 0)
			{
				if(first)
				{
					first=false;
				}
				else
				{
					t=b.r;

					b.r=a.r+r*t-i*b.i;
					b.i=a.i+r*b.i+i*t;
				}

				t=a.r;
				a.r=C[--C_ind].r+r*t-i*a.i;
				a.i=C[C_ind].i+r*a.i+i*t;
			}

			t=b.r;
			b.r=a.r+r*t-i*b.i;
			b.i=a.i+r*b.i+i*t;

			t=a.r;
			a.r=r*t-i*a.i;
			a.i=r*a.i+i*t;

			der=b;

			return a;
		}
	}
}
