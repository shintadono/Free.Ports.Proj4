using System;

namespace Free.Ports.Proj4.Geodesic
{
	/// <summary>
	/// The struct containing information about the ellipsoid. This must be
	/// initialized by geod_init() before use.
	/// </summary>
	public class geod_geodesic : geodesic
	{
		/// <summary>
		/// The equatorial radius.
		/// </summary>
		public double a;

		/// <summary>
		/// The flattening.
		/// </summary>
		public double f;

		public double f1, e2, ep2, n, b, c2, etol2;
		public double[] A3x = new double[6];
		public double[] C3x = new double[15];
		public double[] C4x = new double[21];

		/// <summary>
		/// Initialize a geod_geodesic object.
		/// </summary>
		/// <param name="a">The equatorial radius (meters).</param>
		/// <param name="f">The flattening.</param>
		public geod_geodesic(double a, double f)
		{
			this.a = a;
			this.f = f;
			this.f1 = 1 - this.f;
			this.e2 = this.f * (2 - this.f);
			this.ep2 = this.e2 / sq(this.f1); // e2/(1-e2)
			this.n = this.f / (2 - this.f);
			this.b = this.a * this.f1;
			this.c2 = (sq(this.a) + sq(this.b) * (this.e2 == 0 ? 1 : (this.e2 > 0 ? atanhx(Math.Sqrt(this.e2)) : Math.Atan(Math.Sqrt(-this.e2))) / Math.Sqrt(Math.Abs(this.e2)))) / 2; // authalic radius squared

			// The sig12 threshold for "really short". Using the auxiliary sphere
			// solution with dnm computed at (bet1+bet2)/2, the relative error in the
			// azimuth consistency check is sig12^2*abs(f)*min(1, 1-f/2)/2. (Error
			// measured for 1/100<b/a<100 and abs(f)>=1/1000. For a given f and
			// sig12, the max error occurs for lines near the pole. If the old rule for
			// computing dnm=(dn1+dn2)/2 is used, then the error increases by a
			// factor of 2.) Setting this equal to epsilon gives sig12=etol2. Here
			// 0.1 is a safety factor (error decreased by 100) and max(0.001, abs(f))
			// stops etol2 getting too large in the nearly spherical case.
			this.etol2 = 0.1 * tol2 / Math.Sqrt(maxx(0.001, Math.Abs(this.f)) * minx(1, 1 - this.f / 2) / 2);

			A3coeff();
			C3coeff();
			C4coeff();
		}

		/// <summary>
		/// Solve the direct geodesic problem.
		/// </summary>
		/// <param name="lat1">The atitude of point 1 (degrees).</param>
		/// <param name="lon1">The longitude of point 1 (degrees).</param>
		/// <param name="azi1">The azimuth at point 1 (degrees).</param>
		/// <param name="s12_a12">The distance from point 1 to point 2 (meters); it can be negative.</param>
		/// <param name="lat2">The latitude of point 2 (degrees).</param>
		/// <param name="lon2">The longitude of point 2 (degrees).</param>
		/// <param name="azi2">The (forward) azimuth at point 2 (degrees).</param>
		/// <remarks>
		/// g must have been initialized with a call to geod_init(). lat1
		/// should be in the range [-90deg, 90deg]. The values of lon2
		/// and azi2 returned are in the range [-180deg, 180deg).
		///
		/// If either point is at a pole, the azimuth is defined by keeping the
		/// longitude fixed, writing lat=+/-90deg - epsilon), and
		/// taking the limit epsilon => 0+. An arc length greater that 180deg
		/// signifies a geodesic which is not a shortest path. (For a prolate
		/// ellipsoid, an additional condition is necessary for a shortest path: the
		/// longitudinal extent must not exceed of 180deg.)
		///</remarks>
		/// <returns>a12 arc length of between point 1 and point 2 (degrees).</returns>
		public void geod_direct(double lat1, double lon1, double azi1, double s12, out double lat2, out double lon2, out double azi2)
		{
			double dummy;
			geod_gendirect(lat1, lon1, azi1, GEOD.NOFLAGS, s12, out lat2, out lon2, out azi2, out dummy, out dummy, out dummy, out dummy, out dummy);
		}

		/// <summary>
		/// Solve the direct geodesic problem.
		/// </summary>
		/// <param name="lat1">The atitude of point 1 (degrees).</param>
		/// <param name="lon1">The longitude of point 1 (degrees).</param>
		/// <param name="azi1">The azimuth at point 1 (degrees).</param>
		/// <param name="s12_a12">The distance from point 1 to point 2 (meters); it can be negative.</param>
		/// <param name="lat2">The latitude of point 2 (degrees).</param>
		/// <param name="lon2">The longitude of point 2 (degrees).</param>
		/// <remarks>
		/// g must have been initialized with a call to geod_init(). lat1
		/// should be in the range [-90deg, 90deg]. The values of lon2
		/// and azi2 returned are in the range [-180deg, 180deg).
		///
		/// If either point is at a pole, the azimuth is defined by keeping the
		/// longitude fixed, writing lat=+/-90deg - epsilon), and
		/// taking the limit epsilon => 0+. An arc length greater that 180deg
		/// signifies a geodesic which is not a shortest path. (For a prolate
		/// ellipsoid, an additional condition is necessary for a shortest path: the
		/// longitudinal extent must not exceed of 180deg.)
		///</remarks>
		/// <returns>a12 arc length of between point 1 and point 2 (degrees).</returns>
		public void geod_direct(double lat1, double lon1, double azi1, double s12, out double lat2, out double lon2)
		{
			double dummy;
			geod_gendirect(lat1, lon1, azi1, GEOD.NOFLAGS, s12, out lat2, out lon2, out dummy, out dummy, out dummy, out dummy, out dummy, out dummy);
		}

		/// <summary>
		/// The general direct geodesic problem.
		/// </summary>
		/// <param name="lat1">The atitude of point 1 (degrees).</param>
		/// <param name="lon1">The longitude of point 1 (degrees).</param>
		/// <param name="azi1">The azimuth at point 1 (degrees).</param>
		/// <param name="flags">Bitor'ed combination of <see cref="GEOD"/> flags;
		/// GEOD.ARCMODE determines the meaning of s12_a12 and
		/// GEOD.LONG_UNROLL "unrolls" lon2.</param>
		/// <param name="s12_a12">If flags&amp;GEOD.ARCMODE is 0, this is the distance
		/// from point 1 to point 2 (meters); otherwise it is the arc length
		/// from point 1 to point 2 (degrees); it can be negative.</param>
		/// <param name="plat2">The latitude of point 2 (degrees).</param>
		/// <param name="plon2">The longitude of point 2 (degrees).</param>
		/// <param name="pazi2">The (forward) azimuth at point 2 (degrees).</param>
		/// <param name="ps12">The distance from point 1 to point 2 (meters).</param>
		/// <param name="pm12">The reduced length of geodesic (meters).</param>
		/// <param name="pM12">The geodesic scale of point 2 relative to point 1 (dimensionless).</param>
		/// <param name="pM21">The geodesic scale of point 1 relative to point 2 (dimensionless).</param>
		/// <param name="pS12">The area under the geodesic (square meters).</param>
		/// <remarks>
		/// g must have been initialized with a call to geod_init(). lat1
		/// should be in the range [-90deg, 90deg]. The function
		/// value a12 equals s12_a12 if flags &amp; GEOD.ARCMODE.
		///
		/// With <paramref name="flags"/> &amp; GEOD.LONG_UNROLL bit set, the longitude is "unrolled" so
		/// that the quantity lon2-lon1 indicates how many times and in
		/// what sense the geodesic encircles the ellipsoid.
		///</remarks>
		/// <returns>a12 arc length of from point 1 to point 2 (degrees).</returns>
		public double geod_gendirect(double lat1, double lon1, double azi1, GEOD flags, double s12_a12,
			out double plat2, out double plon2, out double pazi2, out double ps12, out double pm12, out double pM12, out double pM21, out double pS12)
		{
			GEOD outmask = GEOD.LATITUDE | GEOD.LONGITUDE | GEOD.AZIMUTH | GEOD.DISTANCE | GEOD.REDUCEDLENGTH | GEOD.GEODESICSCALE | GEOD.AREA;

			geod_geodesicline l = new geod_geodesicline(this, lat1, lon1, azi1, outmask | ((flags & GEOD.ARCMODE) != 0 ? GEOD.NONE : GEOD.DISTANCE_IN)); // Automatically supply GEOD_DISTANCE_IN if necessary
			return l.geod_genposition(flags, s12_a12, out plat2, out plon2, out pazi2, out ps12, out pm12, out pM12, out pM21, out pS12);
		}

		/// <summary>
		/// Solve the inverse geodesic problem.
		/// </summary>
		/// <param name="lat1">The latitude of point 1 (degrees).</param>
		/// <param name="lon1">The longitude of point 1 (degrees).</param>
		/// <param name="lat2">The latitude of point 2 (degrees).</param>
		/// <param name="lon2">The longitude of point 2 (degrees).</param>
		/// <param name="ps12">The distance from point 1 to point 2 (meters).</param>
		/// <param name="pazi1">The azimuth at point 1 (degrees).</param>
		/// <param name="pazi2">The (forward) azimuth at point 2 (degrees).</param>
		/// <remarks>
		/// g must have been initialized with a call to geod_init. lat1
		/// and lat2 should be in the range [-90deg, 90deg]. The values of
		/// azi1 and azi2 returned are in the range [-180deg, 180deg).
		///
		/// If either point is at a pole, the azimuth is defined by keeping the
		/// longitude fixed, writing lat=+/-(90deg-epsilon), and
		/// taking the limit epsilon => 0+.
		///
		/// The solution to the inverse problem is found using Newton's method. If
		/// this fails to converge (this is very unlikely in geodetic applications
		/// but does occur for very eccentric ellipsoids), then the bisection method
		/// is used to refine the solution.
		/// </remarks>
		public void geod_inverse(double lat1, double lon1, double lat2, double lon2, out double s12, out double azi1, out double azi2)
		{
			double dummy;
			geod_geninverse(lat1, lon1, lat2, lon2, out s12, out azi1, out azi2, out dummy, out dummy, out dummy, out dummy);
		}

		/// <summary>
		/// Solve the inverse geodesic problem.
		/// </summary>
		/// <param name="lat1">The latitude of point 1 (degrees).</param>
		/// <param name="lon1">The longitude of point 1 (degrees).</param>
		/// <param name="lat2">The latitude of point 2 (degrees).</param>
		/// <param name="lon2">The longitude of point 2 (degrees).</param>
		/// <param name="ps12">The distance from point 1 to point 2 (meters).</param>
		/// <remarks>
		/// g must have been initialized with a call to geod_init. lat1
		/// and lat2 should be in the range [-90deg, 90deg].
		///
		/// If either point is at a pole, the azimuth is defined by keeping the
		/// longitude fixed, writing lat=+/-(90deg-epsilon), and
		/// taking the limit epsilon => 0+.
		///
		/// The solution to the inverse problem is found using Newton's method. If
		/// this fails to converge (this is very unlikely in geodetic applications
		/// but does occur for very eccentric ellipsoids), then the bisection method
		/// is used to refine the solution.
		/// </remarks>
		public void geod_inverse(double lat1, double lon1, double lat2, double lon2, out double s12)
		{
			double dummy;
			geod_geninverse(lat1, lon1, lat2, lon2, out s12, out dummy, out dummy, out dummy, out dummy, out dummy, out dummy);
		}

		internal double geod_geninverse_int(double lat1, double lon1, double lat2, double lon2,
			out double ps12,
			out double psalp1, out double pcalp1,
			out double psalp2, out double pcalp2,
			out double pm12, out double pM12, out double pM21, out double pS12)
		{
			double s12 = 0, m12 = 0, M12 = 0, M21 = 0, S12 = 0;
			double lon12, lon12s;
			int latsign, lonsign, swapp;
			double sbet1, cbet1, sbet2, cbet2, s12x = 0, m12x = 0;
			double dn1, dn2, lam12, slam12, clam12;
			double a12 = 0, sig12, calp1 = 0, salp1 = 0, calp2 = 0, salp2 = 0;

			// index zero elements of these arrays are unused
			double[] Ca = new double[nC];
			bool meridian;
			// somg12 > 1 marks that it needs to be calculated
			double omg12 = 0, somg12 = 2, comg12 = 0;

			GEOD outmask = GEOD.DISTANCE | GEOD.REDUCEDLENGTH | GEOD.GEODESICSCALE | GEOD.AREA;

			outmask &= GEOD.OUT_ALL;

			// Compute longitude difference (AngDiff does this carefully). Result is
			// in [-180, 180] but -180 is only for west-going geodesics. 180 is for
			// east-going and meridional geodesics.
			lon12 = AngDiff(lon1, lon2, out lon12s);

			// Make longitude difference positive.
			lonsign = lon12 >= 0 ? 1 : -1;

			// If very close to being on the same half-meridian, then make it so.
			lon12 = lonsign * AngRound(lon12);
			lon12s = AngRound((180 - lon12) - lonsign * lon12s);
			lam12 = lon12 * degree;
			if (lon12 > 90)
			{
				sincosdx(lon12s, out slam12, out clam12);
				clam12 = -clam12;
			}
			else
				sincosdx(lon12, out slam12, out clam12);

			// If really close to the equator, treat as on equator.
			lat1 = AngRound(LatFix(lat1));
			lat2 = AngRound(LatFix(lat2));

			// Swap points so that point with higher (abs) latitude is point 1
			// If one latitude is a nan, then it becomes lat1.
			swapp = Math.Abs(lat1) < Math.Abs(lat2) ? -1 : 1;
			if (swapp < 0)
			{
				lonsign *= -1;
				swapx(ref lat1, ref lat2);
			}

			// Make lat1 <= 0
			latsign = lat1 < 0 ? 1 : -1;
			lat1 *= latsign;
			lat2 *= latsign;

			// Now we have
			//
			//     0 <= lon12 <= 180
			//     -90 <= lat1 <= 0
			//     lat1 <= lat2 <= -lat1
			//
			// longsign, swapp, latsign register the transformation to bring the
			// coordinates to this canonical form.  In all cases, 1 means no change was
			// made. We make these transformations so that there are few cases to
			// check, e.g., on verifying quadrants in atan2. In addition, this
			// enforces some symmetries in the results returned.

			sincosdx(lat1, out sbet1, out cbet1); sbet1 *= f1;

			// Ensure cbet1 = +epsilon at poles
			norm2(ref sbet1, ref cbet1); cbet1 = maxx(tiny, cbet1);

			sincosdx(lat2, out sbet2, out cbet2); sbet2 *= f1;

			// Ensure cbet2 = +epsilon at poles
			norm2(ref sbet2, ref cbet2); cbet2 = maxx(tiny, cbet2);

			// If cbet1 < -sbet1, then cbet2 - cbet1 is a sensitive measure of the
			// |bet1| - |bet2|. Alternatively (cbet1 >= -sbet1), abs(sbet2) + sbet1 is
			// a better measure. This logic is used in assigning calp2 in Lambda12.
			// Sometimes these quantities vanish and in that case we force bet2 = +/-
			// bet1 exactly. An example where is is necessary is the inverse problem
			// 48.522876735459 0 -48.52287673545898293 179.599720456223079643
			// which failed with Visual Studio 10 (Release and Debug)
			if (cbet1 < -sbet1)
			{
				if (cbet2 == cbet1)
					sbet2 = sbet2 < 0 ? sbet1 : -sbet1;
			}
			else
			{
				if (Math.Abs(sbet2) == -sbet1)
					cbet2 = cbet1;
			}

			dn1 = Math.Sqrt(1 + ep2 * sq(sbet1));
			dn2 = Math.Sqrt(1 + ep2 * sq(sbet2));

			meridian = lat1 == -90 || slam12 == 0;

			if (meridian)
			{
				// Endpoints are on a single full meridian, so the geodesic might lie on
				// a meridian.

				double ssig1, csig1, ssig2, csig2;
				calp1 = clam12; salp1 = slam12; // Head to the target longitude
				calp2 = 1; salp2 = 0; // At the target we're heading north

				// tan(bet)=tan(sig)*cos(alp)
				ssig1 = sbet1; csig1 = calp1 * cbet1;
				ssig2 = sbet2; csig2 = calp2 * cbet2;

				// sig12=sig2-sig1
				sig12 = Math.Atan2(maxx(0.0, csig1 * ssig2 - ssig1 * csig2), csig1 * csig2 + ssig1 * ssig2);
				{
					double dummy;
					Lengths(n, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
							cbet1, cbet2, out s12x, out m12x, out dummy,
							out M12, out M21, Ca);
				}

				// Add the check for sig12 since zero length geodesics might yield m12 <0.
				// Test case was
				//
				//    echo 20.001 0 20.001 0 | GeodSolve -i
				//
				// In fact, we will have sig12 > pi/2 for meridional geodesic which is
				// not a shortest path.
				if (sig12 < 1 || m12x >= 0)
				{
					// Need at least 2, to handle 90 0 90 180
					if (sig12 < 3 * tiny) sig12 = m12x = s12x = 0;
					m12x *= b;
					s12x *= b;
					a12 = sig12 / degree;
				}
				else meridian = false; // m12 < 0, i.e., prolate and too close to anti-podal
			}

			if (!meridian &&
				sbet1 == 0 && // and sbet2==0
							  // Mimic the way Lambda12 works with calp1=0
				(f <= 0 || lon12s >= f * 180))
			{
				// Geodesic runs along equator
				calp1 = calp2 = 0; salp1 = salp2 = 1;
				s12x = a * lam12;
				sig12 = omg12 = lam12 / f1;
				m12x = b * Math.Sin(sig12);
				if ((outmask & GEOD.GEODESICSCALE) != 0) M12 = M21 = Math.Cos(sig12);
				a12 = lon12 / f1;
			}
			else if (!meridian)
			{
				// Now point1 and point2 belong within a hemisphere bounded by a
				// meridian and geodesic is neither meridional or equatorial.

				// Figure a starting point for Newton's method
				double dnm = 0;
				sig12 = InverseStart(sbet1, cbet1, dn1, sbet2, cbet2, dn2, lam12, slam12, clam12, out salp1, out calp1, out salp2, out calp2, out dnm, Ca);

				if (sig12 >= 0)
				{
					// Short lines (InverseStart sets salp2, calp2, dnm)
					s12x = sig12 * b * dnm;
					m12x = sq(dnm) * b * Math.Sin(sig12 / dnm);
					if ((outmask & GEOD.GEODESICSCALE) != 0) M12 = M21 = Math.Cos(sig12 / dnm);
					a12 = sig12 / degree;
					omg12 = lam12 / (f1 * dnm);
				}
				else
				{
					// Newton's method. This is a straightforward solution of f(alp1) =
					// lambda12(alp1) - lam12 = 0 with one wrinkle. f(alp) has exactly one
					// root in the interval (0, pi) and its derivative is positive at the
					// root. Thus f(alp) is positive for alp > alp1 and negative for alp <
					// alp1. During the course of the iteration, a range (alp1a, alp1b) is
					// maintained which brackets the root and with each evaluation of
					// f(alp) the range is shrunk, if possible. Newton's method is
					// restarted whenever the derivative of f is negative (because the new
					// value of alp1 is then further from the solution) or if the new
					// estimate of alp1 lies outside (0,pi); in this case, the new starting
					// guess is taken to be (alp1a + alp1b) / 2.
					double ssig1 = 0, csig1 = 0, ssig2 = 0, csig2 = 0, eps = 0;

					uint numit = 0;

					// Bracketing range
					double salp1a = tiny, calp1a = 1, salp1b = tiny, calp1b = -1;
					bool tripn, tripb;
					for (tripn = false, tripb = false; numit < maxit2; ++numit)
					{
						// the WGS84 test set: mean = 1.47, sd = 1.25, max = 16
						// WGS84 and random input: mean = 2.85, sd = 0.60
						double dv = 0;
						double v = Lambda12(sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1,
										slam12, clam12,
										out salp2, out calp2, out sig12, out ssig1, out csig1, out ssig2, out csig2,
										out eps, out somg12, out comg12, numit < maxit1, out dv, Ca) - lam12;

						// 2 * tol0 is approximately 1 ulp for a number in [0, pi].
						// Reversed test to allow escape with NaNs
						if (tripb || !(Math.Abs(v) >= (tripn ? 8 : 1) * tol0)) break;

						// Update bracketing values
						if (v > 0 && (numit > maxit1 || calp1 / salp1 > calp1b / salp1b))
						{ salp1b = salp1; calp1b = calp1; }
						else if (v < 0 && (numit > maxit1 || calp1 / salp1 < calp1a / salp1a))
						{ salp1a = salp1; calp1a = calp1; }

						if (numit < maxit1 && dv > 0)
						{
							double dalp1 = -v / dv;
							double sdalp1 = Math.Sin(dalp1), cdalp1 = Math.Cos(dalp1), nsalp1 = salp1 * cdalp1 + calp1 * sdalp1;
							if (nsalp1 > 0 && Math.Abs(dalp1) < pi)
							{
								calp1 = calp1 * cdalp1 - salp1 * sdalp1;
								salp1 = nsalp1;
								norm2(ref salp1, ref calp1);
								// In some regimes we don't get quadratic convergence because
								// slope -> 0.  So use convergence conditions based on epsilon
								// instead of sqrt(epsilon).
								tripn = Math.Abs(v) <= 16 * tol0;
								continue;
							}
						}

						// Either dv was not postive or updated value was outside legal
						// range.  Use the midpoint of the bracket as the next estimate.
						// This mechanism is not needed for the WGS84 ellipsoid, but it does
						// catch problems with more eccentric ellipsoids.  Its efficacy is
						// such for the WGS84 test set with the starting guess set to alp1 =
						// 90deg:
						// the WGS84 test set: mean = 5.21, sd = 3.93, max = 24
						// WGS84 and random input: mean = 4.74, sd = 0.99
						salp1 = (salp1a + salp1b) / 2;
						calp1 = (calp1a + calp1b) / 2;
						norm2(ref salp1, ref calp1);
						tripn = false;
						tripb = (Math.Abs(salp1a - salp1) + (calp1a - calp1) < tolb || Math.Abs(salp1 - salp1b) + (calp1 - calp1b) < tolb);
					}

					{
						double dummy;
						Lengths(eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
								cbet1, cbet2, out s12x, out m12x, out dummy,
								out M12, out M21, Ca);
					}

					m12x *= b;
					s12x *= b;
					a12 = sig12 / degree;
				}
			}

			if ((outmask & GEOD.DISTANCE) != 0) s12 = 0 + s12x; // Convert -0 to 0

			if ((outmask & GEOD.REDUCEDLENGTH) != 0) m12 = 0 + m12x; // Convert -0 to 0

			if ((outmask & GEOD.AREA) != 0)
			{
				// From Lambda12: sin(alp1) * cos(bet1) = sin(alp0)
				double salp0 = salp1 * cbet1;
				double calp0 = hypotx(calp1, salp1 * sbet1); // calp0 > 0
				double alp12;
				if (calp0 != 0 && salp0 != 0)
				{
					// From Lambda12: tan(bet) = tan(sig) * cos(alp)
					double ssig1 = sbet1;
					double csig1 = calp1 * cbet1;
					double ssig2 = sbet2;
					double csig2 = calp2 * cbet2;
					double k2 = sq(calp0) * ep2;
					double eps = k2 / (2 * (1 + Math.Sqrt(1 + k2)) + k2);

					// Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0).
					double A4 = sq(a) * calp0 * salp0 * e2;

					norm2(ref ssig1, ref csig1);
					norm2(ref ssig2, ref csig2);
					C4f(eps, Ca);

					double B41 = SinCosSeries(false, ssig1, csig1, Ca, nC4);
					double B42 = SinCosSeries(false, ssig2, csig2, Ca, nC4);
					S12 = A4 * (B42 - B41);
				}
				else S12 = 0; // Avoid problems with indeterminate sig1, sig2 on equator

				if (!meridian)
				{
					if (somg12 > 1)
					{
						somg12 = Math.Sin(omg12); comg12 = Math.Cos(omg12);
					}
					else
						norm2(ref somg12, ref comg12);
				}

				if (!meridian &&
					// omg12 < 3/4 * pi
					comg12 > -0.7071 && // Long difference not too big
					sbet2 - sbet1 < 1.75) // Lat difference not too big
				{
					// Use tan(Gamma/2)=tan(omg12/2)*(tan(bet1/2)+tan(bet2/2))/(1+tan(bet1/2)*tan(bet2/2))
					// with tan(x/2)=sin(x)/(1+cos(x))
					double domg12 = 1 + comg12, dbet1 = 1 + cbet1, dbet2 = 1 + cbet2;
					alp12 = 2 * Math.Atan2(somg12 * (sbet1 * dbet2 + sbet2 * dbet1), domg12 * (sbet1 * sbet2 + dbet1 * dbet2));
				}
				else
				{
					// alp12 = alp2 - alp1, used in atan2 so no need to normalize
					double salp12 = salp2 * calp1 - calp2 * salp1;
					double calp12 = calp2 * calp1 + salp2 * salp1;
					// The right thing appears to happen if alp1=+/-180 and alp2=0, viz
					// salp12=-0 and alp12=-180. However this depends on the sign
					// being attached to 0 correctly. The following ensures the correct
					// behavior.
					if (salp12 == 0 && calp12 < 0)
					{
						salp12 = tiny * calp1;
						calp12 = -1;
					}
					alp12 = Math.Atan2(salp12, calp12);
				}
				S12 += c2 * alp12;
				S12 *= swapp * lonsign * latsign;

				// Convert -0 to 0
				S12 += 0;
			}

			// Convert calp, salp to azimuth accounting for lonsign, swapp, latsign.
			if (swapp < 0)
			{
				swapx(ref salp1, ref salp2);
				swapx(ref calp1, ref calp2);
				if ((outmask & GEOD.GEODESICSCALE) != 0)
					swapx(ref M12, ref M21);
			}

			salp1 *= swapp * lonsign; calp1 *= swapp * latsign;
			salp2 *= swapp * lonsign; calp2 *= swapp * latsign;

			psalp1 = salp1;
			pcalp1 = calp1;
			psalp2 = salp2;
			pcalp2 = calp2;

			ps12 = s12;
			pm12 = m12;
			pM12 = M12;
			pM21 = M21;
			pS12 = S12;

			// Returned value in [0, 180]
			return a12;
		}

		/// <summary>
		/// The general inverse geodesic calculation.
		/// </summary>
		/// <param name="lat1">The latitude of point 1 (degrees).</param>
		/// <param name="lon1">The longitude of point 1 (degrees).</param>
		/// <param name="lat2">The latitude of point 2 (degrees).</param>
		/// <param name="lon2">The longitude of point 2 (degrees).</param>
		/// <param name="ps12">The distance between point 1 and point 2 (meters).</param>
		/// <param name="pazi1">The azimuth at point 1 (degrees).</param>
		/// <param name="pazi2">The (forward) azimuth at point 2 (degrees).</param>
		/// <param name="pm12">The reduced length of geodesic (meters).</param>
		/// <param name="pM12">The geodesic scale of point 2 relative to point 1 (dimensionless).</param>
		/// <param name="pM21">The geodesic scale of point 1 relative to point 2 (dimensionless).</param>
		/// <param name="pS12">The area under the geodesic (square meters).</param>
		/// <remarks>
		/// g must have been initialized with a call to geod_init(). lat1
		/// and lat2 should be in the range [-90deg, 90deg].
		/// </remarks>
		/// <returns>a12 arc length of from point 1 to point 2 (degrees).</returns>
		public double geod_geninverse(double lat1, double lon1, double lat2, double lon2,
			out double ps12, out double pazi1, out double pazi2, out double pm12, out double pM12, out double pM21, out double pS12)
		{
			double salp1, calp1, salp2, calp2;
			double a12 = geod_geninverse_int(lat1, lon1, lat2, lon2, out ps12, out salp1, out calp1, out salp2, out calp2, out pm12, out pM12, out pM21, out pS12);
			pazi1 = atan2dx(salp1, calp1);
			pazi2 = atan2dx(salp2, calp2);
			return a12;
		}


		void Lengths(double eps, double sig12, double ssig1, double csig1, double dn1, double ssig2, double csig2, double dn2,
			double cbet1, double cbet2, out double ps12b, out double pm12b, out double pm0, out double pM12, out double pM21,
			// Scratch area of the right size
			double[] Ca)
		{
			double[] Cb = new double[nC];

			// Return m12b = (reduced length)/b; also calculate s12b = distance/b,
			// and m0 = coefficient of secular term in expression for reduced length.
			double A1 = A1m1f(eps);
			C1f(eps, Ca);
			double A2 = A2m1f(eps);
			C2f(eps, Cb);
			double m0 = A1 - A2;
			A2 = 1 + A2;
			A1 = 1 + A1;

			double B1 = SinCosSeries(true, ssig2, csig2, Ca, nC1) - SinCosSeries(true, ssig1, csig1, Ca, nC1);
			// Missing a factor of b
			ps12b = A1 * (sig12 + B1);
			double B2 = SinCosSeries(true, ssig2, csig2, Cb, nC2) - SinCosSeries(true, ssig1, csig1, Cb, nC2);
			double J12 = m0 * sig12 + (A1 * B1 - A2 * B2);

			pm0 = m0;
			// Missing a factor of b.
			// Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure
			// accurate cancellation in the case of coincident points.
			pm12b = dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2) - csig1 * csig2 * J12;

			double csig12 = csig1 * csig2 + ssig1 * ssig2;
			double t = ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2);
			pM12 = csig12 + (t * ssig2 - csig2 * J12) * ssig1 / dn1;
			pM21 = csig12 - (t * ssig1 - csig1 * J12) * ssig2 / dn2;
		}

		static double Astroid(double x, double y)
		{
			// Solve k^4+2*k^3-(x^2+y^2-1)*k^2-2*y^2*k-y^2 = 0 for positive root k.
			// This solution is adapted from Geocentric::Reverse.
			double k;
			double p = sq(x), q = sq(y), r = (p + q - 1) / 6;

			if (!(q == 0 && r <= 0))
			{
				// Avoid possible division by zero when r = 0 by multiplying equations
				// for s and t by r^3 and r, resp.
				double S = p * q / 4; // S = r^3 * s
				double r2 = sq(r);
				double r3 = r * r2;
				// The discriminant of the quadratic equation for T3. This is zero on
				// the evolute curve p^(1/3)+q^(1/3) = 1
				double disc = S * (S + 2 * r3);

				double u = r;
				double v, uv, w;
				if (disc >= 0)
				{
					double T3 = S + r3, T;
					// Pick the sign on the sqrt to maximize abs(T3). This minimizes loss
					// of precision due to cancellation.  The result is unchanged because
					// of the way the T is used in definition of u.
					T3 += T3 < 0 ? -Math.Sqrt(disc) : Math.Sqrt(disc); // T3 = (r * t)^3

					// N.B. cbrtx always returns the real root. cbrtx(-8) = -2.
					T = cbrtx(T3); // T = r * t

					// T can be zero; but then r2 / T -> 0.
					u += T + (T != 0 ? r2 / T : 0);
				}
				else
				{
					// T is complex, but the way u is defined the result is real.
					double ang = Math.Atan2(Math.Sqrt(-disc), -(S + r3));

					// There are three possible cube roots. We choose the root which
					// avoids cancellation. Note that disc < 0 implies that r < 0.
					u += 2 * r * Math.Cos(ang / 3);
				}

				v = Math.Sqrt(sq(u) + q); // guaranteed positive

				// Avoid loss of accuracy when u < 0.
				uv = u < 0 ? q / (v - u) : u + v; // u+v, guaranteed positive
				w = (uv - q) / (2 * v); // positive?

				// Rearrange expression for k to avoid loss of accuracy due to
				// subtraction. Division by 0 not possible because uv > 0, w >= 0.
				k = uv / (Math.Sqrt(uv + sq(w)) + w); // guaranteed positive
			}
			else
			{ // q == 0 && r <= 0
			  // y = 0 with |x| <= 1.  Handle this case directly.
			  // for y small, positive root is k = abs(y)/sqrt(1-x^2)
				k = 0;
			}
			return k;
		}

		double InverseStart(double sbet1, double cbet1, double dn1, double sbet2, double cbet2, double dn2, double lam12, double slam12, double clam12,
			out double psalp1, out double pcalp1, out double psalp2, out double pcalp2, out double pdnm,
			// Scratch area of the right size
			double[] Ca)
		{
			double salp1 = 0, calp1 = 0, salp2 = 0, calp2 = 0, dnm = 0;

			// Return a starting point for Newton's method in salp1 and calp1 (function
			// value is -1). If Newton's method doesn't need to be used, return also
			// salp2 and calp2 and function value is sig12.
			double sig12 = -1; // Return value

			// bet12 = bet2 - bet1 in [0, pi); bet12a = bet2 + bet1 in (-pi, 0]
			double sbet12 = sbet2 * cbet1 - cbet2 * sbet1;
			double cbet12 = cbet2 * cbet1 + sbet2 * sbet1;

			double sbet12a = sbet2 * cbet1 + cbet2 * sbet1;

			bool shortline = cbet12 >= 0 && sbet12 < 0.5 && cbet2 * lam12 < 0.5;
			double somg12, comg12, ssig12, csig12;
			if (shortline)
			{
				double sbetm2 = sq(sbet1 + sbet2);
				// sin((bet1+bet2)/2)^2
				// = (sbet1 + sbet2)^2 / ((sbet1 + sbet2)^2 + (cbet1 + cbet2)^2)
				sbetm2 /= sbetm2 + sq(cbet1 + cbet2);
				dnm = Math.Sqrt(1 + ep2 * sbetm2);
				double omg12 = lam12 / f1 * dnm;
				somg12 = Math.Sin(omg12); comg12 = Math.Cos(omg12);
			}
			else
			{
				somg12 = slam12; comg12 = clam12;
			}

			salp1 = cbet2 * somg12;
			calp1 = comg12 >= 0 ? sbet12 + cbet2 * sbet1 * sq(somg12) / (1 + comg12) : sbet12a - cbet2 * sbet1 * sq(somg12) / (1 - comg12);

			ssig12 = hypotx(salp1, calp1);
			csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12;

			if (shortline && ssig12 < etol2)
			{
				// really short lines
				salp2 = cbet1 * somg12;
				calp2 = sbet12 - cbet1 * sbet2 * (comg12 >= 0 ? sq(somg12) / (1 + comg12) : 1 - comg12);
				norm2(ref salp2, ref calp2);
				// Set return value
				sig12 = Math.Atan2(ssig12, csig12);
			}
			else if (Math.Abs(n) > 0.1 || // No astroid calc if too eccentric
				csig12 >= 0 || ssig12 >= 6 * Math.Abs(n) * pi * sq(cbet1))
			{
				// Nothing to do, zeroth order spherical approximation is OK
			}
			else
			{
				// Scale lam12 and bet2 to x, y coordinate system where antipodal point
				// is at origin and singular point is at y = 0, x = -1.
				double y, lamscale, betscale;
				double x;
				double lam12x = Math.Atan2(-slam12, -clam12); // lam12 - pi
				if (f >= 0)
				{ // In fact f == 0 does not get here
				  // x = dlong, y = dlat
					{
						double k2 = sq(sbet1) * ep2;
						double eps = k2 / (2 * (1 + Math.Sqrt(1 + k2)) + k2);
						lamscale = f * cbet1 * A3f(eps) * pi;
					}
					betscale = lamscale * cbet1;

					x = lam12x / lamscale;
					y = sbet12a / betscale;
				}
				else
				{ // f < 0
				  // x = dlat, y = dlong
					double cbet12a = cbet2 * cbet1 - sbet2 * sbet1;
					double bet12a = Math.Atan2(sbet12a, cbet12a);

					double m12b, m0, dummy;
					// In the case of lon12 = 180, this repeats a calculation made in Inverse.
					Lengths(n, pi + bet12a, sbet1, -cbet1, dn1, sbet2, cbet2, dn2, cbet1, cbet2,
						out dummy, out m12b, out m0, out dummy, out dummy, Ca);

					x = -1 + m12b / (cbet1 * cbet2 * m0 * pi);
					betscale = x < -0.01 ? sbet12a / x : -f * sq(cbet1) * pi;
					lamscale = betscale / cbet1;
					y = lam12x / lamscale;
				}

				if (y > -tol1 && x > -1 - xthresh)
				{
					// strip near cut
					if (f >= 0)
					{
						salp1 = minx(1, -x); calp1 = -Math.Sqrt(1 - sq(salp1));
					}
					else
					{
						calp1 = maxx(x > -tol1 ? 0 : -1, x);
						salp1 = Math.Sqrt(1 - sq(calp1));
					}
				}
				else
				{
					// Estimate alp1, by solving the astroid problem.
					//
					// Could estimate alpha1 = theta + pi/2, directly, i.e.,
					//   calp1 = y/k; salp1 = -x/(1+k);  for f >= 0
					//   calp1 = x/(1+k); salp1 = -y/k;  for f < 0 (need to check)
					//
					// However, it's better to estimate omg12 from astroid and use
					// spherical formula to compute alp1.  This reduces the mean number of
					// Newton iterations for astroid cases from 2.24 (min 0, max 6) to 2.12
					// (min 0 max 5).  The changes in the number of iterations are as
					// follows:
					//
					// change percent
					//    1       5
					//    0      78
					//   -1      16
					//   -2       0.6
					//   -3       0.04
					//   -4       0.002
					//
					// The histogram of iterations is (m = number of iterations estimating
					// alp1 directly, n = number of iterations estimating via omg12, total
					// number of trials = 148605):
					//
					//  iter    m      n
					//    0   148    186
					//    1 13046  13845
					//    2 93315 102225
					//    3 36189  32341
					//    4  5396      7
					//    5   455      1
					//    6    56      0
					//
					// Because omg12 is near pi, estimate work with omg12a = pi - omg12
					double k = Astroid(x, y);
					double omg12a = lamscale * (f >= 0 ? -x * k / (1 + k) : -y * (1 + k) / k);
					somg12 = Math.Sin(omg12a); comg12 = -Math.Cos(omg12a);
					// Update spherical estimate of alp1 using omg12 instead of lam12
					salp1 = cbet2 * somg12;
					calp1 = sbet12a - cbet2 * sbet1 * sq(somg12) / (1 - comg12);
				}
			}

			// Sanity check on starting guess. Backwards check allows NaN through.
			if (!(salp1 <= 0)) norm2(ref salp1, ref calp1);
			else
			{
				salp1 = 1; calp1 = 0;
			}

			psalp1 = salp1;
			pcalp1 = calp1;
			pdnm = dnm;
			psalp2 = salp2;
			pcalp2 = calp2;

			return sig12;
		}

		double Lambda12(double sbet1, double cbet1, double dn1, double sbet2, double cbet2, double dn2, double salp1, double calp1, double slam120, double clam120,
			out double psalp2, out double pcalp2, out double psig12, out double pssig1, out double pcsig1, out double pssig2, out double pcsig2, out double peps,
			out double psomg12, out double pcomg12,
			bool diffp, out double pdlam12,
			// Scratch area of the right size
			double[] Ca)
		{
			double salp2 = 0, calp2 = 0, sig12 = 0,
			ssig1 = 0, csig1 = 0, ssig2 = 0, csig2 = 0, eps = 0,
			somg12 = 0, comg12 = 0, dlam12 = 0;
			double salp0, calp0;
			double somg1, comg1, somg2, comg2, lam12;
			double B312, eta, k2;

			if (sbet1 == 0 && calp1 == 0) calp1 = -tiny; // Break degeneracy of equatorial line. This case has already been handled.

			// sin(alp1) * cos(bet1) = sin(alp0)
			salp0 = salp1 * cbet1;
			calp0 = hypotx(calp1, salp1 * sbet1); // calp0 > 0

			// tan(bet1) = tan(sig1) * cos(alp1)
			// tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1)
			ssig1 = sbet1; somg1 = salp0 * sbet1;
			csig1 = comg1 = calp1 * cbet1;
			norm2(ref ssig1, ref csig1);
			// norm2(ref somg1, ref comg1); -- don't need to normalize!

			// Enforce symmetries in the case abs(bet2) = -bet1.  Need to be careful
			// about this case, since this can yield singularities in the Newton
			// iteration.
			// sin(alp2) * cos(bet2) = sin(alp0)
			salp2 = cbet2 != cbet1 ? salp0 / cbet2 : salp1;

			// calp2 = sqrt(1 - sq(salp2))
			//       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
			// and subst for calp0 and rearrange to give (choose positive sqrt
			// to give alp2 in [0, pi/2]).
			calp2 = cbet2 != cbet1 || Math.Abs(sbet2) != -sbet1 ? Math.Sqrt(sq(calp1 * cbet1) + (cbet1 < -sbet1 ? (cbet2 - cbet1) * (cbet1 + cbet2) : (sbet1 - sbet2) * (sbet1 + sbet2))) / cbet2 : Math.Abs(calp1);

			// tan(bet2) = tan(sig2) * cos(alp2)
			// tan(omg2) = sin(alp0) * tan(sig2).
			ssig2 = sbet2; somg2 = salp0 * sbet2;
			csig2 = comg2 = calp2 * cbet2;
			norm2(ref ssig2, ref csig2);
			// norm2(ref somg2, ref comg2); -- don't need to normalize!

			// sig12 = sig2 - sig1, limit to [0, pi]
			sig12 = Math.Atan2(maxx(0, csig1 * ssig2 - ssig1 * csig2), csig1 * csig2 + ssig1 * ssig2);

			// omg12 = omg2 - omg1, limit to [0, pi]
			somg12 = maxx(0, comg1 * somg2 - somg1 * comg2);
			comg12 = comg1 * comg2 + somg1 * somg2;
			// eta = omg12 - lam120
			eta = Math.Atan2(somg12 * clam120 - comg12 * slam120, comg12 * clam120 + somg12 * slam120);
			k2 = sq(calp0) * ep2;
			eps = k2 / (2 * (1 + Math.Sqrt(1 + k2)) + k2);
			C3f(eps, Ca);
			B312 = (SinCosSeries(true, ssig2, csig2, Ca, nC3 - 1) - SinCosSeries(true, ssig1, csig1, Ca, nC3 - 1));
			lam12 = eta - f * A3f(eps) * salp0 * (sig12 + B312);

			if (diffp)
			{
				if (calp2 == 0) dlam12 = -2 * f1 * dn1 / sbet1;
				else
				{
					double dummy;
					Lengths(eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2, out dummy, out dlam12, out dummy, out dummy, out dummy, Ca);
					dlam12 *= f1 / (calp2 * cbet2);
				}
			}

			psalp2 = salp2;
			pcalp2 = calp2;
			psig12 = sig12;
			pssig1 = ssig1;
			pcsig1 = csig1;
			pssig2 = ssig2;
			pcsig2 = csig2;
			peps = eps;
			psomg12 = somg12;
			pcomg12 = comg12;
			pdlam12 = dlam12;

			return lam12;
		}

		internal double A3f(double eps)
		{
			// Evaluate A3
			return polyval(nA3 - 1, A3x, eps);
		}

		internal void C3f(double eps, double[] c)
		{
			// Evaluate C3 coeffs
			// Elements c[1] thru c[nC3-1] are set
			double mult = 1;
			int o = 0;
			for (int l = 1; l < nC3; ++l)
			{ // l is index of C3[l]
				int m = nC3 - l - 1; // order of polynomial in eps
				mult *= eps;
				c[l] = mult * polyval(m, C3x, o, eps);
				o += m + 1;
			}
		}

		internal void C4f(double eps, double[] c)
		{
			// Evaluate C4 coeffs
			// Elements c[0] thru c[nC4-1] are set
			double mult = 1;
			int o = 0;
			for (int l = 0; l < nC4; ++l)
			{ // l is index of C4[l]
				int m = nC4 - l - 1; // order of polynomial in eps
				c[l] = mult * polyval(m, C4x, o, eps);
				o += m + 1;
				mult *= eps;
			}
		}

		#region Coeffs
		static readonly double[] coeffA1m1f =
			{
				// (1-eps)*A1-1, polynomial in eps2 of order 3
				1, 4, 64, 0, 256,
			};

		// The scale factor A1-1 = mean value of (d/dsigma)I1-1
		internal static double A1m1f(double eps)
		{
			int m = nA1 / 2;
			double t = polyval(m, coeffA1m1f, sq(eps)) / coeffA1m1f[m + 1];
			return (t + eps) / (1 - eps);
		}

		static readonly double[] coeffC1f =
			{
				// C1[1]/eps^1, polynomial in eps2 of order 2
				-1, 6, -16, 32,
				// C1[2]/eps^2, polynomial in eps2 of order 2
				-9, 64, -128, 2048,
				// C1[3]/eps^3, polynomial in eps2 of order 1
				9, -16, 768,
				// C1[4]/eps^4, polynomial in eps2 of order 1
				3, -5, 512,
				// C1[5]/eps^5, polynomial in eps2 of order 0
				-7, 1280,
				// C1[6]/eps^6, polynomial in eps2 of order 0
				-7, 2048,
			};

		// The coefficients C1[l] in the Fourier expansion of B1
		internal static void C1f(double eps, double[] c)
		{
			double eps2 = sq(eps), d = eps;
			int o = 0;
			for (int l = 1; l <= nC1; ++l)
			{ // l is index of C1p[l]
				int m = (nC1 - l) / 2; // order of polynomial in eps^2
				c[l] = d * polyval(m, coeffC1f, o, eps2) / coeffC1f[o + m + 1];
				o += m + 2;
				d *= eps;
			}
		}

		static readonly double[] coeffC1pf =
			{
				// C1p[1]/eps^1, polynomial in eps2 of order 2
				205, -432, 768, 1536,
				// C1p[2]/eps^2, polynomial in eps2 of order 2
				4005, -4736, 3840, 12288,
				// C1p[3]/eps^3, polynomial in eps2 of order 1
				-225, 116, 384,
				// C1p[4]/eps^4, polynomial in eps2 of order 1
				-7173, 2695, 7680,
				// C1p[5]/eps^5, polynomial in eps2 of order 0
				3467, 7680,
				// C1p[6]/eps^6, polynomial in eps2 of order 0
				38081, 61440,
			};

		// The coefficients C1p[l] in the Fourier expansion of B1p
		internal static void C1pf(double eps, double[] c)
		{
			double eps2 = sq(eps), d = eps;
			int o = 0;
			for (int l = 1; l <= nC1p; ++l)
			{ // l is index of C1p[l]
				int m = (nC1p - l) / 2; // order of polynomial in eps^2
				c[l] = d * polyval(m, coeffC1pf, o, eps2) / coeffC1pf[o + m + 1];
				o += m + 2;
				d *= eps;
			}
		}

		static readonly double[] coeffA2m1f =
			{
				// (eps+1)*A2-1, polynomial in eps2 of order 3
				-11, -28, -192, 0, 256,
			};

		// The scale factor A2-1=mean value of (d/dsigma)I2-1
		internal static double A2m1f(double eps)
		{
			int m = nA2 / 2;
			double t = polyval(m, coeffA2m1f, sq(eps)) / coeffA2m1f[m + 1];
			return (t - eps) / (1 + eps);
		}

		static readonly double[] coeffC2f =
			{
				// C2[1]/eps^1, polynomial in eps2 of order 2
				1, 2, 16, 32,
				// C2[2]/eps^2, polynomial in eps2 of order 2
				35, 64, 384, 2048,
				// C2[3]/eps^3, polynomial in eps2 of order 1
				15, 80, 768,
				// C2[4]/eps^4, polynomial in eps2 of order 1
				7, 35, 512,
				// C2[5]/eps^5, polynomial in eps2 of order 0
				63, 1280,
				// C2[6]/eps^6, polynomial in eps2 of order 0
				77, 2048,
			};

		// The coefficients C2[l] in the Fourier expansion of B2
		internal static void C2f(double eps, double[] c)
		{
			double eps2 = sq(eps), d = eps;
			int o = 0;
			for (int l = 1; l <= nC2; ++l)
			{ // l is index of C2[l]
				int m = (nC2 - l) / 2; // order of polynomial in eps^2
				c[l] = d * polyval(m, coeffC2f, o, eps2) / coeffC2f[o + m + 1];
				o += m + 2;
				d *= eps;
			}
		}

		static readonly double[] coeffA3 =
			{
				// A3, coeff of eps^5, polynomial in n of order 0
				-3, 128,
				// A3, coeff of eps^4, polynomial in n of order 1
				-2, -3, 64,
				// A3, coeff of eps^3, polynomial in n of order 2
				-1, -3, -1, 16,
				// A3, coeff of eps^2, polynomial in n of order 2
				3, -1, -2, 8,
				// A3, coeff of eps^1, polynomial in n of order 1
				1, -1, 2,
				// A3, coeff of eps^0, polynomial in n of order 0
				1, 1,
			};

		// The scale factor A3 = mean value of (d/dsigma)I3
		void A3coeff()
		{
			int o = 0, k = 0, j;
			for (j = nA3 - 1; j >= 0; --j)
			{ // coeff of eps^j
				int m = nA3 - j - 1 < j ? nA3 - j - 1 : j; // order of polynomial in n
				A3x[k++] = polyval(m, coeffA3, o, n) / coeffA3[o + m + 1];
				o += m + 2;
			}
		}

		static readonly double[] coeffC3 =
			{
				// C3[1], coeff of eps^5, polynomial in n of order 0
				3, 128,
				// C3[1], coeff of eps^4, polynomial in n of order 1
				2, 5, 128,
				// C3[1], coeff of eps^3, polynomial in n of order 2
				-1, 3, 3, 64,
				// C3[1], coeff of eps^2, polynomial in n of order 2
				-1, 0, 1, 8,
				// C3[1], coeff of eps^1, polynomial in n of order 1
				-1, 1, 4,
				// C3[2], coeff of eps^5, polynomial in n of order 0
				5, 256,
				// C3[2], coeff of eps^4, polynomial in n of order 1
				1, 3, 128,
				// C3[2], coeff of eps^3, polynomial in n of order 2
				-3, -2, 3, 64,
				// C3[2], coeff of eps^2, polynomial in n of order 2
				1, -3, 2, 32,
				// C3[3], coeff of eps^5, polynomial in n of order 0
				7, 512,
				// C3[3], coeff of eps^4, polynomial in n of order 1
				-10, 9, 384,
				// C3[3], coeff of eps^3, polynomial in n of order 2
				5, -9, 5, 192,
				// C3[4], coeff of eps^5, polynomial in n of order 0
				7, 512,
				// C3[4], coeff of eps^4, polynomial in n of order 1
				-14, 7, 512,
				// C3[5], coeff of eps^5, polynomial in n of order 0
				21, 2560,
			};

		// The coefficients C3[l] in the Fourier expansion of B3
		void C3coeff()
		{
			int o = 0, k = 0, l, j;
			for (l = 1; l < nC3; ++l)
			{ // l is index of C3[l]
				for (j = nC3 - 1; j >= l; --j)
				{ // coeff of eps^j
					int m = nC3 - j - 1 < j ? nC3 - j - 1 : j; // order of polynomial in n
					C3x[k++] = polyval(m, coeffC3, o, n) / coeffC3[o + m + 1];
					o += m + 2;
				}
			}
		}

		static readonly double[] coeffC4 =
			{
				// C4[0], coeff of eps^5, polynomial in n of order 0
				97, 15015,
				// C4[0], coeff of eps^4, polynomial in n of order 1
				1088, 156, 45045,
				// C4[0], coeff of eps^3, polynomial in n of order 2
				-224, -4784, 1573, 45045,
				// C4[0], coeff of eps^2, polynomial in n of order 3
				-10656, 14144, -4576, -858, 45045,
				// C4[0], coeff of eps^1, polynomial in n of order 4
				64, 624, -4576, 6864, -3003, 15015,
				// C4[0], coeff of eps^0, polynomial in n of order 5
				100, 208, 572, 3432, -12012, 30030, 45045,
				// C4[1], coeff of eps^5, polynomial in n of order 0
				1, 9009,
				// C4[1], coeff of eps^4, polynomial in n of order 1
				-2944, 468, 135135,
				// C4[1], coeff of eps^3, polynomial in n of order 2
				5792, 1040, -1287, 135135,
				// C4[1], coeff of eps^2, polynomial in n of order 3
				5952, -11648, 9152, -2574, 135135,
				// C4[1], coeff of eps^1, polynomial in n of order 4
				-64, -624, 4576, -6864, 3003, 135135,
				// C4[2], coeff of eps^5, polynomial in n of order 0
				8, 10725,
				// C4[2], coeff of eps^4, polynomial in n of order 1
				1856, -936, 225225,
				// C4[2], coeff of eps^3, polynomial in n of order 2
				-8448, 4992, -1144, 225225,
				// C4[2], coeff of eps^2, polynomial in n of order 3
				-1440, 4160, -4576, 1716, 225225,
				// C4[3], coeff of eps^5, polynomial in n of order 0
				-136, 63063,
				// C4[3], coeff of eps^4, polynomial in n of order 1
				1024, -208, 105105,
				// C4[3], coeff of eps^3, polynomial in n of order 2
				3584, -3328, 1144, 315315,
				// C4[4], coeff of eps^5, polynomial in n of order 0
				-128, 135135,
				// C4[4], coeff of eps^4, polynomial in n of order 1
				-2560, 832, 405405,
				// C4[5], coeff of eps^5, polynomial in n of order 0
				128, 99099,
			};

		// The coefficients C4[l] in the Fourier expansion of I4
		void C4coeff()
		{
			int o = 0, k = 0, l, j;
			for (l = 0; l < nC4; ++l)
			{ // l is index of C4[l]
				for (j = nC4 - 1; j >= l; --j)
				{ // coeff of eps^j
					int m = nC4 - j - 1; // order of polynomial in n
					C4x[k++] = polyval(m, coeffC4, o, n) / coeffC4[o + m + 1];
					o += m + 2;
				}
			}
		}
		#endregion
	}
}
