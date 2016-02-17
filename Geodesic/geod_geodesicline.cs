using System;

namespace Free.Ports.Proj4.Geodesic
{
	/// <summary>
	/// The struct containing information about a single geodesic. This must be
	/// initialized by geod_lineinit() before use.
	/// </summary>
	public class geod_geodesicline : geodesic
	{
		/// <summary>
		/// The starting latitude.
		/// </summary>
		public double lat1;

		/// <summary>
		/// The starting longitude.
		/// </summary>
		public double lon1;

		/// <summary>
		/// The starting azimuth.
		/// </summary>
		public double azi1;

		/// <summary>
		/// The equatorial radius.
		/// </summary>
		public double a;

		/// <summary>
		/// The flattening.
		/// </summary>
		public double f;

		/// <summary>
		/// sine of azi1.
		/// </summary>
		public double salp1;
		/// <summary>
		/// cosine of azi1.
		/// </summary>
		public double calp1;

		/// <summary>
		/// Arc length to reference point.
		/// </summary>
		public double a13;

		/// <summary>
		/// Distance to reference point.
		/// </summary>
		public double s13;

		public double b, c2, f1, salp0, calp0, k2, ssig1, csig1, dn1, stau1, ctau1, somg1, comg1, A1m1, A2m1, A3c, B11, B21, B31, A4, B41;

		public double[] C1a = new double[6 + 1];
		public double[] C1pa = new double[6 + 1];
		public double[] C2a = new double[6 + 1];
		public double[] C3a = new double[6];
		public double[] C4a = new double[6];

		/// <summary>
		/// The capabilities.
		/// </summary>
		public GEOD caps;

		void geod_lineinit_int(geod_geodesic g, double lat1, double lon1, double azi1, double salp1, double calp1, GEOD caps)
		{
			a = g.a;
			f = g.f;
			b = g.b;
			c2 = g.c2;
			f1 = g.f1;

			// If caps is 0 assume the standard direct calculation
			this.caps = (caps != 0 ? caps : GEOD.DISTANCE_IN | GEOD.LONGITUDE) |
				GEOD.LATITUDE | GEOD.AZIMUTH | GEOD.LONG_UNROLL; // always allow latitude and azimuth and unrolling of longitude

			this.lat1 = LatFix(lat1);
			this.lon1 = lon1;
			this.azi1 = azi1;
			this.salp1 = salp1;
			this.calp1 = calp1;

			double cbet1, sbet1;

			sincosdx(AngRound(lat1), out sbet1, out cbet1); sbet1 *= f1;
			// Ensure cbet1 = +epsilon at poles
			norm2(ref sbet1, ref cbet1); cbet1 = maxx(tiny, cbet1);
			dn1 = Math.Sqrt(1 + g.ep2 * sq(sbet1));

			// Evaluate alp0 from sin(alp1) * cos(bet1) = sin(alp0),
			salp0 = salp1 * cbet1;  // alp0 in [0, pi/2 - |bet1|]
									// Alt: calp0 = hypot(sbet1, calp1 * cbet1).  The following
									// is slightly better (consider the case salp1 = 0).
			calp0 = hypotx(calp1, salp1 * sbet1);
			// Evaluate sig with tan(bet1) = tan(sig1) * cos(alp1).
			// sig = 0 is nearest northward crossing of equator.
			// With bet1 = 0, alp1 = pi/2, we have sig1 = 0 (equatorial line).
			// With bet1 =  pi/2, alp1 = -pi, sig1 =  pi/2
			// With bet1 = -pi/2, alp1 =  0 , sig1 = -pi/2
			// Evaluate omg1 with tan(omg1) = sin(alp0) * tan(sig1).
			// With alp0 in (0, pi/2], quadrants for sig and omg coincide.
			// No atan2(0,0) ambiguity at poles since cbet1 = +epsilon.
			// With alp0 = 0, omg1 = 0 for alp1 = 0, omg1 = pi for alp1 = pi.
			ssig1 = sbet1; somg1 = salp0 * sbet1;
			csig1 = comg1 = sbet1 != 0 || calp1 != 0 ? cbet1 * calp1 : 1;
			norm2(ref ssig1, ref csig1); // sig1 in (-pi, pi]
										 // norm2(somg1, comg1); -- don't need to normalize!

			k2 = sq(calp0) * g.ep2;
			double eps = k2 / (2 * (1 + Math.Sqrt(1 + k2)) + k2);

			if ((this.caps & GEOD.CAP_C1) != 0)
			{
				double s, c;
				A1m1 = geod_geodesic.A1m1f(eps);
				geod_geodesic.C1f(eps, C1a);
				B11 = SinCosSeries(true, ssig1, csig1, C1a, nC1);
				s = Math.Sin(B11); c = Math.Cos(B11);

				// tau1=sig1+B11
				stau1 = ssig1 * c + csig1 * s;
				ctau1 = csig1 * c - ssig1 * s;
				// Not necessary because C1pa reverts C1a
				// B11=-SinCosSeries(TRUE, stau1, ctau1, C1pa, nC1p);
			}

			if ((this.caps & GEOD.CAP_C1p) != 0)
				geod_geodesic.C1pf(eps, C1pa);

			if ((this.caps & GEOD.CAP_C2) != 0)
			{
				A2m1 = geod_geodesic.A2m1f(eps);
				geod_geodesic.C2f(eps, C2a);
				B21 = SinCosSeries(true, ssig1, csig1, C2a, nC2);
			}

			if ((this.caps & GEOD.CAP_C3) != 0)
			{
				g.C3f(eps, C3a);
				A3c = -f * salp0 * g.A3f(eps);
				B31 = SinCosSeries(true, ssig1, csig1, C3a, nC3 - 1);
			}

			if ((this.caps & GEOD.CAP_C4) != 0)
			{
				g.C4f(eps, C4a);
				// Multiplier=a^2*e^2*cos(alpha0)*sin(alpha0)
				A4 = sq(a) * calp0 * salp0 * g.e2;
				B41 = SinCosSeries(false, ssig1, csig1, C4a, nC4);
			}

			a13 = s13 = double.NaN;
		}

		/// <summary>
		/// geod_lineinit : Initialize a geod_geodesicline object.
		/// </summary>
		/// <param name="g">The <see cref="geod_geodesic"/> object specifying the ellipsoid.</param>
		/// <param name="lat1">Latitude of point 1 (degrees).</param>
		/// <param name="lon1">Longitude of point 1 (degrees).</param>
		/// <param name="azi1">Azimuth at point 1 (degrees).</param>
		/// <param name="caps">Bitor'ed combination of <see cref="geod_mask"/> values. See remarks for more informations.</param>
		/// <remarks>
		/// <paramref name="caps"/> bitor'ed combination of geod_mask() values specifying
		/// the capabilities the geod_geodesicline object should possess, i.e., which
		/// quantities can be returned in calls to geod_position() and
		/// geod_genposition().
		///
		/// <paramref name="lat1"/> should be in the range [-90deg, 90deg].
		///
		/// The geod_mask values are:
		/// - caps|=GEOD_LATITUDE for the latitude lat2; this is added automatically,
		/// - caps|=GEOD_LONGITUDE for the latitude lon2,
		/// - caps|=GEOD_AZIMUTH for the latitude azi2; this is added automatically,
		/// - caps|=GEOD_DISTANCE for the distance s12,
		/// - caps|=GEOD_REDUCEDLENGTH for the reduced length m12,
		/// - caps|=GEOD_GEODESICSCALE for the geodesic scales M12 and M21,
		/// - caps|=GEOD_AREA for the area S12,
		/// - caps|=GEOD_DISTANCE_IN permits the length of the
		///   geodesic to be given in terms of s12; without this capability the
		///   length can only be specified in terms of arc length.
		///
		/// A value of caps=0 is treated as GEOD.LATITUDE|GEOD.LONGITUDE|GEOD.AZIMUTH|GEOD.DISTANCE_IN
		/// (to support the solution of the "standard" direct problem).
		///
		/// When initialized by this function, point 3 is undefined (s13 = a13 = NaN).
		/// </remarks>
		public geod_geodesicline(geod_geodesic g, double lat1, double lon1, double azi1, GEOD caps)
		{
			azi1 = AngNormalize(azi1);
			double salp1, calp1;

			// Guard against underflow in salp0
			sincosdx(AngRound(azi1), out salp1, out calp1);
			geod_lineinit_int(g, lat1, lon1, azi1, salp1, calp1, caps);
		}

		/// <summary>
		/// geod_directline : Initialize a geod_geodesicline object in terms of the direct geodesic problem.
		/// </summary>
		/// <param name="g">The <see cref="geod_geodesic"/> object specifying the ellipsoid.</param>
		/// <param name="lat1">Latitude of point 1 (degrees).</param>
		/// <param name="lon1">Longitude of point 1 (degrees).</param>
		/// <param name="azi1">Azimuth at point 1 (degrees).</param>
		/// <param name="s12">Distance from point 1 to point 2 (meters); it can be negative.</param>
		/// <param name="caps">Bitor'ed combination of <see cref="geod_mask"/> values. See remarks for more informations.</param>
		/// <remarks>
		/// <paramref name="caps"/> bitor'ed combination of geod_mask() values specifying
		/// the capabilities the geod_geodesicline object should possess, i.e., which
		/// quantities can be returned in calls to geod_position() and
		/// geod_genposition().
		///
		/// This function sets point 3 of the geod_geodesicline to correspond to point
		/// 2 of the direct geodesic problem.  See geod_lineinit() for more
		/// information.
		/// </remarks>
		public geod_geodesicline(geod_geodesic g, double lat1, double lon1, double azi1, double s12, GEOD caps) : this(g, lat1, lon1, azi1, GEOD.NOFLAGS, s12, caps)
		{ }

		/// <summary>
		/// geod_gendirectline : Initialize a geod_geodesicline object in terms of the direct geodesic
		/// problem specified in terms of either distance or arc length.
		/// </summary>
		/// <param name="g">The <see cref="geod_geodesic"/> object specifying the ellipsoid.</param>
		/// <param name="lat1">Latitude of point 1 (degrees).</param>
		/// <param name="lon1">Longitude of point 1 (degrees).</param>
		/// <param name="azi1">Azimuth at point 1 (degrees).</param>
		/// <param name="flags">Must be either <see cref="GEOD.NOFLAGS"/> or <see cref="GEOD.ARCMODE"/> to determining the meaning of the <paramref name="s12_a12"/>.</param>
		/// <param name="s12_a12">if <paramref name="flags"/> = <see cref="GEOD.NOFLAGS"/>, this is the distance from point 1 to point 2 (meters);
		/// if <paramref name="flags"/> = <see cref="GEOD.ARCMODE"/>, it is the arc length from point 1 to point 2 (degrees); it can be negative.</param>
		/// <param name="caps">Bitor'ed combination of <see cref="geod_mask"/> values. See remarks for more informations.</param>
		/// <remarks>
		/// <paramref name="caps"/> bitor'ed combination of geod_mask() values specifying
		/// the capabilities the geod_geodesicline object should possess, i.e., which
		/// quantities can be returned in calls to geod_position() and
		/// geod_genposition().
		///
		/// This function sets point 3 of the geod_geodesicline to correspond to point
		/// 2 of the direct geodesic problem.  See geod_lineinit() for more
		/// information.
		/// </remarks>
		public geod_geodesicline(geod_geodesic g, double lat1, double lon1, double azi1, GEOD flags, double s12_a12, GEOD caps) : this(g, lat1, lon1, azi1, caps)
		{
			geod_gensetdistance(flags, s12_a12);
		}

		/// <summary>
		/// The general position function.
		/// </summary>
		/// <param name="flags">Bitor'ed combination of <see cref="GEOD"/> flags; GEOD.ARCMODE
		/// determines the meaning of s12_a12 and GEOD.LONG_UNROLL "unrolls" lon2; if flags&amp;GEOD.ARCMODE is 0,
		/// then l must have been initialized with caps|=GEOD.DISTANCE_IN.</param>
		/// <param name="s12_a12">If flags&amp;GEOD_ARCMODE is 0, this is the
		/// distance from point 1 to point 2 (meters); otherwise it is the
		/// arc length from point 1 to point 2 (degrees); it can be negative.</param>
		/// <param name="plat2">The latitude of point 2 (degrees).</param>
		/// <param name="plon2">The longitude of point 2 (degrees); requires
		/// that l was initialized with caps|=GEOD.LONGITUDE.</param>
		/// <param name="pazi2">The (forward) azimuth at point 2 (degrees).</param>
		/// <param name="ps12">The distance from point 1 to point 2 (meters); requires that
		/// l was initialized with caps|=GEOD.DISTANCE.</param>
		/// <param name="pm12">The reduced length of geodesic (meters);
		/// requires that l was initialized with caps|=GEOD.REDUCEDLENGTH.</param>
		/// <param name="pM12">The geodesic scale of point 2 relative to point 1 (dimensionless);
		/// requires that l was initialized with caps|=GEOD.GEODESICSCALE.</param>
		/// <param name="pM21">The geodesic scale of point 1 relative to point 2 (dimensionless);
		/// requires that l was initialized with caps|=GEOD.GEODESICSCALE.</param>
		/// <param name="pS12">The area under the geodesic (square meters); requires that l
		/// was initialized with caps|=GEOD.AREA.</param>
		/// <remarks>
		/// l must have been initialized with a call to geod_lineinit() with caps|=GEOD.DISTANCE_IN.
		/// The value azi2 returned is in the range [-180deg, 180deg).
		///
		/// Requesting a value which l is not capable of computing
		/// is not an error; the corresponding argument will not be altered.
		///
		/// With <paramref name="flags"/> &amp; GEOD.LONG_UNROLL bit set, the longitude is "unrolled" so
		/// that the quantity lon2-lon1 indicates how many times and in
		/// what sense the geodesic encircles the ellipsoid.
		///</remarks>
		/// <returns>a12 arc length from point 1 to point 2 (degrees).</returns>
		public double geod_genposition(GEOD flags, double s12_a12,
			out double plat2, out double plon2, out double pazi2, out double ps12, out double pm12, out double pM12, out double pM21, out double pS12)
		{
			double lat2 = 0, lon2 = 0, azi2 = 0, s12 = 0, m12 = 0, M12 = 0, M21 = 0, S12 = 0;

			// Avoid warning about uninitialized B12.
			double sig12, ssig12, csig12, B12 = 0, AB1 = 0;
			double omg12, lam12, lon12;
			double ssig2, csig2, sbet2, cbet2, somg2, comg2, salp2, calp2, dn2;

			GEOD outmask = GEOD.LATITUDE | GEOD.LONGITUDE | GEOD.AZIMUTH | GEOD.DISTANCE | GEOD.REDUCEDLENGTH | GEOD.GEODESICSCALE | GEOD.AREA;
			outmask &= caps & GEOD.OUT_ALL;

			if (!(((flags & GEOD.ARCMODE) != 0 || (caps & (GEOD.DISTANCE_IN & GEOD.OUT_ALL)) != 0)))
			{
				plat2 = plon2 = pazi2 = ps12 = pm12 = pM12 = pM21 = pS12 = 0;

				// Uninitialized or impossible distance calculation requested
				return NaN;
			}

			if ((flags & GEOD.ARCMODE) != 0)
			{
				// Interpret s12_a12 as spherical arc length
				sig12 = s12_a12 * degree;
				sincosdx(s12_a12, out ssig12, out csig12);
			}
			else
			{
				// Interpret s12_a12 as distance
				double tau12 = s12_a12 / (b * (1 + A1m1));
				double s = Math.Sin(tau12);
				double c = Math.Cos(tau12);

				// tau2=tau1+tau12
				B12 = -SinCosSeries(true, stau1 * c + ctau1 * s, ctau1 * c - stau1 * s, C1pa, nC1p);
				sig12 = tau12 - (B12 - B11);
				ssig12 = Math.Sin(sig12); csig12 = Math.Cos(sig12);

				if (Math.Abs(f) > 0.01)
				{
					// Reverted distance series is inaccurate for |f|>1/100, so correct
					// sig12 with 1 Newton iteration. The following table shows the
					// approximate maximum error for a=WGS_a() and various f relative to
					// GeodesicExact.
					//     erri = the error in the inverse solution (nm)
					//     errd = the error in the direct solution (series only) (nm)
					//     errda = the error in the direct solution (series + 1 Newton) (nm)
					//
					//       f     erri  errd errda
					//     -1/5    12e6 1.2e9  69e6
					//     -1/10  123e3  12e6 765e3
					//     -1/20   1110 108e3  7155
					//     -1/50  18.63 200.9 27.12
					//     -1/100 18.63 23.78 23.37
					//     -1/150 18.63 21.05 20.26
					//      1/150 22.35 24.73 25.83
					//      1/100 22.35 25.03 25.31
					//      1/50  29.80 231.9 30.44
					//      1/20   5376 146e3  10e3
					//      1/10  829e3  22e6 1.5e6
					//      1/5   157e6 3.8e9 280e6
					ssig2 = ssig1 * csig12 + csig1 * ssig12;
					csig2 = csig1 * csig12 - ssig1 * ssig12;
					B12 = SinCosSeries(true, ssig2, csig2, C1a, nC1);
					double serr = (1 + A1m1) * (sig12 + (B12 - B11)) - s12_a12 / b;
					sig12 = sig12 - serr / Math.Sqrt(1 + k2 * sq(ssig2));
					ssig12 = Math.Sin(sig12); csig12 = Math.Cos(sig12);
					// Update B12 below
				}
			}

			// sig2=sig1+sig12
			ssig2 = ssig1 * csig12 + csig1 * ssig12;
			csig2 = csig1 * csig12 - ssig1 * ssig12;
			dn2 = Math.Sqrt(1 + k2 * sq(ssig2));

			if ((outmask & (GEOD.DISTANCE | GEOD.REDUCEDLENGTH | GEOD.GEODESICSCALE)) != 0)
			{
				if ((flags & GEOD.ARCMODE) != 0 || Math.Abs(f) > 0.01)
					B12 = SinCosSeries(true, ssig2, csig2, C1a, nC1);
				AB1 = (1 + A1m1) * (B12 - B11);
			}

			// sin(bet2)=cos(alp0)*sin(sig2)
			sbet2 = calp0 * ssig2;

			// Alt: cbet2=hypot(csig2, salp0*ssig2);
			cbet2 = hypotx(salp0, calp0 * csig2);

			if (cbet2 == 0) cbet2 = csig2 = tiny; // I.e., salp0=0, csig2=0. Break the degeneracy in this case

			// tan(alp0) = cos(sig2)*tan(alp2)
			salp2 = salp0; calp2 = calp0 * csig2; // No need to normalize

			if ((outmask & GEOD.DISTANCE) != 0) s12 = (flags & GEOD.ARCMODE) != 0 ? b * ((1 + A1m1) * sig12 + AB1) : s12_a12;

			if ((outmask & GEOD.LONGITUDE) != 0)
			{
				double E = copysignx(1, salp0); // east or west going?

				// tan(omg2)=sin(alp0)*tan(sig2)
				somg2 = salp0 * ssig2; comg2 = csig2; // No need to normalize

				// omg12=omg2-omg1
				omg12 = (flags & GEOD.LONG_UNROLL) != 0 ?
					E * (sig12
						- (Math.Atan2(ssig2, csig2) - Math.Atan2(ssig1, csig1))
						+ (Math.Atan2(E * somg2, comg2) - Math.Atan2(E * somg1, comg1))) :
					Math.Atan2(somg2 * comg1 - comg2 * somg1, comg2 * comg1 + somg2 * somg1);

				lam12 = omg12 + A3c * (sig12 + (SinCosSeries(true, ssig2, csig2, C3a, nC3 - 1) - B31));
				lon12 = lam12 / degree;

				lon2 = (flags & GEOD.LONG_UNROLL) != 0 ? lon1 + lon12 : AngNormalize(AngNormalize(lon1) + AngNormalize(lon12));
			}

			if ((outmask & GEOD.LATITUDE) != 0) lat2 = atan2dx(sbet2, f1 * cbet2);

			if ((outmask & GEOD.AZIMUTH) != 0) azi2 = atan2dx(salp2, calp2);

			if ((outmask & (GEOD.REDUCEDLENGTH | GEOD.GEODESICSCALE)) != 0)
			{
				double B22 = SinCosSeries(true, ssig2, csig2, C2a, nC2);
				double AB2 = (1 + A2m1) * (B22 - B21);
				double J12 = (A1m1 - A2m1) * sig12 + (AB1 - AB2);

				if ((outmask & GEOD.REDUCEDLENGTH) != 0)
					// Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure
					// accurate cancellation in the case of coincident points.
					m12 = b * ((dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2)) - csig1 * csig2 * J12);

				if ((outmask & GEOD.GEODESICSCALE) != 0)
				{
					double t = k2 * (ssig2 - ssig1) * (ssig2 + ssig1) / (dn1 + dn2);
					M12 = csig12 + (t * ssig2 - csig2 * J12) * ssig1 / dn1;
					M21 = csig12 - (t * ssig1 - csig1 * J12) * ssig2 / dn2;
				}
			}

			if ((outmask & GEOD.AREA) != 0)
			{
				double B42 = SinCosSeries(false, ssig2, csig2, C4a, nC4);
				double salp12, calp12;
				if (calp0 == 0 || salp0 == 0)
				{
					// alp12=alp2-alp1, used in atan2 so no need to normalize
					salp12 = salp2 * calp1 - calp2 * salp1;
					calp12 = calp2 * calp1 + salp2 * salp1;
				}
				else
				{
					// tan(alp)=tan(alp0)*sec(sig)
					// tan(alp2-alp1)=(tan(alp2)-tan(alp1))/(tan(alp2)*tan(alp1)+1)
					// =calp0*salp0*(csig1-csig2)/(salp0^2+calp0^2*csig1*csig2)
					// If csig12>0, write
					//   csig1-csig2=ssig12*(csig1*ssig12/(1+csig12)+ssig1)
					// else
					//   csig1-csig2=csig1*(1-csig12)+ssig12*ssig1
					// No need to normalize
					salp12 = calp0 * salp0 * (csig12 <= 0 ? csig1 * (1 - csig12) + ssig12 * ssig1 : ssig12 * (csig1 * ssig12 / (1 + csig12) + ssig1));
					calp12 = sq(salp0) + sq(calp0) * csig1 * csig2;
				}
				S12 = c2 * Math.Atan2(salp12, calp12) + A4 * (B42 - B41);
			}

			plat2 = lat2;
			plon2 = lon2;
			pazi2 = azi2;
			ps12 = s12;
			pm12 = m12;
			pM12 = M12;
			pM21 = M21;
			pS12 = S12;

			return (flags & GEOD.ARCMODE) != 0 ? s12_a12 : sig12 / degree;
		}

		/// <summary>
		/// Specify position of point 3 in terms of distance.
		/// </summary>
		/// <param name="s13">The distance from point 1 to point 3 (meters); it can be negative.</param>
		/// <remarks>This is only useful if the geod_geodesicline object has been constructed with caps |= <see cref="GEOD.DISTANCE_IN"/>.</remarks>
		public void geod_setdistance(double s13)
		{
			this.s13 = s13;
			double dummy;
			a13 = geod_genposition(GEOD.NOFLAGS, s13, out dummy, out dummy, out dummy, out dummy, out dummy, out dummy, out dummy, out dummy);
		}

		void geod_setarc(double a13)
		{
			this.a13 = a13; s13 = NaN;
			double dummy;
			geod_genposition(GEOD.ARCMODE, a13, out dummy, out dummy, out dummy, out s13, out dummy, out dummy, out dummy, out dummy);
		}

		/// <summary>
		/// Specify position of point 3 in terms of either distance or arc length.
		/// </summary>
		/// <param name="flags">Must be either <see cref="GEOD.NOFLAGS"/> or <see cref="GEOD.ARCMODE"/> to determining the meaning of the <paramref name="s13_a13"/>.</param>
		/// <param name="s13_a13">if <paramref name="flags"/> = <see cref="GEOD.NOFLAGS"/>, this is the distance from point 1 to point 3 (meters);
		/// if <paramref name="flags"/> = <see cref="GEOD.ARCMODE"/>, it is the arc length from point 1 to point 3 (degrees); it can be negative.</param>
		/// <remarks>
		/// If <paramref name="flags"/> = <see cref="GEOD.NOFLAGS"/>, this calls <see cref="geod_setdistance"/>.
		/// If <paramref name="flags"/> = <see cref="GEOD.ARCMODE"/>, the s13 is only set if the <see cref="geod_geodesicline"/>
		/// object has been constructed with caps |= <see cref="GEOD.DISTANCE"/>.
		/// </remarks>
		public void geod_gensetdistance(GEOD flags, double s13_a13)
		{
			if ((flags & GEOD.ARCMODE) == GEOD.ARCMODE) geod_setarc(s13_a13);
			else geod_setdistance(s13_a13);
		}

		/// <summary>
		/// Compute the position along a geod_geodesicline.
		/// </summary>
		/// <param name="s12">Distance from point 1 to point 2 (meters); it can be negative.</param>
		/// <param name="lat2">The latitude of point 2 (degrees).</param>
		/// <param name="lon2">The longitude of point 2 (degrees); requires that <paramref name="l"/> was initialized with caps|=GEOD.LONGITUDE.</param>
		/// <param name="azi2">The (forward) azimuth at point 2 (degrees).</param>
		/// <remarks>
		/// <paramref name="l"/> must have been initialized with a call, e.g. to geod_lineinit() with
		/// caps|=GEOD.DISTANCE_IN. The values of <paramref name="lon2"/> and <paramref name="azi2"/> returned are
		/// in the range [-180deg, 180deg).
		/// </remarks>
		public void geod_position(double s12, out double lat2, out double lon2, out double azi2)
		{
			double dummy;
			geod_genposition(GEOD.NONE, s12, out lat2, out lon2, out azi2, out dummy, out dummy, out dummy, out dummy, out dummy);
		}

		/// <summary>
		/// Compute the position along a geod_geodesicline.
		/// </summary>
		/// <param name="s12">Distance from point 1 to point 2 (meters); it can be negative.</param>
		/// <param name="lat2">The latitude of point 2 (degrees).</param>
		/// <param name="lon2">The longitude of point 2 (degrees); requires that the object was initialized with caps|=GEOD.LONGITUDE.</param>
		/// <remarks>
		/// <paramref name="l"/> must have been initialized with a call, e.g. to geod_lineinit() with
		/// caps|=GEOD.DISTANCE_IN. The values of <paramref name="lon2"/> and <paramref name="azi2"/> returned are
		/// in the range [-180deg, 180deg).
		/// </remarks>
		public void geod_position(double s12, out double lat2, out double lon2)
		{
			double dummy;
			geod_genposition(GEOD.NONE, s12, out lat2, out lon2, out dummy, out dummy, out dummy, out dummy, out dummy, out dummy);
		}

		/// <summary>
		/// Initialize a geod_geodesicline object in terms of the inverse geodesic problem.
		/// </summary>
		/// <param name="g">The <see cref="geod_geodesic"/> object specifying the ellipsoid.</param>
		/// <param name="lat1">Latitude of point 1 (degrees).</param>
		/// <param name="lon1">Longitude of point 1 (degrees).</param>
		/// <param name="lat2">Latitude of point 2 (degrees).</param>
		/// <param name="lon2">Longitude of point 2 (degrees).</param>
		/// <param name="caps">Bitor'ed combination of <see cref="geod_mask"/> values. See remarks for more informations.</param>
		/// <remarks>
		/// <paramref name="caps"/> bitor'ed combination of geod_mask() values specifying
		/// the capabilities the geod_geodesicline object should possess, i.e., which
		/// quantities can be returned in calls to geod_position() and
		/// geod_genposition().
		///
		/// This function sets point 3 of the geod_geodesicline to correspond to point
		/// 2 of the direct geodesic problem.  See geod_lineinit() for more
		/// information.
		/// </remarks>
		public void geod_inverseline(geod_geodesic g, double lat1, double lon1, double lat2, double lon2, GEOD caps)
		{
			double salp1, calp1, dummy;
			double a12 = g.geod_geninverse_int(lat1, lon1, lat2, lon2, out dummy, out salp1, out calp1, out dummy, out dummy, out dummy, out dummy, out dummy, out dummy);
			double azi1 = atan2dx(salp1, calp1);

			caps = caps != 0 ? caps : GEOD.DISTANCE_IN | GEOD.LONGITUDE;
			// Ensure that a12 can be converted to a distance
			if ((caps & (GEOD.OUT_ALL & GEOD.DISTANCE_IN)) != 0) caps |= GEOD.DISTANCE;
			geod_lineinit_int(g, lat1, lon1, azi1, salp1, calp1, caps);
			geod_setarc(a12);
		}
	}
}
