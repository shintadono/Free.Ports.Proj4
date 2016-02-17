namespace Free.Ports.Proj4.Geodesic
{
	/// <summary>
	/// The struct for accumulating information about a geodesic polygon. This is
	/// used for computing the perimeter and area of a polygon. This must be
	/// initialized by geod_polygon_init() before use.
	/// </summary>
	public class geod_polygon : geodesic
	{
		/// <summary>
		/// The current latitude.
		/// </summary>
		public double lat;

		/// <summary>
		/// The current longitude.
		/// </summary>
		public double lon;

		public double lat0;
		public double lon0;
		public double[] A = new double[2];
		public double[] P = new double[2];
		public bool polyline;
		public int crossings;

		/// <summary>
		/// The number of points so far.
		/// </summary>
		public int num;

		#region private static methods
		static int transit(double lon1, double lon2)
		{
			// Return 1 or -1 if crossing prime meridian in east or west direction.
			// Otherwise return zero.
			// Compute lon12 the same way as Geodesic::Inverse.
			lon1 = AngNormalize(lon1);
			lon2 = AngNormalize(lon2);
			double dummy;
			double lon12 = AngDiff(lon1, lon2, out dummy);
			return lon1 < 0 && lon2 >= 0 && lon12 > 0 ? 1 : (lon2 < 0 && lon1 >= 0 && lon12 < 0 ? -1 : 0);
		}

		static int transitdirect(double lon1, double lon2)
		{
			lon1 = lon1 % 720;
			lon2 = lon2 % 720;
			return (((lon2 >= 0 && lon2 < 360) || lon2 < -360 ? 0 : 1) - ((lon1 >= 0 && lon1 < 360) || lon1 < -360 ? 0 : 1));
		}

		static void accini(double[] s)
		{
			// Initialize an accumulator; this is an array with two elements.
			s[0] = s[1] = 0;
		}

		static void acccopy(double[] s, double[] t)
		{
			// Copy an accumulator; t = s.
			t[0] = s[0]; t[1] = s[1];
		}

		static void accadd(double[] s, double y)
		{
			// Add y to an accumulator.
			double u, z = sumx(y, s[1], out u);
			s[0] = sumx(z, s[0], out s[1]);
			if (s[0] == 0) s[0] = u;
			else s[1] = s[1] + u;
		}

		static double accsum(double[] s, double y)
		{
			// Return accumulator + y (but don't add to accumulator).
			double[] t = new double[2];
			acccopy(s, t);
			accadd(t, y);
			return t[0];
		}

		static void accneg(double[] s)
		{
			// Negate an accumulator.
			s[0] = -s[0]; s[1] = -s[1];
		}
		#endregion

		/// <summary>
		/// Initialize a geod_polygon object.
		/// </summary>
		/// <param name="polyline">Set true to specify a polyline instead of a polygon.</param>
		/// <remarks>
		/// If polylinep is zero, then the sequence of vertices and edges added by
		/// geod_polygon_addpoint() and geod_polygon_addedge() define a polygon and
		/// the perimeter and area are returned by geod_polygon_compute(). If
		/// polylinep is non-zero, then the vertices and edges define a polyline and
		/// only the perimeter is returned by geod_polygon_compute().
		///
		/// The area and perimeter are accumulated at two times the standard floating
		/// point precision to guard against the loss of accuracy with many-sided
		/// polygons. At any point you can ask for the perimeter and area so far.
		/// </remarks>
		public geod_polygon(bool polyline)
		{
			this.polyline = polyline;
			geod_polygon_clear();
		}

		/// <summary>
		/// Clear the polygon, allowing a new polygon to be started.
		/// </summary>
		void geod_polygon_clear()
		{
			lat0 = lon0 = lat = lon = NaN;
			accini(P);
			accini(A);
			num = crossings = 0;
		}

		/// <summary>
		/// Add a point to the polygon or polyline.
		/// </summary>
		/// <param name="g">The geod_geodesic object specifying the ellipsoid.</param>
		/// <param name="p">The geod_polygon object specifying the polygon.</param>
		/// <param name="lat">The latitude of the point (degrees).</param>
		/// <param name="lon">The longitude of the point (degrees).</param>
		/// <remarks>
		/// g and p must have been initialized with calls to geod_init() and
		/// geod_polygon_init(), respectively. The same g must be used for all the
		/// points and edges in a polygon. lat should be in the range
		/// [-90deg, 90deg].
		/// </remarks>
		public void geod_polygon_addpoint(geod_geodesic g, double lat, double lon)
		{
			lon = AngNormalize(lon);
			if (num == 0)
			{
				lat0 = this.lat = lat;
				lon0 = this.lon = lon;
			}
			else
			{
				double s12, S12 = 0, dummy;
				if (polyline) g.geod_geninverse(this.lat, this.lon, lat, lon, out s12, out dummy, out dummy, out dummy, out dummy, out dummy, out dummy);
				else g.geod_geninverse(this.lat, this.lon, lat, lon, out s12, out dummy, out dummy, out dummy, out dummy, out dummy, out S12);

				accadd(P, s12);
				if (!polyline)
				{
					accadd(A, S12);
					crossings += transit(this.lon, lon);
				}
				this.lat = lat; this.lon = lon;
			}
			++num;
		}

		/// <summary>
		/// Add an edge to the polygon or polyline.
		/// </summary>
		/// <param name="g">The geod_geodesic object specifying the ellipsoid.</param>
		/// <param name="p">The geod_polygon object specifying the polygon.</param>
		/// <param name="azi">The azimuth at current point (degrees).</param>
		/// <param name="s">The distance from current point to next point (meters).</param>
		/// <remarks>
		/// g and p must have been initialized with calls to geod_init() and
		/// geod_polygon_init(), respectively. The same g must be used for all the
		/// points and edges in a polygon. This does nothing if no points have been
		/// added yet. The lat and lon fields of p give the location of
		/// the new vertex.
		/// </remarks>
		public void geod_polygon_addedge(geod_geodesic g, double azi, double s)
		{
			if (num != 0)
			{ // Do nothing is num is zero
				double lat, lon, S12 = 0, dummy;
				if (polyline) g.geod_gendirect(this.lat, this.lon, azi, GEOD.LONG_UNROLL, s, out lat, out lon, out dummy, out dummy, out dummy, out dummy, out dummy, out dummy);
				else g.geod_gendirect(this.lat, this.lon, azi, GEOD.LONG_UNROLL, s, out lat, out lon, out dummy, out dummy, out dummy, out dummy, out dummy, out S12);

				accadd(P, s);
				if (!polyline)
				{
					accadd(A, S12);
					crossings += transitdirect(this.lon, lon);
				}
				this.lat = lat; this.lon = lon;
				++num;
			}
		}

		/// <summary>
		/// Return the results for a polygon.
		/// </summary>
		/// <param name="g">The geod_geodesic object specifying the ellipsoid.</param>
		/// <param name="p">The geod_polygon object specifying the polygon.</param>
		/// <param name="reverse">If set <b>true</b> then clockwise (instead of counter-clockwise)
		/// traversal counts as a positive area.</param>
		/// <param name="sign">If set <b>true</b> then return a signed result for the area if the
		/// polygon is traversed in the "wrong" direction instead of returning the area for the
		/// rest of the earth.</param>
		/// <param name="pA">The area of the polygon (square meters); only set if polyline is set <b>true</b>
		/// in the call to geod_polygon_init().</param>
		/// <param name="pP">The perimeter of the polygon or length of the polyline (meters).</param>
		/// <returns>The number of points.</returns>
		/// <remarks>
		/// The area and perimeter are accumulated at two times the standard floating
		/// point precision to guard against the loss of accuracy with many-sided
		/// polygons. Only simple polygons (which are not self-intersecting) are
		/// allowed. There's no need to "close" the polygon by repeating the first vertex.
		/// 
		/// More points can be added to the polygon after this call.
		/// </remarks>
		public int geod_polygon_compute(geod_geodesic g, bool reverse, bool sign, out double pA, out double pP)
		{
			pP = pA = 0;
			if (num < 2) return num;

			if (polyline)
			{
				pP = P[0];
				return num;
			}

			double s12, S12, dummy;
			g.geod_geninverse(lat, lon, lat0, lon0, out s12, out dummy, out dummy, out dummy, out dummy, out dummy, out S12);

			pP = accsum(P, s12);

			double[] t = new double[2];

			acccopy(A, t);
			accadd(t, S12);
			int crossings = this.crossings + transit(lon, lon0);
			double area0 = 4 * pi * g.c2;

			if ((crossings & 1) != 0) accadd(t, (t[0] < 0 ? 1 : -1) * area0 / 2);

			// area is with the clockwise sense. If !reverse convert to
			// counter-clockwise convention.
			if (!reverse) accneg(t);

			// If sign put area in (-area0/2, area0/2], else put area in [0, area0)
			if (sign)
			{
				if (t[0] > area0 / 2) accadd(t, -area0);
				else if (t[0] <= -area0 / 2) accadd(t, +area0);
			}
			else
			{
				if (t[0] >= area0) accadd(t, -area0);
				else if (t[0] < 0) accadd(t, +area0);
			}

			pA = 0 + t[0];

			return num;
		}

		/// <summary>
		/// Return the results assuming a tentative final test point is added;
		/// however, the data for the test point is not saved. This lets you report a
		/// running result for the perimeter and area as the user moves the mouse
		/// cursor. Ordinary floating point arithmetic is used to accumulate the data
		/// for the test point; thus the area and perimeter returned are less accurate
		/// than if geod_polygon_addpoint() and geod_polygon_compute() are used.
		/// </summary>
		/// <param name="g">The geod_geodesic object specifying the ellipsoid.</param>
		/// <param name="p">The geod_polygon object specifying the polygon.</param>
		/// <param name="lat">The latitude of the test point (degrees).</param>
		/// <param name="lon">The longitude of the test point (degrees).</param>
		/// <param name="reverse">If set <b>true</b> then clockwise (instead of counter-clockwise)
		/// traversal counts as a positive area.</param>
		/// <param name="sign">If set <b>true</b> then return a signed result for the area if the
		/// polygon is traversed in the "wrong" direction instead of returning the area for the
		/// rest of the earth.</param>
		/// <param name="pA">The area of the polygon (square meters); only set if polyline is set <b>true</b>
		/// in the call to geod_polygon_init().</param>
		/// <param name="pP">The perimeter of the polygon or length of the polyline (meters).</param>
		/// <returns>The number of points.</returns>
		/// <remarks>
		/// lat should be in the range [-90deg, 90deg].
		/// </remarks>
		public int geod_polygon_testpoint(geod_geodesic g, double lat, double lon, bool reverse, bool sign, out double pA, out double pP)
		{
			pP = pA = 0;

			int num = this.num + 1;
			if (num == 1) return num;

			double perimeter = P[0];
			double tempsum = polyline ? 0 : A[0];
			int crossings = this.crossings;
			for (int i = 0; i < (polyline ? 1 : 2); ++i)
			{
				double s12, dummy;

				if (!polyline)
				{
					double S12;
					g.geod_geninverse(i == 0 ? this.lat : lat, i == 0 ? this.lon : lon, i != 0 ? lat0 : lat, i != 0 ? lon0 : lon, out s12, out dummy, out dummy, out dummy, out dummy, out dummy, out S12);
					tempsum += S12;
					crossings += transit(i == 0 ? this.lon : lon, i != 0 ? this.lon0 : lon);
				}
				else g.geod_geninverse(i == 0 ? this.lat : lat, i == 0 ? this.lon : lon, i != 0 ? lat0 : lat, i != 0 ? lon0 : lon, out s12, out dummy, out dummy, out dummy, out dummy, out dummy, out dummy);

				perimeter += s12;
			}

			pP = perimeter;
			if (polyline) return num;

			double area0 = 4 * pi * g.c2;
			if ((crossings & 1) != 0) tempsum += (tempsum < 0 ? 1 : -1) * area0 / 2;

			// area is with the clockwise sense. If !reverse convert to counter-clockwise convention.
			if (!reverse) tempsum *= -1;

			// If sign put area in (-area0/2, area0/2], else put area in [0, area0)
			if (sign)
			{
				if (tempsum > area0 / 2) tempsum -= area0;
				else if (tempsum <= -area0 / 2) tempsum += area0;
			}
			else
			{
				if (tempsum >= area0) tempsum -= area0;
				else if (tempsum < 0) tempsum += area0;
			}

			pA = 0 + tempsum;
			return num;
		}

		/// <summary>
		/// Return the results assuming a tentative final test point is added via an
		/// azimuth and distance; however, the data for the test point is not saved.
		/// This lets you report a running result for the perimeter and area as the
		/// user moves the mouse cursor. Ordinary floating point arithmetic is used
		/// to accumulate the data for the test point; thus the area and perimeter
		/// returned are less accurate than if geod_polygon_addedge() and
		/// geod_polygon_compute() are used.
		/// </summary>
		/// <param name="g">The geod_geodesic object specifying the ellipsoid.</param>
		/// <param name="p">The geod_polygon object specifying the polygon.</param>
		/// <param name="azi">The azimuth at current point (degrees).</param>
		/// <param name="s">The distance from current point to final test point (meters).</param>
		/// <param name="reverse">If set <b>true</b> then clockwise (instead of counter-clockwise)
		/// traversal counts as a positive area.</param>
		/// <param name="sign">If set <b>true</b> then return a signed result for the area if the
		/// polygon is traversed in the "wrong" direction instead of returning the area for the
		/// rest of the earth.</param>
		/// <param name="pA">The area of the polygon (square meters); only set if polyline is set <b>true</b>
		/// in the call to geod_polygon_init().</param>
		/// <param name="pP">The perimeter of the polygon or length of the polyline (meters).</param>
		/// <returns>The number of points.</returns>
		public int geod_polygon_testedge(geod_geodesic g, double azi, double s, bool reverse, bool sign, out double pA, out double pP)
		{
			pP = pA = NaN;

			int num = this.num + 1;
			if (num == 1) return 0; // we don't have a starting point!

			double perimeter = P[0] + s;
			if (polyline)
			{
				pP = perimeter;
				return num;
			}

			double tempsum = A[0];
			int crossings = this.crossings;

			{
				double lat, lon, s12, S12, dummy;
				g.geod_gendirect(this.lat, this.lon, azi, GEOD.LONG_UNROLL, s, out lat, out lon, out dummy, out dummy, out dummy, out dummy, out dummy, out S12);

				tempsum += S12;
				crossings += transitdirect(this.lon, lon);
				g.geod_geninverse(lat, lon, lat0, lon0, out s12, out dummy, out dummy, out dummy, out dummy, out dummy, out S12);
				perimeter += s12;
				tempsum += S12;
				crossings += transit(lon, lon0);
			}

			double area0 = 4 * pi * g.c2;
			if ((crossings & 1) != 0) tempsum += (tempsum < 0 ? 1 : -1) * area0 / 2;

			// area is with the clockwise sense. If !reverse convert to counter-clockwise convention.
			if (!reverse) tempsum *= -1;

			// If sign put area in (-area0/2, area0/2], else put area in [0, area0)
			if (sign)
			{
				if (tempsum > area0 / 2) tempsum -= area0;
				else if (tempsum <= -area0 / 2) tempsum += area0;
			}
			else
			{
				if (tempsum >= area0) tempsum -= area0;
				else if (tempsum < 0) tempsum += area0;
			}

			pP = perimeter;
			pA = 0 + tempsum;
			return num;
		}

		/// <summary>
		/// A simple interface for computing the area of a geodesic polygon.
		/// </summary>
		/// <param name="g">The geod_geodesic object specifying the ellipsoid.</param>
		/// <param name="lats">An array of latitudes of the polygon vertices (degrees).</param>
		/// <param name="lons">An array of longitudes of the polygon vertices (degrees).</param>
		/// <param name="n">The number of vertices.</param>
		/// <param name="pA">The area of the polygon (square meters).</param>
		/// <param name="pP">The perimeter of the polygon (meters).</param>
		/// <remarks>
		/// lats should be in the range [-90deg, 90deg].
		///
		/// Only simple polygons (which are not self-intersecting) are allowed.
		/// There's no need to "close" the polygon by repeating the first vertex. The
		/// area returned is signed with counter-clockwise traversal being treated as
		/// positive.
		/// </remarks>
		public static void geod_polygonarea(geod_geodesic g, double[] lats, double[] lons, int n, out double pA, out double pP)
		{
			geod_polygon p = new geod_polygon(false);
			for (int i = 0; i < n; ++i) p.geod_polygon_addpoint(g, lats[i], lons[i]);
			p.geod_polygon_compute(g, false, true, out pA, out pP);
		}
	}
}
