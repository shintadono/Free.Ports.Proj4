﻿// Copyright (c) 2014 Bojan Savric
// Copyright (c) 2016 by the Authors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
//
//
// The Patterson Cylindrical projection was designed by Tom Patterson, US National
// Park Service, in 2014, using Flex Projector. The polynomial equations for the
// projection were developed by Bojan Savric, Oregon State University, in
// collaboration with Tom Patterson and Bernhard Jenny, Oregon State University.
//
// Java reference algorithm implemented by Bojan Savric in Java Map Projection
// Library (a Java port of PROJ.4) in the file PattersonProjection.java.
//
// References:
//    Java Map Projection Library
//       https://github.com/OSUCartography/JMapProjLib
// 
//    Patterson Cylindrical Projection
//       http://shadedrelief.com/patterson/
// 
//    Patterson, T., Savric, B., and Jenny, B. (2015). Cartographic Perspectives
//       (No.78). Describes the projection design and characteristics, and
//       developing the equations.    doi:10.14714/CP78.1270
//       http://dx.doi.org/10.14714/CP78.1270
//
// Port to PROJ.4 by Micah Cochran, 26 March 2016

using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_patterson : PJ
	{
		const double K1 = 1.0148;
		const double K2 = 0.23185;
		const double K3 = -0.14499;
		const double K4 = 0.02406;
		const double C1 = K1;
		const double C2 = 5.0 * K2;
		const double C3 = 7.0 * K3;
		const double C4 = 9.0 * K4;
		const double EPS11 = 1.0e-11;
		const double MAX_Y = 1.790857183;

		public override string Name { get { return "patterson"; } }
		public override string DescriptionName { get { return "Patterson Cylindrical"; } }
		public override string DescriptionType { get { return "Cyl."; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;

			double phi2 = lp.phi * lp.phi;
			xy.x = lp.lam;
			xy.y = lp.phi * (K1 + phi2 * phi2 * (K2 + phi2 * (K3 + K4 * phi2)));

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			double yc = xy.y;

			// make sure y is inside valid range
			if (xy.y > MAX_Y) xy.y = MAX_Y;
			else if (xy.y < -MAX_Y) xy.y = -MAX_Y;

			for (;;)
			{ // Newton-Raphson
				double y2 = yc * yc;
				double f = yc * (K1 + y2 * y2 * (K2 + y2 * (K3 + K4 * y2))) - xy.y;
				double fder = C1 + y2 * y2 * (C2 + y2 * (C3 + C4 * y2));
				double tol = f / fder;
				yc -= tol;

				if (Math.Abs(tol) < EPS11) break;
			}

			LP lp;
			lp.phi = yc;
			lp.lam = xy.x; // longitude

			return lp;
		}

		public override PJ Init()
		{
			es = 0;
			fwd = s_forward;
			inv = s_inverse;

			return this;
		}
	}
}