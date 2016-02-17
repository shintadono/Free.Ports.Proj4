using System;

namespace Free.Ports.Proj4.Geodesic
{
	/// <summary>
	/// Capabilities.
	/// </summary>
	[Flags]
	public enum GEOD
	{
		CAP_NONE = 0,

		CAP_C1 = 0x0001,
		CAP_C1p = 0x0002,
		CAP_C2 = 0x0004,
		CAP_C3 = 0x0008,
		CAP_C4 = 0x0010,

		CAP_ALL = CAP_C1 | CAP_C1p | CAP_C2 | CAP_C3 | CAP_C4,

		OUT_ALL = 0x7F80,

		#region Mask values for the caps argument to geod_lineinit().
		/// <summary>
		/// Calculate nothing.
		/// </summary>
		NONE = 0,

		/// <summary>
		/// Calculate latitude.
		/// </summary>
		LATITUDE = 0x0080 | CAP_NONE,

		/// <summary>
		/// Calculate longitude.
		/// </summary>
		LONGITUDE = 0x0100 | CAP_C3,

		/// <summary>
		/// Calculate azimuth.
		/// </summary>
		AZIMUTH = 0x0200 | CAP_NONE,

		/// <summary>
		/// Calculate distance.
		/// </summary>
		DISTANCE = 0x0400 | CAP_C1,

		/// <summary>
		/// Allow distance as input.
		/// </summary>
		DISTANCE_IN = 0x0800 | CAP_C1 | CAP_C1p,

		/// <summary>
		/// Calculate reduced length.
		/// </summary>
		REDUCEDLENGTH = 0x1000 | CAP_C1 | CAP_C2,

		/// <summary>
		/// Calculate geodesic scale.
		/// </summary>
		GEODESICSCALE = 0x2000 | CAP_C1 | CAP_C2,

		/// <summary>
		/// Calculate reduced length.
		/// </summary>
		AREA = 0x4000 | CAP_C4,

		/// <summary>
		/// Calculate everything.
		/// </summary>
		ALL = OUT_ALL | CAP_ALL,
		#endregion

		#region Flag values for the flags argument to geod_gendirect() and geod_genposition().
		/// <summary>
		/// No flags.
		/// </summary>
		NOFLAGS = 0,

		/// <summary>
		/// Position given in terms of arc distance.
		/// </summary>
		ARCMODE = CAP_C1,

		/// <summary>
		/// Unroll the longitude.
		/// </summary>
		LONG_UNROLL = 0x8000,
		#endregion
	}
}
