#pragma once

#include "earshot/earshot.h"
#include "earshot/testing.h"

/**
	A constant-Q third-octave filter bank at the ISO 266 / IEC 61260-1 centres,
	for needle's VFD spectrum analyser. 64 linear host bins cannot resolve the
	low thirds (the 25 Hz band is 6 Hz wide; one host bin is 344 Hz), which is
	why this needs earshot's capture at all.

	## The bands

	Base ten, as IEC 61260-1:2014 requires: G = 10^(3/10), exact mid-band
	frequencies fm = 1000 * G^(x/3) = 1000 * 10^(x/10) for x = -17 .. 13
	(19.95 Hz .. 19.95 kHz, labelled 20 .. 20k by ISO 266's R10 series), band
	edges fm * G^(-1/6) and fm * G^(1/6). A band is built only when its upper
	edge is below kTopEdge of Nyquist: at 44.1 kHz that drops the 20 kHz band
	(its upper edge, 22.4 kHz, is above Nyquist), at 48 kHz it keeps it.

	## Each filter

	A third-order Butterworth low-pass prototype, transformed to a band-pass
	(six poles, three biquads) with the band edges pre-warped, so the digital
	filter's -3.01 dB points land exactly on the IEC band edges whatever the
	sample rate. Direct form II transposed in double: the 20 Hz band's poles
	sit 0.0004 from the unit circle at 48 kHz, where float coefficients would
	put the band somewhere else. `estest --thirds` drives sines through these
	filters and judges the measured attenuation against the class 1 limits.

	## Level and ballistics

	Mean square through a one-pole of 125 ms (IEC 61672's Fast), primed with
	the first block's mean square so a display does not rise from silence on
	the first frame. Peak hold: the hold value follows the level up, stays for
	thirdHoldSeconds, then falls at thirdFallDbPerSecond, in audio time.

	**Latency**: the filters' group delay at mid-band (about 3/(pi BW): 160 ms
	at 25 Hz, 4 ms at 1 kHz) plus the 125 ms weighting. **Thread**: render
	thread, inside Engine::update, on every sample in sequence: 31 bands x 3
	biquads, about 4.5 M biquad steps a second at 48 kHz. No allocation.
*/
namespace earshot
{

class ThirdOctaveBank
{
public:
	/// IEC 61260-1 band index x (fm = 1000 * 10^(x/10)) for the lowest band.
	static constexpr int kFirstX = -17;
	/// Bands are built only below this fraction of Nyquist (upper edge).
	static constexpr double kTopEdge = 0.95;

	struct Section
	{
		double b0 = 0, b1 = 0, b2 = 0, a1 = 0, a2 = 0;
		double s1 = 0, s2 = 0;
	};

	struct BandFilter
	{
		Section sec[ 3 ];
		int     sections = 3;
		double  exactHz = 0, loHz = 0, hiHz = 0;
		float   nominalHz = 0;

		/// Design for sample rate fs. `order` is the low-pass prototype order (3).
		void design( double fm, double fs, int order, bool prewarp, double G ) noexcept;
		void clear() noexcept;
		double process( double x ) noexcept
		{
			for( int i = 0; i < sections; ++i )
			{
				Section&     s = sec[ i ];
				const double y = s.b0 * x + s.s1;
				s.s1           = s.b1 * x - s.a1 * y + s.s2;
				s.s2           = s.b2 * x - s.a2 * y;
				x              = y;
			}
			return x;
		}
	};

	/// Design every band for `fs`. Clears filter state.
	void setRate( int fs, const testing::Mutations& mut ) noexcept;
	/// A new stream: clear state, prime the level on the next block.
	void reset() noexcept;

	/// Mono samples in sequence.
	void push( const float* mono, int n ) noexcept;
	/// After a frame's samples: ballistics over `dt` audio seconds.
	void publish( Third* out, int& count, double dt, const Config& config ) noexcept;

	int        bandCount() const noexcept { return count_; }
	BandFilter& band( int i ) noexcept { return bands_[ i ]; }

	/// The ISO 266 R10 label for band index x.
	static float nominalFor( int x ) noexcept;

private:
	BandFilter bands_[ kMaxThirds ];
	double     ms_[ kMaxThirds ]    = {};
	double     hold_[ kMaxThirds ]  = {};    ///< mean square at the hold
	double     held_[ kMaxThirds ]  = {};    ///< seconds since the hold was set
	double     blockSq_[ kMaxThirds ] = {};
	int        blockN_ = 0;
	int        count_  = 0;
	int        fs_     = 0;
	double     alpha_  = 0.0;   ///< per-sample one-pole for 125 ms
	bool       primed_ = false;
};

} // namespace earshot
