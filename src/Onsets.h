#pragma once

#include "Fft.h"

#include "earshot/earshot.h"
#include "earshot/testing.h"

#include <complex>
#include <cstdint>
#include <vector>

/**
	Band envelopes and onsets, from a short-time spectrum on a fixed AUDIO hop.

	The hop is the point. Every onset detector in the fleet runs on the video
	frame -- the host hands over one spectrum per frame -- and pattern measured
	what that costs: the same groove read 124.84 BPM at 60 fps and 62.06 at
	30, because a 30 fps host quantises every onset to 33 ms. With samples of
	its own, earshot steps a 21 ms Hann window by a quarter (1024 / 256 at 44.1
	and 48 kHz, scaled by powers of two above), so an onset's time and
	everything downstream are the same at 24 fps and 144, at any block size.
	Onsets carry their audio time.

	## Per hop, per band

	- **Level**: Parseval over the band's bins with the window's energy
	  correction, as the amplitude of the equivalent sine (a sine of amplitude
	  A well inside the band reads A).
	- **Flux**: half-wave-rectified rise of log(1 + g |X|), the log knee at
	  -70 dBFS, summed over the band. Log compression is what lets a soft
	  attack register at all (Boeck et al., DAFx 2012, rank it best online).
	- **Threshold**: the larger of an absolute floor and onsetRatio times a
	  one-second running mean of the flux. Fires on the hop where the flux
	  CROSSES it from below -- one onset per excursion, however many hops the
	  window smears an attack over -- after a refractory period, and only if
	  the masking gate allows.
	- **Masking gate** (pattern's lesson). A hard attack's first milliseconds
	  are broadband: a 150 Hz burst lit all eight of pattern's channels on its
	  first frame. Here a band may fire only if its energy RISE in this hop is
	  within onsetMaskDb of the largest band's rise. A click 30 dB under the
	  attack it belongs to stays in its band; a clap, which really is
	  broadband, lights every band it reaches.

	## Frame one (the fleet's deaf-for-1.5-s bug)

	A detector that starts from zero reads every bin as risen from silence on
	its first frame; if its floor then snaps to that flux it is deaf until the
	floor decays. So the first hop after any (re)start primes and fires
	nothing: the previous spectrum is set to the current one, and each band's
	floor is seeded with an eighth of its log level (pattern v0.1.1's rule: a
	stationary signal's hop-to-hop flux is a few per cent of its level, a real
	onset's is comparable to it, so the seed is above the first and below the
	second). The hop's dt is a constant number of samples, so dt = 0 cannot
	happen on this path; the host path guards it (HostFft.h).

	**Latency**: an onset is reported at the end of the first hop whose window
	sees it, within one window (21 ms) of the attack; the plugin sees it at
	its next update. **Thread**: render thread, inside Engine::update, on the
	audio in sequence. Allocates only in init().
*/
namespace earshot
{

class TempoTracker;

class BandAnalyser
{
public:
	bool init( int maxRate ) noexcept;

	/// A new stream or a gap: the next hop primes.
	void reset() noexcept;
	/// New sample rate: window and hop.
	void setRate( int fs ) noexcept;
	/// Band edges and timings. No allocation; keeps every band's history when
	/// the edges are unchanged, re-primes only the bands whose edges moved.
	void configure( const Config& config ) noexcept;

	/// Mono samples in sequence. `firstSample` is the absolute index of mono[0].
	void push( const float* mono, int n, uint64_t firstSample, const Config& config,
			   const testing::Mutations& mut, TempoTracker* tempo ) noexcept;

	/// Copy results into the frame's bands (and clear the per-frame onset flags).
	void publish( Band* out, int& count ) noexcept;

	int hopSize() const noexcept { return hop_; }
	int windowSize() const noexcept { return size_; }

private:
	void hop( uint64_t endSample, const Config& config, const testing::Mutations& mut, TempoTracker* tempo ) noexcept;
	void assignBins() noexcept;

	struct State
	{
		float  loHz = 0, hiHz = 0;
		int    binLo = 0, binHi = 0;   ///< [binLo, binHi)
		double meanFlux = 0, prevFlux = 0, prevEnergy = 0;
		double level = 0, envelope = 0;
		double lastFire = -1e9;
		bool   primed = false;
		bool   above  = false;   ///< flux was over the threshold on the previous hop
		bool   onset  = false;
		uint32_t count = 0;
		double strength = 0, lastTime = -1;
	};

	int    fs_ = 0, size_ = 0, hop_ = 0;
	int    bandCount_ = 0;
	State  bands_[ kMaxBands ];
	int    regLo_ = 0, regHi_ = 0;     ///< tempo registers: [1, regLo), [regLo, regHi), [regHi, N/2)
	bool   havePrev_ = false;

	Fft                                  fft_;
	std::vector< float >                 window_;
	std::vector< float >                 input_;     ///< circular, size_ samples
	std::vector< std::complex< float > > spec_;
	std::vector< float >                 logMag_, prevLogMag_;
	std::vector< double >                power_;
	uint64_t                             filled_ = 0;    ///< samples since reset
	int                                  writeAt_ = 0;
	int                                  sinceHop_ = 0;
	double                               windowEnergy_ = 1.0, windowSum_ = 1.0;
};

} // namespace earshot
