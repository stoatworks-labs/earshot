#include "Onsets.h"

#include "Dsp.h"
#include "Tempo.h"

#include <algorithm>
#include <cmath>

namespace earshot
{
namespace
{
constexpr double kWindowSeconds  = 0.0213;   ///< 1024 samples at 48 kHz
constexpr double kLogKneeDb      = -70.0;    ///< log(1 + g|X|) turns logarithmic above this
constexpr double kFloorSeconds   = 1.0;      ///< the running mean's time constant (pattern's)
constexpr double kPrimeFraction  = 0.125;    ///< of the band's log level: the floor's seed (pattern v0.1.1)
constexpr double kAbsoluteFloor  = 1.0;      ///< flux never counts below this: a -65 dBFS sine appearing
constexpr double kRegisterLowHz  = 400.0;    ///< tempo registers: pattern's split
constexpr double kRegisterHighHz = 5000.0;

int windowSizeFor( int fs )
{
	const double want = kWindowSeconds * fs;
	int          n    = dsp::nextPow2( (int)want );
	if( (double)n / want > std::sqrt( 2.0 ) )
		n >>= 1;
	return std::max( 256, n );
}
} // namespace

bool BandAnalyser::init( int maxRate ) noexcept
{
	const int n = windowSizeFor( maxRate );
	if( !fft_.init( n ) )
		return false;
	try
	{
		window_.assign( (size_t)n, 0.0f );
		input_.assign( (size_t)n, 0.0f );
		spec_.assign( (size_t)n, std::complex< float >() );
		logMag_.assign( (size_t)n / 2 + 1, 0.0f );
		prevLogMag_.assign( (size_t)n / 2 + 1, 0.0f );
		power_.assign( (size_t)n / 2 + 1, 0.0 );
	}
	catch( ... )
	{
		return false;
	}
	return true;
}

void BandAnalyser::reset() noexcept
{
	havePrev_ = false;
	filled_   = 0;
	writeAt_  = 0;
	sinceHop_ = 0;
	for( State& b : bands_ )
	{
		b.primed   = false;
		b.onset    = false;
		b.lastFire = -1e9;
	}
}

void BandAnalyser::setRate( int fs ) noexcept
{
	fs_   = fs;
	size_ = std::min( windowSizeFor( fs ), fft_.maxSize() );
	hop_  = size_ / 4;
	windowSum_    = 0.0;
	windowEnergy_ = 0.0;
	for( int i = 0; i < size_; ++i )
	{
		window_[ (size_t)i ] = dsp::hann( i, size_ );
		windowSum_ += window_[ (size_t)i ];
		windowEnergy_ += (double)window_[ (size_t)i ] * window_[ (size_t)i ];
	}
	const double binHz = (double)fs_ / size_;
	regLo_ = std::max( 2, (int)std::ceil( kRegisterLowHz / binHz ) );
	regHi_ = std::max( regLo_ + 1, (int)std::ceil( kRegisterHighHz / binHz ) );
	regHi_ = std::min( regHi_, size_ / 2 );
	assignBins();
	reset();
}

void BandAnalyser::configure( const Config& config ) noexcept
{
	const int count = std::max( 1, std::min( config.bandCount, kMaxBands ) );
	for( int b = 0; b < count; ++b )
	{
		const float lo = config.bandEdges[ b ];
		const float hi = config.bandEdges[ b + 1 ];
		State&      s  = bands_[ b ];
		if( b >= bandCount_ || s.loHz != lo || s.hiHz != hi )
		{
			// Only a band whose edges moved starts again; its counter carries on,
			// so a plugin comparing onsetCount never sees it go backwards.
			s.loHz   = lo;
			s.hiHz   = hi;
			s.primed = false;
		}
	}
	bandCount_ = count;
	assignBins();
}

void BandAnalyser::assignBins() noexcept
{
	if( fs_ <= 0 || size_ <= 0 )
		return;
	const double binHz = (double)fs_ / size_;
	const int    top   = size_ / 2;
	for( int b = 0; b < bandCount_; ++b )
	{
		State& s = bands_[ b ];
		// A bin belongs to the band its centre frequency falls in.
		int lo = (int)std::ceil( s.loHz / binHz );
		int hi = (int)std::ceil( s.hiHz / binHz );
		lo     = std::max( 1, std::min( lo, top - 1 ) );
		hi     = std::max( lo + 1, std::min( hi, top ) );
		s.binLo = lo;
		s.binHi = hi;
	}
}

void BandAnalyser::push( const float* mono, int n, uint64_t firstSample, const Config& config,
						 const testing::Mutations& mut, TempoTracker* tempo ) noexcept
{
	if( size_ <= 0 || mono == nullptr )
		return;
	for( int i = 0; i < n; ++i )
	{
		input_[ (size_t)writeAt_ ] = mono[ i ];
		writeAt_ = ( writeAt_ + 1 == size_ ) ? 0 : writeAt_ + 1;
		++filled_;
		++sinceHop_;
		if( filled_ >= (uint64_t)size_ && sinceHop_ >= hop_ )
		{
			sinceHop_ = 0;
			hop( firstSample + (uint64_t)i + 1, config, mut, tempo );
		}
	}
}

void BandAnalyser::hop( uint64_t endSample, const Config& config, const testing::Mutations& mut, TempoTracker* tempo ) noexcept
{
	const int    n     = size_;
	const int    half  = n / 2;
	const double t     = (double)endSample / (double)fs_;
	const double dt    = (double)hop_ / (double)fs_;
	const double gain  = std::pow( 10.0, -kLogKneeDb / 20.0 );
	const double toAmp = 2.0 / windowSum_;
	// Mean square of the band-limited signal from its bins: a sine of
	// amplitude A puts A^2/4 * N * sum(w^2) into the positive half (Parseval).
	const double toMs  = mut.bandNoParseval ? 2.0 / ( (double)n * n ) : 2.0 / ( (double)n * windowEnergy_ );

	for( int i = 0; i < n; ++i )
	{
		const int idx = ( writeAt_ + i ) % n;   // oldest first
		spec_[ (size_t)i ] = std::complex< float >( input_[ (size_t)idx ] * window_[ (size_t)i ], 0.0f );
	}
	fft_.forward( spec_.data(), n );
	for( int k = 0; k <= half; ++k )
	{
		const double p = std::norm( spec_[ (size_t)k ] );
		power_[ (size_t)k ]  = p;
		logMag_[ (size_t)k ] = (float)std::log1p( gain * std::sqrt( p ) * toAmp );
	}

	const bool first = !havePrev_;
	if( first )
	{
		for( int k = 0; k <= half; ++k )
			prevLogMag_[ (size_t)k ] = mut.onsetNoPrime ? 0.0f : logMag_[ (size_t)k ];
	}

	// Tempo registers read every bin, whatever the user's bands.
	double reg[ 3 ] = { 0.0, 0.0, 0.0 };
	for( int k = 1; k < half; ++k )
	{
		const double rise = std::max( 0.0f, logMag_[ (size_t)k ] - prevLogMag_[ (size_t)k ] );
		reg[ k < regLo_ ? 0 : ( k < regHi_ ? 1 : 2 ) ] += rise;
	}

	double energy[ kMaxBands ], flux[ kMaxBands ], logLevel[ kMaxBands ], rise[ kMaxBands ];
	double maxRise = 0.0;
	for( int b = 0; b < bandCount_; ++b )
	{
		const State& s = bands_[ b ];
		double e = 0.0, f = 0.0, l = 0.0;
		for( int k = s.binLo; k < s.binHi; ++k )
		{
			e += power_[ (size_t)k ];
			f += std::max( 0.0f, logMag_[ (size_t)k ] - prevLogMag_[ (size_t)k ] );
			l += logMag_[ (size_t)k ];
		}
		energy[ b ]   = e * toMs;
		flux[ b ]     = f;
		logLevel[ b ] = l;
		const bool comparable = s.primed || ( first && mut.onsetNoPrime );
		const double prevE    = ( first && mut.onsetNoPrime ) ? 0.0 : s.prevEnergy;
		rise[ b ] = comparable ? std::max( 0.0, energy[ b ] - prevE ) : 0.0;
		maxRise   = std::max( maxRise, rise[ b ] );
	}

	const double floorAlpha = dsp::onePoleAlpha( dt, kFloorSeconds );
	const double attack     = dsp::onePoleAlpha( dt, config.envelopeAttack );
	const double release    = dsp::onePoleAlpha( dt, config.envelopeRelease );
	const double maskRatio  = std::pow( 10.0, -config.onsetMaskDb / 10.0 );

	for( int b = 0; b < bandCount_; ++b )
	{
		State& s = bands_[ b ];
		s.level  = std::sqrt( 2.0 * energy[ b ] );

		if( !s.primed && !( first && mut.onsetNoPrime ) )
		{
			// Prime: nothing fires, the floor is seeded, the envelope starts at
			// the level rather than rising from silence.
			s.meanFlux   = ( mut.onsetFullSeed ? 1.0 : kPrimeFraction ) * logLevel[ b ];
			s.prevFlux   = 0.0;
			s.prevEnergy = energy[ b ];
			s.envelope   = s.level;
			s.primed     = true;
			s.above      = false;
			continue;
		}
		if( first && mut.onsetNoPrime )
		{
			s.meanFlux = 0.0;
			s.prevFlux = 0.0;
			s.primed   = true;
			s.above    = false;
		}

		const double threshold = std::max( kAbsoluteFloor, (double)config.onsetRatio * s.meanFlux );
		const bool   over      = flux[ b ] > threshold;
		const bool   rising    = mut.onsetNoRisingEdge || !s.above;
		const bool   armed     = t - s.lastFire >= (double)config.onsetRefractory;
		const bool   masked    = !mut.onsetNoMask && rise[ b ] < maxRise * maskRatio;
		s.above = over;
		if( over && rising && armed && !masked )
		{
			s.onset    = true;
			s.count++;
			s.strength = flux[ b ] / threshold;
			s.lastTime = t;
			s.lastFire = t;
		}

		s.meanFlux += ( flux[ b ] - s.meanFlux ) * floorAlpha;
		s.prevFlux   = flux[ b ];
		s.prevEnergy = energy[ b ];
		s.envelope += ( s.level - s.envelope ) * ( s.level > s.envelope ? attack : release );
	}

	for( int k = 0; k <= half; ++k )
		prevLogMag_[ (size_t)k ] = logMag_[ (size_t)k ];
	havePrev_ = true;

	if( tempo != nullptr )
		tempo->push( t, reg, (double)fs_ / hop_ );
}

void BandAnalyser::publish( Band* out, int& count ) noexcept
{
	count = bandCount_;
	for( int b = 0; b < bandCount_; ++b )
	{
		State& s = bands_[ b ];
		Band&  o = out[ b ];
		o.loHz          = s.loHz;
		o.hiHz          = s.hiHz;
		o.level         = (float)s.level;
		o.db            = (float)dsp::dbFromAmplitude( s.level );
		o.envelope      = (float)s.envelope;
		o.onset         = s.onset;
		o.onsetCount    = s.count;
		o.onsetStrength = (float)s.strength;
		o.lastOnsetTime = s.lastTime;
		s.onset         = false;
	}
}

} // namespace earshot
