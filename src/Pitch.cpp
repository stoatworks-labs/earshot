#include "Pitch.h"

#include "Dsp.h"

#include <algorithm>
#include <cmath>

namespace earshot
{
namespace
{
/// MPM's k: the first key maximum within this fraction of the highest wins.
/// McLeod & Wyvill recommend 0.8 to 1.0; 0.9 keeps a strong second harmonic
/// from winning until the fundamental is 13 dB down (see Pitch.h).
constexpr float kKeyMaxThreshold = 0.9f;

/// Three periods of the lowest pitch.
constexpr double kPeriodsPerWindow = 3.0;

float effectiveMinHz( const Config& c, float floorHz )
{
	return std::max( floorHz, std::min( c.pitchMinHz, c.pitchMaxHz * 0.5f ) );
}
} // namespace

bool PitchDetector::init( int maxRate, float minHzFloor ) noexcept
{
	minHzFloor_ = minHzFloor;
	maxWindow_  = (int)std::ceil( kPeriodsPerWindow * (double)maxRate / (double)minHzFloor );
	const int n = dsp::nextPow2( 2 * maxWindow_ );
	try
	{
		x_.assign( (size_t)maxWindow_, 0.0f );
		nsdf_.assign( (size_t)maxWindow_ + 2, 0.0f );
		buf_.assign( (size_t)n, std::complex< float >() );
	}
	catch( ... )
	{
		return false;
	}
	return true;
}

int PitchDetector::windowFor( int fs, const Config& config ) const noexcept
{
	const double minHz = effectiveMinHz( config, minHzFloor_ );
	return std::min( maxWindow_, (int)std::ceil( kPeriodsPerWindow * (double)fs / minHz ) );
}

void PitchDetector::analyse( const float* mono, int n, int fs, const Config& config, const testing::Mutations& mut,
							 const Fft& fft, Pitch& out ) noexcept
{
	out        = Pitch();
	out.valid  = true;
	if( mono == nullptr || n < 8 || n > maxWindow_ || fs <= 0 )
		return;

	// Remove DC first: a cheap interface's offset is a zero-frequency "period"
	// that keeps n(tau) positive and hides every zero crossing.
	double mean = 0.0;
	for( int i = 0; i < n; ++i )
		mean += mono[ i ];
	mean /= n;
	double energy = 0.0;
	for( int i = 0; i < n; ++i )
	{
		x_[ (size_t)i ] = (float)( mono[ i ] - mean );
		energy += (double)x_[ (size_t)i ] * x_[ (size_t)i ];
	}
	const double rms = std::sqrt( energy / n );
	out.levelDb      = (float)dsp::dbFromAmplitude( rms * std::sqrt( 2.0 ) );

	if( mut.pitchSpectralPeak )
	{
		spectralPeak( n, fs, config, fft, out );
		return;
	}
	if( !( energy > 0.0 ) )
		return;

	// r(tau) by FFT: zero-pad to at least 2n so the circular correlation is the
	// linear one, |X|^2, inverse.
	const int nfft = dsp::nextPow2( 2 * n );
	if( nfft > (int)buf_.size() || nfft > fft.maxSize() )
		return;
	for( int i = 0; i < nfft; ++i )
		buf_[ (size_t)i ] = std::complex< float >( i < n ? x_[ (size_t)i ] : 0.0f, 0.0f );
	fft.forward( buf_.data(), nfft );
	for( int i = 0; i < nfft; ++i )
		buf_[ (size_t)i ] = std::complex< float >( std::norm( buf_[ (size_t)i ] ), 0.0f );
	fft.inverse( buf_.data(), nfft );

	const double minHz  = effectiveMinHz( config, minHzFloor_ );
	const double maxHz  = std::max( minHz * 2.0, (double)config.pitchMaxHz );
	const int    tauMax = std::min( n - 2, (int)std::ceil( (double)fs / minHz ) + 1 );
	const double tauMin = (double)fs / maxHz;

	// m(tau), incrementally, in double: m(0) = 2 sum x^2, and each step drops
	// the sample that leaves each end of the overlap.
	double m = 2.0 * energy;
	nsdf_[ 0 ] = 1.0f;
	const float scale = 1.0f / (float)nfft;
	for( int tau = 1; tau <= tauMax + 1 && tau < n; ++tau )
	{
		const double a = x_[ (size_t)( tau - 1 ) ];
		const double b = x_[ (size_t)( n - tau ) ];
		m -= a * a + b * b;
		const double r = (double)buf_[ (size_t)tau ].real() * scale;
		nsdf_[ (size_t)tau ] = ( m > 1e-30 ) ? (float)( 2.0 * r / m ) : 0.0f;
	}

	// Key maxima: one per positive lobe, starting after the first
	// negative-going crossing so the lobe around tau = 0 never counts.
	constexpr int kMaxKeys = 512;
	int           keyLag[ kMaxKeys ];
	float         keyVal[ kMaxKeys ];
	int           keys    = 0;
	bool          started = false;   // passed the first negative-going crossing
	bool          inLobe  = false;
	int           bestLag = 0;
	float         bestVal = 0.0f;
	for( int tau = 1; tau <= tauMax; ++tau )
	{
		const float prev = nsdf_[ (size_t)( tau - 1 ) ];
		const float cur  = nsdf_[ (size_t)tau ];
		if( !started )
		{
			if( prev > 0.0f && cur <= 0.0f )
				started = true;
			continue;
		}
		if( !inLobe && prev <= 0.0f && cur > 0.0f )
		{
			inLobe  = true;
			bestLag = tau;
			bestVal = cur;
		}
		else if( inLobe && cur <= 0.0f )
		{
			inLobe = false;
			if( keys < kMaxKeys )
			{
				keyLag[ keys ] = bestLag;
				keyVal[ keys ] = bestVal;
				++keys;
			}
		}
		else if( inLobe && cur > bestVal )
		{
			bestLag = tau;
			bestVal = cur;
		}
	}
	if( inLobe && keys < kMaxKeys && bestLag < tauMax )
	{
		keyLag[ keys ] = bestLag;
		keyVal[ keys ] = bestVal;
		++keys;
	}
	if( keys == 0 )
		return;

	float highest = 0.0f;
	int   highestIndex = 0;
	for( int i = 0; i < keys; ++i )
		if( keyVal[ i ] > highest )
		{
			highest      = keyVal[ i ];
			highestIndex = i;
		}
	if( !( highest > 0.0f ) )
		return;

	int chosen = highestIndex;
	if( !mut.pitchGlobalMax )
	{
		for( int i = 0; i < keys; ++i )
			if( keyVal[ i ] >= kKeyMaxThreshold * highest )
			{
				chosen = i;
				break;
			}
	}

	const int tau = keyLag[ chosen ];
	double    lag = tau;
	double    peak = keyVal[ chosen ];
	if( !mut.pitchNoInterp && tau >= 1 && tau + 1 <= tauMax + 1 )
	{
		const double a = nsdf_[ (size_t)( tau - 1 ) ], b = nsdf_[ (size_t)tau ], c = nsdf_[ (size_t)( tau + 1 ) ];
		const double d = a - 2.0 * b + c;
		if( d < 0.0 )
		{
			const double delta = 0.5 * ( a - c ) / d;
			if( std::fabs( delta ) <= 1.0 )
			{
				lag += delta;
				peak = b - 0.25 * ( a - c ) * delta;
			}
		}
	}

	out.confidence = (float)dsp::clampd( peak, 0.0, 1.0 );
	out.hz         = (float)( (double)fs / lag );
	out.midi       = (float)dsp::midiFromHz( out.hz );
	const bool inRange = lag >= tauMin * 0.999 && out.hz >= minHz * 0.999;
	out.voiced = inRange && out.confidence >= config.pitchClarity && out.levelDb >= config.pitchGateDb;
}

/// The negative control for missing fundamentals and strong harmonics: the
/// loudest spectral peak in range, which is what a naive "pitch from the
/// spectrum" reads. Not reachable from Config.
void PitchDetector::spectralPeak( int n, int fs, const Config& config, const Fft& fft, Pitch& out ) noexcept
{
	const int nfft = dsp::nextPow2( 2 * n );
	if( nfft > (int)buf_.size() || nfft > fft.maxSize() )
		return;
	for( int i = 0; i < nfft; ++i )
		buf_[ (size_t)i ] = std::complex< float >( i < n ? x_[ (size_t)i ] * dsp::hann( i, n ) : 0.0f, 0.0f );
	fft.forward( buf_.data(), nfft );
	const double binHz = (double)fs / nfft;
	const int    lo    = std::max( 1, (int)( config.pitchMinHz / binHz ) );
	const int    hi    = std::min( nfft / 2 - 2, (int)( config.pitchMaxHz / binHz ) + 1 );
	int          best  = lo;
	for( int k = lo; k <= hi; ++k )
		if( std::norm( buf_[ (size_t)k ] ) > std::norm( buf_[ (size_t)best ] ) )
			best = k;
	const double a = std::log( std::abs( buf_[ (size_t)( best - 1 ) ] ) + 1e-30 );
	const double b = std::log( std::abs( buf_[ (size_t)best ] ) + 1e-30 );
	const double c = std::log( std::abs( buf_[ (size_t)( best + 1 ) ] ) + 1e-30 );
	const double d = a - 2.0 * b + c;
	const double delta = d < 0.0 ? 0.5 * ( a - c ) / d : 0.0;
	out.hz         = (float)( ( best + delta ) * binHz );
	out.midi       = (float)dsp::midiFromHz( out.hz );
	out.confidence = 1.0f;
	out.voiced     = out.levelDb >= config.pitchGateDb;
}

} // namespace earshot
