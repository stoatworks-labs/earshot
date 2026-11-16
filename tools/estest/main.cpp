// estest: earshot's headless harness.
//
//   estest                 every group, every cell of the matrix
//   estest --pitch --tap   just those groups
//   estest --rate 48000 --block 128      restrict the matrix
//   estest --negative      every mutation must make the checks it targets fail
//   estest --file x.wav [--expect-bpm 125] [--from 6]
//   estest --perf          analysis cost at 48 kHz
//   estest --devices       list capture devices (never opens one)
//   estest --live 10 [--device NAME]     open a real device for 10 s (opt-in;
//                          may raise the microphone-permission prompt)

#include "Harness.h"
#include "Signals.h"
#include "Wav.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <set>
#include <string>
#include <thread>
#include <vector>

using namespace estest;

namespace
{
struct Options
{
	std::set< std::string > groups;
	int                     rate = 0, block = 0;
	bool                    negative = false, perf = false, devices = false;
	double                  live = 0.0;
	std::string             device = earshot::kDefaultDevice;
	std::string             file;
	double                  expectBpm = 0.0, from = 6.0, fps = 60.0;
	std::string             writeWav;
	std::string             only;   ///< --negative: run just the mutations whose names contain this
	std::string             with;   ///< run the groups with one negative case's mutation set
};

std::vector< Run > cells( const Options& o )
{
	std::vector< Run > out;
	for( const Run& r : matrix() )
		if( ( o.rate == 0 || r.fs == o.rate ) && ( o.block == 0 || r.block == o.block ) )
			out.push_back( r );
	return out;
}

void runGroups( Report& rep, const Options& o, const std::set< std::string >& names, const Mutations& mut )
{
	const std::vector< Run > cs = cells( o );
	for( const Group& g : groups() )
	{
		if( !names.empty() && names.count( g.name ) == 0 )
			continue;
		if( !rep.quiet )
			std::printf( "\n[%s] %s\n", g.name.c_str(), g.what.c_str() );
		if( g.perRun )
			for( const Run& r : cs )
				g.fn( rep, r, mut );
		else
			g.fn( rep, cs.empty() ? matrix()[ 0 ] : cs[ 0 ], mut );
	}
}

//---------------------------------------------------------------------------
// Negative controls
//---------------------------------------------------------------------------

struct NegativeCase
{
	const char*                         name;
	std::function< void( Mutations& ) > set;
	std::vector< std::string >          targets;   ///< checks that must fail
	std::set< std::string >             groups;    ///< groups to run
	/// A statistical defect that shows only at some rates (the maximum
	/// filter's benefit depends on where a band's bins fall): it must fail in
	/// at least one cell, and the report says it is a partial catch.
	bool                                anyCell = false;
};

const std::vector< NegativeCase >& negativeCases()
{
	static const std::vector< NegativeCase > c = {
		{ "pitchGlobalMax", []( Mutations& m ) { m.pitchGlobalMax = true; }, { "pitch.octave" }, { "pitch" } },
		{ "pitchSpectralPeak", []( Mutations& m ) { m.pitchSpectralPeak = true; }, { "pitch.mf", "pitch.h2" }, { "pitch" } },
		{ "pitchNoInterp", []( Mutations& m ) { m.pitchNoInterp = true; }, { "pitch.sine" }, { "pitch" } },
		{ "pitchLagDelay=0.1", []( Mutations& m ) { m.pitchLagDelay = 0.1; }, { "pitch.vibrato" }, { "pitch" } },
		{ "partialsNoInterp", []( Mutations& m ) { m.partialsNoInterp = true; }, { "partials.sines" }, { "partials" } },
		{ "onsetNoPrime", []( Mutations& m ) { m.onsetNoPrime = true; }, { "onset.prime.quiet", "fallback.prime.quiet" },
		  { "onset-prime", "fallback" } },
		{ "onsetFullSeed", []( Mutations& m ) { m.onsetFullSeed = true; }, { "onset.prime.hears", "fallback.prime.hears" },
		  { "onset-prime", "fallback" } },
		{ "onsetSnapZeroDt", []( Mutations& m ) { m.onsetSnapZeroDt = true; }, { "fallback.dt0" }, { "fallback" } },
		{ "onsetNoMask", []( Mutations& m ) { m.onsetNoMask = true; }, { "onset.attack.hard" }, { "onset-attack" } },
		{ "onsetNoRisingEdge", []( Mutations& m ) { m.onsetNoRisingEdge = true; }, { "onset.attack.once" }, { "onset-attack" } },
		{ "onsetPlainFlux", []( Mutations& m ) { m.onsetPlainFlux = true; }, { "onset.steady" }, { "onset-prime" }, true },
		{ "bandNoParseval", []( Mutations& m ) { m.bandNoParseval = true; }, { "bands.level", "bands.noise" }, { "bands" } },
		{ "tempoLegacy", []( Mutations& m ) { m.tempoLegacy = true; }, { "tempo.4otf-125" }, { "tempo" } },
		{ "tempoNoOnsetGate", []( Mutations& m ) { m.tempoNoOnsetGate = true; }, { "tempo.noise" }, { "tempo" } },
		{ "thirdOrder2", []( Mutations& m ) { m.thirdOrder2 = true; }, { "thirds.mask" }, { "thirds" } },
		{ "thirdFixedOrder", []( Mutations& m ) { m.thirdFixedOrder = true; }, { "thirds.mask" }, { "thirds" } },
		{ "thirdBase2", []( Mutations& m ) { m.thirdBase2 = true; }, { "thirds.centre" }, { "thirds" } },
		{ "thirdNoPrewarp", []( Mutations& m ) { m.thirdNoPrewarp = true; }, { "thirds.mask", "thirds.edge" }, { "thirds" } },
		{ "stereoSwap", []( Mutations& m ) { m.stereoSwap = true; }, { "stereo.balance", "stereo.pair" }, { "stereo" } },
		{ "stereoNoNormalise", []( Mutations& m ) { m.stereoNoNormalise = true; }, { "stereo.correlation" }, { "stereo" } },
		{ "tapPeekLatest", []( Mutations& m ) { m.tapPeekLatest = true; }, { "tap.contiguous" }, { "tap" } },
		{ "tapNoAntiAlias", []( Mutations& m ) { m.tapNoAntiAlias = true; }, { "tap.alias" }, { "tap" } },
		{ "silentByLevel", []( Mutations& m ) { m.silentByLevel = true; }, { "silent.quiet" }, { "silent" } },
		{ "silentNever", []( Mutations& m ) { m.silentNever = true; }, { "silent.zero" }, { "silent" } },
		{ "configureClears", []( Mutations& m ) { m.configureClears = true; }, { "params.keep" }, { "params" } },
		{ "ringNoValidate", []( Mutations& m ) { m.ringNoValidate = true; }, { "ring.lapped" }, { "ring" } },
		{ "fallbackInventsPitch", []( Mutations& m ) { m.fallbackInventsPitch = true; }, { "fallback.empty" }, { "fallback" } },
	};
	return c;
}

int runNegative( const Options& o )
{
	std::printf( "Negative controls: each mutation must make every check it targets fail in every cell.\n\n" );
	int caught = 0, total = 0;
	for( const NegativeCase& c : negativeCases() )
	{
		if( !o.only.empty() && std::string( c.name ).find( o.only ) == std::string::npos )
			continue;
		Mutations m;
		c.set( m );
		Report rep;
		rep.quiet = true;
		runGroups( rep, o, c.groups, m );
		bool        all = true;
		std::string detail;
		for( const std::string& t : c.targets )
		{
			const int ran = rep.cellsRun( t ), failed = rep.cellsFailed( t );
			all = all && ran > 0 && ( c.anyCell ? failed > 0 : failed == ran );
			detail += format( "  %s %d/%d", t.c_str(), failed, ran );
		}
		++total;
		if( all )
			++caught;
		std::printf( "  %s  %-22s%s%s\n", all ? "caught" : "MISSED", c.name, detail.c_str(),
					 c.anyCell ? "  (statistical: one cell is enough)" : "" );
		if( !all )
			for( const Outcome& out : rep.outcomes() )
				if( std::find( c.targets.begin(), c.targets.end(), out.check ) != c.targets.end() && out.pass )
					std::printf( "          still passes: %s %s %s\n", out.check.c_str(), out.run.c_str(), out.detail.c_str() );
	}
	std::printf( "\nnegative: %d of %d mutations caught\n", caught, total );
	return caught == total ? 0 : 1;
}

//---------------------------------------------------------------------------
// A file through the engine
//---------------------------------------------------------------------------

int runFile( const Options& o )
{
	Signal      s;
	std::string err;
	if( !readWav( o.file, s, err ) )
	{
		std::fprintf( stderr, "estest: %s\n", err.c_str() );
		return 2;
	}
	std::printf( "%s: %d Hz, %d ch, %.2f s, update at %.0f fps\n", o.file.c_str(), s.fs, s.channels, s.seconds(), o.fps );
	Run run{ s.fs, 512, o.fps };
	Driver d( run );
	double settledAt = -1.0, lastBpm = 0.0, worst = 0.0;
	int    onsets[ earshot::kMaxBands ] = {};
	int    readings = 0, inside = 0;
	uint32_t last[ earshot::kMaxBands ] = {};
	d.play( s, [ & ]( const Frame& f, double vt ) {
		for( int b = 0; b < f.bandCount; ++b )
		{
			onsets[ b ] += (int)( f.bands[ b ].onsetCount - last[ b ] );
			last[ b ] = f.bands[ b ].onsetCount;
		}
		if( f.tempo.settled )
		{
			if( settledAt < 0.0 )
				settledAt = vt;
			lastBpm = f.tempo.bpm;
			if( vt >= o.from && o.expectBpm > 0.0 )
			{
				++readings;
				const double e = std::fabs( f.tempo.bpm - o.expectBpm );
				worst = std::max( worst, e );
				inside += e <= 1.0 ? 1 : 0;
			}
		}
	} );
	std::printf( "tempo: %.2f BPM (settled at %.2f s)\n", lastBpm, settledAt );
	std::printf( "onsets per band:" );
	for( int b = 0; b < earshot::kMaxBands && b < d.engine().frame().bandCount; ++b )
		std::printf( " %d", onsets[ b ] );
	std::printf( "\n" );
	if( o.expectBpm > 0.0 )
	{
		const bool ok = readings > 0 && inside == readings;
		std::printf( "%s  file.tempo: every settled frame from %.1f s within 1 BPM of %.1f: %d of %d, worst %.2f BPM off\n",
					 ok ? "ok  " : "FAIL", o.from, o.expectBpm, inside, readings, worst );
		return ok ? 0 : 1;
	}
	return 0;
}

int runWavRoundTrip( const Options& o )
{
	Signal s = silence( 48000, 0.5, 2 );
	addSine( s, 1000.0, 0.5, 0.0, 0 );
	addNoise( s, 0.1, 3, 1 );
	std::string err;
	Signal      back;
	if( !writeWav( o.writeWav, s, err ) || !readWav( o.writeWav, back, err ) )
	{
		std::printf( "FAIL  wav.roundtrip %s\n", err.c_str() );
		return 1;
	}
	const bool ok = back.fs == s.fs && back.channels == s.channels && back.x == s.x;
	std::printf( "%s  wav.roundtrip  float32 stereo written and read back bit-exactly\n", ok ? "ok  " : "FAIL" );
	return ok ? 0 : 1;
}

//---------------------------------------------------------------------------

int runPerf()
{
	Run    run{ 48000, 256, 60.0 };
	Driver d( run );
	Signal s = render( grooves()[ 0 ], 48000, 20.0 );
	addSine( s, 220.0, 0.2 );
	const Signal st = widen( s, 2 );
	const double t0 = nowMs();
	d.play( st, nullptr );
	const double ms = nowMs() - t0;
	std::printf( "perf: 20 s of 48 kHz stereo in %.1f ms (%.0fx real time); update mean %.3f ms, max %.3f ms over %d frames\n",
				 ms, 20000.0 / ms, d.totalUpdateMs / std::max( 1, d.updates ), d.maxUpdateMs, d.updates );
	return 0;
}

int runDevices()
{
	const auto list = earshot::listDevices();
	std::printf( "%zu capture device(s):\n", list.size() );
	for( const auto& dev : list )
		std::printf( "  %s%s  (%d ch)  key %.0f\n", dev.name.c_str(), dev.isDefault ? "  [default]" : "", dev.channels,
					 earshot::deviceKey( dev.name.c_str() ) );
	return 0;
}

int runLive( const Options& o )
{
	earshot::Engine e;
	e.selectDevice( o.device.c_str() );
	const double t0 = nowMs();
	double       next = 0.0;
	while( nowMs() - t0 < o.live * 1000.0 )
	{
		const Frame& f = e.update( ( nowMs() - t0 ) / 1000.0, nullptr, 0 );
		if( nowMs() - t0 >= next )
		{
			earshot::Status st;
			e.status( st );
			std::printf( "%6.2f s  %-8s %-30s %d Hz %d/%d ch  pitch %s %.1f Hz  tempo %.1f  band0 %.1f dB\n",
						 ( nowMs() - t0 ) / 1000.0, earshot::toString( st.state ), st.opened, st.sampleRate, st.channels,
						 st.nativeChannels, f.pitch.voiced ? "on " : "off", f.pitch.hz, f.tempo.bpm,
						 f.bandCount > 0 ? f.bands[ 0 ].db : -200.0f );
			next += 500.0;
		}
		std::this_thread::sleep_for( std::chrono::milliseconds( 16 ) );
	}
	return 0;
}

void usage()
{
	std::printf( "usage: estest [--GROUP ...] [--rate HZ] [--block N] [--negative [--only NAME]] [--with MUTATION] [--perf]\n"
				 "              [--file WAV [--expect-bpm BPM] [--from S] [--fps N]] [--write-wav PATH]\n"
				 "              [--devices] [--live SECONDS [--device NAME]]\n\ngroups:\n" );
	for( const Group& g : groups() )
		std::printf( "  --%-14s %s\n", g.name.c_str(), g.what.c_str() );
}
} // namespace

int main( int argc, char** argv )
{
	Options o;
	for( int i = 1; i < argc; ++i )
	{
		const std::string a    = argv[ i ];
		auto              next = [ & ]() -> const char* { return i + 1 < argc ? argv[ ++i ] : ""; };
		if( a == "--rate" )
			o.rate = std::atoi( next() );
		else if( a == "--block" )
			o.block = std::atoi( next() );
		else if( a == "--negative" )
			o.negative = true;
		else if( a == "--with" )
			o.with = next();
		else if( a == "--only" )
			o.only = next();
		else if( a == "--perf" )
			o.perf = true;
		else if( a == "--devices" )
			o.devices = true;
		else if( a == "--live" )
			o.live = std::atof( next() );
		else if( a == "--device" )
			o.device = next();
		else if( a == "--file" )
			o.file = next();
		else if( a == "--expect-bpm" )
			o.expectBpm = std::atof( next() );
		else if( a == "--fps" )
			o.fps = std::atof( next() );
		else if( a == "--from" )
			o.from = std::atof( next() );
		else if( a == "--write-wav" )
			o.writeWav = next();
		else if( a == "--help" || a == "-h" )
		{
			usage();
			return 0;
		}
		else if( a.rfind( "--", 0 ) == 0 )
		{
			const std::string g = a.substr( 2 );
			const bool known = std::any_of( groups().begin(), groups().end(), [ & ]( const Group& x ) { return x.name == g; } );
			if( !known )
			{
				std::fprintf( stderr, "estest: unknown option %s\n", a.c_str() );
				usage();
				return 2;
			}
			o.groups.insert( g );
		}
	}
	if( o.devices )
		return runDevices();
	if( o.live > 0.0 )
		return runLive( o );
	if( !o.file.empty() )
		return runFile( o );
	if( !o.writeWav.empty() )
		return runWavRoundTrip( o );
	if( o.perf )
		return runPerf();
	if( o.negative )
		return runNegative( o );

	Mutations mut;
	if( !o.with.empty() )
	{
		const auto& nc = negativeCases();
		auto it = std::find_if( nc.begin(), nc.end(), [ & ]( const NegativeCase& c ) { return o.with == c.name; } );
		if( it == nc.end() )
		{
			std::fprintf( stderr, "estest: no mutation called %s\n", o.with.c_str() );
			return 2;
		}
		it->set( mut );
		std::printf( "with mutation %s: the checks it targets should FAIL\n", it->name );
	}

	Report     rep;
	const double t0 = nowMs();
	runGroups( rep, o, o.groups, mut );
	rep.printFigures();
	std::printf( "\nestest: %d passed, %d failed (%.1f s)\n", rep.passed(), rep.failed(), ( nowMs() - t0 ) / 1000.0 );
	if( rep.failed() > 0 )
	{
		std::printf( "failures:\n" );
		for( const Outcome& out : rep.outcomes() )
			if( !out.pass )
				std::printf( "  %s %s: %s\n", out.check.c_str(), out.run.c_str(), out.detail.c_str() );
	}
	return rep.failed() == 0 ? 0 : 1;
}
