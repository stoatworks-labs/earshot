# earshot

> **AI-assisted project.** This codebase was created with [Claude](https://claude.com/claude-code)
> (Anthropic), directed and reviewed by a human author. What its offline
> harness (`estest`) proves, on synthetic audio, at 44.1 and 48 kHz and at two
> block sizes each: pitch within **0.1 cent** on sines from 40 Hz to 2 kHz and
> within **2.7 cents** at 10 dB SNR, with the confidence equal to S/(S+N) within
> 0.0004; onsets that hear a kick **8-13 ms** after the stream starts (no
> deaf first second) and fire **nothing** on held chords, vibrato or a 70 dB
> swell; a third-octave bank whose every band sits **inside the IEC 61260
> class 1 mask** (0.89 dB to spare on the skirts); a tempo within **0.09 BPM**
> on five grooves by 4 s, including the four-on-the-floor at 125 that pattern
> v0.1.0 read as 62, and within 0.22 BPM on pattern's real release-video
> soundtrack; a waveform tap that is the input sample for sample; and **no
> allocation** on the render or audio thread. Every one of its 65 checks
> has a negative control, a planted defect it must fail, and all 51 do (one
> only at one sample rate, a statistical effect; one, for tempo, fails four
> of the five grooves but never the shuffle); fifteen further planted bugs
> are all caught.
> **It has never opened a live audio device**, from a host or anywhere else,
> and no released plugin ships it yet. The same harness passes in CI on
> Windows x64 under MSVC and on Intel macOS, as on Apple Silicon.

A shared live-audio front end for the Stoatworks FFGL plugins. One C++17
library, pulled into a plugin with one `add_subdirectory`, that opens its own
capture device and hands the plugin, once a frame, what the fleet's
audio-reactive plugins keep needing and Resolume's FFT cannot give them.

## The one idea

An FFGL plugin gets one audio input: `FF_USAGE_FFT`, 64 bins of a mono sum,
no phase, no channels, no time domain, and a bin law nobody has measured.
That is enough for a bar that bounces. It is not enough for a tuner, a
third-octave analyser whose 25 Hz band is 6 Hz wide, a stereo scope, a
tracker that hears the first beat, or a loudspeaker model that needs the
waveform as a stream. Each plugin that wanted one of those reinvented its
own capture. earshot is that work done once, with the fleet's lessons kept:

- **Devices by name, never by index**, and a missing device is Missing, not
  silently swapped for the laptop microphone (spasis).
- **Silent is a state**: a device that opens and returns bit-exact zeros
  (no microphone permission, an unpatched Dante receiver) is shown as dead,
  not papered over (spasis).
- **Opening a device never blocks the render thread** (spasis's measured
  bug, fixed here with a worker thread).
- **Onsets are primed on the first frame** (pattern's deaf-for-1.5-s bug),
  and a hard attack's click lights only the bands it really reaches
  (pattern's broadband lesson).
- **A tempo at the right metre** (pattern v0.1.0 read 62 for 125).
- With no device chosen, it falls back to the host's FFT for band envelopes
  and onsets, and leaves everything else **empty** rather than guessed.

## What a plugin gets

Every frame, from `engine.update( hostSeconds, fftBins, 64 )`:

| | What | Latency |
| --- | --- | --- |
| Pitch | MPM, 40 Hz-2 kHz, a confidence that is S/(S+N), octave-safe | ~38 ms |
| Partials | the strongest sines, QIFFT-refined, labelled by harmonic | 85 ms |
| Onsets | per band, primed on frame one, masked against clicks | one 21 ms window |
| Band envelopes | Parseval-correct levels, attack/release | one hop |
| Thirds | ISO 266 / IEC 61260 base-ten bank, 20 Hz-20 kHz, Fast weighting, peak hold | filter delay + 125 ms |
| Tempo | pattern v0.1.1's three-register detector, gated on onsets | settles in ~4 s |
| Stereo | correlation, width, balance, levels of a chosen channel pair | 21 ms |
| Waveform tap | the pair resampled to a fixed rate, as a stream | 1.6 ms |

The calls, the thread rules, how to save the device choice in a composition,
and what to show for each capture state are in
[docs/INTEGRATING.md](docs/INTEGRATING.md). The API is one header,
[include/earshot/earshot.h](include/earshot/earshot.h), documented in place.

```cmake
add_subdirectory( external/earshot )
target_link_libraries( MyPlugin PRIVATE earshot::earshot )
```

## Used by

Each of these vendors earshot as an unpatched snapshot of one commit in
`external/earshot/`, with a `SNAPSHOT` file naming it
([docs/INTEGRATING.md](docs/INTEGRATING.md), "Taking a snapshot"). All of
them took `8124191`.

- [needle](https://github.com/stoatworks-labs/needle): its VFD spectrum
  display, a third-octave bank of earshot's. On needle's `main` since
  2026-11-16, not yet in a release.
- woofer, eidophone, rubens, marbling and snowglobe: plugins built on it and
  about to be released, each to its own public repo.

There are no tags and no releases here: a plugin takes a commit, so `main`'s
history is never rewritten.

## Status

**2026-11-18: public, a library with a harness, vendored by six plugins and
shipped in none yet.**

Verified by measurement on one Apple Silicon Mac (macOS 26.4), every check in
four cells (44.1 and 48 kHz, blocks of 128 and 1031 samples, 60 fps), with
`tools/verify.sh` green (universal build, estest native and under Rosetta,
and again built without fused multiply-add):

- **Pitch.** Sines 40 Hz-2 kHz: worst 0.10 cents (limit 1, from MPM's
  parabola bias). Missing fundamental (harmonics 2-6): 0.02 cents. Residue
  tones (harmonics 9-12, 15-19): 0.01 cents. Sub-octave and strong second
  harmonic traps: 2.6 cents (limit 5). A sine in white noise at 10 dB SNR:
  2.7 cents, confidence within 0.0004 of S/(S+N). Vibrato: always inside the
  window's instantaneous range.
- **Partials.** Frequency within 0.0159 bin (Hann QIFFT theory: 0.0160).
- **Onsets.** A kick 0.25 s after the stream starts on a playing mix: heard
  in 8-13 ms; the mix itself fires nothing. 1 ms attacks on 90 Hz, 300 Hz,
  1.5 kHz, 8 kHz: each heard within 6.3 ms, nothing beyond the click's
  physical reach. 30 and 80 ms rises: one onset each. Claps: every band they
  reach. 12 s each of noisy pads, a 41 Hz bass, sung vibrato, an organ
  chord: no onsets. White noise: one per 12 s in the lowest band (a known
  limit, AGENTS.md).
- **Thirds.** Mid-bands exact, ISO labels, 30 bands at 44.1 kHz and 31 at
  48. Every band's response inside IEC 61260:1995 class 1 on 4,000 points
  per band: 0.89 dB to spare on the skirts; edge attenuation 3.0099-3.0103 dB
  (class 1 allows 2.0-5.0). The running filters match that response end to
  end within 0.0004 dB. The Fast level decays at 34.744 dB/s (the 125 ms
  time constant's 34.74), and the peak hold holds for its second and falls
  at its 20 dB/s.
- **Tempo.** Five synthetic grooves (four-on-the-floor 125, kick 128, break
  96, drum and bass 174, triplet shuffle 110): settled at 4.0 s, every frame
  from 6 s within 0.09 BPM; the same at 24, 30 and 60 fps. White noise and a
  held chord never settle. pattern's release-video soundtrack: 124.96 BPM,
  every frame from 6 s within 0.22 BPM, at 30 and 60 fps.
- **Stereo and tap.** Correlation, width and balance exact to 1e-5 on tones;
  the tap equal to the analytic signal within 2.3e-5 of full scale, 87 dB of
  alias rejection, flat to 0.0004 dB at 2 kHz.
- **Capture logic**, against a scripted device: a 300 ms open never blocks
  `update()` (worst 0.7 ms unloaded); a missing device is retried and opens at 2.1 s,
  the default never asked for; unplugged, Missing in 7 ms and back as a new
  stream; Silent at exactly 3.000 s of zeros, never on -120 dBFS noise;
  210 million frames through the ring from a racing thread, no torn read; a
  parameter change keeps the audio history and the settled tempo.
- **No allocation** in `feed()`, `update()`, `readTap()` or the device
  callback.

And in CI (`.github/workflows/ci.yml`), on every push and pull request, first
green on 2026-11-18 ([run 44225830323](https://github.com/stoatworks-labs/earshot/actions/runs/44225830323)):

- **Windows x64, MSVC** (windows-2022): the library at `/W4 /permissive-`
  and the harness build with no warnings; estest 238 of 238 in every cell,
  the negative controls 51 of 51 with coverage 65 of 65, the WAV round trip,
  and the consumer built as a DLL, loaded by a host and pitching 440 Hz.
- **macOS x86_64** on an Intel runner (macos-15-intel), natively rather than
  under Rosetta: estest 238 of 238 and the consumer.
- **macOS arm64** (macos-15): the universal build with both slices, estest,
  the negative controls, `EARSHOT_CAPTURE=OFF`, and the consumer exporting
  its one symbol and nothing of earshot's.

The runners open no audio device. The mutation test, the build without fused
multiply-add and the soundtrack tempo run only in `tools/verify.sh`.

Not verified:

- **It has never opened a live device.** The miniaudio backend has listed
  this Mac's 11 capture devices and nothing more; `estest --live` exists and
  has not been run.
- **It has never run inside Resolume.** Six plugins vendor it, and none of
  them had been loaded into Resolume with earshot in it as of 2026-11-18
  (each plugin's README says what its own release found).
- **Capture on Windows.** MSVC builds the WASAPI backend, but no Windows
  device has ever been listed or opened. On the CI runner `estest --devices`
  lists none, and an empty list does not say whether the backend started
  (enumeration is silent about failing).
- Everything is synthetic audio except one real track (the tempo).

## Building and testing

```
cmake -S . -B build -DCMAKE_OSX_ARCHITECTURES=arm64
cmake --build build -j4
build/tools/estest/estest          # every check, every cell, about a minute
tools/verify.sh                    # everything, about ten minutes
```

On Windows, with Visual Studio's CMake:

```
cmake -S . -B build -A x64
cmake --build build --config Release
build\tools\estest\Release\estest.exe
```

`CLAUDE.md` is the command reference; `AGENTS.md` is why everything is the
way it is, and the traps.

## Licence

MIT, copyright Stoatworks Labs. miniaudio (vendored) is public domain or MIT-0.
spasis and pattern, which much of this comes from, are MIT. See
[ATTRIBUTIONS.md](ATTRIBUTIONS.md).
