# Validation

This page records how accurate SKYQUAKE is, and where it fails. Results are reported
whether they look good or bad.

## Level 1: the physics is implemented correctly (verification)

| Check | Result |
|---|---|
| Closed-form arrival vs brute-force minimisation (25 stations) | agree to < 0.1 ms |
| Geodetic → ECEF → geodetic round trip | < 1e-9 deg, < 0.1 mm |
| Tangent-plane curvature drop at 300 km | 7.0 km (expected d²/2R ≈ 7.06 km) |

## Level 2: the inversion recovers known answers (synthetic)

Network: 60 stations over about 400 × 330 km. Pick noise: 0.3 s (1σ).

| Test | Result |
|---|---|
| Noise-free data | all 6 parameters exact to 1e-3 |
| Blind heading search (no prior) | heading within 2° (the demo case: 0.01°) |
| 3 corrupted picks (+5 to +20 s) | all rejected, rms 0.25 s |
| **Uncertainty calibration**, 40 realisations | ≥ 85 % of errors within 2σ for every parameter (nominal 95 %) |
| Geographic heading independent of fitting frame (frame shifted ~1 deg) | agrees within 4 sigma; frame-relative azimuth would be off by ~0.5 deg |
| Sound speed fitted from a wrong start (300 vs 320 m/s) | recovered within 3σ, along-track bias removed |

### Measured bias: wrong sound speed (synthetic, seed 1)

Measured with v0.1's parameterisation (reference point = where the line is at 40 km).
Since the switch to the closest-approach parameterisation the bias shows up as height and
cross-track error instead; the lesson is unchanged.

| Assumed c | Along-track shift of the 40 km point | Heading | rms |
|---|---|---|---|
| −5 % (304 m/s) | ~147 km | −0.34° | 0.57 s |
| exact (320 m/s) | 0.8 km | +0.01° | 0.25 s |
| +5 % (336 m/s) | ~86 km | +0.25° | 0.48 s |
| **fitted** | 0.8 km (σ 0.9 km) | +0.01° | 0.25 s |

**Lesson:** the heading is robust, but the position along the track is degenerate with the
sound speed. Estimating `c` from the data (the default in `skyquake fit`) removes the bias
at the cost of roughly doubling the position uncertainty. The synthetic test uses a uniform
atmosphere, so in real air some model error remains. That error is not yet quantified.

## Level 3: real events

| Event | Ground truth | SKYQUAKE result | Verdict |
|---|---|---|---|
| Shenzhou-15 orbital module, S. California, 2024-04-02 | published seismic reconstruction (Fernando & Charalambous, *Science*, doi:10.1126/science.adz4676) | see below | **partial**: geometry and time agree; speed poorly constrained (best value above orbital speed) |
| Falcon 9 S2 (Starlink 11-4), Europe, 2025-02-19 | optical trajectory (Vierinen et al. 2026, arXiv:2605.29124) + recovered debris in W. Poland | see below | **good**: speed and time match optics; track extrapolated 650 km passes 9–42 km from all debris sites (systematic northward offset) |

### Shenzhou-15 orbital module (2024-04-02), first real run

Data: 297 public HHZ seismometers (CI network via SCEDC, NN via EarthScope). The boom shows
up as a clean V-shaped front in `examples/shenzhou15/section_cross.png`. One M1.2
earthquake (USGS ci40703496, 08:38:28 UTC near Corona) also triggers many stations; it is
outside the boom window and does not enter the fit.

Automatic iterative picking: 137 picks, **62 used** after 4σ outlier rejection,
rms 0.84 s, reduced χ² 4.7.

Uncertainties below are **jackknife** (10 station groups), which also capture unmodelled
systematics. The formal covariance errors were 2–6× smaller and are in `fit.txt`.

| Quantity | SKYQUAKE (jackknife 1σ) | Independent information | Agreement |
|---|---|---|---|
| Time over the network | 08:43:51 UTC (±6 s) | eyewitness/press reports ~08:43 UTC | ✅ |
| Heading | 65.5° ± 0.24° (true) | Santa Barbara → Las Vegas ≈ 63–64° (published path) | ✅ within ~2° (path end points are approximate) |
| Track | passes over Santa Barbara, ~8 km south of the SB–LV line near 118° W | "over Santa Barbara … and Las Vegas" | ✅ qualitatively |
| Height over the network | 56.7 ± 1.6 km | not published in the sources read | — |
| Descent angle | 0.05° ± 0.65° (near level) | — | plausible for a grazing re-entry |
| Speed | **9.9 ± 2.4 km/s** | Mach 25–30 published (≈ 8.0–9.6 km/s at c ≈ 320 m/s); physical maximum ≈ 7.8 km/s | ⚠️ consistent with the publication; best value above orbital speed |

**Correction.** The first version of this page called the speed "biased high by ≥ 3.7σ".
That used the formal ±0.55 km/s, which turned out to be far too optimistic. With jackknife
errors the best value is about 1σ (all picks) to 1.8σ (inliers only, ±1.0 km/s) above
orbital speed: worth investigating, but not a demonstrated bias.

### Why is speed so weakly constrained? (diagnostics, `examples/shenzhou15/diagnose.py`)

- **The data barely separate speeds.** Refitting with speed fixed: rms misfit 0.859 s at
  7.8 km/s vs 0.828 s at the best fit (9.6 km/s), a difference of 0.03 s.
- **Geometry couples speed to distance.** Along-track position and distance from the track
  are correlated (r = +0.47) in this network, so any error in the distance term leaks into
  the speed.

### Hypotheses tested for the high best-fit speed

| Hypothesis | Test | Result |
|---|---|---|
| Uniform-atmosphere approximation | New layered atmosphere (US Standard 1976, exact ray tracing) on the real picks; synthetic layered data fitted with the uniform model | ❌ not the cause: real data still give 9.6–9.8 km/s; on synthetic data the uniform model is biased by only +0.0 to +0.3 km/s |
| Deceleration (drag at ~57 km, roughly 5 g) | Synthetic decelerating trajectories fitted with the constant-speed model | ❌ not the cause: deceleration biases the fitted speed **down** (5 g: 7.6 → 7.1 km/s) |
| Weak, late or early picks | Residual vs SNR; fit with SNR ≥ 10 picks only | ❌ not the cause: weak correlation (r = 0.18); high-SNR picks alone give 9.55 km/s |
| Winds (jet stream at 10–15 km, stratospheric winds) | — | open: needs ERA5 wind profiles |
| Fragmentation (several boom sources) | — | open: the module broke up gradually (published) |

Layered-atmosphere result on the same picks: speed 9.76 km/s, height 57.6 km, azimuth 65.7°.
Heading and height are robust across all models; speed is not.

### Falcon 9 upper stage, Starlink 11-4 (2025-02-19), second real run

**Data.** 575 public stations (EIDA routing: NL, BE, GR, BW, … plus 31 more from BGR).
Most are the dense Groningen borehole network in the northern Netherlands, which lies
almost directly under the track. Polish (PL), Danish and Raspberry Shake stations could not
be retrieved this time (GEOFON outage; Raspberry Shake refused bulk requests), so there are
no stations north of the track (North Sea) and none near the debris field.

**Boom.** Clear half-V front: apex ~03:49:15 UTC near the track, arriving later with
distance to the south (`examples/falcon9_poland/image.png`).

**Fit.** 219 picks, **213 used**, rms 0.88 s, reduced χ² 1.67, a much cleaner result
than Shenzhou-15.

| Quantity | SKYQUAKE (jackknife 1σ) | Independent reference | Agreement |
|---|---|---|---|
| Time | 65 km over Groningen at 03:45:56 UTC (±3 s) | optics: fragment families cross 60 km at 03:45:51 and 03:46:22 UTC | ✅ |
| Speed | **6.25 ± 0.51 km/s** | optics: 6.9–7.1 km/s at 70 km, decelerating to 2.6–4.1 km/s at 40 km | ✅ consistent (boom emitted near 65 km) |
| Height | 65.2 ± 0.7 km | optics: fragments between 85 and 36 km, peak detections near 60 km | ✅ |
| Heading | 94.8° ± 0.4° (true) | optics/press: east-south-east over N. Netherlands and N. Germany to W. Poland | ✅ qualitatively |
| Track 650 km downrange | see table below | recovered debris (town centres) | ⚠️ all sites north of the track |

Debris sites vs the SKYQUAKE ground track extrapolated ~650 km east
(`compare_debris.py`; covariance inflated ×3.9 to match the jackknife; town-centre
coordinates are approximate to a few km):

| Debris site | Distance from track (+ = north) | 1σ |
|---|---|---|
| Komorniki (tank, near Poznań) | +9.4 km | 3.4 km |
| Wiry | +9.1 km | 3.4 km |
| Międzychód | +24.7 km | 3.0 km |
| Śliwno | +29.8 km | 3.2 km |
| Sędziny (Szamotuły) | +41.5 km | 3.2 km |

**Reading.** Extrapolated 650 km beyond the network, the seismic track passes within
~9 km of the main recovered tank and within 42 km of every debris site. For comparison, the
pre-event forecast of an uncontrolled re-entry is uncertain by hours, i.e. thousands of km.
But every site lies on the north side, beyond the statistical error. That points to a
systematic effect not yet modelled. Candidates: the real path is not straight over 650 km
(gravity, drag, Earth rotation), a ~1–3° heading bias from one-sided station geometry, and
fragment dispersion (the optical study finds two fragment families, the recovered pieces
most likely from the F2 family).

#### Direct comparison with the optical fragment positions

`compare_optical.py` compares the seismic line with the 1114 camera-triangulated positions
(30 fragments, 03:44:00–03:47:25 UTC) published by Vierinen et al. (data not
redistributed; the script says how to get them). Figure: `optical_vs_seismic.png`.

| Quantity, within ±150 km along-track of the network | Optical minus SKYQUAKE |
|---|---|
| Cross-track position | **−0.4 km** (16–84 %: −2.4 … +1.6 km) ✅ |
| Height | +5.1 km (optics higher) |
| Time at the same place | **−35 s** (optics earlier) ⚠️ |
| Speed | 6.99 km/s optical vs 5.9–6.25 km/s seismic (−11 … −16 %) |

**The cross-track offset grows linearly along the track**: about +15 km at 580 km uprange,
zero at the network, about −13 km at 700 km downrange. That is a **heading bias of
~1.2°**, not path curvature (which would be quadratic). It explains most of the
northward offset of the debris sites; fragment dispersion explains the rest. Most stations
used (Groningen) lie on a north–south line, almost perpendicular to the track, so the
along-track baseline that sets the heading is short. The jackknife ±0.36° does not capture
this geometric weakness: the bias is ~3.5 jackknife σ.

**Timing and height are biased together** (35 s late, 5 km low). This is the known
trade-off between reference time, height and sound speed. The fitted path-averaged sound
speed (322 m/s) is probably too high; a real atmosphere is slower along most of the
path. Next steps: winds, a calibrated atmosphere, and stations along the track (Germany,
Poland).

#### Heading fix: more stations along the track, declustering, spatial jackknife

Three changes, then the same comparison with the camera data (`heading_study.py`):

1. **21 more stations** via the IRIS federator (DK, GE, CZ, PL, temporary networks …),
   8 of them east of 8° E. This gives the network an along-track baseline.
2. **Declustering** (`decluster_sigma`): pick errors are multiplied by √(stations within
   20 km), so ~200 Groningen stations no longer dominate as if they were independent.
   Outliers are rejected with the true errors first.
3. **Spatial-block jackknife** (`jackknife(blocks="spatial")`): drops whole regions (k-means
   clusters) instead of random stations, which exposes geometric weakness.

| Variant | Heading error vs cameras | Spatial-jackknife 1σ | error / σ |
|---|---|---|---|
| Before: Netherlands only, plain | +1.24° | 0.35° | 3.6 ❌ |
| Before: Netherlands only, declustered | +0.68° | 1.00° | 0.7 |
| **Now: + eastern stations, plain** | **+0.12°** | 0.94° | 0.1 ✅ |
| Now: + eastern stations, declustered | +0.17° | 0.48° | 0.4 ✅ |
| Now: Groningen removed (73 picks) | +0.14° | 0.84° | 0.2 ✅ |

Against the 1114 camera positions over the **whole 1300 km path**, the seismic line is now
offset by a median **0.0 km** (16–84 %: −1.4 … +1.3 km); before, by 7 km and up to 15 km.

Debris sites vs the new track (650 km downrange; spatial-block jackknife of the whole
fit-and-extrapolate chain):

| Debris site | Distance from track (+ = north) | 1σ |
|---|---|---|
| Komorniki (tank) | −3.2 km | 10.9 km |
| Wiry | −3.6 km | 11.0 km |
| Międzychód | +13.4 km | 9.8 km |
| Śliwno | +17.9 km | 10.3 km |
| Sędziny (Szamotuły) | +29.5 km | 10.4 km |

The main tank is now consistent with zero offset. The other sites lie 1.3–2.8σ north,
which is consistent with the ~33 km cross-track spread of the strewn field (several
fragment families).

**What got worse: speed and timing.** With stations spanning ~600 km along-track, the
constant-speed model returns 11.1 km/s (cameras: 7.06 km/s near the network) and the time
is ~41 s late. The cameras show the stage slowing from ~7 to ~3.5 km/s across this stretch.
A straight constant-speed line cannot represent that, and the longer baseline now exposes
it. **Next model step: a decelerating trajectory.**

**Layered atmosphere did worse here:** speed 3.5 km/s, inconsistent with the optical
7 km/s at 70 km. Many southern stations lie 200–300 km from the track, beyond the
geometric shadow-zone limit, where our layered travel times are only approximate. The
homogeneous model with a fitted sound speed agrees better with the optical reference for
this event.

Reproduce: `python examples/falcon9_poland/fetch.py`, `fetch_more.py`, then
`pick_and_fit.py` and `compare_debris.py`.

Reproduce: `python examples/shenzhou15/fetch.py`, then `python examples/shenzhou15/pick_and_fit.py`.
Outputs: `picks.csv`, `fit.txt`, `map.png`, `picks_on_section.png`.
