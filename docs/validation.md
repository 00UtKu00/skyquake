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
| Shenzhou-15 orbital module, S. California, 2024-04-02 | published seismic reconstruction (Fernando & Charalambous, *Science*, doi:10.1126/science.adz4676) | see below | **partial**: geometry and time agree; speed is biased high |
| Falcon 9 S2 (Starlink 11-4), W. Poland, 2025-02-19 | recovered debris locations | — | not yet run |

### Shenzhou-15 orbital module (2024-04-02), first real run

Data: 297 public HHZ seismometers (CI network via SCEDC, NN via EarthScope). The boom shows
up as a clean V-shaped front in `examples/shenzhou15/section_cross.png`. One M1.2
earthquake (USGS ci40703496, 08:38:28 UTC near Corona) also triggers many stations; it is
outside the boom window and does not enter the fit.

Automatic iterative picking: 137 picks, **62 used** after 4σ outlier rejection,
rms 0.84 s, reduced χ² 4.7.

| Quantity | SKYQUAKE | Independent information | Agreement |
|---|---|---|---|
| Time over the network | 08:43:51 UTC (±2 s) | eyewitness/press reports ~08:43 UTC | ✅ |
| Heading | 65.5° ± 0.04° (true) | Santa Barbara → Las Vegas ≈ 63–64° (published path) | ✅ within ~2° (path end points are approximate) |
| Track | passes over Santa Barbara, ~8 km south of the SB–LV line near 118° W | "over Santa Barbara … and Las Vegas" | ✅ qualitatively |
| Height over the network | 56.7 ± 0.6 km | not published in the sources read | — |
| Descent angle | 0.05° ± 0.16° (near level) | — | plausible for a grazing re-entry |
| Speed | **9.9 ± 0.55 km/s** | Mach 25–30 published; **physical maximum ≈ 7.8 km/s** (orbital speed) | ❌ biased high by ≥ 3.7σ |

**What the speed failure tells us.** The speed comes from the small along-track moveout
(about 30 s across 300 km), which is the part most sensitive to the uniform-atmosphere
assumption, winds, and pick quality. The formal ±0.55 km/s therefore understates the real
error for this parameter. Next steps: (1) a stratified atmosphere with winds (ERA5), (2) a
better picker (fewer rejected picks), (3) optionally a physical prior `v ≤ 7.9 km/s`. The
χ² of 4.7 already flags that the model does not explain the data down to the pick noise.

Reproduce: `python examples/shenzhou15/fetch.py`, then `python examples/shenzhou15/pick_and_fit.py`.
Outputs: `picks.csv`, `fit.txt`, `map.png`, `picks_on_section.png`.
