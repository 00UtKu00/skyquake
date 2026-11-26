# Validation

This page records how accurate FALLTRACE is, and where it fails. Results are reported
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

| Assumed c | Along-track shift of the 40 km point | Heading | rms |
|---|---|---|---|
| −5 % (304 m/s) | ~147 km | −0.34° | 0.57 s |
| exact (320 m/s) | 0.8 km | +0.01° | 0.25 s |
| +5 % (336 m/s) | ~86 km | +0.25° | 0.48 s |
| **fitted** | 0.8 km (σ 0.9 km) | +0.01° | 0.25 s |

**Lesson:** the heading is robust, but the position along the track is degenerate with the
sound speed. Estimating `c` from the data (the default in `falltrace fit`) removes the bias
at the cost of roughly doubling the position uncertainty. The synthetic test uses a uniform
atmosphere, so in real air some model error remains. That error is not yet quantified.

## Level 3: real events (not started)

| Event | Ground truth | FALLTRACE result | Verdict |
|---|---|---|---|
| Shenzhou-15 orbital module, S. California, 2024-04-02 | published seismic reconstruction | — | not yet run |
| Falcon 9 S2 (Starlink 11-4), W. Poland, 2025-02-19 | recovered debris locations | — | not yet run |

Until this table has rows with results, FALLTRACE is **not validated on real data**.
