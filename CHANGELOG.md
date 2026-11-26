# Changelog

## 0.1.0.dev0 (2026-11-26)
- Initial physics core: WGS84 geodesy, Mach-cone forward model, robust multi-start inversion
- Sound speed estimated from data (on by default in `skyquake fit`)
- Calibrated uncertainties, Monte-Carlo corridor, GeoJSON export
- ObsPy-based seismic download and STA/LTA sonic-boom picker (not yet validated on real data)
- Validation event catalogue (`skyquake events`)
- Geographic heading reported at the object's own location (frame-independent)
- Closest-approach parameterisation (cross-track offset + height): works for level flight
- First real event: Shenzhou-15 orbital module (examples/shenzhou15)
- Layered atmosphere: exact ray tracing through the US Standard Atmosphere 1976 (`--layered`)
- Jackknife uncertainties (`--jackknife`) and a physical speed bound (`--max-speed`)
- Shenzhou-15 speed diagnostics: layering, deceleration and pick SNR ruled out
- skyquake.pipeline: reusable load / iterative pick-and-fit
- Second real event: Falcon 9 upper stage over Europe (examples/falcon9_poland), compared with optical trajectory and recovered debris
- 29 tests
