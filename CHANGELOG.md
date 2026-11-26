# Changelog

## 0.1.0.dev0 (2026-11-26)
- Initial physics core: WGS84 geodesy, Mach-cone forward model, robust multi-start inversion
- Sound speed estimated from data (on by default in `skyquake fit`)
- Calibrated uncertainties, Monte-Carlo corridor, GeoJSON export
- ObsPy-based seismic download and STA/LTA sonic-boom picker (not yet validated on real data)
- Validation event catalogue (`skyquake events`)
- Geographic heading reported at the object's own location (frame-independent)
- 18 tests
