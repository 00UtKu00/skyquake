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
- Declustering (decluster_sigma) and spatial-block jackknife for clustered networks
- Falcon 9: +21 stations via the IRIS federator; heading error vs cameras 1.24 -> 0.12 deg
- Decelerating trajectory model (fit_deceleration, --decelerating), exact vs brute force
- Falcon 9: deceleration does not remove the speed/timing bias; diagnosed as the single-sound-speed approximation (timing needs ~295 m/s, moveout ~324 m/s)
- Robust layered fit: grid search over height, time and profile scale before optimising
- Metadata guard: impossible station elevations (e.g. 5R placeholder 123456 m) replaced by 0 m; Falcon 9 timing bias vs cameras 40 s -> 23 s, height 6.2 -> 1-2 km
- Speed biased high in both events with every model; fragment-cloud hypothesis documented
- Fragment-cloud test with the 30 camera-tracked fragments: hypothesis refuted; real fragments through the standard layered atmosphere reproduce observed arrivals to 1.24 s plus a constant 23.5 s delay; the single-line speed is a model artefact
- 41 tests
