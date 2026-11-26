# Roadmap

## v0.1: physics core ✅
Mach-cone model, robust inversion, calibrated uncertainties, sound-speed estimation,
GeoJSON output, seismic picker.

## v0.2: first real event
- [ ] Run Shenzhou-15 (2024-04-02, S. California) end to end and compare with the published reconstruction
- [ ] Interactive pick-review plot (spectrogram, STA/LTA and pick per station)
- [ ] Publish the result, good or bad, in `docs/validation.md`

## v0.3: second event and a better atmosphere
- [ ] Falcon 9 Starlink 11-4 (2025-02-19, Poland) against the recovered debris sites
- [ ] Ray tracing through ERA5 temperature and wind profiles
- [ ] Decelerating trajectory model

## v0.4: more sensors
- [ ] NEXRAD radar debris detection
- [ ] GOES GLM fireball flashes
- [ ] Infrasound stations
- [ ] Joint inversion of all sensors

## v0.5: dark flight and debris footprint
- [ ] Ballistic fragment propagation with winds, as a footprint probability map

## v1.0: automatic
- [ ] Watch re-entry predictions, fetch data after each event, publish reports automatically
- [ ] Browser UI
