# Data sources

All inputs come from public archives. FALLTRACE never needs private data.

| Source | What | Access | Used in |
|---|---|---|---|
| FDSN web services: EarthScope (formerly IRIS), ORFEUS/EIDA, GFZ, SCEDC, and others | Seismic waveforms and station metadata | `obspy.clients.fdsn`, no key needed | `falltrace pick` |
| CelesTrak / Space-Track | TLE orbital elements (heading prior and cross-check) | Space-Track requires a free account | `falltrace.orbit` |
| Agency re-entry notices (e.g. ESA, POLSA) | Ground truth for validation | Public reports | `data/events.toml` |

Planned: NEXRAD weather radar (debris clouds), GOES GLM (fireball flashes), infrasound
stations, and ERA5 reanalysis (winds and temperature for the atmosphere model).

## Usage notes

- Respect each data centre's usage policy and rate limits. `falltrace pick` downloads one
  station at a time and caps the total at `max_stations`.
- Cite the seismic networks you use. Every network has a DOI listed at
  [fdsn.org/networks](https://www.fdsn.org/networks/).
