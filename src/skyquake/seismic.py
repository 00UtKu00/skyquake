# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 SKYQUAKE contributors
"""Fetch public seismic data (FDSN web services) and pick sonic-boom arrivals.

Requires the optional dependency ``obspy`` (``pip install skyquake[seismic]``).

Picking is the weakest link of the chain, so it is deliberately conservative:
for each station we take the onset of the *strongest* STA/LTA trigger in the window,
report its signal-to-noise ratio, and assign a pick uncertainty that grows as SNR falls.
Every pick should still be reviewed by eye before a result is published.
"""

from __future__ import annotations

from datetime import datetime

import numpy as np


def _obspy():
    try:
        import obspy  # noqa: F401
        from obspy import UTCDateTime
        from obspy.clients.fdsn import Client
        from obspy.signal.trigger import classic_sta_lta, trigger_onset
    except ImportError as e:  # pragma: no cover
        raise ImportError("seismic features need obspy: pip install 'skyquake[seismic]'") from e
    return UTCDateTime, Client, classic_sta_lta, trigger_onset


def pick_sonic_booms(
    lat: float,
    lon: float,
    radius_km: float,
    start: datetime,
    end: datetime,
    *,
    provider: str = "IRIS",
    channels: str = "HHZ,BHZ,EHZ,HNZ",
    freqmin: float = 2.0,
    freqmax: float = 15.0,
    sta_s: float = 0.5,
    lta_s: float = 15.0,
    trigger_on: float = 4.0,
    min_snr: float = 5.0,
    max_stations: int = 300,
    log=print,
) -> list[dict]:
    """Return one pick per station (dicts ready for :func:`skyquake.io.write_picks`)."""
    UTCDateTime, Client, classic_sta_lta, trigger_onset = _obspy()
    client = Client(provider)
    t1, t2 = UTCDateTime(start), UTCDateTime(end)
    inv = client.get_stations(
        latitude=lat, longitude=lon, maxradius=radius_km / 111.195,
        starttime=t1, endtime=t2, channel=channels, level="channel",
    )
    picks: list[dict] = []
    seen = 0
    for net in inv:
        for sta in net:
            if seen >= max_stations:
                break
            seen += 1
            code = f"{net.code}.{sta.code}"
            try:
                st = client.get_waveforms(net.code, sta.code, "*", channels,
                                          t1 - lta_s, t2)
            except Exception as e:  # noqa: BLE001 - network services fail in many ways
                log(f"  {code}: no data ({type(e).__name__})")
                continue
            st.merge(fill_value=0)
            tr = st[0]
            df = tr.stats.sampling_rate
            if df < 2.5 * freqmax:
                log(f"  {code}: sampling rate {df} Hz too low, skipped")
                continue
            tr.detrend("demean")
            tr.taper(0.02)
            tr.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=True)
            cft = classic_sta_lta(tr.data, int(sta_s * df), int(lta_s * df))
            on = trigger_onset(cft, trigger_on, 1.0)
            if len(on) == 0:
                continue
            peaks = [cft[a:b + 1].max() for a, b in on]
            a, b = on[int(np.argmax(peaks))]
            noise = np.std(tr.data[max(0, a - int(lta_s * df)):a]) or 1e-12
            snr = float(np.abs(tr.data[a:b + 1]).max() / noise)
            if snr < min_snr:
                continue
            onset = tr.stats.starttime + a / df
            picks.append({
                "station": code,
                "lat": sta.latitude,
                "lon": sta.longitude,
                "elev_m": sta.elevation,
                "arrival_utc": onset.datetime.isoformat() + "Z",
                # 0.1 s floor (onset sharpness), growing as SNR drops
                "sigma_s": round(max(0.1, 3.0 / snr), 3),
                "snr": round(snr, 1),
            })
            log(f"  {code}: pick {onset}  snr={snr:.1f}")
    return picks
