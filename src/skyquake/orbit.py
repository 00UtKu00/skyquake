"""Catalogue cross-check: heading and sub-satellite point from a TLE.

Requires the optional dependency ``skyfield`` (``pip install skyquake[orbit]``).

Caveat: TLEs in the last hours before re-entry can be off by hundreds of kilometres along
track. We use them only as a *prior on heading* and as an independent sanity check, never
as ground truth.
"""

from __future__ import annotations

from datetime import UTC, datetime, timedelta

import numpy as np


def ground_heading(tle1: str, tle2: str, when: datetime) -> dict:
    try:
        from skyfield.api import EarthSatellite, load, wgs84
    except ImportError as e:  # pragma: no cover
        raise ImportError("orbit features need skyfield: pip install 'skyquake[orbit]'") from e
    ts = load.timescale()
    sat = EarthSatellite(tle1, tle2, ts=ts)
    when = when.astimezone(UTC)

    def sub(t):
        p = wgs84.subpoint_of(sat.at(ts.from_datetime(t)))
        return p.latitude.degrees, p.longitude.degrees

    lat1, lon1 = sub(when)
    lat2, lon2 = sub(when + timedelta(seconds=1))
    phi1, phi2 = np.radians(lat1), np.radians(lat2)
    dlon = np.radians(lon2 - lon1)
    y = np.sin(dlon) * np.cos(phi2)
    x = np.cos(phi1) * np.sin(phi2) - np.sin(phi1) * np.cos(phi2) * np.cos(dlon)
    az = np.degrees(np.arctan2(y, x))
    return {"lat": float(lat1), "lon": float(lon1), "azimuth_deg": float(az % 360.0)}
