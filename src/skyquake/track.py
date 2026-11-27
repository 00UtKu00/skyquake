# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 SKYQUAKE contributors
"""Turn a fitted trajectory into geographic products: ground track, corridor, GeoJSON."""

from __future__ import annotations

import numpy as np

from skyquake.geo import LocalFrame
from skyquake.invert import FitResult
from skyquake.physics import Trajectory


def geographic_state(traj: Trajectory, frame: LocalFrame) -> dict:
    """Reference point and direction expressed at the object's own location.

    Azimuth and descent angle in a fitting frame are relative to *that frame's* north and
    vertical. Hundreds of km from the frame origin these differ from true geographic
    values by meridian convergence and vertical tilt (~0.25 deg at 0.45 deg of longitude).
    Always report these frame-independent values.
    """
    lat, lon, h = frame.from_enu(traj.origin)
    here = LocalFrame(float(lat), float(lon), 0.0)
    u_ecef = traj.direction @ frame._rotation
    e, n, up = here._rotation @ u_ecef
    return {
        "lat": float(lat),
        "lon": float(lon),
        "height_m": float(h),
        "azimuth_deg": float(np.degrees(np.arctan2(e, n)) % 360.0),
        "fpa_deg": float(-np.degrees(np.arcsin(np.clip(up, -1.0, 1.0)))),
        "speed_mps": traj.speed_mps,
    }


def track_points(traj: Trajectory, frame: LocalFrame, before_km=300.0, after_km=300.0,
                 n=61, min_height_m=0.0):
    """Geodetic points (lat, lon, h) along the line around its reference point.

    Points below ``min_height_m`` (true ellipsoidal height) are dropped.
    """
    t = traj.t0_s + np.linspace(-before_km, after_km, n) * 1e3 / traj.speed_mps
    lat, lon, h = frame.from_enu(traj.position(t))
    pts = np.column_stack([lat, lon, h])
    return pts[pts[:, 2] >= min_height_m]


def corridor(fit: FitResult, frame: LocalFrame, n=1000, seed=0) -> dict:
    """Monte-Carlo uncertainty of where the line passes the network (its reference point).

    Returns the mean point plus 1-sigma along-track, cross-track and height spread. Debris
    that survives continues downrange in ballistic "dark flight", which is not modelled.
    """
    rng = np.random.default_rng(seed)
    o = np.array([tr.canonical().origin for tr in fit.sample(n, rng)
                  if tr.speed_mps > fit.sound_speed])
    mean = o.mean(axis=0)
    az = np.radians(fit.trajectory.azimuth_deg)
    along = np.array([np.sin(az), np.cos(az)])
    cross = np.array([np.cos(az), -np.sin(az)])
    d = o - mean
    lat, lon, h = frame.from_enu(mean)
    return {
        "lat": float(lat),
        "lon": float(lon),
        "height_m": float(h),
        "sigma_along_km": float(np.std(d[:, :2] @ along) / 1e3),
        "sigma_cross_km": float(np.std(d[:, :2] @ cross) / 1e3),
        "sigma_height_km": float(np.std(d[:, 2]) / 1e3),
        "samples": int(len(o)),
    }


def to_geojson(fit: FitResult, frame: LocalFrame, **track_kw) -> dict:
    pts = track_points(fit.trajectory, frame, **track_kw)
    cor = corridor(fit, frame)
    geo = {f"geo_{k}": v for k, v in geographic_state(fit.trajectory, frame).items()}
    return {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": {"type": "LineString",
                             "coordinates": [[p[1], p[0], p[2]] for p in pts]},
                "properties": {"kind": "trajectory", **geo,
                               **{f"frame_{k}": v for k, v in fit.params.items()},
                               **{f"sigma_{k}": v for k, v in fit.sigma.items()},
                               "rms_s": fit.rms_s, "stations": fit.n_used},
            },
            {
                "type": "Feature",
                "geometry": {"type": "Point", "coordinates": [cor["lon"], cor["lat"]]},
                "properties": {"kind": "closest_approach_to_network", **cor},
            },
        ],
    }
