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


def track_points(traj: Trajectory, frame: LocalFrame, top_m=80e3, bottom_m=15e3, n=60):
    """Geodetic points (lat, lon, h) along the trajectory between two local heights."""
    t = np.linspace(traj.time_at_altitude(top_m), traj.time_at_altitude(bottom_m), n)
    lat, lon, h = frame.from_enu(traj.position(t))
    return np.column_stack([lat, lon, h])


def corridor(fit: FitResult, frame: LocalFrame, height_m=15e3, n=1000, seed=0) -> dict:
    """Monte-Carlo position uncertainty where the trajectory reaches ``height_m``.

    Returns the mean point and the 1-sigma along-track / cross-track spread (km). This is
    where the hypersonic segment ends; fragments then follow a ballistic "dark flight"
    that is *not* modelled yet, so the true debris field extends downrange of this point.
    """
    rng = np.random.default_rng(seed)
    pts = []
    for tr in fit.sample(n, rng):
        if tr.speed_mps <= fit.sound_speed or tr.fpa_deg <= 0:
            continue
        pts.append(tr.position(tr.time_at_altitude(height_m))[:2])
    pts = np.array(pts)
    mean = pts.mean(axis=0)
    az = np.radians(fit.trajectory.azimuth_deg)
    along = np.array([np.sin(az), np.cos(az)])
    cross = np.array([np.cos(az), -np.sin(az)])
    d = pts - mean
    lat, lon, _ = frame.from_enu(np.array([mean[0], mean[1], height_m]))
    return {
        "lat": float(lat),
        "lon": float(lon),
        "height_m": height_m,
        "sigma_along_km": float(np.std(d @ along) / 1e3),
        "sigma_cross_km": float(np.std(d @ cross) / 1e3),
        "samples": int(len(pts)),
    }


def to_geojson(fit: FitResult, frame: LocalFrame, **corridor_kw) -> dict:
    pts = track_points(fit.trajectory, frame)
    cor = corridor(fit, frame, **corridor_kw)
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
                "properties": {"kind": "end_of_hypersonic_segment", **cor},
            },
        ],
    }
