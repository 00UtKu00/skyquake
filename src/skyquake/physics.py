# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 SKYQUAKE contributors
"""Forward model: sonic-boom (Mach cone) arrival times from a hypersonic line source.

An object moving in a straight line at constant speed ``v`` through air with sound speed
``c < v`` drags a conical shock behind it. For a station at along-track coordinate ``x``
and perpendicular distance ``r`` from the line (measured from the reference point that
the object passes at time ``t0``), the first shock arrival is

    t = t0 + x / v + r * sqrt(1/c**2 - 1/v**2)

Derivation: a signal emitted at along-track position ``s`` arrives at
``s/v + sqrt((x-s)**2 + r**2)/c``; minimising over ``s`` (Fermat / stationary phase)
gives the expression above. See ``docs/method.md``.

Assumptions (version 0.1): straight line, constant speed, homogeneous sound speed.
Each is a deliberate first-order approximation, listed with its error budget in the docs.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

DEFAULT_SOUND_SPEED = 320.0  # m/s, path-averaged troposphere/stratosphere value


@dataclass(frozen=True)
class Trajectory:
    """Straight-line hypersonic trajectory in a local ENU frame.

    The reference point (east_m, north_m, up_m) is where the object is at time ``t0_s``.
    ``azimuth_deg`` is the heading clockwise from north; ``fpa_deg`` is the descent angle
    below the horizontal (positive = descending).
    """

    east_m: float
    north_m: float
    up_m: float
    azimuth_deg: float
    fpa_deg: float
    speed_mps: float
    t0_s: float

    @property
    def direction(self) -> np.ndarray:
        az, g = np.radians(self.azimuth_deg), np.radians(self.fpa_deg)
        return np.array([np.sin(az) * np.cos(g), np.cos(az) * np.cos(g), -np.sin(g)])

    @property
    def origin(self) -> np.ndarray:
        return np.array([self.east_m, self.north_m, self.up_m])

    def position(self, t_s) -> np.ndarray:
        t = np.asarray(t_s, dtype=float)[..., None]
        return self.origin + self.speed_mps * (t - self.t0_s) * self.direction

    def time_at_altitude(self, up_m: float) -> float:
        """Time at which the line crosses local height ``up_m`` (requires fpa != 0)."""
        vz = -self.speed_mps * np.sin(np.radians(self.fpa_deg))
        if abs(vz) < 1e-9:
            raise ValueError("level trajectory never changes altitude")
        return self.t0_s + (up_m - self.up_m) / vz

    def normalised(self) -> Trajectory:
        return replace(self, azimuth_deg=self.azimuth_deg % 360.0)

    def canonical(self) -> Trajectory:
        """Same line, reference point moved to where it passes the frame origin horizontally.

        This is the point the inversion is parameterised on (cross-track offset + height),
        which stays well defined for level flight, unlike "the point at height h".
        """
        az, g = np.radians(self.azimuth_deg), np.radians(self.fpa_deg)
        heading = np.array([np.sin(az), np.cos(az)])
        s = -(self.origin[:2] @ heading) / np.cos(g)
        o = self.origin + s * self.direction
        return replace(self, east_m=o[0], north_m=o[1], up_m=o[2],
                       t0_s=self.t0_s + s / self.speed_mps).normalised()


def arrival_times(traj: Trajectory, stations_enu, sound_speed: float = DEFAULT_SOUND_SPEED):
    """First Mach-cone arrival time (s, same clock as ``t0_s``) at each station (N, 3)."""
    v, c = traj.speed_mps, sound_speed
    if v <= c:
        raise ValueError(f"object must be supersonic (v={v:.1f} m/s <= c={c:.1f} m/s)")
    d = np.asarray(stations_enu, dtype=float) - traj.origin
    u = traj.direction
    x = d @ u
    r = np.linalg.norm(d - np.outer(x, u), axis=1)
    return traj.t0_s + x / v + r * np.sqrt(1.0 / c**2 - 1.0 / v**2)


def emission_points(traj: Trajectory, stations_enu, sound_speed: float = DEFAULT_SOUND_SPEED):
    """Where on the trajectory each station's first arrival was emitted, (N, 3) in ENU.

    Used to check physical validity: an emission point below ground (or below the height
    where the object went subsonic) means the straight-line model does not apply there.
    """
    v, c = traj.speed_mps, sound_speed
    d = np.asarray(stations_enu, dtype=float) - traj.origin
    u = traj.direction
    x = d @ u
    r = np.linalg.norm(d - np.outer(x, u), axis=1)
    s = x - r * c / np.sqrt(v**2 - c**2)
    return traj.origin + np.outer(s, u)


def mach_angle_deg(speed_mps: float, sound_speed: float = DEFAULT_SOUND_SPEED) -> float:
    return float(np.degrees(np.arcsin(sound_speed / speed_mps)))


# --- decelerating trajectories ---------------------------------------------------------
#
# Re-entering objects slow down strongly below ~70 km (Falcon 9, 2025: 7 -> 3.5 km/s within
# a few hundred km). Model: speed changes linearly with path length,
#     v(s) = v0 + g s          (g = dv/ds, units 1/s; negative = decelerating)
# so the time to reach s is  tau(s) = ln(v(s)/v0) / g  (-> s/v0 as g -> 0). The boom
# arrival is the minimum over emission position s of  t0 + tau(s) + |q - P(s)|/c, taken only
# where the object is still supersonic. g is reported as dv per 100 km of path (km/s).


def emission_time(s, v0: float, g: float) -> np.ndarray:
    """Time after t0 to travel path length ``s`` (m) with v(s) = v0 + g s."""
    s = np.asarray(s, dtype=float)
    if abs(g) < 1e-12:
        return s / v0
    v = v0 + g * s
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(v > 0, np.log(np.maximum(v, 1e-9) / v0) / g, np.inf)


def decelerating_arrival_times(traj: Trajectory, stations_enu, sound_speed: float,
                               dv_kms_per_100km: float, return_emission: bool = False):
    """First boom arrival for a straight line with speed varying linearly along the path."""
    st = np.asarray(stations_enu, dtype=float)
    u, o, v0, c = traj.direction, traj.origin, traj.speed_mps, sound_speed
    g = dv_kms_per_100km * 1e3 / 100e3
    if v0 <= c:
        raise ValueError("object must be supersonic at the reference point")
    d = st - o
    xs = d @ u
    r = np.linalg.norm(d - np.outer(xs, u), axis=1)
    s0 = xs - r * c / np.sqrt(v0**2 - c**2)       # constant-speed stationary point

    def total(s):  # s: (N, M)
        p = o + s[..., None] * u
        dist = np.linalg.norm(st[:, None, :] - p, axis=-1)
        tt = traj.t0_s + emission_time(s, v0, g) + dist / c
        supersonic = (v0 + g * s) > 1.02 * c
        return np.where(supersonic & np.isfinite(tt), tt, 1e6 + np.abs(s) / 1e3)

    best = s0
    rows = np.arange(len(st))
    for half, n in ((400e3, 161), (6e3, 25), (0.6e3, 13)):
        grid = best[:, None] + np.linspace(-half, half, n)[None, :]
        tt = total(grid)
        i = np.clip(np.argmin(tt, axis=1), 1, n - 2)
        y0, y1, y2 = tt[rows, i - 1], tt[rows, i], tt[rows, i + 1]
        step = grid[0, 1] - grid[0, 0]
        denom = y0 - 2 * y1 + y2
        ok = (np.abs(denom) > 1e-12) & (np.maximum(y0, y2) < 1e5)
        shift = np.where(ok, 0.5 * (y0 - y2) / np.where(ok, denom, 1.0), 0.0)
        best = grid[rows, i] + np.clip(shift, -1, 1) * step
    t_arr = total(best[:, None])[:, 0]
    if return_emission:
        return t_arr, o + np.outer(best, u)
    return t_arr
