"""Layered (1-D) atmosphere: acoustic travel times by exact ray tracing.

The homogeneous model uses one path-averaged sound speed. The real atmosphere is not
uniform: sound speed falls from ~340 m/s at the ground to ~295 m/s in the tropopause,
rises to ~330 m/s at the stratopause (~50 km) and falls again above. Rays bend, and the
travel-time curve T(distance) has a different shape from the straight-ray one. Because
station distance from the track is correlated with along-track position in real networks,
that shape error leaks into the fitted speed (see docs/validation.md, Shenzhou-15).

Method
------
The profile c(z) is piecewise linear between nodes. In a linear-gradient layer a ray is a
circular arc, so horizontal distance and travel time across each layer are analytic for a
given ray parameter p = sin(theta)/c:

    dx = (cos_a - cos_b) / (p g)
    dt = (1/g) ln( c_b (1 + cos_a) / (c_a (1 + cos_b)) ),   cos = sqrt(1 - p^2 c^2)

Summing layers gives X(p), T(p) for a source at each node height and a receiver at the
ground; T(X, h) is tabulated on a regular grid and interpolated.

Shadow zone: downward rays only reach horizontal distances up to X_max(h) (set by the
largest sound speed below the source). Beyond it the direct ray does not exist; we
continue T linearly with slope 1/c(0) (a ground-guided/diffracted arrival). This is an
approximation and is reported via :meth:`LayeredAtmosphere.in_shadow`.

Winds are not modelled yet.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import RegularGridInterpolator

GAMMA_R = 1.4 * 287.053  # J/(kg K) * gamma, dry air

# US Standard Atmosphere 1976: (base height km, base temperature K, lapse K/km)
_STD1976 = [
    (0.0, 288.15, -6.5),
    (11.0, 216.65, 0.0),
    (20.0, 216.65, 1.0),
    (32.0, 228.65, 2.8),
    (47.0, 270.65, 0.0),
    (51.0, 270.65, -2.8),
    (71.0, 214.65, -2.0),
    (84.852, 186.946, 0.0),
]


def temperature_1976(z_m) -> np.ndarray:
    """US Standard Atmosphere 1976 temperature (K); geometric ~ geopotential height here."""
    z = np.asarray(z_m, dtype=float) / 1e3
    t = np.empty_like(z)
    for i, (zb, tb, lapse) in enumerate(_STD1976):
        top = _STD1976[i + 1][0] if i + 1 < len(_STD1976) else np.inf
        m = (z >= zb) & (z < top) if i else (z < top)
        t[m] = tb + lapse * (z[m] - zb)
    return t


def sound_speed_1976(z_m) -> np.ndarray:
    return np.sqrt(GAMMA_R * temperature_1976(z_m))


class LayeredAtmosphere:
    """Tabulated acoustic travel times T(horizontal distance, source height above receiver)."""

    def __init__(self, z_nodes_m=None, c_nodes=None, *, dz_m=500.0, zmax_m=110e3,
                 xmax_m=800e3, dx_m=500.0, n_rays=600):
        if z_nodes_m is None:
            z_nodes_m = np.arange(0.0, zmax_m + dz_m / 2, dz_m)
            c_nodes = sound_speed_1976(z_nodes_m)
        self.z = np.asarray(z_nodes_m, dtype=float)
        self.c = np.asarray(c_nodes, dtype=float)
        if self.z[0] != 0.0 or np.any(np.diff(self.z) <= 0):
            raise ValueError("z nodes must start at 0 and increase")
        self.x_grid = np.arange(0.0, xmax_m + dx_m / 2, dx_m)
        self.x_max = np.zeros(self.z.size)
        table = np.zeros((self.z.size, self.x_grid.size))
        table[0] = self.x_grid / self.c[0]  # source on the ground
        q = np.concatenate([[1.0], np.logspace(-0.0001, -9, n_rays)])  # 1 - p/pmax
        for k in range(1, self.z.size):
            pmax = 1.0 / self.c[: k + 1].max()
            p = pmax * (1.0 - q)
            x, t = self._trace(p, k)
            order = np.argsort(x)
            x, t = x[order], t[order]
            self.x_max[k] = x[-1]
            inside = self.x_grid <= x[-1]
            table[k, inside] = np.interp(self.x_grid[inside], x, t)
            table[k, ~inside] = t[-1] + (self.x_grid[~inside] - x[-1]) / self.c[0]
        self.table = table
        self._interp = RegularGridInterpolator((self.z, self.x_grid), table,
                                               bounds_error=False, fill_value=None)

    def _trace(self, p: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
        x = np.zeros_like(p)
        t = np.zeros_like(p)
        for i in range(k):
            ca, cb = self.c[i], self.c[i + 1]
            h = self.z[i + 1] - self.z[i]
            cos_a = np.sqrt(np.clip(1.0 - (p * ca) ** 2, 0.0, None))
            cos_b = np.sqrt(np.clip(1.0 - (p * cb) ** 2, 0.0, None))
            g = (cb - ca) / h
            if abs(g) < 1e-12:
                with np.errstate(divide="ignore", invalid="ignore"):
                    x += np.where(cos_a > 0, h * p * ca / cos_a, np.inf)
                    t += np.where(cos_a > 0, h / (ca * cos_a), np.inf)
            else:
                with np.errstate(divide="ignore", invalid="ignore"):
                    dx = np.where(p > 0, (cos_a - cos_b) / (p * g), 0.0)
                x += dx
                t += np.log(cb * (1.0 + cos_a) / (ca * (1.0 + cos_b))) / g
        return x, t

    def travel_time(self, horizontal_m, height_m, scale: float = 1.0) -> np.ndarray:
        """Travel time (s) from a source ``height_m`` above the receiver.

        ``scale`` multiplies the whole sound-speed profile (a fitted temperature offset).
        """
        h = np.clip(np.asarray(height_m, dtype=float), 0.0, self.z[-1])
        d = np.clip(np.asarray(horizontal_m, dtype=float), 0.0, self.x_grid[-1])
        pts = np.stack(np.broadcast_arrays(h, d), axis=-1)
        return self._interp(pts) / scale

    def in_shadow(self, horizontal_m, height_m) -> np.ndarray:
        h = np.clip(np.asarray(height_m, dtype=float), 0.0, self.z[-1])
        return np.asarray(horizontal_m) > np.interp(h, self.z, self.x_max)

    @property
    def min_sound_speed(self) -> float:
        return float(self.c.min())


def msis_atmosphere(when, lat_deg: float, lon_deg: float, *, dz_m: float = 500.0,
                    zmax_m: float = 110e3, **kw) -> LayeredAtmosphere:
    """Event-specific sound-speed profile from the NRLMSIS 2.x empirical model.

    MSIS gives temperature from the ground to the thermosphere for a given date, time and
    place, including season, latitude and solar activity. It is a climatology, not a
    weather analysis: it knows that a February night at 53 N has a cold stratosphere and
    mesosphere, but not that day's weather or winds. Requires ``pymsis`` (no account;
    it downloads public space-weather indices from CelesTrak on first use).
    """
    try:
        import pymsis
    except ImportError as e:  # pragma: no cover
        raise ImportError("msis_atmosphere needs pymsis: pip install pymsis") from e
    z = np.arange(0.0, zmax_m + dz_m / 2, dz_m)
    out = pymsis.calculate(np.datetime64(when, "m"), lon_deg, lat_deg, z / 1e3)
    temp = np.squeeze(out)[..., pymsis.Variable.TEMPERATURE]
    return LayeredAtmosphere(z, np.sqrt(GAMMA_R * temp), **kw)


def layered_arrival_times(traj, stations_enu, atmosphere: LayeredAtmosphere,
                          scale: float = 1.0, return_emission: bool = False):
    """First Mach-cone arrival through a layered atmosphere.

    Arrival = min over emission position s of  t0 + s/v + T(D(s), h(s)), found by a coarse
    grid search around the homogeneous stationary point and a parabolic refinement.
    """
    st = np.asarray(stations_enu, dtype=float)
    u, o, v = traj.direction, traj.origin, traj.speed_mps
    c0 = 320.0 * scale
    d = st - o
    xs = d @ u
    r = np.linalg.norm(d - np.outer(xs, u), axis=1)
    s0 = xs - r * c0 / np.sqrt(max(v**2 - c0**2, 1.0))

    def total(s):  # s: (N, M)
        p = o + s[..., None] * u
        horiz = np.hypot(p[..., 0] - st[:, None, 0], p[..., 1] - st[:, None, 1])
        height = p[..., 2] - st[:, None, 2]
        return traj.t0_s + s / v + atmosphere.travel_time(horiz, height, scale)

    best_s = s0
    for half, n in ((150e3, 61), (6e3, 25), (0.6e3, 13)):
        grid = best_s[:, None] + np.linspace(-half, half, n)[None, :]
        tt = total(grid)
        i = np.clip(np.argmin(tt, axis=1), 1, n - 2)
        rows = np.arange(len(st))
        y0, y1, y2 = tt[rows, i - 1], tt[rows, i], tt[rows, i + 1]
        step = grid[0, 1] - grid[0, 0]
        denom = y0 - 2 * y1 + y2
        shift = np.where(np.abs(denom) > 1e-12, 0.5 * (y0 - y2) / denom, 0.0)
        best_s = grid[rows, i] + np.clip(shift, -1, 1) * step
    t_arr = total(best_s[:, None])[:, 0]
    if return_emission:
        return t_arr, o + np.outer(best_s, u)
    return t_arr
