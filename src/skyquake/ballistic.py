"""Physical re-entry trajectories: gravity, drag, Earth curvature and rotation.

The straight constant-speed line is a good model for *where* an object went (heading,
cross-track position, height), but its speed parameter absorbs everything the line cannot
represent. Real fragments follow curved, decelerating paths (Falcon 9, 2025: the
single-line fit returned 10-12 km/s for a fragment the cameras measured at 7 km/s).

This module integrates the equations of motion in the Earth-fixed (ECEF) frame:

    a = -mu r/|r|^3  -  2 Omega x v  -  Omega x (Omega x r)  -  rho(h) |v| v / (2 beta)

with the WGS84 gravitational parameter, Earth rotation Omega, air density rho(h), and the
ballistic coefficient beta = m / (Cd A) in kg/m^2. Winds are neglected (v relative to the
ground ~ relative to the air). Booms are emitted along the integrated path; each station
hears the first arrival, the minimum over path samples of  t + T(horizontal, height)
through a layered atmosphere.

Parameters (same reference point as the straight-line model, plus beta):
    cross_km, height_km, azimuth_deg, fpa_deg, speed_kms (at the reference point), t0_s,
    log10_beta, [atm_scale]
"""

from __future__ import annotations

import numpy as np

from skyquake.atmosphere import LayeredAtmosphere
from skyquake.geo import LocalFrame, geodetic_to_ecef

MU = 3.986004418e14       # m^3/s^2
OMEGA = 7.2921159e-5      # rad/s
H_SCALE = 7.2e3           # m, fallback exponential atmosphere
RHO0 = 1.225              # kg/m^3

# US Standard Atmosphere 1976 density (kg/m^3) at 10 km steps, 0..100 km (fallback)
_RHO_1976 = np.array([1.225, 4.135e-1, 8.891e-2, 1.841e-2, 3.996e-3, 1.027e-3, 3.097e-4,
                      8.283e-5, 1.846e-5, 3.416e-6, 5.604e-7])


def density_1976(h_m) -> np.ndarray:
    h = np.clip(np.asarray(h_m, dtype=float), 0.0, 100e3)
    return np.exp(np.interp(h, np.arange(0.0, 100e3 + 1, 10e3), np.log(_RHO_1976)))


def msis_density(when, lat_deg, lon_deg, zmax_m=120e3, dz_m=1000.0):
    """Density profile callable from NRLMSIS (requires pymsis)."""
    import pymsis

    z = np.arange(0.0, zmax_m + dz_m / 2, dz_m)
    out = pymsis.calculate(np.datetime64(when, "m"), lon_deg, lat_deg, z / 1e3)
    rho = np.squeeze(out)[..., pymsis.Variable.MASS_DENSITY]
    lr = np.log(rho)
    return lambda h: np.exp(np.interp(np.clip(h, 0.0, zmax_m), z, lr))


class BallisticModel:
    """Propagates a re-entry state and predicts first boom arrivals at stations."""

    def __init__(self, frame: LocalFrame, atmosphere: LayeredAtmosphere, density=None,
                 t_before=150.0, t_after=200.0, dt=1.0):
        self.frame = frame
        self.atm = atmosphere
        self.density = density or density_1976
        self.t_before, self.t_after, self.dt = t_before, t_after, dt
        self._R = frame._rotation              # ENU <- ECEF
        self._o = frame._origin
        lat = np.radians(frame.lat0)
        a, b = 6378137.0, 6356752.314245
        self._r_earth = np.sqrt(((a**2 * np.cos(lat))**2 + (b**2 * np.sin(lat))**2)
                                / ((a * np.cos(lat))**2 + (b * np.sin(lat))**2))
        self._omega = np.array([0.0, 0.0, OMEGA])

    def initial_state(self, p):
        """ECEF position and velocity at the reference point from the parameter vector."""
        cross, h, az, fpa, v0 = p[0] * 1e3, p[1] * 1e3, np.radians(p[2]), np.radians(p[3]), \
            p[4] * 1e3
        right = np.array([np.cos(az), -np.sin(az), 0.0])
        pos_enu = cross * right
        pos_enu[2] = h
        vel_enu = v0 * np.array([np.sin(az) * np.cos(fpa), np.cos(az) * np.cos(fpa),
                                 -np.sin(fpa)])
        return self._o + pos_enu @ self._R, vel_enu @ self._R

    def _make_accel(self, inv2beta):
        """Acceleration in plain floats (numpy on 3-vectors is ~10x slower here)."""
        import math

        w, w2, R = OMEGA, OMEGA * OMEGA, self._r_earth
        zg, lr = self._rho_table

        def accel(x, y, z, vx, vy, vz):
            rn = math.sqrt(x * x + y * y + z * z)
            hh = rn - R
            if hh <= 0.0:
                rho = math.exp(lr[0])
            elif hh >= zg[-1]:
                rho = math.exp(lr[-1])
            else:
                i = int(hh / self._rho_dz)
                f = hh / self._rho_dz - i
                rho = math.exp(lr[i] + f * (lr[i + 1] - lr[i]))
            g = -MU / (rn * rn * rn)
            dr = -rho * math.sqrt(vx * vx + vy * vy + vz * vz) * inv2beta
            return (g * x + 2 * w * vy + w2 * x + dr * vx,
                    g * y - 2 * w * vx + w2 * y + dr * vy,
                    g * z + dr * vz)

        return accel

    @property
    def _rho_table(self):
        if not hasattr(self, "_rho_cache"):
            self._rho_dz = 500.0
            zg = np.arange(0.0, 150e3 + 1, self._rho_dz)
            self._rho_cache = (zg, np.log(np.maximum(self.density(zg), 1e-15)).tolist())
        return self._rho_cache

    def propagate(self, p):
        """Path samples: times (s, relative to the reference), ENU positions (m), speeds."""
        import math

        r0, v0 = self.initial_state(p)
        accel = self._make_accel(0.5 / 10.0 ** p[6])
        R = self._r_earth
        out = {}
        for sign, span in ((1.0, self.t_after), (-1.0, self.t_before)):
            x, y, z = map(float, r0)
            vx, vy, vz = map(float, v0)
            h = sign * self.dt
            ts, rs, vs = [0.0], [(x, y, z)], [math.sqrt(vx * vx + vy * vy + vz * vz)]
            for k in range(1, int(span / self.dt) + 1):
                a1 = accel(x, y, z, vx, vy, vz)
                b = 0.5 * h
                a2 = accel(x + b * vx, y + b * vy, z + b * vz,
                           vx + b * a1[0], vy + b * a1[1], vz + b * a1[2])
                v2 = (vx + b * a1[0], vy + b * a1[1], vz + b * a1[2])
                a3 = accel(x + b * v2[0], y + b * v2[1], z + b * v2[2],
                           vx + b * a2[0], vy + b * a2[1], vz + b * a2[2])
                v3 = (vx + b * a2[0], vy + b * a2[1], vz + b * a2[2])
                a4 = accel(x + h * v3[0], y + h * v3[1], z + h * v3[2],
                           vx + h * a3[0], vy + h * a3[1], vz + h * a3[2])
                v4 = (vx + h * a3[0], vy + h * a3[1], vz + h * a3[2])
                c = h / 6.0
                x += c * (vx + 2 * v2[0] + 2 * v3[0] + v4[0])
                y += c * (vy + 2 * v2[1] + 2 * v3[1] + v4[1])
                z += c * (vz + 2 * v2[2] + 2 * v3[2] + v4[2])
                vx += c * (a1[0] + 2 * a2[0] + 2 * a3[0] + a4[0])
                vy += c * (a1[1] + 2 * a2[1] + 2 * a3[1] + a4[1])
                vz += c * (a1[2] + 2 * a2[2] + 2 * a3[2] + a4[2])
                speed = math.sqrt(vx * vx + vy * vy + vz * vz)
                alt = math.sqrt(x * x + y * y + z * z) - R
                # stop: hit the ground, subsonic, or (integrating backwards) drag running away
                if not math.isfinite(speed) or alt < 5e3 or alt > 200e3 or speed < 200.0 \
                        or speed > 12e3:
                    break
                ts.append(sign * k * self.dt)
                rs.append((x, y, z))
                vs.append(speed)
            out[sign] = (np.array(ts), np.array(rs), np.array(vs))
        t = np.concatenate([out[-1.0][0][::-1], out[1.0][0][1:]])
        r = np.concatenate([out[-1.0][1][::-1], out[1.0][1][1:]])
        v = np.concatenate([out[-1.0][2][::-1], out[1.0][2][1:]])
        enu = (r - self._o) @ self._R.T
        return t, enu, v

    def arrivals(self, p, stations_enu, return_detail=False):
        st = np.asarray(stations_enu, dtype=float)
        scale = p[7] if len(p) > 7 else 1.0
        t, enu, speed = self.propagate(p)
        horiz = np.hypot(enu[None, :, 0] - st[:, None, 0], enu[None, :, 1] - st[:, None, 1])
        height = enu[None, :, 2] - st[:, None, 2]
        tt = p[5] + t[None, :] + self.atm.travel_time(horiz, height, scale)
        cmin = self.atm.min_sound_speed * scale
        tt = np.where((speed > 1.05 * cmin)[None, :], tt, np.inf)
        n = t.size
        k = np.clip(np.argmin(tt, axis=1), 1, n - 2)
        rows = np.arange(len(st))
        y0, y1, y2 = tt[rows, k - 1], tt[rows, k], tt[rows, k + 1]
        denom = y0 - 2 * y1 + y2
        ok = np.isfinite(y0) & np.isfinite(y2) & (np.abs(denom) > 1e-12)
        y = np.where(ok, y1 - 0.125 * (y0 - y2) ** 2 / np.where(ok, denom, 1.0), y1)
        y = np.where(np.isfinite(y), y, 1e6)
        if return_detail:
            return y, dict(t=t, enu=enu, speed=speed, k=k)
        return y


def reference_speed_profile(model: BallisticModel, p, s_km):
    """Speed (km/s) vs along-track distance (km) from the reference point."""
    t, enu, v = model.propagate(p)
    ds = np.linalg.norm(np.diff(enu, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(ds)])
    s -= s[np.argmin(np.abs(t))]
    return np.interp(s_km * 1e3, s, v) / 1e3


BALLISTIC_NAMES = ("cross_km", "height_km", "azimuth_deg", "fpa_deg", "speed_kms", "t0_s",
                   "log10_beta")


def fit_ballistic(stations_enu, arrival_s, sigma_s, model: BallisticModel, *,
                  start=None, fit_scale=False, outlier_sigma=4.0, max_speed_kms=8.5,
                  log10_beta_range=(1.0, 5.0), log=None):
    """Fit a physical re-entry trajectory to boom arrival times.

    ``start``: a straight-line FitResult (same frame) used as the starting geometry;
    computed if omitted. The speed bound defaults to 8.5 km/s (orbital speed plus margin).
    ``log10_beta_range``: allowed ballistic coefficient (kg/m^2, log10). Speed and drag
    trade off ("fast with strong drag" vs "slower with normal drag" produce similar boom
    patterns), so a physical range for the object type is a meaningful prior: e.g.
    (1.7, 4.0) = 50 ... 10 000 kg/m^2 for dense rocket-stage fragments.
    """
    from scipy.optimize import least_squares

    from skyquake.invert import FitResult, _covariance, _to_traj, fit_trajectory

    st = np.asarray(stations_enu, dtype=float)
    t = np.asarray(arrival_s, dtype=float)
    sig = np.broadcast_to(np.asarray(sigma_s, dtype=float), t.shape).copy()
    if start is None:
        start = fit_trajectory(st, t, sig, fit_sound_speed=True, outlier_sigma=outlier_sigma)
    names = BALLISTIC_NAMES + (("atm_scale",) if fit_scale else ())
    lo = np.array([-2e3, 20.0, -1e4, -20.0, 1.0, -np.inf, log10_beta_range[0]]
                  + ([0.85] if fit_scale else []))
    hi = np.array([2e3, 130.0, 1e4, 30.0, max_speed_kms, np.inf, log10_beta_range[1]]
                  + ([1.15] if fit_scale else []))

    def resid(p, mask):
        return (model.arrivals(p, st[mask]) - t[mask]) / sig[mask]

    def cost(p):
        r = resid(p, np.ones(t.size, bool))
        return float(np.sum(4.0 * (np.sqrt(1.0 + (r / 2.0) ** 2) - 1.0)))

    base = np.array([start.params[k] for k in BALLISTIC_NAMES[:6]] + [2.5]
                    + ([1.0] if fit_scale else []))
    # Grid over height, speed and drag. The reference time is not gridded: for each
    # candidate it is set to the median residual, which removes the strong time-height-
    # speed trade-off that trapped the optimiser (Falcon 9: the straight-line start was
    # ~25 s and ~6 km off the true fragment).
    grid = []
    for dh in np.arange(-6.0, 12.1, 2.0):
        for v0 in np.arange(5.5, min(max_speed_kms, 8.0) + 1e-9, 0.5):
            for lb in np.linspace(max(log10_beta_range[0], 1.5), min(log10_beta_range[1], 3.5),
                                  5):
                q = base.copy()
                q[1], q[4], q[6], q[5] = q[1] + dh, v0, lb, 0.0
                q = np.clip(q, lo + 1e-6, hi - 1e-6)
                q[5] = float(np.median(t - model.arrivals(q, st)))
                grid.append((cost(q), q))
    grid.sort(key=lambda g: g[0])
    if log:
        log(f"ballistic grid: best cost {grid[0][0]:.0f} at v0={grid[0][1][4]:.1f}, "
            f"log10 beta={grid[0][1][6]:.1f}")
    all_mask = np.ones(t.size, bool)
    best = None
    for _, q in grid[:3]:
        r = least_squares(resid, q, args=(all_mask,), bounds=(lo, hi), loss="soft_l1",
                          f_scale=2.0, x_scale="jac", diff_step=1e-4, max_nfev=200)
        if best is None or r.cost < best.cost:
            best = r
    notes = ["physical trajectory: gravity, drag, Earth curvature and rotation",
             "robust (soft-L1) final fit: a few stations beyond the direct-ray shadow-zone "
             "limit, where layered travel times are approximate, must not steer the speed"]
    mask = np.abs(resid(best.x, all_mask)) < outlier_sigma
    if (~mask).any():
        notes.append(f"rejected {int((~mask).sum())} outlier pick(s) > {outlier_sigma} sigma")
    final = least_squares(resid, best.x, args=(mask,), bounds=(lo, hi), loss="soft_l1",
                          f_scale=2.0, x_scale="jac", diff_step=1e-4, max_nfev=400)
    p = final.x
    cov, chi2 = _covariance(final, int(mask.sum()), notes)
    std = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    for i, name in enumerate(names):
        if np.isclose(p[i], lo[i], atol=1e-3) or np.isclose(p[i], hi[i], atol=1e-3):
            notes.append(f"{name} hit its bound; result is not a true optimum for it")
    p[2] = p[2] % 360.0
    pred = model.arrivals(p, st)
    return FitResult(
        trajectory=_to_traj(p[:6]),   # tangent line at the reference point
        params=dict(zip(names, map(float, p), strict=True)),
        sigma=dict(zip(names, map(float, std), strict=True)),
        covariance=cov,
        residuals_s=pred - t,
        inliers=mask,
        rms_s=float(np.sqrt(np.mean((pred - t)[mask] ** 2))),
        reduced_chi2=chi2,
        sound_speed=model.atm.min_sound_speed * (p[7] if fit_scale else 1.0),
        n_used=int(mask.sum()),
        notes=notes,
    )


__all__ = ["BallisticModel", "BALLISTIC_NAMES", "density_1976", "fit_ballistic",
           "msis_density", "reference_speed_profile", "geodetic_to_ecef"]
