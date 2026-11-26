"""Inverse problem: arrival times -> trajectory, with honest uncertainties.

Parameterisation
----------------
The reference point is where the line passes the frame origin horizontally (its closest
horizontal approach to the network centre). Parameters:

    cross_km    signed horizontal offset of that point, + = to the right of the heading
    height_km   height of the line at that point
    azimuth_deg heading, clockwise from the frame's north
    fpa_deg     descent angle below the frame's horizontal (+ = descending)
    speed_kms   speed
    t0_s        time at that point
    [sound_speed_kms]  optional: path-averaged sound speed

This stays well conditioned for level flight (real re-entries are often within a degree
of horizontal over a regional network), unlike "the point where the line is at 40 km".

Pipeline
--------
1. Multi-start robust least squares (soft-L1 loss) to survive bad picks and avoid local
   minima in heading.
2. Outlier rejection at ``outlier_sigma`` normalised residuals.
3. Final ordinary least-squares refit on the inliers; covariance from the Gauss-Newton
   Hessian, scaled by the reduced chi-square when the fit is worse than the stated pick
   errors (never scaled *down*, so we do not over-claim precision).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

from skyquake.physics import DEFAULT_SOUND_SPEED, Trajectory, arrival_times

PARAM_NAMES = ("cross_km", "height_km", "azimuth_deg", "fpa_deg", "speed_kms", "t0_s")
SOUND_PARAM = "sound_speed_kms"


def _to_traj(p: np.ndarray) -> Trajectory:
    az = np.radians(p[2])
    right = np.array([np.cos(az), -np.sin(az)])
    return Trajectory(
        east_m=p[0] * 1e3 * right[0],
        north_m=p[0] * 1e3 * right[1],
        up_m=p[1] * 1e3,
        azimuth_deg=p[2],
        fpa_deg=p[3],
        speed_mps=p[4] * 1e3,
        t0_s=p[5],
    )


def traj_to_params(traj: Trajectory) -> np.ndarray:
    """Inverse of the parameterisation (moves the reference point to the canonical one)."""
    c = traj.canonical()
    az = np.radians(c.azimuth_deg)
    right = np.array([np.cos(az), -np.sin(az)])
    return np.array([c.origin[:2] @ right / 1e3, c.up_m / 1e3, c.azimuth_deg, c.fpa_deg,
                     c.speed_mps / 1e3, c.t0_s])


@dataclass
class FitResult:
    trajectory: Trajectory
    params: dict[str, float]
    sigma: dict[str, float]
    covariance: np.ndarray
    residuals_s: np.ndarray
    inliers: np.ndarray
    rms_s: float
    reduced_chi2: float
    sound_speed: float
    n_used: int
    notes: list[str] = field(default_factory=list)

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(self.params)

    def sample(self, n: int, rng: np.random.Generator | None = None) -> list[Trajectory]:
        """Draw trajectories from the (Gaussian) posterior approximation.

        When the sound speed was fitted, its uncertainty and its correlation with the
        geometry are part of the covariance, so the samples include that error source.
        """
        rng = rng or np.random.default_rng(0)
        mean = np.array([self.params[k] for k in self.names])
        draws = rng.multivariate_normal(mean, self.covariance, size=n)
        return [_to_traj(d) for d in draws]

    def summary(self) -> str:
        lines = [f"stations used: {self.n_used}   rms residual: {self.rms_s:.3f} s   "
                 f"reduced chi2: {self.reduced_chi2:.2f}"]
        for k in self.names:
            lines.append(f"  {k:<16} {self.params[k]:>12.3f}  +/- {self.sigma[k]:.3f}")
        lines += [f"  note: {n}" for n in self.notes]
        return "\n".join(lines)


def fit_trajectory(
    stations_enu,
    arrival_s,
    sigma_s=0.5,
    *,
    sound_speed: float = DEFAULT_SOUND_SPEED,
    outlier_sigma: float = 4.0,
    azimuth_guess_deg: float | None = None,
    fit_sound_speed: bool = False,
) -> FitResult:
    """Fit a straight-line trajectory to sonic-boom arrival times.

    Parameters
    ----------
    stations_enu : (N, 3) station positions in a local ENU frame (m). Put the frame near
        the network centre: the trajectory is referenced to its closest approach to it.
    arrival_s    : (N,) arrival times in seconds on any common clock.
    sigma_s      : pick uncertainty, scalar or (N,) (s, 1-sigma).
    azimuth_guess_deg : e.g. the catalogued orbit's heading, if known. Narrows the search.
    fit_sound_speed : also estimate the path-averaged sound speed (7th parameter).
        Strongly recommended for real data: a 5 % sound-speed error moves the trajectory
        by tens of km in the synthetic benchmark (see docs/validation.md).
    """
    st = np.asarray(stations_enu, dtype=float)
    t = np.asarray(arrival_s, dtype=float)
    sig = np.broadcast_to(np.asarray(sigma_s, dtype=float), t.shape).copy()
    if st.shape != (t.size, 3):
        raise ValueError("stations_enu must have shape (N, 3) matching arrival_s")
    if t.size < 8:
        raise ValueError("need at least 8 arrivals for a 6-parameter fit with outlier control")

    names = PARAM_NAMES + ((SOUND_PARAM,) if fit_sound_speed else ())
    lo = np.array([-2e3, 5.0, -1e4, -15.0, 0.45, -np.inf] + ([0.25] if fit_sound_speed else []))
    hi = np.array([2e3, 150.0, 1e4, 75.0, 12.0, np.inf] + ([0.40] if fit_sound_speed else []))

    def c_of(p):
        return p[6] * 1e3 if fit_sound_speed else sound_speed

    def resid(p, mask):
        pred = arrival_times(_to_traj(p), st[mask], c_of(p))
        return (pred - t[mask]) / sig[mask]

    first = int(np.argmin(t))
    azimuths = (
        np.arange(0.0, 360.0, 30.0)
        if azimuth_guess_deg is None
        else azimuth_guess_deg + np.array([-20.0, 0.0, 20.0])
    )

    all_mask = np.ones(t.size, bool)
    best = None
    for az in azimuths:
        right = np.array([np.cos(np.radians(az)), -np.sin(np.radians(az))])
        cross0 = st[first, :2] @ right / 1e3  # the first station is closest to the line
        for fpa in (1.0, 15.0):
            for v in (2.5, 7.0):
                for h in (40.0, 80.0):
                    start = [cross0, h, az, fpa, v, t[first] - h * 1e3 / sound_speed]
                    if fit_sound_speed:
                        start.append(sound_speed / 1e3)
                    p0 = np.clip(start, lo + 1e-6, hi - 1e-6)
                    r = least_squares(resid, p0, args=(all_mask,), bounds=(lo, hi),
                                      loss="soft_l1", f_scale=2.0, x_scale="jac",
                                      max_nfev=400)
                    if best is None or r.cost < best.cost:
                        best = r

    notes: list[str] = []
    full = resid(best.x, all_mask)
    mask = np.abs(full) < outlier_sigma
    if (~mask).any():
        notes.append(f"rejected {int((~mask).sum())} outlier pick(s) > {outlier_sigma} sigma")
    if mask.sum() < 8:
        raise RuntimeError("too few inliers after outlier rejection; check picks")

    final = least_squares(resid, best.x, args=(mask,), bounds=(lo, hi), loss="linear",
                          x_scale="jac", max_nfev=2000)
    p = final.x
    r_in = final.fun
    dof = max(int(mask.sum()) - p.size, 1)
    chi2 = float(r_in @ r_in / dof)
    jac = final.jac
    try:
        cov = np.linalg.inv(jac.T @ jac) * max(chi2, 1.0)
    except np.linalg.LinAlgError:
        cov = np.full((p.size, p.size), np.nan)
        notes.append("singular Jacobian: geometry cannot constrain all parameters")
    std = np.sqrt(np.clip(np.diag(cov), 0.0, None))

    for i, name in enumerate(names):
        if np.isclose(p[i], lo[i], atol=1e-3) or np.isclose(p[i], hi[i], atol=1e-3):
            notes.append(f"{name} hit its bound; result is not a true optimum for it")
    if std[3] > 10:
        notes.append("descent angle poorly constrained (typical when stations lie on one side)")
    if not fit_sound_speed:
        notes.append(f"sound speed fixed at {sound_speed:.0f} m/s; its error is NOT in the "
                     "uncertainties (use fit_sound_speed=True)")

    p[2] = p[2] % 360.0
    traj = _to_traj(p)
    c_final = c_of(p)
    pred_all = arrival_times(traj, st, c_final)
    return FitResult(
        trajectory=traj,
        params=dict(zip(names, map(float, p), strict=True)),
        sigma=dict(zip(names, map(float, std), strict=True)),
        covariance=cov,
        residuals_s=pred_all - t,
        inliers=mask,
        rms_s=float(np.sqrt(np.mean((pred_all - t)[mask] ** 2))),
        reduced_chi2=chi2,
        sound_speed=c_final,
        n_used=int(mask.sum()),
        notes=notes,
    )
