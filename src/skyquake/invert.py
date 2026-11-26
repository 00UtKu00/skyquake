"""Inverse problem: arrival times -> trajectory, with honest uncertainties.

Pipeline
--------
1. Multi-start robust least squares (soft-L1 loss) to survive bad picks and avoid local
   minima in heading.
2. Outlier rejection at ``outlier_sigma`` normalised residuals.
3. Final ordinary least-squares refit on the inliers; covariance from the Gauss-Newton
   Hessian, scaled by the reduced chi-square when the fit is worse than the stated pick
   errors (never scaled *down*, so we do not over-claim precision).
4. Optional Monte-Carlo propagation of that covariance to a ground-track corridor.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

from skyquake.physics import DEFAULT_SOUND_SPEED, Trajectory, arrival_times

PARAM_NAMES = ("east_km", "north_km", "azimuth_deg", "fpa_deg", "speed_kms", "t0_s")
SOUND_PARAM = "sound_speed_kms"


def _to_traj(p: np.ndarray, ref_up_m: float) -> Trajectory:
    return Trajectory(
        east_m=p[0] * 1e3,
        north_m=p[1] * 1e3,
        up_m=ref_up_m,
        azimuth_deg=p[2],
        fpa_deg=p[3],
        speed_mps=p[4] * 1e3,
        t0_s=p[5],
    )


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
    ref_up_m: float
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
        return [_to_traj(d, self.ref_up_m) for d in draws]

    def summary(self) -> str:
        lines = [f"stations used: {self.n_used}   rms residual: {self.rms_s:.3f} s   "
                 f"reduced chi2: {self.reduced_chi2:.2f}"]
        for k in self.names:
            lines.append(f"  {k:<12} {self.params[k]:>12.3f}  +/- {self.sigma[k]:.3f}")
        lines += [f"  note: {n}" for n in self.notes]
        return "\n".join(lines)


def fit_trajectory(
    stations_enu,
    arrival_s,
    sigma_s=0.5,
    *,
    sound_speed: float = DEFAULT_SOUND_SPEED,
    ref_up_m: float = 40_000.0,
    outlier_sigma: float = 4.0,
    azimuth_guess_deg: float | None = None,
    fit_sound_speed: bool = False,
) -> FitResult:
    """Fit a straight-line trajectory to sonic-boom arrival times.

    Parameters
    ----------
    stations_enu : (N, 3) station positions in a local ENU frame (m).
    arrival_s    : (N,) arrival times in seconds on any common clock.
    sigma_s      : pick uncertainty, scalar or (N,) (s, 1-sigma).
    ref_up_m     : the height at which the reference point (and t0) is defined.
    azimuth_guess_deg : e.g. the catalogued orbit's heading, if known. Narrows the search.
    fit_sound_speed : also estimate the path-averaged sound speed (7th parameter).
        Strongly recommended for real data: a 5 % sound-speed error moves the along-track
        position by ~100 km in the synthetic benchmark (see docs/validation.md).
    """
    st = np.asarray(stations_enu, dtype=float)
    t = np.asarray(arrival_s, dtype=float)
    sig = np.broadcast_to(np.asarray(sigma_s, dtype=float), t.shape).copy()
    if st.shape != (t.size, 3):
        raise ValueError("stations_enu must have shape (N, 3) matching arrival_s")
    if t.size < 8:
        raise ValueError("need at least 8 arrivals for a 6-parameter fit with outlier control")

    names = PARAM_NAMES + ((SOUND_PARAM,) if fit_sound_speed else ())
    lo = np.array([-5e3, -5e3, -1e4, 0.0, 0.45, -np.inf] + ([0.25] if fit_sound_speed else []))
    hi = np.array([5e3, 5e3, 1e4, 60.0, 12.0, np.inf] + ([0.40] if fit_sound_speed else []))

    def c_of(p):
        return p[6] * 1e3 if fit_sound_speed else sound_speed

    def resid(p, mask):
        pred = arrival_times(_to_traj(p, ref_up_m), st[mask], c_of(p))
        return (pred - t[mask]) / sig[mask]

    first = int(np.argmin(t))
    e0, n0 = st[first, 0] / 1e3, st[first, 1] / 1e3
    t0_guess = t[first] - ref_up_m / sound_speed
    azimuths = (
        np.arange(0.0, 360.0, 30.0)
        if azimuth_guess_deg is None
        else azimuth_guess_deg + np.array([-20.0, 0.0, 20.0])
    )

    all_mask = np.ones(t.size, bool)
    best = None
    for az in azimuths:
        for fpa in (3.0, 15.0):
            for v in (2.5, 6.0):
                start = [e0, n0, az, fpa, v, t0_guess]
                if fit_sound_speed:
                    start.append(sound_speed / 1e3)
                p0 = np.clip(start, lo + 1e-6, hi - 1e-6)
                r = least_squares(resid, p0, args=(all_mask,), bounds=(lo, hi),
                                  loss="soft_l1", f_scale=2.0, x_scale="jac", max_nfev=400)
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

    traj = _to_traj(p, ref_up_m).normalised()
    p[2] = traj.azimuth_deg
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
        ref_up_m=ref_up_m,
        n_used=int(mask.sum()),
        notes=notes,
    )
