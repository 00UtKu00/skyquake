# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 SKYQUAKE contributors
"""Why is the fitted speed too high? Quick diagnostics on picks.csv."""

from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from skyquake import LocalFrame, fit_trajectory
from skyquake.invert import PARAM_NAMES, _to_traj
from skyquake.io import read_picks
from skyquake.physics import arrival_times

HERE = Path(__file__).parent
p = read_picks(HERE / "picks.csv")
frame = LocalFrame.centred_on(p.lat, p.lon)
st = frame.to_enu(p.lat, p.lon, p.elev_m)
fit = fit_trajectory(st, p.t_s, p.sigma_s, fit_sound_speed=True, azimuth_guess_deg=65)
print(fit.summary())

tr = fit.trajectory
u = tr.direction
d = st - tr.origin
x = d @ u / 1e3
right = np.array([np.cos(np.radians(tr.azimuth_deg)), -np.sin(np.radians(tr.azimuth_deg)), 0])
cross = d @ right / 1e3
res = fit.residuals_s
print("\n1) residual trend (s per 100 km):")
for name, q in (("along", x), ("cross", cross), ("|cross|", np.abs(cross))):
    k = np.polyfit(q, res, 1)[0] * 100
    print(f"   vs {name:<8} {k:+.3f}   (range {q.min():.0f} .. {q.max():.0f} km)")
print(f"   corr(along, |cross|) = {np.corrcoef(x, np.abs(cross))[0, 1]:+.2f}")

print("\n2) chi2 with speed fixed:")
p0 = np.array([fit.params[k] for k in PARAM_NAMES] + [fit.params["sound_speed_kms"]])
for v in (6.5, 7.0, 7.5, 7.8, fit.params["speed_kms"]):
    def r(q, v=v):
        qq = q.copy()
        full = np.insert(qq, 4, v)
        return (arrival_times(_to_traj(full[:6]), st, full[6] * 1e3) - p.t_s) / p.sigma_s
    q0 = np.delete(p0, 4)
    sol = least_squares(r, q0, loss="linear", x_scale="jac")
    chi2 = (sol.fun @ sol.fun) / (len(p.t_s) - 6)
    print(f"   v = {v:5.2f} km/s  reduced chi2 = {chi2:6.2f}  rms = "
          f"{np.sqrt(np.mean((sol.fun * p.sigma_s) ** 2)):.3f} s")

print("\n3) jackknife (drop 10% of stations at a time, 10 groups):")
rng = np.random.default_rng(0)
groups = np.array_split(rng.permutation(len(p.t_s)), 10)
vs = []
for g in groups:
    keep = np.setdiff1d(np.arange(len(p.t_s)), g)
    f = fit_trajectory(st[keep], p.t_s[keep], p.sigma_s[keep], fit_sound_speed=True,
                       azimuth_guess_deg=65)
    vs.append(f.params["speed_kms"])
vs = np.array(vs)
n = len(vs)
jk_sigma = np.sqrt((n - 1) / n * np.sum((vs - vs.mean()) ** 2))
print(f"   speeds: {np.round(vs, 2)}")
print(f"   jackknife sigma = {jk_sigma:.2f} km/s  (formal {fit.sigma['speed_kms']:.2f})")
