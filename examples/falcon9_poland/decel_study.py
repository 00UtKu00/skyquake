# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 SKYQUAKE contributors
"""Constant-speed vs decelerating model on the Falcon 9 picks, checked against the camera
speed profile and timing (Vierinen et al. 2026). See compare_optical.py for the data.

    python examples/falcon9_poland/decel_study.py path/to/ref/fragments
"""

import sys
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from skyquake.atmosphere import LayeredAtmosphere
from skyquake.geo import LocalFrame
from skyquake.invert import decluster_sigma, fit_trajectory, jackknife
from skyquake.io import read_picks
from skyquake.physics import emission_time
from skyquake.track import geographic_state

HERE = Path(__file__).parent


def load_optical(folder, frame):
    frag, t, pts = [], [], []
    for f in sorted(Path(folder).glob("*.h5")):
        with h5py.File(f, "r") as h:
            lat, lon, alt = h["geoloc"][()]
            t.append(float(h["time"][()]))
        frag.append(f.name.split("_")[0] + "_" + "_".join(f.name.split("_")[1:5]))
        pts.append(frame.to_enu(lat, lon, alt))
    return np.array(frag), np.array(t), np.array(pts)


def optical_speed_profile(frag, t, pts, traj):
    """Per-track finite-difference speeds, located by along-track position s (km)."""
    s_all, v_all = [], []
    for f in sorted(set(frag)):
        k = np.flatnonzero(frag == f)
        if k.size < 4:
            continue
        o = k[np.argsort(t[k])]
        dp = np.linalg.norm(np.diff(pts[o], axis=0), axis=1)
        dt = np.diff(t[o])
        good = dt > 0
        mid = 0.5 * (pts[o][1:] + pts[o][:-1])
        s = (mid - traj.origin) @ traj.direction / 1e3
        s_all.append(s[good])
        v_all.append(dp[good] / dt[good] / 1e3)
    return np.concatenate(s_all), np.concatenate(v_all)


def main(folder):
    p = read_picks(HERE / "picks.csv")
    frame = LocalFrame.centred_on(p.lat, p.lon)
    st = frame.to_enu(p.lat, p.lon, p.elev_m)
    sig = decluster_sigma(st, p.sigma_s, 20e3)
    kw = dict(azimuth_guess_deg=94.0, outlier_sigma=1e9)
    const = fit_trajectory(st, p.t_s, sig, fit_sound_speed=True, **kw)
    dec = fit_trajectory(st, p.t_s, sig, fit_deceleration=True, **kw)
    lay = fit_trajectory(st, p.t_s, sig, fit_sound_speed=True, atmosphere=LayeredAtmosphere(),
                         **kw)
    jk = jackknife(st, p.t_s, sig, blocks="spatial", fit_sound_speed=True,
                   atmosphere=LayeredAtmosphere(), **kw)
    epoch = p.epoch.timestamp()

    frag, t_opt, pts = load_optical(folder, frame)
    s_opt, v_opt = optical_speed_profile(frag, t_opt, pts, dec.trajectory)
    lo, hi = np.percentile(s_opt, [2, 98])
    bins = np.arange(np.floor(lo / 50) * 50, hi + 50, 50)
    centres, med = [], []
    for a, b in zip(bins[:-1], bins[1:], strict=True):
        m = (s_opt >= a) & (s_opt < b)
        if m.sum() >= 8:
            centres.append(0.5 * (a + b))
            med.append(np.median(v_opt[m]))
    centres, med = np.array(centres), np.array(med)

    def v_model(fit, s_km):
        g = fit.params.get("dv_kms_per_100km", 0.0) / 100.0
        return fit.params["speed_kms"] + g * s_km

    # timing: optical time minus model time at the same along-track position
    def height_offset(fit):
        tr = fit.trajectory
        d = pts - tr.origin
        perp = d - np.outer(d @ tr.direction, tr.direction)
        return perp[:, 2] / 1e3

    def dt_model(fit):
        tr = fit.trajectory
        s = (pts - tr.origin) @ tr.direction
        g = fit.params.get("dv_kms_per_100km", 0.0) * 1e3 / 100e3
        t_line = epoch + tr.t0_s + emission_time(s, tr.speed_mps, g)
        return t_opt - t_line, s / 1e3

    lines = [f"{len(p)} picks (declustered weights); spatial-block jackknife for the "
             "layered model", ""]
    models = (("constant speed", const), ("decelerating", dec), ("layered atmosphere", lay))
    for name, fit in models:
        g = geographic_state(fit.trajectory, frame)
        dt, s = dt_model(fit)
        near = np.abs(s) < 150
        speed_err = np.median(med - v_model(fit, centres))
        lines += [f"{name}:",
                  f"  speed at reference {fit.params['speed_kms']:.2f} km/s, "
                  f"dv per 100 km {fit.params.get('dv_kms_per_100km', 0.0):+.2f} km/s, "
                  f"heading {g['azimuth_deg']:.2f}, height {g['height_m'] / 1e3:.1f} km",
                  f"  rms {fit.rms_s:.3f} s, reduced chi2 {fit.reduced_chi2:.2f}"
                  + (f", profile scale {fit.params['atm_scale']:.3f}"
                     if "atm_scale" in fit.params else ""),
                  f"  height vs cameras near network: "
                  f"{np.median(height_offset(fit)[np.abs(s) < 150]):+.1f} km (cameras - model)",
                  f"  vs cameras: median speed error {speed_err:+.2f} km/s over "
                  f"{len(centres)} 50-km bins; timing (optics - model) near network "
                  f"{np.median(dt[near]):+.1f} s", ""]
    lines += ["jackknife 1-sigma (layered):",
              *[f"  {k:<18} +/- {v:.3f}   (formal {lay.sigma[k]:.3f})" for k, v in jk.items()],
              "", "camera speed profile (50 km bins, along-track km: km/s):",
              "  " + ", ".join(f"{c:+.0f}: {v:.2f}" for c, v in zip(centres, med, strict=True))]
    txt = "\n".join(lines)
    print(txt)
    (HERE / "decel_study.txt").write_text(txt + "\n", encoding="utf-8")

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(s_opt, v_opt, ".", color="0.8", ms=2, label="camera, per fragment pair")
    ax.plot(centres, med, "o", color="k", label="camera, 50 km median")
    ss = np.linspace(bins[0], bins[-1], 100)
    ax.plot(ss, v_model(const, ss), "--", color="C0", label="SKYQUAKE constant speed")
    ax.plot(ss, v_model(dec, ss), "-", color="C3", lw=1, label="SKYQUAKE decelerating")
    ax.plot(ss, v_model(lay, ss), "-", color="C2", lw=2, label="SKYQUAKE layered atmosphere")
    ax.set_ylim(0, 13)
    ax.set_xlabel("along-track distance from the seismic reference point (km)")
    ax.set_ylabel("speed (km/s)")
    ax.legend()
    ax.set_title("Falcon 9 upper stage: seismic speed model vs camera speeds")
    fig.tight_layout()
    fig.savefig(HERE / "speed_profile.png", dpi=110)


if __name__ == "__main__":
    main(sys.argv[1])
