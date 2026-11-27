# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 SKYQUAKE contributors
"""Physical trajectory fit (gravity, drag, curvature, rotation; MSIS atmosphere) on the
Falcon 9 picks, compared with camera fragment F1, whose booms the picks were shown to be
(fragment_test.py: observed - predicted = -0.5 s, spread 0.81 s).

    python examples/falcon9_poland/ballistic_study.py path/to/ref/fragments
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from skyquake.atmosphere import msis_atmosphere
from skyquake.ballistic import BallisticModel, fit_ballistic, msis_density
from skyquake.geo import LocalFrame
from skyquake.invert import decluster_sigma, fit_trajectory, spatial_blocks
from skyquake.io import read_picks
from skyquake.track import geographic_state

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from fragment_test import load_fragments  # noqa: E402

WHEN, LAT, LON = "2025-02-19T03:45", 53.0, 7.0


def main(folder):
    p = read_picks(HERE / "picks.csv")
    frame = LocalFrame.centred_on(p.lat, p.lon)
    st = frame.to_enu(p.lat, p.lon, p.elev_m)
    sig = decluster_sigma(st, p.sigma_s, 20e3)
    epoch = p.epoch.timestamp()
    model = BallisticModel(frame, msis_atmosphere(WHEN, LAT, LON),
                           density=msis_density(WHEN, LAT, LON))
    line = fit_trajectory(st, p.t_s, sig, fit_sound_speed=True, azimuth_guess_deg=94.0,
                          outlier_sigma=1e9)
    # speed and drag trade off in the boom data; a physical prior for a dense rocket-stage
    # fragment (beta >= 50 kg/m^2) is used for the reported result, and the sensitivity to
    # that prior is printed below.
    prior = (1.7, 4.0)
    fit = fit_ballistic(st, p.t_s, sig, model, start=line, outlier_sigma=1e9,
                        log10_beta_range=prior, log=print)
    sens = {r: fit_ballistic(st, p.t_s, sig, model, start=line, outlier_sigma=1e9,
                             log10_beta_range=r) for r in ((1.0, 5.0), (2.0, 4.0))}

    # spatial-block jackknife of the physical fit
    est = []
    for g in spatial_blocks(st, 8):
        keep = np.setdiff1d(np.arange(len(p)), g)
        f = fit_ballistic(st[keep], p.t_s[keep], sig[keep], model, start=line,
                          outlier_sigma=1e9, log10_beta_range=prior)
        est.append(f.params)
    nb = len(est)
    jk = {k: float(np.sqrt((nb - 1) / nb * np.sum(
        (np.array([e[k] for e in est]) - np.mean([e[k] for e in est])) ** 2)))
        for k in fit.params}

    # camera fragment F1
    t1, q1 = load_fragments(folder, frame)["1"]
    tt, enu, spd = model.propagate(np.array([fit.params[k] for k in fit.names]))
    t_abs = epoch + fit.params["t0_s"] + tt
    # time / height / cross offsets of F1 relative to the fitted path, matched by position
    d = np.linalg.norm(q1[:, None, :2] - enu[None, :, :2], axis=2)
    j = np.argmin(d, axis=1)
    near = np.hypot(q1[:, 0], q1[:, 1]) < 150e3
    dt = t1 - t_abs[j]
    dh = (q1[:, 2] - enu[j, 2]) / 1e3
    v1 = np.linalg.norm(np.diff(q1, axis=0), axis=1) / np.diff(t1)
    s1 = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(q1, axis=0), axis=1))])
    sm = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(enu, axis=0), axis=1))])
    v_model_at_f1 = np.interp(np.arange(len(q1) - 1), np.arange(len(q1) - 1), spd[j[:-1]])
    k_ref = int(np.argmin(np.hypot(q1[:, 0], q1[:, 1])))
    g = geographic_state(fit.trajectory, frame)

    lines = [
        "SKYQUAKE physical trajectory fit - Falcon 9 upper stage, 2025-02-19",
        f"atmosphere and density: NRLMSIS for {WHEN} UTC, {LAT}N {LON}E",
        fit.summary(), "",
        f"heading {g['azimuth_deg']:.2f} deg true at the reference point",
        "spatial-block jackknife 1-sigma (8 regions):",
        *[f"  {k:<14} +/- {v:.3f}" for k, v in jk.items()], "",
        "compared with camera fragment F1 within 150 km of the network:",
        f"  speed: SKYQUAKE {fit.params['speed_kms']:.2f} km/s at the reference point, "
        f"F1 {np.median(v1[np.hypot(q1[:-1, 0], q1[:-1, 1]) < 150e3]) / 1e3:.2f} km/s",
        f"  time (F1 - model, same place): median {np.median(dt[near]):+.1f} s",
        f"  height (F1 - model): median {np.median(dh[near]):+.1f} km",
        f"  straight-line fit for comparison: speed {line.params['speed_kms']:.2f} km/s", "",
        f"sensitivity to the drag prior (reported: log10 beta in {prior}):",
        *[f"  log10 beta in {r}: speed {f.params['speed_kms']:.2f} km/s, log10 beta "
          f"{f.params['log10_beta']:.2f}, height {f.params['height_km']:.1f} km, "
          f"rms {f.rms_s:.3f} s" for r, f in sens.items()],
    ]
    txt = "\n".join(lines)
    print(txt)
    (HERE / "ballistic_study.txt").write_text(txt + "\n", encoding="utf-8")

    fig, ax = plt.subplots(1, 2, figsize=(14, 5))
    ax[0].plot((s1[:-1] - s1[k_ref]) / 1e3, v1 / 1e3, ".", color="0.6", ms=3,
               label="camera F1")
    ax[0].plot((sm - sm[np.argmin(np.abs(tt))]) / 1e3, spd / 1e3, "-", color="C3", lw=2,
               label="SKYQUAKE physical fit (seismic only)")
    ax[0].axhline(line.params["speed_kms"], color="C0", ls="--",
                  label="SKYQUAKE straight line (seismic only)")
    ax[0].set_xlim(-500, 700)
    ax[0].set_ylim(0, 13)
    ax[0].set_xlabel("along-track distance from the network (km)")
    ax[0].set_ylabel("speed (km/s)")
    ax[0].legend(fontsize=8)
    ax[1].plot((s1 - s1[k_ref]) / 1e3, q1[:, 2] / 1e3, ".", color="0.6", ms=3,
               label="camera F1")
    ax[1].plot((sm - sm[np.argmin(np.abs(tt))]) / 1e3, enu[:, 2] / 1e3, "-", color="C3",
               lw=2, label="SKYQUAKE physical fit")
    ax[1].set_xlim(-500, 700)
    ax[1].set_ylim(20, 110)
    ax[1].set_xlabel("along-track distance from the network (km)")
    ax[1].set_ylabel("height (km)")
    ax[1].legend(fontsize=8)
    fig.suptitle("Falcon 9 upper stage: physical seismic fit vs camera fragment F1")
    fig.tight_layout()
    fig.savefig(HERE / "ballistic_vs_camera.png", dpi=110)
    _ = v_model_at_f1


if __name__ == "__main__":
    main(sys.argv[1])
