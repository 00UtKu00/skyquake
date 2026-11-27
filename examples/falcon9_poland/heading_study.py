# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 SKYQUAKE contributors
"""Heading accuracy study: SKYQUAKE fit variants vs the camera-derived heading.

The optical reference line is a straight 3-D line fitted (least squares) to the camera-
triangulated fragment positions within +/-300 km along-track of the seismic network, in
the same local frame as the seismic fit. See compare_optical.py for how to get the data.

    python examples/falcon9_poland/heading_study.py path/to/ref/fragments
"""

import sys
from pathlib import Path

import h5py
import numpy as np

from skyquake.geo import LocalFrame
from skyquake.invert import decluster_sigma, fit_trajectory, jackknife
from skyquake.io import read_picks

HERE = Path(__file__).parent


def optical_line(folder, frame, centre_enu, half_km=300.0):
    pts = []
    for f in sorted(Path(folder).glob("*.h5")):
        with h5py.File(f, "r") as h:
            lat, lon, alt = h["geoloc"][()]
        pts.append(frame.to_enu(lat, lon, alt))
    pts = np.array(pts)
    d = pts - centre_enu
    near = np.hypot(d[:, 0], d[:, 1]) < half_km * 1e3
    p = pts[near]
    c = p.mean(axis=0)
    u = np.linalg.svd(p - c)[2][0]
    if u[0] < 0:                      # point eastward (direction of travel)
        u = -u
    az = np.degrees(np.arctan2(u[0], u[1])) % 360
    return az, c, near.sum()


def heading_at(fit, frame, ref_enu):
    """Fit heading expressed in the same frame (no conversion needed: same frame)."""
    return fit.params["azimuth_deg"]


def main(folder):
    p = read_picks(HERE / "picks.csv")
    frame = LocalFrame.centred_on(p.lat, p.lon)
    st = frame.to_enu(p.lat, p.lon, p.elev_m)
    base = dict(fit_sound_speed=True, azimuth_guess_deg=95.0)
    az_opt, c_opt, n_opt = optical_line(folder, frame, np.zeros(3))
    print(f"{len(p)} picks; optical reference heading {az_opt:.2f} deg "
          f"(line through {n_opt} camera positions within 300 km)\n")

    lon_e = p.lon > 8.0
    variants = {
        "plain (all picks)": (np.ones(len(p), bool), p.sigma_s),
        "declustered (20 km)": (np.ones(len(p), bool), decluster_sigma(st, p.sigma_s, 20e3)),
        "Groningen removed": (~((p.lat > 52.9) & (p.lon > 6.0) & (p.lon < 7.4)), p.sigma_s),
    }
    print(f"stations east of 8 E: {lon_e.sum()}")
    print(f"{'variant':<24} {'n':>4} {'heading':>8} {'error':>7} {'formal':>7} "
          f"{'jk-rand':>8} {'jk-spat':>8} {'err/spat':>9}")
    rows = []
    for name, (m, sig) in variants.items():
        if m.sum() < 15:
            print(f"{name:<24} {m.sum():>4}  too few picks")
            continue
        f = fit_trajectory(st[m], p.t_s[m], sig[m], **base)
        jr = jackknife(st[m], p.t_s[m], sig[m], blocks="random", **base)
        js = jackknife(st[m], p.t_s[m], sig[m], blocks="spatial", **base)
        err = (f.params["azimuth_deg"] - az_opt + 180) % 360 - 180
        rows.append((name, m.sum(), f.params["azimuth_deg"], err, f.sigma["azimuth_deg"],
                     jr["azimuth_deg"], js["azimuth_deg"]))
        print(f"{name:<24} {m.sum():>4} {f.params['azimuth_deg']:>8.2f} {err:>+7.2f} "
              f"{f.sigma['azimuth_deg']:>7.2f} {jr['azimuth_deg']:>8.2f} "
              f"{js['azimuth_deg']:>8.2f} {abs(err) / js['azimuth_deg']:>9.1f}")
    (HERE / "heading_study.txt").write_text(
        f"optical heading {az_opt:.2f}\n" + "\n".join(
            f"{r[0]}\tn={r[1]}\taz={r[2]:.2f}\terr={r[3]:+.2f}\tformal={r[4]:.2f}"
            f"\tjk_random={r[5]:.2f}\tjk_spatial={r[6]:.2f}" for r in rows) + "\n",
        encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv[1])
