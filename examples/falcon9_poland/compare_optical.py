# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 SKYQUAKE contributors
"""Compare the SKYQUAKE track with the camera-triangulated fragment positions of
Vierinen et al. 2026 (github.com/jvierine/falcon9, folder fragments/, doi:10.5281/zenodo.20070800).

The optical data are not redistributed here. Get them with:
    git clone --depth 1 --filter=blob:none --sparse https://github.com/jvierine/falcon9.git ref
    cd ref && git sparse-checkout set fragments
then run:
    python examples/falcon9_poland/compare_optical.py path/to/ref/fragments
"""

import sys
from datetime import UTC, datetime
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from skyquake.geo import LocalFrame
from skyquake.invert import fit_trajectory
from skyquake.io import read_picks
from skyquake.track import track_points

HERE = Path(__file__).parent


def load_optical(folder):
    rows = []
    for f in sorted(Path(folder).glob("*.h5")):
        with h5py.File(f, "r") as h:
            lat, lon, alt = h["geoloc"][()]
            rows.append((f.name.split("_")[0], float(h["time"][()]), lat, lon, alt))
    frag = np.array([r[0] for r in rows])
    t, lat, lon, alt = (np.array([r[i] for r in rows]) for i in range(1, 5))
    return frag, t, lat, lon, alt


def main(folder):
    frag, tu, lat, lon, alt = load_optical(folder)
    p = read_picks(HERE / "picks.csv")
    frame = LocalFrame.centred_on(p.lat, p.lon)
    st = frame.to_enu(p.lat, p.lon, p.elev_m)
    fit = fit_trajectory(st, p.t_s, p.sigma_s, fit_sound_speed=True, azimuth_guess_deg=95.0)
    tr = fit.trajectory
    epoch = p.epoch.timestamp()            # picks.csv times are relative to the first pick

    enu = frame.to_enu(lat, lon, alt)
    d = enu - tr.origin
    s = d @ tr.direction                   # along-track metres from the reference point
    perp = d - np.outer(s, tr.direction)
    right = np.array([np.cos(np.radians(tr.azimuth_deg)), -np.sin(np.radians(tr.azimuth_deg)),
                      0.0])
    cross = perp @ right / 1e3             # + = right of travel (south for an eastward track)
    vert = perp[:, 2] / 1e3                # + = optical point above our line
    t_line = epoch + tr.t0_s + s / tr.speed_mps
    dt = tu - t_line                       # + = optical fragment later than our line

    near = np.abs(s) < 150e3               # within 150 km along-track of the network centre
    print(f"{len(tu)} optical positions, {len(set(frag))} fragments; "
          f"{near.sum()} within +/-150 km along-track of the seismic reference point")
    for name, m in (("near network (+/-150 km)", near), ("all", np.ones_like(near))):
        print(f"\n{name}:")
        lo, hi = np.percentile(cross[m], [16, 84])
        print(f"  cross-track offset of optics from our line: median {np.median(cross[m]):+.1f} km"
              f"  (16-84%: {lo:+.1f} .. {hi:+.1f})")
        print(f"  height of optics above our line:            median {np.median(vert[m]):+.1f} km")
        print(f"  time of optics minus our line:              median {np.median(dt[m]):+.1f} s")
    # optical speed near the network: per fragment, finite differences
    v = []
    for f in sorted(set(frag)):
        k = np.flatnonzero((frag == f) & near)
        if k.size > 3:
            o = np.argsort(tu[k])
            e = enu[k][o]
            v.append(np.median(np.linalg.norm(np.diff(e, axis=0), axis=1) / np.diff(tu[k][o])))
    if v:
        print(f"\noptical speed near the network: median {np.median(v) / 1e3:.2f} km/s over "
              f"{len(v)} fragments;  SKYQUAKE: {tr.speed_mps / 1e3:.2f} km/s")

    t0 = datetime.fromtimestamp(tu.min(), UTC)
    t1 = datetime.fromtimestamp(tu.max(), UTC)
    print(f"optical coverage {t0:%H:%M:%S} .. {t1:%H:%M:%S} UTC, "
          f"lon {lon.min():.1f}..{lon.max():.1f}")

    fig, ax = plt.subplots(1, 2, figsize=(15, 6))
    sc = ax[0].scatter(lon, lat, c=alt / 1e3, s=4, cmap="plasma", label="optical fragments")
    plt.colorbar(sc, ax=ax[0], label="optical altitude (km)")
    tp = track_points(tr, frame, before_km=900, after_km=900, min_height_m=0)
    ax[0].plot(tp[:, 1], tp[:, 0], "-", color="C2", lw=1.5, label="SKYQUAKE (seismic)")
    ax[0].plot(p.lon, p.lat, "^", color="0.6", ms=2, label="seismic stations used")
    ax[0].set_xlabel("longitude")
    ax[0].set_ylabel("latitude")
    ax[0].legend(fontsize=8)
    ax[0].set_title("ground tracks: optical (Vierinen et al. 2026) vs seismic")
    ax[1].scatter(s / 1e3, cross, c=alt / 1e3, s=4, cmap="plasma")
    ax[1].axhline(0, color="C2")
    ax[1].set_xlabel("along-track distance from seismic reference point (km)")
    ax[1].set_ylabel("cross-track offset of optics (km, + = south)")
    ax[1].set_title("offset of optical fragments from the seismic line")
    fig.tight_layout()
    fig.savefig(HERE / "optical_vs_seismic.png", dpi=110)


if __name__ == "__main__":
    main(sys.argv[1])
