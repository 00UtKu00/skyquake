# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 SKYQUAKE contributors
"""Record section sorted by signed cross-track distance (expect a V around the track)."""

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from obspy import Inventory, UTCDateTime, read, read_inventory

HERE = Path(__file__).parent
RAW = HERE / "raw"
T0 = UTCDateTime("2024-04-02T08:35:00")
SB, LV = (34.42, -119.70), (36.17, -115.14)
FMIN, FMAX = float(sys.argv[1]) if len(sys.argv) > 1 else 2.0, 20.0
WIN = (500.0, 1300.0)


def local_xy(lat, lon):
    k = 111.195
    return (lon - SB[1]) * k * np.cos(np.radians(SB[0])), (lat - SB[0]) * k


def main():
    inv = Inventory()
    for f in RAW.glob("*_inventory.xml"):
        inv += read_inventory(str(f))
    ax_e, ax_n = local_xy(*LV)
    axis = np.array([ax_e, ax_n]) / np.hypot(ax_e, ax_n)
    rows = []
    for f in sorted(RAW.glob("*.mseed")):
        st = read(str(f))
        st.merge(fill_value=0)
        tr = st[0]
        try:
            c = inv.select(network=tr.stats.network, station=tr.stats.station)[0][0]
        except IndexError:
            continue
        tr.trim(T0 + WIN[0] - 60, T0 + WIN[1])
        if tr.stats.npts < 100:
            continue
        tr.detrend("demean")
        tr.taper(0.02)
        tr.filter("bandpass", freqmin=FMIN, freqmax=FMAX, corners=4, zerophase=True)
        env = np.abs(tr.data)
        df = tr.stats.sampling_rate
        k = int(df)  # 1 s smoothing
        env = np.convolve(env, np.ones(k) / k, mode="same")
        noise = np.median(env) or 1.0
        t = (tr.stats.starttime - T0) + np.arange(tr.stats.npts) / df
        x, y = local_xy(c.latitude, c.longitude)
        rows.append(dict(id=tr.id, lat=c.latitude, lon=c.longitude, elev=c.elevation,
                         along=x * axis[0] + y * axis[1], cross=-x * axis[1] + y * axis[0],
                         t=t[::5], env=(env / noise)[::5]))
    rows.sort(key=lambda r: r["cross"])
    fig, ax = plt.subplots(figsize=(11, 14))
    for r in rows:
        m = (r["t"] >= WIN[0]) & (r["t"] <= WIN[1])
        ax.plot(r["t"][m], r["cross"] + 1.5 * np.log10(np.clip(r["env"][m], 1, 1e4)),
                lw=0.4, color="k")
    ax.set_xlabel("seconds after 2024-04-02 08:35:00 UTC")
    ax.set_ylabel("signed distance from Santa Barbara -> Las Vegas line (km, + = north-west)")
    ax.set_title(f"log envelope {FMIN}-{FMAX} Hz, {len(rows)} stations")
    fig.tight_layout()
    fig.savefig(HERE / "section_cross.png", dpi=110)
    with open(HERE / "stations_xy.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "lat", "lon", "elev_m", "along_km", "cross_km"])
        for r in rows:
            w.writerow([r["id"], r["lat"], r["lon"], r["elev"], round(r["along"], 1),
                        round(r["cross"], 1)])


if __name__ == "__main__":
    main()
