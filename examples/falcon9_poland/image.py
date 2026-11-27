# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 SKYQUAKE contributors
"""Time x cross-track image of log envelopes (handles dense networks like Groningen)."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime, read, read_inventory

HERE = Path(__file__).parent
RAW = HERE / "raw"
T0 = UTCDateTime("2025-02-19T03:42:00")
A, B = (54.0, -4.0), (52.4, 16.8)
FMIN, FMAX = 2.0, 15.0
TW = (300.0, 900.0)
DT = 1.0
XBIN = 4.0


def local_xy(lat, lon):
    k = 111.195
    return (lon - A[1]) * k * np.cos(np.radians(53.2)), (lat - A[0]) * k


def main():
    inv = read_inventory(str(RAW / "inventory.xml"))
    ax = np.array(local_xy(*B))
    ax /= np.hypot(*ax)
    tgrid = np.arange(TW[0], TW[1], DT)
    xs, envs, alongs = [], [], []
    for f in sorted(RAW.glob("*.mseed")):
        st = read(str(f))
        st.merge(fill_value=0)
        tr = st[0]
        try:
            c = inv.select(network=tr.stats.network, station=tr.stats.station)[0][0]
        except IndexError:
            continue
        tr.trim(T0 + TW[0] - 60, T0 + TW[1] + 10)
        if tr.stats.npts < 100 or tr.stats.sampling_rate < 2.5 * FMAX:
            continue
        tr.detrend("demean")
        tr.taper(0.02)
        tr.filter("bandpass", freqmin=FMIN, freqmax=FMAX, corners=4, zerophase=True)
        df = tr.stats.sampling_rate
        env = np.convolve(np.abs(tr.data), np.ones(int(df)) / int(df), mode="same")
        t = (tr.stats.starttime - T0) + np.arange(tr.stats.npts) / df
        e = np.interp(tgrid, t, env)
        e = np.log10(np.clip(e / (np.median(e) or 1.0), 1e-3, None))
        x, y = local_xy(c.latitude, c.longitude)
        xs.append(-x * ax[1] + y * ax[0])
        alongs.append(x * ax[0] + y * ax[1])
        envs.append(e)
    xs, envs = np.array(xs), np.array(envs)
    bins = np.arange(-320, 440, XBIN)
    img = np.full((bins.size - 1, tgrid.size), np.nan)
    count = np.zeros(bins.size - 1, int)
    for i in range(bins.size - 1):
        m = (xs >= bins[i]) & (xs < bins[i + 1])
        count[i] = m.sum()
        if m.any():
            img[i] = np.median(envs[m], axis=0)
    fig, axp = plt.subplots(figsize=(13, 9))
    im = axp.imshow(img, aspect="auto", origin="lower", cmap="magma", vmin=0, vmax=1.2,
                    extent=[tgrid[0], tgrid[-1], bins[0], bins[-1]], interpolation="nearest")
    plt.colorbar(im, label="log10 envelope / median (per-bin median over stations)")
    axp.set_xlabel("seconds after 2025-02-19 03:42:00 UTC")
    axp.set_ylabel("signed cross-track distance (km, + = north of Irish Sea -> Poznan line)")
    axp.set_title(f"{len(xs)} stations, {FMIN}-{FMAX} Hz, {XBIN:.0f} km bins")
    axp.set_ylim(-320, 160)
    fig.tight_layout()
    fig.savefig(HERE / "image.png", dpi=100)
    print(f"{len(xs)} stations; bins with data: {(count > 0).sum()}")


if __name__ == "__main__":
    main()
