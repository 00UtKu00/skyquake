"""Record section sorted by signed distance from an approximate track (expect a V).

    python examples/falcon9_poland/section.py [fmin] [t_start t_end]
"""

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime, read, read_inventory

HERE = Path(__file__).parent
RAW = HERE / "raw"
T0 = UTCDateTime("2025-02-19T03:42:00")
A, B = (54.0, -4.0), (52.4, 16.8)      # Irish Sea -> Poznan (rough, for sorting only)
FMIN = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
FMAX = 15.0
WIN = (float(sys.argv[2]), float(sys.argv[3])) if len(sys.argv) > 3 else (60.0, 1500.0)


def local_xy(lat, lon):
    k = 111.195
    return (lon - A[1]) * k * np.cos(np.radians(53.2)), (lat - A[0]) * k


def main():
    inv = read_inventory(str(RAW / "inventory.xml"))
    ax = np.array(local_xy(*B))
    ax /= np.hypot(*ax)
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
        if tr.stats.npts < 100 or tr.stats.sampling_rate < 2.5 * FMAX:
            continue
        tr.detrend("demean")
        tr.taper(0.02)
        tr.filter("bandpass", freqmin=FMIN, freqmax=FMAX, corners=4, zerophase=True)
        df = tr.stats.sampling_rate
        env = np.convolve(np.abs(tr.data), np.ones(int(df)) / int(df), mode="same")
        env /= np.median(env) or 1.0
        t = (tr.stats.starttime - T0) + np.arange(tr.stats.npts) / df
        x, y = local_xy(c.latitude, c.longitude)
        rows.append(dict(id=tr.id, lat=c.latitude, lon=c.longitude, elev=c.elevation,
                         along=x * ax[0] + y * ax[1], cross=-x * ax[1] + y * ax[0],
                         t=t[::5], env=env[::5]))
    rows.sort(key=lambda r: r["cross"])
    fig, axp = plt.subplots(1, 2, figsize=(16, 14), gridspec_kw={"width_ratios": [3, 1]})
    for r in rows:
        m = (r["t"] >= WIN[0]) & (r["t"] <= WIN[1])
        axp[0].plot(r["t"][m], r["cross"] + 1.2 * np.log10(np.clip(r["env"][m], 1, 1e4)),
                    lw=0.3, color="k")
    axp[0].set_xlabel("seconds after 2025-02-19 03:42:00 UTC")
    axp[0].set_ylabel("signed distance from Irish Sea -> Poznan line (km, + = north)")
    axp[0].set_title(f"log envelope {FMIN}-{FMAX} Hz, {len(rows)} stations")
    axp[1].scatter([r["lon"] for r in rows], [r["lat"] for r in rows], s=4,
                   c=[r["cross"] for r in rows], cmap="coolwarm")
    axp[1].plot([A[1], B[1]], [A[0], B[0]], "k--")
    axp[1].set_xlim(2, 20)
    axp[1].set_ylim(50, 57)
    axp[1].set_title("stations (colour = cross-track)")
    fig.tight_layout()
    fig.savefig(HERE / "section.png", dpi=100)
    with open(HERE / "stations_xy.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "lat", "lon", "elev_m", "along_km", "cross_km"])
        for r in rows:
            w.writerow([r["id"], r["lat"], r["lon"], r["elev"], round(r["along"], 1),
                        round(r["cross"], 1)])
    print(f"{len(rows)} stations; cross range {rows[0]['cross']:.0f}..{rows[-1]['cross']:.0f} km")


if __name__ == "__main__":
    main()
