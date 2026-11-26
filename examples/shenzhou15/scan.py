"""Exploratory scan: STA/LTA of every station, drawn as a record section.

Stations are sorted by their position projected on the Santa Barbara -> Las Vegas line
(the published direction of travel). A sonic boom shows up as a coherent, sloping front.

    python examples/shenzhou15/scan.py   -> scan.png, scan_triggers.csv
"""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from obspy import Inventory, UTCDateTime, read, read_inventory
from obspy.signal.trigger import classic_sta_lta

HERE = Path(__file__).parent
RAW = HERE / "raw"
T0 = UTCDateTime("2024-04-02T08:35:00")
SB = (34.42, -119.70)          # Santa Barbara
LV = (36.17, -115.14)          # Las Vegas
FMIN, FMAX = 4.0, 20.0
STA, LTA = 0.5, 20.0


def local_xy(lat, lon):
    """Equirectangular km around Santa Barbara (only used to sort and plot)."""
    k = 111.195
    return (lon - SB[1]) * k * np.cos(np.radians(SB[0])), (lat - SB[0]) * k


def main(win=(0.0, 1500.0), out="scan") -> None:
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
        if tr.stats.npts < 600 * tr.stats.sampling_rate:
            continue
        net, sta = tr.stats.network, tr.stats.station
        try:
            c = inv.select(network=net, station=sta)[0][0]
        except IndexError:
            continue
        tr.detrend("demean")
        tr.taper(0.01)
        tr.filter("bandpass", freqmin=FMIN, freqmax=FMAX, corners=4, zerophase=True)
        df = tr.stats.sampling_rate
        cft = classic_sta_lta(tr.data, int(STA * df), int(LTA * df))
        cft[: int(LTA * df)] = 0
        t = (tr.stats.starttime - T0) + np.arange(tr.stats.npts) / df
        x, y = local_xy(c.latitude, c.longitude)
        along = x * axis[0] + y * axis[1]
        cross = -x * axis[1] + y * axis[0]
        inwin = (t >= win[0]) & (t <= win[1])
        i = int(np.flatnonzero(inwin)[np.argmax(cft[inwin])])
        rows.append(dict(id=f"{net}.{sta}", lat=c.latitude, lon=c.longitude,
                         elev=c.elevation, along_km=along, cross_km=cross,
                         tmax_s=t[i], cftmax=float(cft[i]), t=t[::10], cft=cft[::10]))

    rows.sort(key=lambda r: r["along_km"])
    fig, ax = plt.subplots(figsize=(11, 14))
    for r in rows:
        y = r["along_km"]
        m = (r["t"] >= win[0]) & (r["t"] <= win[1])
        ax.plot(r["t"][m], y + 6 * np.clip(r["cft"][m], 0, 12) / 12, lw=0.4, color="k")
        if r["cftmax"] > 6:
            ax.plot(r["tmax_s"], y, ".", color="r", ms=4)
    ax.set_xlabel("seconds after 2024-04-02 08:35:00 UTC")
    ax.set_ylabel("station position along Santa Barbara -> Las Vegas (km)")
    ax.set_title(f"STA/LTA {FMIN}-{FMAX} Hz, {len(rows)} stations (red: max STA/LTA > 6)")
    fig.tight_layout()
    fig.savefig(HERE / f"{out}.png", dpi=110)

    with open(HERE / f"{out}_triggers.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "lat", "lon", "elev_m", "along_km", "cross_km", "tmax_s", "cftmax"])
        for r in rows:
            w.writerow([r["id"], r["lat"], r["lon"], r["elev"], round(r["along_km"], 1),
                        round(r["cross_km"], 1), round(r["tmax_s"], 2), round(r["cftmax"], 1)])
    print(f"{len(rows)} stations scanned; {sum(r['cftmax'] > 6 for r in rows)} with STA/LTA > 6")


if __name__ == "__main__":
    import sys

    if len(sys.argv) == 3:
        main((float(sys.argv[1]), float(sys.argv[2])), out="zoom")
    else:
        main()
