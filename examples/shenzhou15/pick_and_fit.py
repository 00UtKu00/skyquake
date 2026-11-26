"""Shenzhou-15 orbital module: pick the sonic boom and fit a trajectory with SKYQUAKE.

Iterative predict -> pick -> fit:
  0. a coarse V-shaped moveout read off section_cross.png (apex ~690 s, ~2.1 s/km)
  1. pick the STA/LTA onset inside a window around the prediction
  2. fit the Mach-cone model (sound speed estimated)
  3. re-pick in a narrower window around the fitted prediction, refit

Outputs: picks.csv, fit.txt, track.geojson, picks_on_section.png, map.png
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from obspy import Inventory, UTCDateTime, read, read_inventory
from obspy.signal.trigger import classic_sta_lta

from skyquake import LocalFrame, arrival_times, fit_trajectory
from skyquake.io import write_picks
from skyquake.track import corridor, geographic_state, to_geojson, track_points

HERE = Path(__file__).parent
RAW = HERE / "raw"
T0 = UTCDateTime("2024-04-02T08:35:00")
SB, LV = (34.42, -119.70), (36.17, -115.14)
FMIN, FMAX = 2.0, 20.0
STA_S, LTA_S, ON = 0.3, 10.0, 3.0
QUAKE = ("2024-04-02T08:38:28", 33.836, -117.657)  # USGS ci40703496, M1.2; excluded by window


def load():
    inv = Inventory()
    for f in RAW.glob("*_inventory.xml"):
        inv += read_inventory(str(f))
    out = []
    for f in sorted(RAW.glob("*.mseed")):
        st = read(str(f))
        st.merge(fill_value=0)
        tr = st[0]
        try:
            c = inv.select(network=tr.stats.network, station=tr.stats.station)[0][0]
        except IndexError:
            continue
        tr.trim(T0 + 450, T0 + 1400)
        if tr.stats.npts < 1000:
            continue
        tr.detrend("demean")
        tr.taper(0.02)
        tr.filter("bandpass", freqmin=FMIN, freqmax=FMAX, corners=4, zerophase=True)
        df = tr.stats.sampling_rate
        cft = classic_sta_lta(tr.data, int(STA_S * df), int(LTA_S * df))
        cft[: int(LTA_S * df)] = 0
        t = (tr.stats.starttime - T0) + np.arange(tr.stats.npts) / df
        out.append(dict(id=f"{tr.stats.network}.{tr.stats.station}", lat=c.latitude,
                        lon=c.longitude, elev=c.elevation, t=t, cft=cft,
                        data=tr.data.astype(float)))
    return out


def pick(s, t_pred, half):
    """Onset of the strongest STA/LTA peak within t_pred +/- half. Returns (t, snr) or None."""
    m = (s["t"] >= t_pred - half) & (s["t"] <= t_pred + half)
    if m.sum() < 10:
        return None
    idx = np.flatnonzero(m)
    i_max = idx[np.argmax(s["cft"][idx])]
    if s["cft"][i_max] < ON:
        return None
    i = i_max
    while i > idx[0] and s["cft"][i] > ON:   # walk back to the threshold crossing
        i -= 1
    df = 1.0 / (s["t"][1] - s["t"][0])
    pre = s["data"][max(0, i - int(10 * df)):i]
    post = s["data"][i:i + int(3 * df)]
    if len(pre) < 10 or len(post) < 10:
        return None
    snr = float(np.abs(post).max() / (np.std(pre) * 3 + 1e-12))
    return float(s["t"][i]), snr


def main():
    stations = load()
    frame = LocalFrame.centred_on([s["lat"] for s in stations], [s["lon"] for s in stations])
    enu = frame.to_enu([s["lat"] for s in stations], [s["lon"] for s in stations],
                       [s["elev"] for s in stations])

    # step 0: coarse V model in the Santa Barbara -> Las Vegas coordinates
    k = 111.195
    def cross_km(lat, lon):
        x = (lon - SB[1]) * k * np.cos(np.radians(SB[0]))
        y = (lat - SB[0]) * k
        ax = np.array([(LV[1] - SB[1]) * k * np.cos(np.radians(SB[0])), (LV[0] - SB[0]) * k])
        ax /= np.hypot(*ax)
        return -x * ax[1] + y * ax[0]
    pred = np.array([690.0 + 2.1 * abs(cross_km(s["lat"], s["lon"]) + 5.0) for s in stations])

    fit = None
    for it, half in enumerate((30.0, 10.0, 5.0)):
        picks = [pick(s, p, half) for s, p in zip(stations, pred, strict=True)]
        ok = np.array([p is not None and p[1] >= 2.0 for p in picks])
        t = np.array([p[0] if p else np.nan for p in picks])
        snr = np.array([p[1] if p else np.nan for p in picks])
        sig = np.where(ok, np.clip(2.0 / np.nan_to_num(snr, nan=1.0), 0.15, 2.0), 1.0)
        fit = fit_trajectory(enu[ok], t[ok], sig[ok], fit_sound_speed=True,
                             azimuth_guess_deg=None if it == 0 else fit.trajectory.azimuth_deg)
        pred = arrival_times(fit.trajectory, enu, fit.sound_speed)
        print(f"iteration {it}: window +/-{half:.0f} s, {ok.sum()} picks, "
              f"{fit.n_used} used, rms {fit.rms_s:.2f} s")

    used_ids = np.flatnonzero(ok)[fit.inliers]

    # model comparison and honest (jackknife) speed uncertainty on the final picks
    from skyquake.atmosphere import LayeredAtmosphere
    from skyquake.invert import jackknife

    kw = dict(fit_sound_speed=True, azimuth_guess_deg=fit.trajectory.azimuth_deg)
    layered = fit_trajectory(enu[ok], t[ok], sig[ok], atmosphere=LayeredAtmosphere(), **kw)
    jk = jackknife(enu[ok], t[ok], sig[ok], **kw)
    rows = []
    for j in used_ids:
        s = stations[j]
        rows.append({"station": s["id"], "lat": s["lat"], "lon": s["lon"], "elev_m": s["elev"],
                     "arrival_utc": (T0 + t[j]).isoformat() + "Z",
                     "sigma_s": round(float(sig[j]), 3), "snr": round(float(snr[j]), 1)})
    write_picks(HERE / "picks.csv", rows)

    g = geographic_state(fit.trajectory, frame)
    c = corridor(fit, frame)
    t_ref = T0 + fit.params["t0_s"]
    lines = [
        "SKYQUAKE - Shenzhou-15 orbital module, 2024-04-02",
        f"frame centre {frame.lat0:.4f}, {frame.lon0:.4f}; picks within window, "
        f"{fit.n_used} stations used",
        fit.summary(),
        "",
        f"closest approach to the network centre at {t_ref.isoformat()}Z",
        f"  position lat {g['lat']:.4f}, lon {g['lon']:.4f}, height {g['height_m'] / 1e3:.1f} km",
        f"  heading {g['azimuth_deg']:.2f} deg true, descent {g['fpa_deg']:.2f} deg, "
        f"speed {g['speed_mps'] / 1e3:.2f} km/s "
        f"(Mach {g['speed_mps'] / fit.sound_speed:.1f} vs fitted c)",
        f"  fitted path-averaged sound speed {fit.sound_speed:.1f} m/s",
        f"  1-sigma: {c['sigma_cross_km']:.2f} km cross-track, "
        f"{c['sigma_height_km']:.2f} km height",
        "",
        "jackknife 1-sigma (10 groups) vs formal:",
        *[f"  {k:<16} +/- {v:.3f}   (formal {fit.sigma[k]:.3f})" for k, v in jk.items()],
        "",
        "layered atmosphere (US Standard 1976, no winds):",
        f"  speed {layered.params['speed_kms']:.2f} km/s, height "
        f"{layered.params['height_km']:.1f} km, azimuth {layered.params['azimuth_deg']:.2f}, "
        f"rms {layered.rms_s:.3f} s",
    ]
    txt = "\n".join(lines)
    print(txt)
    (HERE / "fit.txt").write_text(txt + "\n", encoding="utf-8")
    (HERE / "track.geojson").write_text(json.dumps(to_geojson(fit, frame), indent=1))

    # figures
    res = fit.residuals_s[fit.inliers]
    cr = np.array([cross_km(stations[j]["lat"], stations[j]["lon"]) for j in used_ids])
    fig, ax = plt.subplots(1, 2, figsize=(14, 7))
    allc = np.array([cross_km(s["lat"], s["lon"]) for s in stations])
    ax[0].plot(t[ok], allc[ok], ".", color="0.7", ms=4, label="picked, rejected")
    ax[0].plot(t[used_ids], cr, ".", color="C3", ms=5, label="picked, used")
    order = np.argsort(allc)
    ax[0].plot(pred[order], allc[order], "-", color="C0", lw=0.6, label="model")
    ax[0].set_xlabel("s after 08:35:00 UTC")
    ax[0].set_ylabel("cross-track km (Santa Barbara -> Las Vegas line)")
    ax[0].legend()
    ax[0].set_title(f"sonic-boom picks vs model (rms {fit.rms_s:.2f} s)")
    ax[1].hist(res, bins=30, color="C3")
    ax[1].set_xlabel("residual (s)")
    ax[1].set_title("residuals of used picks")
    fig.tight_layout()
    fig.savefig(HERE / "picks_on_section.png", dpi=110)

    fig, ax = plt.subplots(figsize=(10, 8))
    lat_all = [s["lat"] for s in stations]
    lon_all = [s["lon"] for s in stations]
    ax.plot(lon_all, lat_all, "^", color="0.75", ms=4, label="stations")
    sc = ax.scatter([stations[j]["lon"] for j in used_ids], [stations[j]["lat"] for j in used_ids],
                    c=t[used_ids], cmap="viridis", s=18, zorder=3, label="boom picks")
    plt.colorbar(sc, label="arrival (s after 08:35:00 UTC)")
    tp = track_points(fit.trajectory, frame, before_km=350, after_km=350)
    ax.plot(tp[:, 1], tp[:, 0], "-", color="C3", lw=2, label="SKYQUAKE track")
    for draw in fit.sample(60):
        tq = track_points(draw, frame, before_km=350, after_km=350, n=10)
        ax.plot(tq[:, 1], tq[:, 0], "-", color="C3", lw=0.3, alpha=0.3)
    ax.plot([SB[1], LV[1]], [SB[0], LV[0]], "--", color="k", lw=1,
            label="Santa Barbara -> Las Vegas (published direction)")
    ax.plot(QUAKE[2], QUAKE[1], "*", color="orange", ms=12, label="M1.2 quake (excluded)")
    places = {"Santa Barbara": SB, "Las Vegas": LV, "Los Angeles": (34.05, -118.24)}
    for name, (la, lo) in places.items():
        ax.annotate(name, (lo, la), fontsize=9)
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")
    ax.set_aspect(1 / np.cos(np.radians(35)))
    ax.legend(loc="lower left", fontsize=8)
    ax.set_title("Shenzhou-15 orbital module re-entry, 2024-04-02, from seismic sonic booms")
    fig.tight_layout()
    fig.savefig(HERE / "map.png", dpi=110)


if __name__ == "__main__":
    main()
