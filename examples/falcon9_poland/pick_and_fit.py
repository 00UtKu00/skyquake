"""Falcon 9 upper stage (Starlink 11-4), 2025-02-19: pick the boom and fit with SKYQUAKE,
then compare with the optical trajectory (Vierinen et al. 2026) and the debris sites.

Outputs: picks.csv, fit.txt, map.png, picks_on_section.png, track.geojson
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime

from skyquake.atmosphere import LayeredAtmosphere
from skyquake.invert import fit_trajectory, jackknife
from skyquake.io import write_picks
from skyquake.pipeline import iterative_pick_fit, load_traces
from skyquake.track import corridor, geographic_state, to_geojson, track_points

HERE = Path(__file__).parent
T0 = UTCDateTime("2025-02-19T03:42:00")
A, B = (54.0, -4.0), (52.4, 16.8)   # rough Irish Sea -> Poznan line, for the first guess
# debris recovery towns (approximate town-centre coordinates; POLSA/press reports)
DEBRIS = {
    "Komorniki": (52.335, 16.810), "Wiry": (52.325, 16.860), "Sliwno": (52.590, 16.320),
    "Sedziny (Szamotuly)": (52.680, 16.450), "Miedzychod": (52.600, 15.890),
}


def cross_km(lat, lon):
    k = 111.195
    x = (lon - A[1]) * k * np.cos(np.radians(53.2))
    y = (lat - A[0]) * k
    ax = np.array([(B[1] - A[1]) * k * np.cos(np.radians(53.2)), (B[0] - A[0]) * k])
    ax /= np.hypot(*ax)
    return -x * ax[1] + y * ax[0]


def first_guess(lat, lon):
    # V read off image.png: apex ~438 s near the line, hyperbola for h = 65 km, c = 0.32 km/s
    r = np.abs(cross_km(lat, lon) - 5.0)
    return 438.0 + (np.sqrt(65.0**2 + r**2) - 65.0) / 0.32


def main():
    stations = load_traces(HERE / "raw", T0, t_range=(250.0, 1200.0), band=(2.0, 15.0))
    print(f"{len(stations)} stations loaded")
    pf = iterative_pick_fit(stations, first_guess)
    fit, frame = pf.fit, pf.frame
    ok = pf.picked
    kw = dict(fit_sound_speed=True, azimuth_guess_deg=fit.trajectory.azimuth_deg)
    layered = fit_trajectory(pf.enu[ok], pf.t[ok], pf.sigma[ok], atmosphere=LayeredAtmosphere(),
                             **kw)
    jk = jackknife(pf.enu[ok], pf.t[ok], pf.sigma[ok], **kw)

    rows = [{"station": stations[j].id, "lat": stations[j].lat, "lon": stations[j].lon,
             "elev_m": stations[j].elev, "arrival_utc": (T0 + pf.t[j]).isoformat() + "Z",
             "sigma_s": round(float(pf.sigma[j]), 3), "snr": round(float(pf.snr[j]), 1)}
            for j in pf.used]
    write_picks(HERE / "picks.csv", rows)

    g = geographic_state(fit.trajectory, frame)
    c = corridor(fit, frame)
    lines = [
        "SKYQUAKE - Falcon 9 upper stage (Starlink 11-4), 2025-02-19",
        f"{len(stations)} stations, {ok.sum()} picks, {fit.n_used} used",
        fit.summary(), "",
        f"closest approach to the network centre at {(T0 + fit.params['t0_s']).isoformat()}Z",
        f"  position lat {g['lat']:.4f}, lon {g['lon']:.4f}, height {g['height_m'] / 1e3:.1f} km",
        f"  heading {g['azimuth_deg']:.2f} deg true, descent {g['fpa_deg']:.2f} deg, "
        f"speed {g['speed_mps'] / 1e3:.2f} km/s",
        f"  1-sigma: {c['sigma_cross_km']:.2f} km cross-track, "
        f"{c['sigma_height_km']:.2f} km height", "",
        "jackknife 1-sigma (10 groups) vs formal:",
        *[f"  {k:<16} +/- {v:.3f}   (formal {fit.sigma[k]:.3f})" for k, v in jk.items()], "",
        "layered atmosphere (US Standard 1976, no winds):",
        f"  speed {layered.params['speed_kms']:.2f} km/s, height "
        f"{layered.params['height_km']:.1f} km, azimuth {layered.params['azimuth_deg']:.2f}, "
        f"rms {layered.rms_s:.3f} s",
    ]
    txt = "\n".join(lines)
    print(txt)
    (HERE / "fit.txt").write_text(txt + "\n", encoding="utf-8")
    (HERE / "track.geojson").write_text(json.dumps(to_geojson(fit, frame), indent=1))

    # section with picks and model
    allc = np.array([cross_km(s.lat, s.lon) for s in stations])
    fig, ax = plt.subplots(1, 2, figsize=(14, 7))
    ax[0].plot(pf.t[ok], allc[ok], ".", color="0.7", ms=3, label="picked, rejected")
    ax[0].plot(pf.t[pf.used], allc[pf.used], ".", color="C3", ms=4, label="picked, used")
    o = np.argsort(allc)
    ax[0].plot(pf.pred[o], allc[o], "-", color="C0", lw=0.6, label="model")
    ax[0].set_xlabel("s after 03:42:00 UTC")
    ax[0].set_ylabel("cross-track km (Irish Sea -> Poznan line)")
    ax[0].legend()
    ax[0].set_title(f"boom picks vs model (rms {fit.rms_s:.2f} s)")
    ax[1].hist(fit.residuals_s[fit.inliers], bins=30, color="C3")
    ax[1].set_xlabel("residual (s)")
    fig.tight_layout()
    fig.savefig(HERE / "picks_on_section.png", dpi=110)

    # map
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot([s.lon for s in stations], [s.lat for s in stations], "^", color="0.75", ms=3,
            label="stations")
    sc = ax.scatter([stations[j].lon for j in pf.used], [stations[j].lat for j in pf.used],
                    c=pf.t[pf.used], cmap="viridis", s=14, zorder=3, label="boom picks")
    plt.colorbar(sc, label="arrival (s after 03:42:00 UTC)")
    tp = track_points(fit.trajectory, frame, before_km=500, after_km=900, min_height_m=0)
    ax.plot(tp[:, 1], tp[:, 0], "-", color="C3", lw=2, label="SKYQUAKE track (extrapolated)")
    for draw in fit.sample(60):
        tq = track_points(draw, frame, before_km=500, after_km=900, n=12, min_height_m=0)
        ax.plot(tq[:, 1], tq[:, 0], "-", color="C3", lw=0.3, alpha=0.3)
    for name, (la, lo) in DEBRIS.items():
        ax.plot(lo, la, "x", color="k", ms=8)
        ax.annotate(name, (lo, la), fontsize=8, xytext=(3, 3), textcoords="offset points")
    ax.plot([], [], "x", color="k", label="debris recovered (town centres)")
    ax.set_xlim(2, 19)
    ax.set_ylim(50, 56.5)
    ax.set_aspect(1 / np.cos(np.radians(53)))
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")
    ax.legend(loc="lower left", fontsize=8)
    ax.set_title("Falcon 9 upper-stage re-entry, 2025-02-19, from seismic sonic booms")
    fig.tight_layout()
    fig.savefig(HERE / "map.png", dpi=110)


if __name__ == "__main__":
    main()
