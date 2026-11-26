"""Test the fragment-cloud hypothesis with the real camera-tracked fragments.

Idea: every camera-tracked fragment (Vierinen et al. 2026) emits a boom along its own path.
A station hears the earliest of all of them. If our single-line fit, applied to arrival times
predicted from this *real fragment cloud*, returns a speed well above the fragments' true
~7 km/s, then the cloud alone explains the high seismic speeds.

    python examples/falcon9_poland/fragment_test.py path/to/ref/fragments
"""

import sys
from pathlib import Path

import h5py
import numpy as np

from skyquake.atmosphere import LayeredAtmosphere, msis_atmosphere
from skyquake.geo import LocalFrame
from skyquake.invert import decluster_sigma, fit_trajectory
from skyquake.io import read_picks

HERE = Path(__file__).parent
DT = 0.2  # s, resampling of each fragment track


def load_fragments(folder, frame):
    """{fragment id: (t, enu)} with all camera pairs merged and resampled to DT."""
    raw = {}
    for f in sorted(Path(folder).glob("*.h5")):
        fid = f.name.split("_")[0]
        with h5py.File(f, "r") as h:
            lat, lon, alt = h["geoloc"][()]
            raw.setdefault(fid, []).append((float(h["time"][()]), frame.to_enu(lat, lon, alt)))
    out = {}
    for fid, rows in raw.items():
        rows.sort(key=lambda r: r[0])
        t = np.array([r[0] for r in rows])
        p = np.array([r[1] for r in rows])
        t, idx = np.unique(t, return_index=True)   # duplicate epochs from camera pairs
        p = p[idx]
        if t.size < 4:
            continue
        tt = np.arange(t[0], t[-1], DT)
        out[fid] = (tt, np.column_stack([np.interp(tt, t, p[:, i]) for i in range(3)]))
    return out


def cloud_arrivals(st, frags, travel, ids=None):
    """Earliest boom arrival over the chosen fragments; also whether the minimum is interior
    to a track (an end-point minimum means the true emission lies outside camera coverage)."""
    best = np.full(len(st), np.inf)
    interior = np.zeros(len(st), bool)
    for fid, (t, p) in frags.items():
        if ids is not None and fid not in ids:
            continue
        horiz = np.hypot(p[None, :, 0] - st[:, None, 0], p[None, :, 1] - st[:, None, 1])
        height = p[None, :, 2] - st[:, None, 2]
        arr = t[None, :] + travel(horiz, height)
        k = np.argmin(arr, axis=1)
        a = arr[np.arange(len(st)), k]
        better = a < best
        best[better] = a[better]
        interior[better] = (k[better] > 0) & (k[better] < t.size - 1)
    return best, interior


def main(folder):
    p = read_picks(HERE / "picks.csv")
    frame = LocalFrame.centred_on(p.lat, p.lon)
    st = frame.to_enu(p.lat, p.lon, p.elev_m)
    obs = p.epoch.timestamp() + p.t_s
    frags = load_fragments(folder, frame)
    print(f"{len(frags)} camera-tracked fragments, {len(p)} stations with boom picks\n")

    atm = LayeredAtmosphere()
    msis = msis_atmosphere("2025-02-19T03:45", 53.0, 7.0)
    travels = {
        "homogeneous 320 m/s": lambda h, z: np.hypot(h, z) / 320.0,
        "layered 1976 x0.94": lambda h, z: atm.travel_time(h, z, 0.94),
        "layered 1976 x1.00": lambda h, z: atm.travel_time(h, z, 1.0),
        "MSIS 2025-02-19 03:45 UTC, 53N 7E": lambda h, z: msis.travel_time(h, z, 1.0),
    }
    sig = decluster_sigma(st, p.sigma_s, 20e3)
    kw = dict(fit_sound_speed=True, azimuth_guess_deg=94.0, outlier_sigma=1e9)
    real = fit_trajectory(st, p.t_s, sig, **kw)
    lines = [f"observed picks, single-line fit: speed {real.params['speed_kms']:.2f} km/s, "
             f"heading {real.params['azimuth_deg']:.2f}", ""]

    # the fragment with the longest camera track serves as the single-source control
    longest = max(frags, key=lambda k: frags[k][0].size)
    for name, travel in travels.items():
        pred, interior = cloud_arrivals(st, frags, travel)
        ok = interior & np.isfinite(pred)
        res = obs[ok] - pred[ok]
        lines.append(f"{name}:")
        lines.append(f"  stations with an interior emission point: {ok.sum()} of {len(st)}")
        lines.append(f"  observed - predicted: median {np.median(res):+.1f} s, "
                     f"robust spread {1.4826 * np.median(np.abs(res - np.median(res))):.2f} s")
        for label, ids in (("whole fragment cloud", None), (f"single fragment '{longest}'",
                                                             {longest})):
            pr, inte = cloud_arrivals(st, frags, travel, ids)
            m = inte & np.isfinite(pr)
            if m.sum() < 20:
                lines.append(f"  {label}: only {m.sum()} usable stations")
                continue
            f = fit_trajectory(st[m], pr[m] - p.epoch.timestamp(), sig[m], **kw)
            lines.append(f"  single-line fit to {label} ({m.sum()} stations): "
                         f"speed {f.params['speed_kms']:.2f} km/s, heading "
                         f"{f.params['azimuth_deg']:.2f}, rms {f.rms_s:.2f} s")
        lines.append("")
    # true fragment speeds near the network
    v = []
    for _t, q in frags.values():
        near = np.hypot(q[:, 0], q[:, 1]) < 200e3
        if near.sum() > 10:
            v.append(np.median(np.linalg.norm(np.diff(q[near], axis=0), axis=1) / DT))
    lines.append(f"true camera speeds within 200 km of the network: median "
                 f"{np.median(v) / 1e3:.2f} km/s over {len(v)} fragments")
    txt = "\n".join(lines)
    print(txt)
    (HERE / "fragment_test.txt").write_text(txt + "\n", encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv[1])
