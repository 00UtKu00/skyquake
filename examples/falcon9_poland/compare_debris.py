"""Distance of each recovered-debris town from the SKYQUAKE ground track (extrapolated
~650 km east of the seismic network), with Monte-Carlo uncertainty.

The ground track of a straight line is a great circle, so we measure the perpendicular
(cross-track) distance on the sphere from each town to the track.
"""

from pathlib import Path

import numpy as np

from skyquake.geo import LocalFrame
from skyquake.invert import fit_trajectory, jackknife
from skyquake.io import read_picks
from skyquake.track import track_points

HERE = Path(__file__).parent
R = 6371.0
DEBRIS = {
    "Komorniki": (52.335, 16.810), "Wiry": (52.325, 16.860), "Sliwno": (52.590, 16.320),
    "Sedziny (Szamotuly)": (52.680, 16.450), "Miedzychod": (52.600, 15.890),
}


def unit(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)


def cross_track_km(track_latlon, lat, lon):
    """Signed distance (km, + = left/north of travel) from a point to the track's great circle."""
    a, b = unit(*track_latlon[0]), unit(*track_latlon[-1])
    n = np.cross(a, b)
    n /= np.linalg.norm(n)
    return R * np.arcsin(np.clip(unit(lat, lon) @ n, -1, 1))


def main():
    p = read_picks(HERE / "picks.csv")
    frame = LocalFrame.centred_on(p.lat, p.lon)
    st = frame.to_enu(p.lat, p.lon, p.elev_m)
    kw = dict(fit_sound_speed=True, azimuth_guess_deg=95.0)
    fit = fit_trajectory(st, p.t_s, p.sigma_s, **kw)
    jk = jackknife(st, p.t_s, p.sigma_s, **kw)
    inflate = max(jk["azimuth_deg"] / fit.sigma["azimuth_deg"],
                  jk["cross_km"] / fit.sigma["cross_km"], 1.0)
    nominal = track_points(fit.trajectory, frame, before_km=0, after_km=800, n=2)[:, :2]
    rng = np.random.default_rng(1)
    mean = np.array([fit.params[k] for k in fit.names])
    draws = rng.multivariate_normal(mean, fit.covariance * inflate**2, size=2000)
    from skyquake.invert import _to_traj

    tracks = [track_points(_to_traj(d), frame, before_km=0, after_km=800, n=2)[:, :2]
              for d in draws]
    print(f"track from the network centre ({nominal[0][0]:.3f}, {nominal[0][1]:.3f}) "
          f"to 800 km downrange ({nominal[1][0]:.3f}, {nominal[1][1]:.3f})")
    print(f"uncertainty: covariance inflated x{inflate:.1f} to match the jackknife\n")
    print(f"{'debris town':<22} {'miss (km)':>10} {'1-sigma':>8}   (+ = north of track)")
    out = []
    for name, (la, lo) in DEBRIS.items():
        d0 = cross_track_km(nominal, la, lo)
        ds = np.array([cross_track_km(t, la, lo) for t in tracks])
        out.append((name, d0, ds.std()))
        print(f"{name:<22} {d0:>+10.1f} {ds.std():>8.1f}")
    miss = np.array([o[1] for o in out])
    print(f"\nall {len(out)} sites within {np.abs(miss).max():.1f} km of the track; "
          f"mean offset {miss.mean():+.1f} km; strewn field spans "
          f"{miss.max() - miss.min():.1f} km across track")
    (HERE / "debris_comparison.txt").write_text(
        "\n".join(f"{n}\t{d:+.1f}\t{s:.1f}" for n, d, s in out) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
