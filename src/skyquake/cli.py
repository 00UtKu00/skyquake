"""Command line interface.

    skyquake demo                         synthetic end-to-end run (no network)
    skyquake fit picks.csv -o track.geojson
    skyquake pick --lat 34 --lon -118 --radius-km 300 --start ... --end ... -o picks.csv
    skyquake events                       list catalogued validation events
"""

from __future__ import annotations

import argparse
import json
import sys
import tomllib
from datetime import datetime
from importlib import resources

import numpy as np

from skyquake.geo import LocalFrame
from skyquake.invert import PARAM_NAMES, fit_trajectory, traj_to_params
from skyquake.physics import Trajectory, arrival_times, emission_points


def synthetic_scenario(seed: int = 1, n_stations: int = 60, noise_s: float = 0.3,
                       n_outliers: int = 3):
    """A Southern-California-like network and a known trajectory. Returns (frame, st, t, truth)."""
    rng = np.random.default_rng(seed)
    frame = LocalFrame(34.0, -118.0)
    truth = Trajectory(east_m=-60e3, north_m=40e3, up_m=40e3, azimuth_deg=118.0,
                       fpa_deg=6.0, speed_mps=4200.0, t0_s=0.0)
    st_all = []
    while len(st_all) < n_stations:
        lat = 34.0 + rng.uniform(-1.8, 1.8)
        lon = -118.0 + rng.uniform(-2.6, 2.6)
        enu = frame.to_enu(lat, lon, rng.uniform(0, 1500))
        # keep only stations whose boom was emitted while the object was high & hypersonic
        if emission_points(truth, enu[None, :])[0, 2] > 20e3:
            st_all.append(enu)
    st = np.array(st_all)
    t = arrival_times(truth, st) + rng.normal(0, noise_s, len(st))
    bad = rng.choice(len(st), n_outliers, replace=False)
    t[bad] += rng.uniform(5, 20, n_outliers)  # e.g. a passing truck triggered the picker
    return frame, st, t, truth


def cmd_demo(a) -> int:
    frame, st, t, truth = synthetic_scenario(a.seed)
    fit = fit_trajectory(st, t, sigma_s=0.3)
    true = dict(zip(PARAM_NAMES, traj_to_params(truth), strict=True))
    print("SKYQUAKE synthetic demo: 60 stations, 0.3 s pick noise, 3 bad picks\n")
    print(f"{'parameter':<12} {'truth':>10} {'recovered':>10} {'1-sigma':>9} {'error/sigma':>12}")
    for k in PARAM_NAMES:
        err = fit.params[k] - true[k]
        print(f"{k:<12} {true[k]:>10.3f} {fit.params[k]:>10.3f} {fit.sigma[k]:>9.3f} "
              f"{err / fit.sigma[k]:>12.2f}")
    print(f"\nrms residual {fit.rms_s:.3f} s, {fit.n_used} stations used")
    for n in fit.notes:
        print("note:", n)
    return 0


def cmd_fit(a) -> int:
    from skyquake.io import read_picks
    from skyquake.track import to_geojson

    p = read_picks(a.picks)
    frame = LocalFrame.centred_on(p.lat, p.lon)
    st = frame.to_enu(p.lat, p.lon, p.elev_m)
    kw = dict(sound_speed=a.sound_speed, azimuth_guess_deg=a.azimuth,
              fit_sound_speed=not a.fixed_sound_speed, max_speed_kms=a.max_speed)
    if a.decelerating:
        kw["fit_deceleration"] = True
    if a.layered:
        from skyquake.atmosphere import LayeredAtmosphere

        kw["atmosphere"] = LayeredAtmosphere()
    fit = fit_trajectory(st, p.t_s, p.sigma_s, **kw)
    from skyquake.track import corridor, geographic_state

    print(f"epoch (t=0): {p.epoch.isoformat()}")
    print(f"fit in local frame centred at {frame.lat0:.4f}, {frame.lon0:.4f}:")
    print(fit.summary())
    if a.jackknife:
        from skyquake.invert import jackknife

        print("\njackknife 1-sigma (10 groups; more honest when the model is imperfect):")
        if kw.get("atmosphere") is None:
            kw["azimuth_guess_deg"] = fit.params["azimuth_deg"]
        jk = jackknife(st, p.t_s, p.sigma_s, **kw)
        for k, v in jk.items():
            print(f"  {k:<16} +/- {v:.3f}   (formal {fit.sigma[k]:.3f})")
    g = geographic_state(fit.trajectory, frame)
    c = corridor(fit, frame)
    print(f"\nclosest approach to the network (t0 = {fit.params['t0_s']:.2f} s):")
    print(f"  position     lat {g['lat']:.4f}, lon {g['lon']:.4f}, height "
          f"{g['height_m'] / 1e3:.1f} km")
    print(f"  heading      {g['azimuth_deg']:.2f} deg true   descent {g['fpa_deg']:.2f} deg   "
          f"speed {g['speed_mps'] / 1e3:.3f} km/s")
    print(f"  1-sigma      {c['sigma_cross_km']:.2f} km cross-track, "
          f"{c['sigma_height_km']:.2f} km height")
    if a.output:
        gj = to_geojson(fit, frame)
        gj["properties"] = {"epoch_utc": p.epoch.isoformat(), "frame": vars(frame)}
        with open(a.output, "w", encoding="utf-8") as f:
            json.dump(gj, f, indent=2)
        print(f"wrote {a.output}")
    return 0


def cmd_pick(a) -> int:
    from skyquake.io import write_picks
    from skyquake.seismic import pick_sonic_booms

    rows = pick_sonic_booms(a.lat, a.lon, a.radius_km, datetime.fromisoformat(a.start),
                            datetime.fromisoformat(a.end), provider=a.provider)
    write_picks(a.output, rows)
    print(f"{len(rows)} picks -> {a.output}  (review them before fitting!)")
    return 0


def cmd_events(_a) -> int:
    text = resources.files("skyquake").joinpath("data/events.toml").read_text("utf-8")
    data = tomllib.loads(text)
    for ev in data["event"]:
        print(f"{ev['id']:<28} {ev['date']}  {ev['region']:<28} [{ev['status']}]")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="skyquake", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("demo", help="synthetic end-to-end run")
    d.add_argument("--seed", type=int, default=1)
    d.set_defaults(func=cmd_demo)

    f = sub.add_parser("fit", help="fit a trajectory to a picks CSV")
    f.add_argument("picks")
    f.add_argument("-o", "--output", help="GeoJSON output path")
    f.add_argument("--sound-speed", type=float, default=320.0,
                   help="starting (or, with --fixed-sound-speed, fixed) sound speed in m/s")
    f.add_argument("--fixed-sound-speed", action="store_true",
                   help="do not estimate the sound speed from the data (not recommended)")
    f.add_argument("--layered", action="store_true",
                   help="trace rays through a layered standard atmosphere (slower)")
    f.add_argument("--max-speed", type=float, default=12.0,
                   help="upper speed bound in km/s (7.9 = orbital speed limit)")
    f.add_argument("--decelerating", action="store_true",
                   help="fit a speed that changes linearly along the path")
    f.add_argument("--jackknife", action="store_true",
                   help="also report jackknife uncertainties (10 refits)")
    f.add_argument("--azimuth", type=float, help="heading prior (deg), e.g. from the TLE")
    f.set_defaults(func=cmd_fit)

    p = sub.add_parser("pick", help="download seismic data and pick sonic booms")
    p.add_argument("--lat", type=float, required=True)
    p.add_argument("--lon", type=float, required=True)
    p.add_argument("--radius-km", type=float, default=300)
    p.add_argument("--start", required=True, help="UTC, ISO 8601")
    p.add_argument("--end", required=True, help="UTC, ISO 8601")
    p.add_argument("--provider", default="IRIS")
    p.add_argument("-o", "--output", default="picks.csv")
    p.set_defaults(func=cmd_pick)

    e = sub.add_parser("events", help="list validation events")
    e.set_defaults(func=cmd_events)

    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
