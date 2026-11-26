"""Reported geographic values must not depend on which local frame the fit ran in."""

from dataclasses import replace

import numpy as np

from skyquake.cli import synthetic_scenario
from skyquake.geo import LocalFrame
from skyquake.invert import fit_trajectory
from skyquake.track import geographic_state


def test_geographic_state_is_frame_independent():
    frame, st, t, truth = synthetic_scenario(seed=1)
    lat, lon, h = frame.from_enu(st)
    true_geo = geographic_state(truth, frame)

    # refit in a frame shifted ~1 degree away, as `skyquake fit` does (centred on stations)
    other = LocalFrame(frame.lat0 + 0.6, frame.lon0 + 0.9)
    fit = fit_trajectory(other.to_enu(lat, lon, h), t, sigma_s=0.3)
    got = geographic_state(fit.trajectory, other)

    # frame-relative azimuth differs by meridian convergence (~0.5 deg here) ...
    assert abs(fit.params["azimuth_deg"] - truth.azimuth_deg) > 0.3

    # The two frames put their "40 km" reference point at different places on the line, and
    # a straight line's geographic heading changes along it. Compare at the same point:
    p_true = frame.to_enu(true_geo["lat"], true_geo["lon"], true_geo["height_m"])
    p_got = frame.to_enu(got["lat"], got["lon"], got["height_m"])
    off = p_got - p_true
    s = off @ truth.direction
    assert np.linalg.norm(off - s * truth.direction) < 3e3  # same line, cross-track
    o = truth.origin + s * truth.direction
    same_point = replace(truth, east_m=o[0], north_m=o[1], up_m=o[2])
    ref = geographic_state(same_point, frame)

    # ... but the geographic heading and descent angle agree within their uncertainty
    assert abs(got["azimuth_deg"] - ref["azimuth_deg"]) < 4 * fit.sigma["azimuth_deg"]
    assert abs(got["fpa_deg"] - ref["fpa_deg"]) < 4 * fit.sigma["fpa_deg"]
