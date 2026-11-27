# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 SKYQUAKE contributors
"""Layered-atmosphere inversion: recovers truth on layered synthetic data, and shows why it
matters (a homogeneous fit to the same data is biased)."""

import numpy as np
import pytest

from skyquake.atmosphere import LayeredAtmosphere, layered_arrival_times
from skyquake.geo import LocalFrame
from skyquake.invert import PARAM_NAMES, fit_trajectory, traj_to_params
from skyquake.physics import Trajectory


@pytest.fixture(scope="module")
def scenario():
    atm = LayeredAtmosphere()
    rng = np.random.default_rng(21)
    frame = LocalFrame(35.0, -117.5)
    truth = Trajectory(east_m=0, north_m=0, up_m=57e3, azimuth_deg=65.0, fpa_deg=0.4,
                       speed_mps=7600.0, t0_s=0.0)
    lat = 35.0 + rng.uniform(-1.0, 1.0, 70)
    lon = -117.5 + rng.uniform(-1.3, 1.3, 70)
    st = frame.to_enu(lat, lon, rng.uniform(0, 1500, 70))
    t = layered_arrival_times(truth, st, atm) + rng.normal(0, 0.3, 70)
    return atm, st, t, truth


def test_layered_fit_recovers_truth(scenario):
    atm, st, t, truth = scenario
    fit = fit_trajectory(st, t, 0.3, atmosphere=atm, fit_sound_speed=True,
                         azimuth_guess_deg=65.0)
    want = traj_to_params(truth)
    for i, k in enumerate(PARAM_NAMES):
        assert abs(fit.params[k] - want[i]) < 4 * fit.sigma[k] + 1e-6, (k, fit.params[k], want[i])
    assert abs(fit.params["atm_scale"] - 1.0) < 4 * fit.sigma["atm_scale"] + 1e-4


def test_layered_fit_handles_a_slower_than_standard_atmosphere():
    """Regression for Falcon 9: true air 6 % slower than standard. The homogeneous start is
    tens of seconds off in absolute time; the layered fit must still find the truth."""
    atm = LayeredAtmosphere()
    rng = np.random.default_rng(8)
    frame = LocalFrame(53.0, 7.0)
    truth = Trajectory(east_m=0, north_m=0, up_m=68e3, azimuth_deg=94.0, fpa_deg=1.0,
                       speed_mps=7000.0, t0_s=0.0)
    lat = np.concatenate([53.2 + rng.normal(0, 0.1, 60), rng.uniform(51.3, 52.9, 40)])
    lon = np.concatenate([6.8 + rng.normal(0, 0.15, 60), rng.uniform(4.5, 10.0, 40)])
    st = frame.to_enu(lat, lon, 0.0)
    t = layered_arrival_times(truth, st, atm, scale=0.94) + rng.normal(0, 0.3, len(st))
    fit = fit_trajectory(st, t, 0.3, atmosphere=atm, fit_sound_speed=True,
                         azimuth_guess_deg=94.0)
    want = traj_to_params(truth)
    assert fit.rms_s < 0.6
    assert abs(fit.params["t0_s"] - want[5]) < 4 * fit.sigma["t0_s"] + 1.0
    assert abs(fit.params["height_km"] - want[1]) < 4 * fit.sigma["height_km"] + 0.5
    assert abs(fit.params["atm_scale"] - 0.94) < 4 * fit.sigma["atm_scale"] + 0.005


def test_max_speed_bound_is_respected(scenario):
    _, st, t, _ = scenario
    fit = fit_trajectory(st, t, 0.3, fit_sound_speed=True, max_speed_kms=7.0,
                         azimuth_guess_deg=65.0)
    assert fit.params["speed_kms"] <= 7.0 + 1e-9
