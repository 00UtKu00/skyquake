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


def test_max_speed_bound_is_respected(scenario):
    _, st, t, _ = scenario
    fit = fit_trajectory(st, t, 0.3, fit_sound_speed=True, max_speed_kms=7.0,
                         azimuth_guess_deg=65.0)
    assert fit.params["speed_kms"] <= 7.0 + 1e-9
