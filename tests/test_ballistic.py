# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 SKYQUAKE contributors
"""Physical trajectory model: sanity of the dynamics and recovery of speed and drag."""

import numpy as np
import pytest

from skyquake.atmosphere import LayeredAtmosphere
from skyquake.ballistic import BallisticModel, density_1976, fit_ballistic
from skyquake.geo import LocalFrame


@pytest.fixture(scope="module")
def model():
    return BallisticModel(LocalFrame(53.0, 7.0), LayeredAtmosphere())


def test_density_reference_values():
    assert density_1976(0.0) == pytest.approx(1.225)
    assert density_1976(50e3) == pytest.approx(1.027e-3, rel=1e-3)


def test_no_drag_keeps_orbital_energy(model):
    # beta = 1e5 kg/m^2 at 120 km: essentially drag-free over the window
    p = np.array([0.0, 120.0, 90.0, 0.0, 7.8, 0.0, 5.0])
    t, enu, v = model.propagate(p)
    assert np.ptp(v) < 120.0          # only gravity / curvature exchange, < 1.5 %


def test_drag_decelerates_and_path_curves_down(model):
    p = np.array([0.0, 70.0, 94.0, 1.0, 7.0, 0.0, 2.5])
    t, enu, v = model.propagate(p)
    k0 = np.argmin(np.abs(t))
    assert v[-1] < v[k0] < v[0]       # slower later, faster earlier
    later = t > 60
    assert np.all(np.diff(enu[later, 2]) < 0)   # keeps descending


def test_fit_recovers_speed_and_beta(model):
    rng = np.random.default_rng(3)
    lat = np.concatenate([53.2 + rng.normal(0, 0.12, 50), rng.uniform(51.5, 52.9, 40)])
    lon = np.concatenate([6.9 + rng.normal(0, 0.2, 50), rng.uniform(4.5, 11.0, 40)])
    st = model.frame.to_enu(lat, lon, 0.0)
    truth = np.array([-10.0, 72.0, 94.0, 1.0, 7.0, 0.0, 2.6])
    t = model.arrivals(truth, st) + rng.normal(0, 0.3, len(st))
    fit = fit_ballistic(st, t, 0.3, model)
    assert abs(fit.params["speed_kms"] - 7.0) < 4 * fit.sigma["speed_kms"] + 0.1
    assert abs(fit.params["height_km"] - 72.0) < 4 * fit.sigma["height_km"] + 0.5
    assert abs(fit.params["t0_s"]) < 4 * fit.sigma["t0_s"] + 1.0
    assert fit.rms_s < 0.6
