# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 SKYQUAKE contributors
import numpy as np
import pytest
from scipy.optimize import minimize_scalar

from skyquake.physics import (
    Trajectory,
    arrival_times,
    decelerating_arrival_times,
    emission_time,
)

TR = Trajectory(east_m=0, north_m=0, up_m=65e3, azimuth_deg=94.0, fpa_deg=1.0,
                speed_mps=7000.0, t0_s=0.0)
C = 320.0
rng = np.random.default_rng(0)
ST = np.column_stack([rng.uniform(-3e5, 3e5, 30), rng.uniform(-2e5, 1e5, 30),
                      rng.uniform(0, 500, 30)])


def test_zero_deceleration_matches_closed_form():
    np.testing.assert_allclose(decelerating_arrival_times(TR, ST, C, 0.0),
                               arrival_times(TR, ST, C), atol=0.02)


def test_emission_time_limits():
    assert emission_time(7000.0, 7000.0, 0.0) == pytest.approx(1.0)
    # v halves after 350 km when g = -0.01 /s  ->  tau = ln(0.5)/(-0.01)
    assert emission_time(350e3, 7000.0, -0.01) == pytest.approx(np.log(0.5) / -0.01)


def brute(tr, q, c, dv):
    g = dv * 1e3 / 100e3

    def f(s):
        if tr.speed_mps + g * s <= 1.02 * c:
            return 1e9
        p = tr.origin + s * tr.direction
        return tr.t0_s + emission_time(s, tr.speed_mps, g) + np.linalg.norm(q - p) / c

    s = np.linspace(-6e5, 6e5, 24001)
    k = int(np.argmin([f(x) for x in s]))
    return minimize_scalar(f, bounds=(s[max(k - 1, 0)], s[min(k + 1, len(s) - 1)]),
                           method="bounded", options={"xatol": 0.1}).fun


@pytest.mark.parametrize("dv", [-0.5, -1.2])
def test_decelerating_matches_brute_force(dv):
    got = decelerating_arrival_times(TR, ST[:12], C, dv)
    want = np.array([brute(TR, q, C, dv) for q in ST[:12]])
    np.testing.assert_allclose(got, want, atol=0.02)


def test_fit_recovers_speed_and_deceleration():
    from skyquake.invert import fit_trajectory, traj_to_params

    r = np.random.default_rng(4)
    st = np.column_stack([r.uniform(-1e5, 5e5, 120), r.uniform(-2.5e5, 0.3e5, 120),
                          r.uniform(0, 500, 120)])
    truth = Trajectory(east_m=0, north_m=0, up_m=65e3, azimuth_deg=94.0, fpa_deg=1.0,
                       speed_mps=7000.0, t0_s=0.0)
    t = decelerating_arrival_times(truth, st, C, -0.9) + r.normal(0, 0.3, len(st))
    const = fit_trajectory(st, t, 0.3, fit_sound_speed=True, azimuth_guess_deg=94.0)
    dec = fit_trajectory(st, t, 0.3, fit_deceleration=True, azimuth_guess_deg=94.0)
    want = traj_to_params(truth)
    assert abs(dec.params["speed_kms"] - 7.0) < 4 * dec.sigma["speed_kms"] + 0.05
    assert abs(dec.params["dv_kms_per_100km"] + 0.9) < 4 * dec.sigma["dv_kms_per_100km"] + 0.05
    assert abs(dec.params["azimuth_deg"] - want[2]) < 4 * dec.sigma["azimuth_deg"] + 0.05
    # the constant-speed model on the same data is worse and biased
    assert dec.reduced_chi2 < const.reduced_chi2
    assert abs(const.params["speed_kms"] - 7.0) > abs(dec.params["speed_kms"] - 7.0)


def test_deceleration_delays_downrange_stations():
    east = np.array([[250e3, -40e3, 0.0]])
    assert decelerating_arrival_times(TR, east, C, -1.0)[0] > arrival_times(TR, east, C)[0]
