# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 SKYQUAKE contributors
"""Inversion tests. The calibration test is the important one: it checks that the reported
1-sigma uncertainties are honest (neither over- nor under-confident)."""

import numpy as np

from skyquake.cli import synthetic_scenario
from skyquake.geo import LocalFrame
from skyquake.invert import PARAM_NAMES, fit_trajectory, traj_to_params
from skyquake.physics import Trajectory, arrival_times
from skyquake.track import corridor


def test_parameterisation_roundtrip():
    from skyquake.invert import _to_traj

    _, _, _, truth = synthetic_scenario(seed=1)
    p = traj_to_params(truth)
    tr = _to_traj(p)
    st = np.array([[1e4, 2e4, 0.0], [-5e4, 3e4, 500.0]])
    np.testing.assert_allclose(arrival_times(tr, st), arrival_times(truth, st), atol=1e-9)


def test_noise_free_recovery_is_exact():
    _, st, t, truth = synthetic_scenario(seed=3, noise_s=0.0, n_outliers=0)
    t = arrival_times(truth, st)
    fit = fit_trajectory(st, t, sigma_s=0.1)
    got = np.array([fit.params[k] for k in PARAM_NAMES])
    np.testing.assert_allclose(got, traj_to_params(truth), atol=1e-3)


def test_blind_heading_search_finds_the_right_direction():
    _, st, t, truth = synthetic_scenario(seed=5)
    fit = fit_trajectory(st, t, sigma_s=0.3)  # no azimuth prior at all
    assert abs((fit.params["azimuth_deg"] - truth.azimuth_deg + 180) % 360 - 180) < 2.0


def test_outliers_are_rejected():
    _, st, t, _ = synthetic_scenario(seed=2, n_outliers=3)
    fit = fit_trajectory(st, t, sigma_s=0.3)
    assert (~fit.inliers).sum() >= 3
    assert fit.rms_s < 0.5


def test_uncertainties_are_calibrated():
    """Over many noisy realisations, |error| < 2 sigma should hold ~95% of the time."""
    z = []
    for seed in range(40):
        _, st, t, truth = synthetic_scenario(seed=100 + seed, n_outliers=0)
        fit = fit_trajectory(st, t, sigma_s=0.3, azimuth_guess_deg=118.0)
        got = np.array([fit.params[k] for k in PARAM_NAMES])
        sig = np.array([fit.sigma[k] for k in PARAM_NAMES])
        z.append((got - traj_to_params(truth)) / sig)
    z = np.abs(np.array(z))
    coverage = (z < 2.0).mean(axis=0)
    # 40 trials: allow sampling scatter around the nominal 95.4%
    assert np.all(coverage >= 0.85), dict(zip(PARAM_NAMES, coverage, strict=True))


def test_level_hypersonic_flight():
    """Regression: real re-entries can be almost level (Shenzhou-15: ~0.3 deg). The old
    'point at 40 km height' parameterisation broke down here; this one must not."""
    rng = np.random.default_rng(11)
    frame = LocalFrame(35.0, -117.5)
    truth = Trajectory(east_m=0, north_m=0, up_m=65e3, azimuth_deg=64.0, fpa_deg=0.3,
                       speed_mps=9000.0, t0_s=0.0)
    lat = 35.0 + rng.uniform(-1.5, 1.5, 80)
    lon = -117.5 + rng.uniform(-2.0, 2.0, 80)
    st = frame.to_enu(lat, lon, rng.uniform(0, 1500, 80))
    t = arrival_times(truth, st) + rng.normal(0, 0.3, 80)
    fit = fit_trajectory(st, t, sigma_s=0.3, fit_sound_speed=True)
    want = traj_to_params(truth)
    for i, k in enumerate(PARAM_NAMES):
        assert abs(fit.params[k] - want[i]) < 4 * fit.sigma[k] + 1e-6, k
    assert fit.sigma["height_km"] < 5.0


def test_sound_speed_is_recovered_from_a_wrong_start():
    _, st, t, truth = synthetic_scenario(seed=7, n_outliers=0)
    fit = fit_trajectory(st, t, sigma_s=0.3, sound_speed=300.0, fit_sound_speed=True)
    c, sc = fit.params["sound_speed_kms"], fit.sigma["sound_speed_kms"]
    assert abs(c - 0.320) < 3 * sc
    want = traj_to_params(truth)
    assert abs(fit.params["height_km"] - want[1]) < 3 * fit.sigma["height_km"] + 0.5


def test_wrong_fixed_sound_speed_is_flagged_by_residuals():
    _, st, t, _ = synthetic_scenario(seed=1, n_outliers=0)
    good = fit_trajectory(st, t, sigma_s=0.3, sound_speed=320.0, azimuth_guess_deg=118)
    bad = fit_trajectory(st, t, sigma_s=0.3, sound_speed=304.0, azimuth_guess_deg=118)
    assert bad.reduced_chi2 > 2 * good.reduced_chi2


def test_corridor_is_finite():
    frame, st, t, _ = synthetic_scenario(seed=4)
    fit = fit_trajectory(st, t, sigma_s=0.3)
    c = corridor(fit, frame, n=300)
    assert c["samples"] > 250
    assert 0 < c["sigma_cross_km"] < 20 and 0 < c["sigma_height_km"] < 20
