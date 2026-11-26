"""Inversion tests. The calibration test is the important one: it checks that the reported
1-sigma uncertainties are honest (neither over- nor under-confident)."""

import numpy as np

from skyquake.cli import synthetic_scenario
from skyquake.invert import PARAM_NAMES, fit_trajectory
from skyquake.track import corridor


def truth_vector(tr):
    return np.array([tr.east_m / 1e3, tr.north_m / 1e3, tr.azimuth_deg, tr.fpa_deg,
                     tr.speed_mps / 1e3, tr.t0_s])


def test_noise_free_recovery_is_exact():
    _, st, t, truth = synthetic_scenario(seed=3, noise_s=0.0, n_outliers=0)
    from skyquake.physics import arrival_times
    t = arrival_times(truth, st)
    fit = fit_trajectory(st, t, sigma_s=0.1)
    got = np.array([fit.params[k] for k in PARAM_NAMES])
    np.testing.assert_allclose(got, truth_vector(truth), atol=1e-3)


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
        z.append((got - truth_vector(truth)) / sig)
    z = np.abs(np.array(z))
    coverage = (z < 2.0).mean(axis=0)
    # 40 trials: allow sampling scatter around the nominal 95.4%
    assert np.all(coverage >= 0.85), dict(zip(PARAM_NAMES, coverage, strict=True))
    assert np.all(coverage <= 1.0)


def test_sound_speed_is_recovered_from_a_wrong_start():
    _, st, t, truth = synthetic_scenario(seed=7, n_outliers=0)
    fit = fit_trajectory(st, t, sigma_s=0.3, sound_speed=300.0, fit_sound_speed=True)
    c, sc = fit.params["sound_speed_kms"], fit.sigma["sound_speed_kms"]
    assert abs(c - 0.320) < 3 * sc
    # along-track position no longer biased by ~100 km
    assert abs(fit.params["east_km"] - truth.east_m / 1e3) < 3 * fit.sigma["east_km"] + 1.0


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
    assert 0 < c["sigma_cross_km"] < 20 and 0 < c["sigma_along_km"] < 50
