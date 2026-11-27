# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 SKYQUAKE contributors
"""Spatial-block jackknife and declustering: honest errors for clustered networks."""

import numpy as np

from skyquake.geo import LocalFrame
from skyquake.invert import decluster_sigma, fit_trajectory, jackknife, spatial_blocks
from skyquake.physics import Trajectory, arrival_times


def test_spatial_blocks_partition_all_stations():
    rng = np.random.default_rng(0)
    st = np.column_stack([rng.uniform(-2e5, 2e5, 300), rng.uniform(-2e5, 2e5, 300),
                          np.zeros(300)])
    parts = spatial_blocks(st, 8)
    assert len(parts) == 8
    allidx = np.sort(np.concatenate(parts))
    np.testing.assert_array_equal(allidx, np.arange(300))


def test_decluster_inflates_only_dense_stations():
    st = np.array([[0, 0, 0], [1e3, 0, 0], [2e3, 0, 0], [3e3, 0, 0], [3e5, 0, 0]], float)
    s = decluster_sigma(st, 0.3, radius_m=20e3)
    assert s[-1] == 0.3                       # isolated station unchanged
    np.testing.assert_allclose(s[:4], 0.3 * 2)  # 4 neighbours (incl. itself) -> x2


def clustered_scenario(seed=0):
    """A dense cluster under the track plus a sparse, one-sided spread: Groningen-like."""
    rng = np.random.default_rng(seed)
    frame = LocalFrame(53.0, 7.0)
    truth = Trajectory(east_m=0, north_m=0, up_m=65e3, azimuth_deg=95.0, fpa_deg=1.0,
                       speed_mps=6500.0, t0_s=0.0)
    dense = np.column_stack([7.0 + rng.normal(0, 0.15, 200), 53.2 + rng.normal(0, 0.1, 200)])
    sparse = np.column_stack([rng.uniform(4, 10, 40), rng.uniform(50.6, 52.8, 40)])
    ll = np.vstack([dense, sparse])
    st = frame.to_enu(ll[:, 1], ll[:, 0], 0.0)
    # correlated error: one shared offset per 30 km cell (local atmosphere), plus pick noise
    cell = np.floor(st[:, :2] / 30e3).astype(int)
    keys = {tuple(c): rng.normal(0, 0.6) for c in map(tuple, cell)}
    t = arrival_times(truth, st) + np.array([keys[tuple(c)] for c in cell])
    t += rng.normal(0, 0.2, len(t))
    return st, t, truth


def test_spatial_jackknife_is_more_honest_than_random_on_clustered_networks():
    st, t, truth = clustered_scenario()
    kw = dict(fit_sound_speed=True, azimuth_guess_deg=95.0)
    rand = jackknife(st, t, 0.2, blocks="random", **kw)
    spat = jackknife(st, t, 0.2, blocks="spatial", **kw)
    assert spat["azimuth_deg"] > rand["azimuth_deg"]


def test_declustered_fit_does_not_shrink_errors():
    st, t, _ = clustered_scenario(1)
    plain = fit_trajectory(st, t, 0.2, fit_sound_speed=True, azimuth_guess_deg=95.0)
    dec = fit_trajectory(st, t, decluster_sigma(st, 0.2), fit_sound_speed=True,
                         azimuth_guess_deg=95.0)
    assert dec.sigma["azimuth_deg"] >= plain.sigma["azimuth_deg"]
