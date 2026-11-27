# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 SKYQUAKE contributors
import numpy as np
import pytest

pytest.importorskip("obspy")

from obspy.signal.trigger import classic_sta_lta  # noqa: E402

from skyquake.pipeline import StationTrace, pick_onset  # noqa: E402


def make_trace(onset=120.0, df=100.0, snr=20.0, seed=0):
    rng = np.random.default_rng(seed)
    t = np.arange(0, 240, 1 / df)
    x = rng.normal(0, 1, t.size)
    k = t >= onset
    x[k] += snr * np.exp(-(t[k] - onset) / 0.8) * np.sin(2 * np.pi * 6 * (t[k] - onset))
    cft = classic_sta_lta(x, int(0.3 * df), int(10 * df))
    cft[: int(10 * df)] = 0
    return StationTrace("XX.TEST", 0, 0, 0, t, cft, x)


def test_pick_finds_onset_within_a_few_samples():
    s = make_trace(onset=120.0)
    t, snr = pick_onset(s, t_pred=118.0, half=10.0)
    assert abs(t - 120.0) < 0.1
    assert snr > 3


def test_no_pick_when_nothing_in_window():
    s = make_trace(onset=200.0)
    assert pick_onset(s, t_pred=100.0, half=5.0) is None
