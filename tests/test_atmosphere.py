import numpy as np
import pytest

from skyquake.atmosphere import (
    LayeredAtmosphere,
    layered_arrival_times,
    sound_speed_1976,
    temperature_1976,
)
from skyquake.physics import Trajectory, arrival_times


def test_standard_atmosphere_reference_values():
    # textbook values of the 1976 standard atmosphere
    t = temperature_1976(np.array([0.0, 11e3, 20e3, 47e3, 51e3, 71e3]))
    np.testing.assert_allclose(t, [288.15, 216.65, 216.65, 270.65, 270.65, 214.65], atol=0.01)
    assert sound_speed_1976(0.0) == pytest.approx(340.29, abs=0.05)


@pytest.fixture(scope="module")
def uniform():
    z = np.arange(0.0, 100e3 + 1, 1000.0)
    return LayeredAtmosphere(z, np.full(z.size, 320.0), xmax_m=600e3, dx_m=250.0)


def test_uniform_profile_gives_straight_rays(uniform):
    d = np.array([0.0, 10e3, 50e3, 120e3, 300e3])
    h = np.array([60e3, 60e3, 40e3, 57e3, 80e3])
    np.testing.assert_allclose(uniform.travel_time(d, h), np.hypot(d, h) / 320.0, rtol=2e-4)


def test_layered_mach_matches_closed_form_for_uniform_air(uniform):
    tr = Trajectory(east_m=0, north_m=0, up_m=57e3, azimuth_deg=65.0, fpa_deg=0.5,
                    speed_mps=8000.0, t0_s=0.0)
    rng = np.random.default_rng(0)
    st = np.column_stack([rng.uniform(-2e5, 2e5, 40), rng.uniform(-2e5, 2e5, 40),
                          rng.uniform(0, 2e3, 40)])
    np.testing.assert_allclose(layered_arrival_times(tr, st, uniform),
                               arrival_times(tr, st, 320.0), atol=0.05)


def test_scale_speeds_up_sound(uniform):
    t1 = uniform.travel_time(100e3, 50e3)
    t2 = uniform.travel_time(100e3, 50e3, scale=1.1)
    assert t2 == pytest.approx(t1 / 1.1)


def test_standard_atmosphere_has_a_shadow_zone():
    atm = LayeredAtmosphere()
    # the ground is the fastest layer below 57 km, so direct rays only reach a finite range
    x_max_57 = np.interp(57e3, atm.z, atm.x_max)
    assert 50e3 < x_max_57 < 800e3
    assert not atm.in_shadow(20e3, 57e3)
    assert atm.in_shadow(x_max_57 + 10e3, 57e3)
