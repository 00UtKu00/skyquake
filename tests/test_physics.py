import numpy as np
import pytest
from scipy.optimize import minimize_scalar

from falltrace.physics import Trajectory, arrival_times, emission_points, mach_angle_deg

TRAJ = Trajectory(east_m=-50e3, north_m=20e3, up_m=40e3, azimuth_deg=110.0,
                  fpa_deg=7.0, speed_mps=4000.0, t0_s=12.0)
C = 320.0


def brute_force_arrival(traj, station, c):
    """Minimise emission-time + travel-time over the source position (no closed form used)."""
    u, o, v = traj.direction, traj.origin, traj.speed_mps

    def t_of(s):
        return traj.t0_s + s / v + np.linalg.norm(station - (o + s * u)) / c

    return minimize_scalar(t_of, bounds=(-2e6, 2e6), method="bounded",
                           options={"xatol": 1e-6}).fun


def test_closed_form_matches_brute_force():
    rng = np.random.default_rng(0)
    st = np.column_stack([rng.uniform(-2e5, 2e5, 25), rng.uniform(-2e5, 2e5, 25),
                          rng.uniform(-5e3, 2e3, 25)])
    closed = arrival_times(TRAJ, st, C)
    brute = np.array([brute_force_arrival(TRAJ, s, C) for s in st])
    np.testing.assert_allclose(closed, brute, atol=1e-4)


def test_station_on_track_hears_it_as_it_passes():
    # directly under the line (r = height): t = t0 + x/v + h*sqrt(1/c^2 - 1/v^2)
    tr = Trajectory(0, 0, 30e3, 0.0, 0.0, 3000.0, 0.0)
    t = arrival_times(tr, np.array([[0.0, 10e3, 0.0]]), C)[0]
    assert t == pytest.approx(10e3 / 3000 + 30e3 * np.sqrt(1 / C**2 - 1 / 3000**2))


def test_emission_point_is_consistent():
    st = np.array([[30e3, -40e3, 0.0]])
    e = emission_points(TRAJ, st, C)[0]
    s = (e - TRAJ.origin) @ TRAJ.direction
    t = TRAJ.t0_s + s / TRAJ.speed_mps + np.linalg.norm(st[0] - e) / C
    assert t == pytest.approx(arrival_times(TRAJ, st, C)[0], abs=1e-9)


def test_subsonic_rejected():
    with pytest.raises(ValueError):
        arrival_times(Trajectory(0, 0, 1e4, 0, 0, 300.0, 0), np.zeros((1, 3)), C)


def test_mach_angle():
    assert mach_angle_deg(640.0, 320.0) == pytest.approx(30.0)
