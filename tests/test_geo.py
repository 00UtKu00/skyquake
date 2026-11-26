import numpy as np

from falltrace.geo import LocalFrame, ecef_to_geodetic, geodetic_to_ecef


def test_ecef_roundtrip():
    lat = np.array([0.0, 34.05, -33.9, 52.4, 70.0])
    lon = np.array([0.0, -118.25, 151.2, 16.9, -150.0])
    h = np.array([0.0, 120.0, 30.0, 80_000.0, -50.0])
    lat2, lon2, h2 = ecef_to_geodetic(geodetic_to_ecef(lat, lon, h))
    np.testing.assert_allclose(lat2, lat, atol=1e-9)
    np.testing.assert_allclose(lon2, lon, atol=1e-9)
    np.testing.assert_allclose(h2, h, atol=1e-4)


def test_equator_radius():
    np.testing.assert_allclose(geodetic_to_ecef(0, 0, 0), [6_378_137.0, 0, 0])


def test_enu_origin_and_axes():
    f = LocalFrame(34.0, -118.0)
    np.testing.assert_allclose(f.to_enu(34.0, -118.0, 0.0), 0, atol=1e-6)
    up = f.to_enu(34.0, -118.0, 1000.0)
    np.testing.assert_allclose(up, [0, 0, 1000.0], atol=1e-6)
    north = f.to_enu(34.01, -118.0, 0.0)
    assert north[1] > 1100 and abs(north[0]) < 1e-6


def test_curvature_drop():
    # a point on the ellipsoid 300 km east sits ~7 km below the tangent plane (d^2 / 2R)
    f = LocalFrame(0.0, 0.0)
    p = f.to_enu(0.0, np.degrees(300e3 / 6_378_137.0), 0.0)
    assert 6_900 < -p[2] < 7_200


def test_enu_roundtrip():
    f = LocalFrame(52.4, 16.9, 80.0)
    enu = np.array([[1e5, -2e5, 3e4], [0, 0, 0], [-4e4, 1e3, 1e5]])
    lat, lon, h = f.from_enu(enu)
    np.testing.assert_allclose(f.to_enu(lat, lon, h), enu, atol=1e-5)
