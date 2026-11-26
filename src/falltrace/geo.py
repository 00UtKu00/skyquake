"""WGS84 geodesy: geodetic <-> ECEF <-> local East-North-Up (ENU) frame.

All inversions run in a Cartesian ENU frame centred on a reference point. Because the
frame is built from exact ECEF coordinates, Earth curvature is handled correctly for the
stations: a station 300 km away sits about 7 km *below* the frame's horizontal plane.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

WGS84_A = 6_378_137.0
WGS84_F = 1.0 / 298.257223563
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)


def geodetic_to_ecef(lat_deg, lon_deg, h_m) -> np.ndarray:
    """Geodetic latitude/longitude (deg) and ellipsoidal height (m) to ECEF (m), shape (..., 3)."""
    lat = np.radians(np.asarray(lat_deg, dtype=float))
    lon = np.radians(np.asarray(lon_deg, dtype=float))
    h = np.asarray(h_m, dtype=float)
    n = WGS84_A / np.sqrt(1.0 - WGS84_E2 * np.sin(lat) ** 2)
    x = (n + h) * np.cos(lat) * np.cos(lon)
    y = (n + h) * np.cos(lat) * np.sin(lon)
    z = (n * (1.0 - WGS84_E2) + h) * np.sin(lat)
    return np.stack(np.broadcast_arrays(x, y, z), axis=-1)


def ecef_to_geodetic(xyz) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """ECEF (m), shape (..., 3), to geodetic latitude/longitude (deg) and height (m).

    Fixed-point iteration; converges to sub-millimetre accuracy away from the poles.
    """
    xyz = np.asarray(xyz, dtype=float)
    x, y, z = xyz[..., 0], xyz[..., 1], xyz[..., 2]
    p = np.hypot(x, y)
    lon = np.arctan2(y, x)
    lat = np.arctan2(z, p * (1.0 - WGS84_E2))
    h = np.zeros_like(p)
    for _ in range(12):
        n = WGS84_A / np.sqrt(1.0 - WGS84_E2 * np.sin(lat) ** 2)
        h = p / np.cos(lat) - n
        lat = np.arctan2(z, p * (1.0 - WGS84_E2 * n / (n + h)))
    return np.degrees(lat), np.degrees(lon), h


@dataclass(frozen=True)
class LocalFrame:
    """Local tangent-plane ENU frame anchored at (lat0, lon0, h0)."""

    lat0: float
    lon0: float
    h0: float = 0.0

    @property
    def _rotation(self) -> np.ndarray:
        lat, lon = np.radians(self.lat0), np.radians(self.lon0)
        sl, cl, so, co = np.sin(lat), np.cos(lat), np.sin(lon), np.cos(lon)
        return np.array(
            [
                [-so, co, 0.0],
                [-sl * co, -sl * so, cl],
                [cl * co, cl * so, sl],
            ]
        )

    @property
    def _origin(self) -> np.ndarray:
        return geodetic_to_ecef(self.lat0, self.lon0, self.h0)

    def to_enu(self, lat_deg, lon_deg, h_m=0.0) -> np.ndarray:
        """Geodetic -> ENU metres, shape (..., 3)."""
        d = geodetic_to_ecef(lat_deg, lon_deg, h_m) - self._origin
        return d @ self._rotation.T

    def from_enu(self, enu) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """ENU metres, shape (..., 3) -> geodetic (lat deg, lon deg, h m)."""
        xyz = np.asarray(enu, dtype=float) @ self._rotation + self._origin
        return ecef_to_geodetic(xyz)

    @classmethod
    def centred_on(cls, lat_deg, lon_deg) -> LocalFrame:
        """Frame at the mean position of a set of points (good default for a station network)."""
        lat = np.asarray(lat_deg, dtype=float)
        lon = np.asarray(lon_deg, dtype=float)
        # circular mean for longitude so networks straddling +/-180 deg work
        rad = np.radians(lon)
        lon_mean = np.degrees(np.arctan2(np.sin(rad).mean(), np.cos(rad).mean()))
        return cls(float(lat.mean()), float(lon_mean), 0.0)
