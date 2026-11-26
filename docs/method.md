# Method

## 1. Forward model

The object moves on a straight line `p(t) = p0 + v·u·(t − t0)` in a local East-North-Up
frame. Here `u` is a unit vector set by the heading (azimuth, clockwise from north) and the
descent angle (flight-path angle, positive downward), and `p0` is the point at reference
height `h_ref` (default 40 km).

A station at position `q` receives sound emitted at along-track coordinate `s` at time

    T(s) = t0 + s/v + |q − (p0 + s·u)| / c

The first arrival is the stationary point `dT/ds = 0` (Fermat's principle). Writing `x` for
the station's along-track coordinate and `r` for its perpendicular distance from the line:

    t_arrival = t0 + x/v + r·sqrt(1/c² − 1/v²)        (v > c)

The emission point is at `s* = x − r·c / sqrt(v² − c²)`. The geometry is the Mach cone,
with half-angle `asin(c/v)`. `tests/test_physics.py` checks the closed form against direct
numerical minimisation of `T(s)` to 0.1 ms.

## 2. Geodesy

Stations are converted from WGS84 geodetic coordinates to ECEF and then to the ENU frame
exactly, with no flat-Earth approximation. Earth curvature is therefore represented
correctly: a station 300 km away sits about 7 km below the tangent plane.

## 3. Inverse problem

Parameters: `east, north` (of the reference point), `azimuth`, `fpa`, `v`, `t0`, and
optionally the path-averaged sound speed `c`.

1. **Multi-start robust fit.** 12 headings × 2 descent angles × 2 speeds, using a soft-L1
   loss (Huber-like), so bad picks cannot drag the solution. A heading prior from a TLE
   narrows the search when one is available.
2. **Outlier rejection.** Picks with normalised residuals above 4σ are dropped.
3. **Final fit.** Ordinary least squares on the inliers. The covariance is
   `(JᵀJ)⁻¹ · max(1, χ²_red)`, so it is inflated when the data fit worse than the stated
   pick errors and never shrunk.
4. **Corridor.** Covariance samples are pushed to the end of the hypersonic segment
   (default 15 km height) and summarised as along-track and cross-track 1σ.

## 4. Assumptions and their error budget

| Assumption | Effect when wrong | Status |
|---|---|---|
| Straight line | Curvature over a few hundred km is small for steep, fast segments | Acceptable for v0.1 |
| Constant speed | Deceleration below about 50 km biases `v` and `fpa` | Roadmap: piecewise or decelerating model |
| Homogeneous sound speed | Mitigated by fitting `c`. Optional layered atmosphere (`--layered`) traces exact rays through the US Standard Atmosphere 1976; on Shenzhou-15 it changes speed by < 0.2 km/s | Winds: roadmap (ERA5) |
| Single line source | Fragmentation produces several boom trains | Roadmap: multi-source picking |
| Seismometer = pressure sensor | The air-to-ground coupling delay is negligible for the timing | Acceptable |

## 5. What the result means

The fitted line describes the **hypersonic, boom-producing segment**. Fragments continue
in ballistic "dark flight" below that segment, drifting with the wind. Any debris-recovery
area therefore lies **downrange** of the corridor point and is not yet computed by
SKYQUAKE.
