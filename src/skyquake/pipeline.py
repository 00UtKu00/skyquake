"""Reusable event pipeline: load seismograms, iterative predict -> pick -> fit.

Requires ``obspy``. Used by the examples in ``examples/``.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from skyquake.geo import LocalFrame
from skyquake.invert import FitResult, fit_trajectory
from skyquake.physics import arrival_times


def sane_elevation(elev_m: float, station: str = "", log=print) -> float:
    """Replace impossible station elevations (metadata placeholders such as 123456) by 0 m.

    Real stations lie between about -12 km (deep boreholes) and +9 km. Network 5R, for
    example, publishes 123456.0 as a dummy value; trusted blindly, that puts a station
    123 km up, above the re-entering object.
    """
    if elev_m is None or not np.isfinite(elev_m) or not (-12_000.0 <= elev_m <= 9_000.0):
        log(f"  warning: {station} has impossible elevation {elev_m!r} m; using 0 m")
        return 0.0
    return float(elev_m)


@dataclass
class StationTrace:
    id: str
    lat: float
    lon: float
    elev: float
    t: np.ndarray      # seconds after the event reference time
    cft: np.ndarray    # STA/LTA characteristic function
    data: np.ndarray   # band-passed samples


def load_traces(raw_dir: str | Path, t_ref, t_range=(0.0, 1500.0), band=(2.0, 20.0),
                sta_s=0.3, lta_s=10.0) -> list[StationTrace]:
    """Read every ``*.mseed`` in ``raw_dir`` with the StationXML files found there."""
    from obspy import Inventory, UTCDateTime, read, read_inventory
    from obspy.signal.trigger import classic_sta_lta

    raw = Path(raw_dir)
    t_ref = UTCDateTime(t_ref)
    inv = Inventory()
    for f in raw.glob("*.xml"):
        inv += read_inventory(str(f))
    out = []
    for f in sorted(raw.glob("*.mseed")):
        st = read(str(f))
        st.merge(fill_value=0)
        tr = st[0]
        try:
            c = inv.select(network=tr.stats.network, station=tr.stats.station)[0][0]
        except IndexError:
            continue
        tr.trim(t_ref + t_range[0] - lta_s - 30, t_ref + t_range[1])
        df = tr.stats.sampling_rate
        if tr.stats.npts < 1000 or df < 2.5 * band[1]:
            continue
        tr.detrend("demean")
        tr.taper(0.02)
        tr.filter("bandpass", freqmin=band[0], freqmax=band[1], corners=4, zerophase=True)
        cft = classic_sta_lta(tr.data, int(sta_s * df), int(lta_s * df))
        cft[: int(lta_s * df)] = 0
        t = (tr.stats.starttime - t_ref) + np.arange(tr.stats.npts) / df
        sid = f"{tr.stats.network}.{tr.stats.station}"
        out.append(StationTrace(sid, c.latitude, c.longitude, sane_elevation(c.elevation, sid),
                                t, cft, tr.data.astype(float)))
    return out


def pick_onset(s: StationTrace, t_pred: float, half: float, on: float = 3.0):
    """Onset of the strongest STA/LTA peak within t_pred +/- half -> (t, snr) or None."""
    m = (s.t >= t_pred - half) & (s.t <= t_pred + half)
    if m.sum() < 10:
        return None
    idx = np.flatnonzero(m)
    i_max = idx[np.argmax(s.cft[idx])]
    if s.cft[i_max] < on:
        return None
    i = i_max
    while i > idx[0] and s.cft[i] > on:
        i -= 1
    df = 1.0 / (s.t[1] - s.t[0])
    pre = s.data[max(0, i - int(10 * df)):i]
    post = s.data[i:i + int(3 * df)]
    if len(pre) < 10 or len(post) < 10:
        return None
    return float(s.t[i]), float(np.abs(post).max() / (np.std(pre) * 3 + 1e-12))


@dataclass
class PickFit:
    fit: FitResult
    frame: LocalFrame
    enu: np.ndarray
    t: np.ndarray
    snr: np.ndarray
    sigma: np.ndarray
    picked: np.ndarray      # bool, station has a usable pick
    used: np.ndarray        # indices into stations of picks used by the final fit
    pred: np.ndarray        # final model prediction for every station


def iterative_pick_fit(stations: list[StationTrace], initial_prediction: Callable,
                       windows=(30.0, 10.0, 5.0), min_snr=2.0, log=print,
                       **fit_kw) -> PickFit:
    """Pick around a prediction, fit, re-predict, narrow the window, repeat.

    ``initial_prediction(lat, lon)`` returns the first guess arrival (s) per station.
    """
    lat = np.array([s.lat for s in stations])
    lon = np.array([s.lon for s in stations])
    frame = LocalFrame.centred_on(lat, lon)
    enu = frame.to_enu(lat, lon, np.array([s.elev for s in stations]))
    pred = np.asarray(initial_prediction(lat, lon), dtype=float)
    fit = None
    kw = dict(fit_sound_speed=True) | fit_kw
    for it, half in enumerate(windows):
        picks = [pick_onset(s, p, half) for s, p in zip(stations, pred, strict=True)]
        ok = np.array([p is not None and p[1] >= min_snr for p in picks])
        t = np.array([p[0] if p else np.nan for p in picks])
        snr = np.array([p[1] if p else np.nan for p in picks])
        sig = np.where(ok, np.clip(2.0 / np.nan_to_num(snr, nan=1.0), 0.15, 2.0), 1.0)
        if it and "azimuth_guess_deg" not in fit_kw:
            kw["azimuth_guess_deg"] = fit.trajectory.azimuth_deg
        fit = fit_trajectory(enu[ok], t[ok], sig[ok], **kw)
        pred = arrival_times(fit.trajectory, enu, fit.sound_speed)
        log(f"iteration {it}: window +/-{half:.0f} s, {ok.sum()} picks, {fit.n_used} used, "
            f"rms {fit.rms_s:.2f} s")
    used = np.flatnonzero(ok)[fit.inliers]
    return PickFit(fit, frame, enu, t, snr, sig, ok, used, pred)
