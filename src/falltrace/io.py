"""Arrival-pick CSV format.

Columns: station, lat, lon, elev_m, arrival_utc (ISO 8601), sigma_s (optional, default 0.5)
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import numpy as np


@dataclass
class Picks:
    station: list[str]
    lat: np.ndarray
    lon: np.ndarray
    elev_m: np.ndarray
    epoch: datetime
    t_s: np.ndarray  # seconds after epoch
    sigma_s: np.ndarray

    def __len__(self) -> int:
        return len(self.station)


def _parse_utc(s: str) -> datetime:
    dt = datetime.fromisoformat(s.strip().replace("Z", "+00:00"))
    return dt.replace(tzinfo=UTC) if dt.tzinfo is None else dt.astimezone(UTC)


def read_picks(path: str | Path) -> Picks:
    with open(path, newline="", encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(f) if r.get("station")]
    if not rows:
        raise ValueError(f"no picks in {path}")
    times = [_parse_utc(r["arrival_utc"]) for r in rows]
    epoch = min(times)
    return Picks(
        station=[r["station"] for r in rows],
        lat=np.array([float(r["lat"]) for r in rows]),
        lon=np.array([float(r["lon"]) for r in rows]),
        elev_m=np.array([float(r.get("elev_m") or 0.0) for r in rows]),
        epoch=epoch,
        t_s=np.array([(t - epoch).total_seconds() for t in times]),
        sigma_s=np.array([float(r.get("sigma_s") or 0.5) for r in rows]),
    )


def write_picks(path: str | Path, rows: list[dict]) -> None:
    cols = ["station", "lat", "lon", "elev_m", "arrival_utc", "sigma_s", "snr"]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
