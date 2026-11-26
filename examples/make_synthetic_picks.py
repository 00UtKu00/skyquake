"""Regenerate examples/synthetic_picks.csv from the built-in synthetic scenario."""

from datetime import UTC, datetime, timedelta
from pathlib import Path

from skyquake.cli import synthetic_scenario
from skyquake.io import write_picks

frame, st, t, truth = synthetic_scenario(seed=1)
lat, lon, h = frame.from_enu(st)
epoch = datetime(2030, 1, 1, 12, 0, 0, tzinfo=UTC)  # fictional date: this is NOT a real event
rows = [
    {
        "station": f"SY.S{i:03d}",
        "lat": round(float(la), 5),
        "lon": round(float(lo), 5),
        "elev_m": round(float(hh), 1),
        "arrival_utc": (epoch + timedelta(seconds=float(ti))).isoformat().replace("+00:00", "Z"),
        "sigma_s": 0.3,
        "snr": "",
    }
    for i, (la, lo, hh, ti) in enumerate(zip(lat, lon, h, t, strict=True))
]
out = Path(__file__).with_name("synthetic_picks.csv")
write_picks(out, rows)
print(f"wrote {out} ({len(rows)} picks); truth: {truth}")
