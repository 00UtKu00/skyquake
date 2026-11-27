# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 SKYQUAKE contributors
"""Download raw seismograms for the Shenzhou-15 orbital module re-entry (2024-04-02).

Saves one MiniSEED file per station plus a StationXML inventory to ./raw (git-ignored),
so the analysis can be re-run offline and exactly reproduced.

    python examples/shenzhou15/fetch.py
"""

from pathlib import Path

from obspy import UTCDateTime
from obspy.clients.fdsn import Client

OUT = Path(__file__).with_name("raw")
T1 = UTCDateTime("2024-04-02T08:35:00")
T2 = UTCDateTime("2024-04-02T09:00:00")
CENTER = (35.0, -117.5)   # between Santa Barbara and Las Vegas
RADIUS_DEG = 3.5
PROVIDERS = [("SCEDC", "CI"), ("IRIS", "NN")]  # S. California, Nevada


def main() -> None:
    OUT.mkdir(exist_ok=True)
    for provider, net in PROVIDERS:
        client = Client(provider, timeout=120)
        inv = client.get_stations(network=net, channel="HHZ", starttime=T1, endtime=T2,
                                  latitude=CENTER[0], longitude=CENTER[1],
                                  maxradius=RADIUS_DEG, level="channel")
        inv.write(str(OUT / f"{net}_inventory.xml"), format="STATIONXML")
        stations = [(n.code, s.code) for n in inv for s in n]
        print(f"{provider} {net}: {len(stations)} stations with HHZ")
        bulk = [(n, s, "*", "HHZ", T1, T2) for n, s in stations]
        for i in range(0, len(bulk), 40):
            try:
                st = client.get_waveforms_bulk(bulk[i:i + 40])
            except Exception as e:  # noqa: BLE001
                print(f"  batch {i}: {type(e).__name__}: {e}")
                continue
            # several segments per station when there are gaps: keep them all in one file
            for sid in sorted({tr.id for tr in st}):
                st.select(id=sid).write(str(OUT / f"{sid}.mseed"), format="MSEED")
            print(f"  batch {i}: {len(st)} traces")


if __name__ == "__main__":
    main()
