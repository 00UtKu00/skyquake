# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 SKYQUAKE contributors
"""Download raw seismograms for the Falcon 9 upper-stage re-entry (2025-02-19, Europe).

Optical reference (Vierinen et al. 2026, arXiv:2605.29124): first seen 03:43:00 UTC over
the UK; fragment families crossed 60 km at ~03:45:51 and ~03:46:22 UTC; last fragment
at 40 km at 03:47:25 UTC; debris recovered in western Poland.

Saves one MiniSEED file per station plus StationXML to ./raw (git-ignored).

    python examples/falcon9_poland/fetch.py
"""

from pathlib import Path

from obspy import Inventory, UTCDateTime
from obspy.clients.fdsn import RoutingClient

OUT = Path(__file__).with_name("raw")
T1 = UTCDateTime("2025-02-19T03:42:00")
T2 = UTCDateTime("2025-02-19T04:08:00")
BOX = dict(minlatitude=50.5, maxlatitude=56.5, minlongitude=3.0, maxlongitude=19.5)


def main() -> None:
    OUT.mkdir(exist_ok=True)
    client = RoutingClient("eida-routing", timeout=180)
    inv = Inventory()
    for cha in ("HHZ", "BHZ", "EHZ"):
        try:
            inv += client.get_stations(channel=cha, starttime=T1, endtime=T2, level="channel",
                                       **BOX)
        except Exception as e:  # noqa: BLE001
            print(f"stations {cha}: {type(e).__name__}: {e}")
    # one vertical channel per station, preferring HH > BH > EH
    best: dict[tuple[str, str], str] = {}
    rank = {"HHZ": 0, "BHZ": 1, "EHZ": 2}
    for net in inv:
        for sta in net:
            for ch in sta:
                key = (net.code, sta.code)
                if key not in best or rank[ch.code] < rank[best[key]]:
                    best[key] = ch.code
    inv.write(str(OUT / "inventory.xml"), format="STATIONXML")
    print(f"{len(best)} stations: " + ", ".join(sorted({k[0] for k in best})))

    bulk = [(n, s, "*", c, T1, T2) for (n, s), c in sorted(best.items())]
    got = 0
    for i in range(0, len(bulk), 50):
        try:
            st = client.get_waveforms_bulk(bulk[i:i + 50])
        except Exception as e:  # noqa: BLE001
            print(f"  batch {i}: {type(e).__name__}: {e}")
            continue
        for sid in sorted({tr.id for tr in st}):
            st.select(id=sid).write(str(OUT / f"{sid}.mseed"), format="MSEED")
            got += 1
        print(f"  batch {i}: {len(st)} traces")
    print(f"saved {got} stations")


if __name__ == "__main__":
    main()
