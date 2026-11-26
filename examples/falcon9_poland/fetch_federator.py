"""Add stations found through the IRIS federator (DK, GE, YD, CZ, PL, ...) that the
EIDA routing and BGR downloads did not return. Saves into ./raw."""

from pathlib import Path

from obspy import UTCDateTime, read_inventory
from obspy.clients.fdsn import RoutingClient

OUT = Path(__file__).with_name("raw")
T1 = UTCDateTime("2025-02-19T03:42:00")
T2 = UTCDateTime("2025-02-19T04:08:00")
BOX = dict(minlatitude=50.5, maxlatitude=56.5, minlongitude=3.0, maxlongitude=19.5)
RANK = {"HHZ": 0, "BHZ": 1, "EHZ": 2, "SHZ": 3}


def main() -> None:
    client = RoutingClient("iris-federator", timeout=180)
    inv = client.get_stations(starttime=T1, endtime=T2, channel=",".join(RANK), level="channel",
                              **BOX)
    have = {".".join(f.name.split(".")[:2]) for f in OUT.glob("*.mseed")}
    todo = {}
    for n in inv:
        for s in n:
            key = f"{n.code}.{s.code}"
            if key in have:
                continue
            for ch in s:
                if ch.code in RANK and (key not in todo or RANK[ch.code] < RANK[todo[key][3]]):
                    todo[key] = (n.code, s.code, "*", ch.code, T1, T2)
    print(f"{len(todo)} new stations: "
          + ", ".join(sorted({k.split('.')[0] for k in todo})))
    inv.write(str(OUT / "inventory_federator.xml"), format="STATIONXML")
    bulk = sorted(todo.values())
    got = 0
    for i in range(0, len(bulk), 30):
        try:
            st = client.get_waveforms_bulk(bulk[i:i + 30])
        except Exception as e:  # noqa: BLE001
            print(f"  batch {i}: {type(e).__name__}")
            continue
        for sid in sorted({tr.id for tr in st}):
            st.select(id=sid).write(str(OUT / f"{sid}.mseed"), format="MSEED")
            got += 1
    print(f"saved {got}")
    _ = read_inventory  # keep import for users extending this script


if __name__ == "__main__":
    main()
