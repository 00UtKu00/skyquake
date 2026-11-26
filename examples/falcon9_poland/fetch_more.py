"""Fetch networks the EIDA routing did not return: GEOFON (GE, DK...), BGR (GR),
Raspberry Shake citizen seismometers (AM). Saves into ./raw alongside fetch.py output."""

from pathlib import Path

from obspy import Inventory, UTCDateTime, read_inventory
from obspy.clients.fdsn import Client

OUT = Path(__file__).with_name("raw")
T1 = UTCDateTime("2025-02-19T03:42:00")
T2 = UTCDateTime("2025-02-19T04:08:00")
BOX = dict(minlatitude=50.5, maxlatitude=56.5, minlongitude=3.0, maxlongitude=19.5)
SOURCES = [("GFZ", "*", "HHZ,BHZ"), ("BGR", "*", "HHZ,BHZ"), ("RASPISHAKE", "AM", "EHZ,SHZ")]


def main() -> None:
    inv_path = OUT / "inventory.xml"
    inv = read_inventory(str(inv_path)) if inv_path.exists() else Inventory()
    have = {f.name.split(".")[0] + "." + f.name.split(".")[1] for f in OUT.glob("*.mseed")}
    for provider, net, cha in SOURCES:
        try:
            c = Client(provider, timeout=120)
            new = c.get_stations(network=net, channel=cha, starttime=T1, endtime=T2,
                                 level="channel", **BOX)
        except Exception as e:  # noqa: BLE001
            print(f"{provider}: stations failed: {type(e).__name__}: {e}")
            continue
        todo = []
        for n in new:
            for s in n:
                key = f"{n.code}.{s.code}"
                if key in have:
                    continue
                chans = sorted({ch.code for ch in s}, key=lambda x: ("HBES".index(x[0]), x))
                todo.append((n.code, s.code, "*", chans[0], T1, T2))
        print(f"{provider}: {len(todo)} new stations")
        inv += new
        got = 0
        for i in range(0, len(todo), 40):
            try:
                st = c.get_waveforms_bulk(todo[i:i + 40])
            except Exception as e:  # noqa: BLE001
                print(f"  batch {i}: {type(e).__name__}")
                continue
            for sid in sorted({tr.id for tr in st}):
                st.select(id=sid).write(str(OUT / f"{sid}.mseed"), format="MSEED")
                have.add(".".join(sid.split(".")[:2]))
                got += 1
        print(f"  saved {got}")
    inv.write(str(inv_path), format="STATIONXML")


if __name__ == "__main__":
    main()
