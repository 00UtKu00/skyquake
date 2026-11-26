import warnings

import numpy as np

from skyquake.io import read_picks


def test_placeholder_elevations_are_replaced(tmp_path):
    f = tmp_path / "p.csv"
    f.write_text(
        "station,lat,lon,elev_m,arrival_utc,sigma_s\n"
        "XX.GOOD,52.0,8.0,45.0,2025-02-19T03:49:00Z,0.3\n"
        "5R.NIKEB,52.1,8.3,123456.0,2025-02-19T03:49:05Z,0.3\n"
        "XX.BORE,53.2,6.8,-200.0,2025-02-19T03:49:10Z,0.3\n",
        encoding="utf-8",
    )
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        p = read_picks(f)
    np.testing.assert_allclose(p.elev_m, [45.0, 0.0, -200.0])
    assert any("123456" in str(x.message) for x in w)
