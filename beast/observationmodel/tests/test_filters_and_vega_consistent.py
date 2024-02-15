from astropy.table import QTable

from beast.config import __ROOT__
from beast.observationmodel.vega import Vega
from beast.observationmodel import phot

def test_filters_and_vega_consistent():
    """
    Test to ensure that the filters.hd5 and vega.hd5 are consistent.
    In other words, both have the same filters.
    """
    ftab = QTable.read(__ROOT__ + "filters.hd5", path="content")
    vtab = QTable.read(__ROOT__ + "vega.hd5", path="sed")

    otxt = ""
    for cfilt in ftab["TABLENAME"].data:
        if cfilt not in vtab["FNAME"].data:
            otxt = f"{otxt} {cfilt}"
    assert otxt == "", "filters in filters.hd5 missing from vega.hd5:" + otxt

    otxt = ""
    for cfilt in vtab["FNAME"].data:
        if cfilt not in ftab["TABLENAME"].data:
            otxt = f"{otxt} {cfilt}"
    assert otxt == "", "filters in vega.hd5 missing from filters.hd5:" + otxt

if __name__ == "__main__":  # pragma: no cover
    test_filters_and_vega_consistent()