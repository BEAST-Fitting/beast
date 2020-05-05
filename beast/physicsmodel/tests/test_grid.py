import numpy as np
from tempfile import NamedTemporaryFile
from astropy.table import Table

from beast.physicsmodel.grid import SEDGrid
from beast.tests.helpers import compare_tables


def test_sedgrid():
    """
    Simple tests of the SpectralGrid class
    """
    n_bands = 3
    filter_names = ["BAND1", "BAND2", "BAND3"]
    n_models = 100
    lamb = [1.0, 2.0, 3.0]
    seds = np.zeros((n_models, n_bands))
    cols = {"Av": [1.0, 1.1, 1.3], "Rv": [2.0, 3.0, 4.0]}
    header = {"Origin": "test_code"}
    gtable = Table(cols)
    gtable.meta = header
    tgrid = SEDGrid(lamb, seds=seds, grid=gtable, header=header, backend="memory")
    tgrid.header["filters"] = " ".join(filter_names)

    # check that the grid has the expected properties
    assert hasattr(tgrid, "lamb"), "grid missing lambda property"
    assert hasattr(tgrid, "seds"), "grid missing seds property"
    assert hasattr(tgrid, "header"), "grid missing header property"
    assert hasattr(tgrid, "grid"), "grid missing grid property"
    assert hasattr(tgrid, "nbytes"), "grid missing nbytes property"

    np.testing.assert_allclose(tgrid.lamb, lamb, err_msg="grid lambdas not equal")
    np.testing.assert_allclose(tgrid.seds, seds, err_msg="grid seds not equal")
    assert isinstance(tgrid.nbytes, int), "grid nbytes property not int"
    compare_tables(tgrid.grid, gtable)
    assert tgrid.grid.keys() == list(cols.keys()), "colnames of grid not equal"

    # test writing and reading to disk in different formats
    fileformats = [".fits", ".hdf5"]
    for cformat in fileformats:
        print(f"testing {cformat}")
        tfile = NamedTemporaryFile(suffix=cformat)
        tgrid.write(tfile.name)

        # read in the HD5 file
        dgrid = SEDGrid(tfile.name, backend="memory")

        # check that the grid has the expected values
        np.testing.assert_allclose(dgrid.lamb, lamb, err_msg=f"{cformat} file grid lambdas not equal")
        np.testing.assert_allclose(dgrid.seds, seds, err_msg=f"{cformat} file grid seds not equal")
        assert isinstance(dgrid.nbytes, int), f"{cformat} file grid nbytes property not int"
        compare_tables(dgrid.grid, gtable, otag=f"{cformat} file")
        assert dgrid.grid.keys() == list(cols.keys()), f"{cformat} file colnames of grid not equal"


if __name__ == "__main__":
    test_sedgrid()
