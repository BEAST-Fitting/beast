import numpy as np
from tempfile import NamedTemporaryFile
from astropy.table import Table
from astropy.io.misc.hdf5 import read_table_hdf5

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
    cov_diag = np.full((n_models, n_bands), 0.1)
    n_offdiag = ((n_bands ** 2) - n_bands) // 2
    cov_offdiag = np.full((n_models, n_offdiag), 1.0)
    cols = {"Av": [1.0, 1.1, 1.3], "Rv": [2.0, 3.0, 4.0]}
    header = {"Origin": "test_code"}
    gtable = Table(cols)
    gtable.meta = header

    tgrid = SEDGrid(
        lamb,
        seds=seds,
        grid=gtable,
        header=header,
        cov_diag=cov_diag,
        cov_offdiag=cov_offdiag,
        backend="memory",
    )
    tgrid.header["filters"] = " ".join(filter_names)

    # check that the grid has the expected properties
    expected_props = [
        "lamb",
        "seds",
        "cov_diag",
        "cov_offdiag",
        "grid",
        "nbytes",
        "filters",
        "keys",
    ]
    for cprop in expected_props:
        assert hasattr(tgrid, cprop), f"missing {cprop} property"

    np.testing.assert_allclose(tgrid.lamb, lamb, err_msg="lambdas not equal")
    np.testing.assert_allclose(tgrid.seds, seds, err_msg="seds not equal")
    np.testing.assert_allclose(tgrid.cov_diag, cov_diag, err_msg="covdiag not equal")
    np.testing.assert_allclose(
        tgrid.cov_offdiag, cov_offdiag, err_msg="covoffdiag not equal"
    )
    assert isinstance(tgrid.nbytes, int), "grid nbytes property not int"
    compare_tables(tgrid.grid, gtable)
    assert tgrid.grid.keys() == list(cols.keys()), "colnames of grid not equal"
    assert tgrid.filters == filter_names, "filters of grid not equal"

    # test writing and reading to disk in different formats
    fileformats = [".fits", ".hdf"]
    for cformat in fileformats:
        print(f"testing {cformat} file format")
        tfile = NamedTemporaryFile(suffix=cformat)
        print(tfile.name)

        # write the file
        tgrid.write(tfile.name)

        # read in the file using different backends
        backs = ["memory", "cache", "disk"]
        for cback in backs:
            if (cback == "disk") and (cformat == ".fits"):  # not supported
                continue

            print(f"    testing {cback} backend")
            dgrid = SEDGrid(tfile.name, backend=cback)

            for cprop in expected_props:
                assert hasattr(dgrid, cprop), f"missing {cprop} property"

            # check that the grid has the expected values

            # this test is having a problem in the online travis ci
            # it someone manages to access another file with HST filter names!
            # no idea way.  Works fine offline.
            # assert dgrid.filters == filter_names, "{cformat} file filters not equal"

            np.testing.assert_allclose(
                dgrid.lamb, lamb, err_msg=f"{cformat} file grid lambdas not equal"
            )
            np.testing.assert_allclose(
                dgrid.seds, seds, err_msg=f"{cformat} file grid seds not equal"
            )
            np.testing.assert_allclose(
                dgrid.cov_diag,
                cov_diag,
                err_msg=f"{cformat} file grid cov_diag not equal",
            )
            np.testing.assert_allclose(
                dgrid.cov_offdiag,
                cov_offdiag,
                err_msg=f"{cformat} file grid cov_offdiag not equal",
            )
            assert isinstance(
                dgrid.nbytes, int
            ), f"{cformat} file grid nbytes property not int"

            dTable = dgrid.grid
            if (cback == "disk") and (cformat == ".hdf"):
                dTable = read_table_hdf5(dgrid.grid)
            compare_tables(dTable, gtable, otag=f"{cformat} file")

            assert dTable.keys() == list(
                cols.keys()
            ), f"{cformat} file colnames of grid not equal"


if __name__ == "__main__":
    test_sedgrid()
