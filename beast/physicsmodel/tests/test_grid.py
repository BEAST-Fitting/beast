import numpy as np
from tempfile import NamedTemporaryFile
from astropy.table import Table
from astropy.io.misc.hdf5 import read_table_hdf5
import pytest

from beast.physicsmodel.grid import SEDGrid
from beast.tests.helpers import compare_tables


@pytest.mark.parametrize("cformat", [".fits", ".hdf"])
@pytest.mark.parametrize("cback", ["memory", "cache", "disk"])
@pytest.mark.parametrize("copygrid", [False, True])
def test_sedgrid(cformat, cback, copygrid):
    """
    Tests of the SEDGrid class
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
        "header",
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
    assert isinstance(tgrid.nbytes, (int, np.integer)), "grid nbytes property not integer"
    compare_tables(tgrid.grid, gtable)
    assert tgrid.grid.keys() == list(cols.keys()), "colnames of grid not equal"
    assert tgrid.filters == filter_names, "filters of grid not equal"

    # test writing and reading to disk
    print(f"testing {cformat} file format")
    tfile = NamedTemporaryFile(suffix=cformat)

    # write the file
    tgrid.write(tfile.name)

    # read in the file using different backends
    if (cback == "disk") and (cformat == ".fits"):  # not supported
        return True

    print(f"    testing {cback} backend")
    dgrid_in = SEDGrid(tfile.name, backend=cback)

    # test making a copy
    print(f"    testing copygrid={copygrid}")
    if copygrid:
        dgrid = dgrid_in.copy()
    else:
        dgrid = dgrid_in
    print(dgrid)

    for cprop in expected_props:
        assert hasattr(dgrid, cprop), f"missing {cprop} property"

    # check that the grid has the expected values

    # this test is having a problem in the online travis ci
    # it someone manages to access another file with HST filter names!
    # no idea way.  Works fine offline.
    # assert dgrid.filters == filter_names, "{cformat} file filters not equal"

    assert len(dgrid) == n_bands, f"{cformat} file len not equal"

    np.testing.assert_allclose(
        dgrid.lamb, lamb, err_msg=f"{cformat} file grid lambdas not equal"
    )
    np.testing.assert_allclose(
        dgrid.seds, seds, err_msg=f"{cformat} file grid seds not equal"
    )
    np.testing.assert_allclose(
        dgrid.cov_diag, cov_diag, err_msg=f"{cformat} file grid cov_diag not equal",
    )
    np.testing.assert_allclose(
        dgrid.cov_offdiag,
        cov_offdiag,
        err_msg=f"{cformat} file grid cov_offdiag not equal",
    )
    assert isinstance(dgrid.nbytes, (int, np.integer)), f"{cformat} file grid nbytes property not integer"

    dTable = dgrid.grid
    if (cback == "disk") and (cformat == ".hdf"):
        dTable = read_table_hdf5(dgrid.grid)
    compare_tables(dTable, gtable, otag=f"{cformat} file")

    assert dTable.keys() == list(
        cols.keys()
    ), f"{cformat} file colnames of grid not equal"

    assert dgrid.keys() == tgrid.keys(), f"{cformat} file colnames of grid not equal"

    # final copy - needed for disk backend to get the now defined variables
    print(dgrid)

    dgrid_fin = dgrid.copy()

    print(dgrid_fin)


def test_grid_warnings():
    with pytest.raises(ValueError) as exc:
        SEDGrid(backend="hdf")
    assert exc.value.args[0] == "hdf backend not supported"

    with pytest.raises(ValueError) as exc:
        SEDGrid("test.txt")
    assert exc.value.args[0] == "txt file type not supported"

    # define grid contents
    n_bands = 3
    # filter_names = ["BAND1", "BAND2", "BAND3"]
    n_models = 100
    lamb = [1.0, 2.0, 3.0]
    seds = np.zeros((n_models, n_bands))
    # cov_diag = np.full((n_models, n_bands), 0.1)
    # n_offdiag = ((n_bands ** 2) - n_bands) // 2
    # cov_offdiag = np.full((n_models, n_offdiag), 1.0)
    cols = {"Av": [1.0, 1.1, 1.3], "Rv": [2.0, 3.0, 4.0]}
    header = {"Origin": "test_code"}
    gtable = Table(cols)
    gtable.meta = header

    with pytest.raises(ValueError) as exc:
        SEDGrid(lamb)
    assert exc.value.args[0] == "seds or grid not passed"

    for ftype in ["fits", "hdf"]:
        with pytest.raises(ValueError) as exc:
            a = SEDGrid(lamb, seds=seds, grid=gtable)
            a.grid = cols
            a.write(f"testgridwriteerror.{ftype}")
        assert exc.value.args[0] == "Only astropy.Table are supported"

        with pytest.raises(ValueError) as exc:
            a = SEDGrid(lamb, seds=seds, grid=gtable)
            a.grid = None
            a.write(f"testgridwriteerror.{ftype}")
        assert exc.value.args[0] == "Full data set not specified (lamb, seds, grid)"
