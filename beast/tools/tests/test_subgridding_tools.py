import numpy as np

from astropy.tests.helper import remote_data

from beast.tools import subgridding_tools
from beast.tests.helpers import download_rename
from beast.physicsmodel.grid import FileSEDGrid

def split_and_check(grid_fname, num_subgrids):
    sub_fnames = subgridding_tools.split_grid(grid_fname, num_subgrids)

    complete_g = FileSEDGrid(grid_fname)

    # count the number of grid cells
    sub_seds = []
    sub_grids = []

    for sub_fname in sub_fnames:
        sub_g = FileSEDGrid(sub_fname)

        sub_seds.append(sub_g.seds)
        sub_grids.append(sub_g.grid.data)

        np.testing.assert_equal(complete_g.lamb, sub_g.lamb)
        np.testing.assert_equal(complete_g.grid.columns, sub_g.grid.columns)

    sub_seds_reconstructed = np.concatenate(sub_seds)
    np.testing.assert_equal(sub_seds_reconstructed, complete_g.seds)

    sub_grids_reconstructed = np.concatenate(sub_grids)
    np.testing.assert_equal(sub_grids_reconstructed, complete_g.grid.data)

@remote_data
def test_split_grid():
    seds_trim_fname = download_rename('beast_example_phat_seds_trim.grid.hd5')
    split_and_check(seds_trim_fname, 4) # an even number
    split_and_check(seds_trim_fname, 3) # an odd numer
    split_and_check(seds_trim_fname, 1) # an edge case

@remote_data
def test_reduct_grid_info():
    pass
