import os

import numpy as np
from astropy.tests.helper import remote_data

from beast.tools import subgridding_tools
from beast.tests.helpers import download_rename
from beast.physicsmodel.grid import FileSEDGrid


def split_and_check(grid_fname, num_subgrids):
    complete_g = FileSEDGrid(grid_fname)
    sub_fnames = subgridding_tools.split_grid(grid_fname, num_subgrids)

    # count the number of grid cells
    sub_seds = []
    sub_grids = []

    for sub_fname in sub_fnames:
        sub_g = FileSEDGrid(sub_fname)

        sub_seds.append(sub_g.seds)
        sub_grids.append(sub_g.grid.data)

        np.testing.assert_equal(complete_g.lamb, sub_g.lamb)
        assert complete_g.grid.columns.items() == sub_g.grid.columns.items()

    sub_seds_reconstructed = np.concatenate(sub_seds)
    np.testing.assert_equal(sub_seds_reconstructed, complete_g.seds)

    sub_grids_reconstructed = np.concatenate(sub_grids)
    np.testing.assert_equal(sub_grids_reconstructed, complete_g.grid.data)

    # the split method skips anything that already exists, so if we
    # want to use this function multiple times for the same test
    # grid, we need to do this.
    for f in sub_fnames:
        os.remove(f)


@remote_data
def test_split_grid():
    seds_trim_fname = download_rename('beast_example_phat_seds_trim.grid.hd5')
    split_and_check(seds_trim_fname, 1)  # an edge case
    split_and_check(seds_trim_fname, 3)  # an odd numer
    split_and_check(seds_trim_fname, 4)  # an even number


@remote_data
def test_reduce_grid_info():
    seds_trim_fname = download_rename('beast_example_phat_seds_trim.grid.hd5')
    sub_fnames = subgridding_tools.split_grid(seds_trim_fname, 3)

    complete_g_info = subgridding_tools.subgrid_info(seds_trim_fname)
    cap_unique = 50
    sub_g_info = subgridding_tools.reduce_grid_info(
        sub_fnames, nprocs=3, cap_unique=cap_unique)

    for q in complete_g_info:
        assert q in sub_g_info
        assert complete_g_info[q]['min'] == sub_g_info[q]['min']
        assert complete_g_info[q]['max'] == sub_g_info[q]['max']
        num_unique = len(complete_g_info[q]['unique'])
        if num_unique > cap_unique:
            # Cpan still be larger if one of the sub results during the
            # reduction is larger. This is as intended.
            assert sub_g_info[q]['num_unique'] >= cap_unique
        else:
            assert sub_g_info[q]['num_unique'] == num_unique
