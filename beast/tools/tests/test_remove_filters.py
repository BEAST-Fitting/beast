import pytest
from beast.tools import remove_filters, read_beast_data
from beast.tests.helpers import download_rename
import os


@pytest.mark.remote_data
def test_remove_filters():
    """
    Test for remove_filters.py

    The SED grid for this test has two entries for F475W: HST_ACS_WFC_F475W and
    HST_WFC3_F475W (the actual F475W observations are with ACS). This tests four
    combinations of running remove_filters:

    1. Use the catalog to choose which filters in the SED grid to keep
        A. without using beast_filt keyword: the code doesn't know what F475W in
           the catalog means, so it won't delete either F475W entry from the
           grid
        B. using beast_filt='HST_ACS_WFC_F475W': the code knows that when it
           sees F475W in the catalog, it means ACS, so it should remove
           HST_WFC3_F475W from the grid

    2. Use the rm_filters='F475W' keyword to choose filter removal
        A. without using beast_filt keyword: any time the code sees F475W in the
           grid, it will delete it
        B. using beast_filt='HST_WFC3_F475W': the code knows that F475W in
           rm_filters means WFC3, so it will only delete the F475W WFC3 entry in
           the grid
    """

    # download the needed files
    obs_fname = download_rename("phat_small/b15_4band_det_27_A.fits")
    seds_fname = download_rename("phat_small/beast_example_phat_seds_extrafilter.hd5")

    # name to use for the output grid
    temp_physgrid_file = "temp_newgrid.hd5"

    # ==== case 1A ====

    # run filter removal
    remove_filters.remove_filters_from_files(
        obs_fname,
        physgrid=seds_fname,
        physgrid_outfile=temp_physgrid_file,
        # beast_filt=['HST_ACS_WFC_F475W'],
    )

    # check that the proper filters are retained
    expected_filters = [
        "HST_WFC3_F275W",
        "HST_WFC3_F336W",
        "HST_ACS_WFC_F475W",
        "HST_WFC3_F475W",
        "HST_ACS_WFC_F814W",
        "HST_WFC3_F110W",
        "HST_WFC3_F160W",
    ]
    temp = read_beast_data.read_sed_data(temp_physgrid_file, param_list=["filters"])
    assert set(temp["filters"]) == set(
        expected_filters
    ), "remove_filters case 1A doesn't match"

    # remove temp file
    os.remove(temp_physgrid_file)

    # ==== case 1B ====

    # run filter removal
    remove_filters.remove_filters_from_files(
        obs_fname,
        physgrid=seds_fname,
        physgrid_outfile=temp_physgrid_file,
        beast_filt=["HST_ACS_WFC_F475W"],
    )

    # check that the proper filters are retained
    expected_filters = [
        "HST_WFC3_F275W",
        "HST_WFC3_F336W",
        "HST_ACS_WFC_F475W",
        # 'HST_WFC3_F475W',
        "HST_ACS_WFC_F814W",
        "HST_WFC3_F110W",
        "HST_WFC3_F160W",
    ]
    temp = read_beast_data.read_sed_data(temp_physgrid_file, param_list=["filters"])
    assert set(temp["filters"]) == set(
        expected_filters
    ), "remove_filters case 1B doesn't match"

    # remove temp file
    os.remove(temp_physgrid_file)

    # ==== case 2A ====

    # run filter removal
    remove_filters.remove_filters_from_files(
        obs_fname,
        physgrid=seds_fname,
        physgrid_outfile=temp_physgrid_file,
        rm_filters=["F475W"],
        # beast_filt=['HST_WFC3_F475W'],
    )

    # check that the proper filters are retained
    expected_filters = [
        "HST_WFC3_F275W",
        "HST_WFC3_F336W",
        # 'HST_ACS_WFC_F475W',
        # 'HST_WFC3_F475W',
        "HST_ACS_WFC_F814W",
        "HST_WFC3_F110W",
        "HST_WFC3_F160W",
    ]
    temp = read_beast_data.read_sed_data(temp_physgrid_file, param_list=["filters"])
    assert set(temp["filters"]) == set(
        expected_filters
    ), "remove_filters case 2A doesn't match"

    # remove temp file
    os.remove(temp_physgrid_file)

    # ==== case 2B ====

    # run filter removal
    remove_filters.remove_filters_from_files(
        obs_fname,
        physgrid=seds_fname,
        physgrid_outfile=temp_physgrid_file,
        rm_filters=["F475W"],
        beast_filt=["HST_WFC3_F475W"],
    )

    # check that the proper filters are retained
    expected_filters = [
        "HST_WFC3_F275W",
        "HST_WFC3_F336W",
        "HST_ACS_WFC_F475W",
        # 'HST_WFC3_F475W',
        "HST_ACS_WFC_F814W",
        "HST_WFC3_F110W",
        "HST_WFC3_F160W",
    ]
    temp = read_beast_data.read_sed_data(temp_physgrid_file, param_list=["filters"])
    assert set(temp["filters"]) == set(
        expected_filters
    ), "remove_filters case 2B doesn't match"

    # remove temp file
    os.remove(temp_physgrid_file)
