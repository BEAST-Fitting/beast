from beast.tests.helpers import download_rename, compare_hdf5, compare_tables
from beast.tools.tests import beast_settings_mockup
from beast.tools.run import create_physicsmodel
from astropy.tests.helper import remote_data
import os
from astropy.table import Table


@remote_data
def test_create_physicsmodel_no_subgrid():
    """
    Test create_physicsmodel, assuming no subgrids
    """

    # get beast settings mockup
    settings = beast_settings_mockup.beast_settings_mockup()

    # download files
    # - input files
    download_rename("filters.hd5")
    # - intermediate files
    iso_fname_cache = download_rename("beast_example_phat_iso.csv")
    priors_fname_cache = download_rename("beast_example_phat_spec_w_priors.grid.hd5")
    # - anticipated output files
    seds_fname_cache = download_rename("beast_example_phat_seds.grid.hd5")

    # run create_physicsmodel
    create_physicsmodel.create_physicsmodel(
        settings, nsubs=settings.n_subgrid, nprocs=1
    )

    # check that files match
    # - isochrones
    table_cache = Table.read(
        iso_fname_cache, format="ascii.csv", comment="#", delimiter=",",
    )
    table_new = Table.read(
        "./beast_no_subgrid/beast_no_subgrid_iso.csv",
        format="ascii.csv",
        comment="#",
        delimiter=",",
    )
    compare_tables(table_cache, table_new)
    # - spectra with priors
    compare_hdf5(
        priors_fname_cache, "./beast_no_subgrid/beast_no_subgrid_spec_w_priors.grid.hd5"
    )
    # - SEDs grid
    compare_hdf5(seds_fname_cache, "./beast_no_subgrid/beast_no_subgrid_seds.grid.hd5")


@remote_data
def test_create_physicsmodel_with_subgrid():
    """
    Test create_physicsmodel, assuming two subgrids
    """

    # get beast settings mockup
    settings = beast_settings_mockup.beast_settings_mockup()
    # edit to have 2 subgrids
    settings.n_subgrid = 2

    # download files
    # - input files
    download_rename("filters.hd5")
    # - intermediate files
    iso_fname_cache = download_rename("beast_example_phat_iso.csv")
    priors_fname_cache = download_rename("beast_example_phat_spec_w_priors.grid.hd5")
    priors_sub0_fname_cache = download_rename(
        "beast_example_phat_subgrids_spec_w_priors.gridsub0.hd5"
    )
    priors_sub1_fname_cache = download_rename(
        "beast_example_phat_subgrids_spec_w_priors.gridsub1.hd5"
    )
    # - anticipated output files
    seds_sub0_fname_cache = download_rename(
        "beast_example_phat_subgrids_seds.gridsub0.hd5"
    )
    seds_sub1_fname_cache = download_rename(
        "beast_example_phat_subgrids_seds.gridsub1.hd5"
    )

    # run create_physicsmodel
    create_physicsmodel.create_physicsmodel(
        settings, nsubs=settings.n_subgrid, nprocs=1
    )

    # check that files match

    # - isochrones
    table_cache = Table.read(
        iso_fname_cache, format="ascii.csv", comment="#", delimiter=",",
    )
    table_new = Table.read(
        "beast_example_phat_subgrids/beast_example_phat_subgrids_iso.csv",
        format="ascii.csv",
        comment="#",
        delimiter=",",
    )
    compare_tables(table_cache, table_new)

    # - spectra with priors
    compare_hdf5(
        priors_fname_cache,
        "./beast_example_phat_subgrids/beast_example_phat_subgrids_spec_w_priors.grid.hd5",
    )
    compare_hdf5(
        priors_sub0_fname_cache,
        "beast_example_phat_subgrids/beast_example_phat_subgrids_spec_w_priors.gridsub0.hd5",
    )
    compare_hdf5(
        priors_sub1_fname_cache,
        "beast_example_phat_subgrids/beast_example_phat_subgrids_spec_w_priors.gridsub1.hd5",
    )

    # - SEDs grid
    compare_hdf5(
        seds_sub0_fname_cache,
        "beast_example_phat_subgrids/beast_example_phat_subgrids_seds.gridsub0.hd5",
    )
    compare_hdf5(
        seds_sub1_fname_cache,
        "beast_example_phat_subgrids/beast_example_phat_subgrids_seds.gridsub1.hd5",
    )

    # - list of subgrids
    with open("./beast_example_phat_subgrids/subgrid_fnames.txt") as f:
        temp = f.read()
    subgrid_list = [x for x in temp.split("\n") if x != ""]
    expected_list = [
        "beast_example_phat_subgrids/beast_example_phat_subgrids_seds.gridsub0.hd5",
        "beast_example_phat_subgrids/beast_example_phat_subgrids_seds.gridsub1.hd5",
    ]
    assert subgrid_list == expected_list, "subgrid_fnames.txt has incorrect content"
