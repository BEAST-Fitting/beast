from beast.tests.helpers import download_rename, compare_hdf5, compare_tables
from astropy.tests.helper import remote_data
import os
import importlib
from astropy.table import Table

import shutil

@remote_data
def test_create_physicsmodel_no_subgrid():
    """
    Test create_physicsmodel, assuming no subgrids
    """

    # download files
    # - datamodel settings
    #datamodel_fname = download_rename("datamodel_no_subgrid.py")
    #os.rename(datamodel_fname, "datamodel.py")
    shutil.copy2('/astro/dust_kg3/lhagen/stsci/beast-examples_lea-hagen/testing/datamodel_no_subgrid.py', './datamodel.py')
    # - input files (referenced in datamodel)
    filter_fname = download_rename("filters.hd5")
    # - intermediate files
    iso_fname_cache = download_rename("beast_example_phat_iso.csv")
    priors_fname_cache = download_rename("beast_example_phat_spec_w_priors.grid.hd5")
    # - anticipated output files
    seds_fname_cache = download_rename("beast_example_phat_seds.grid.hd5")
    # seds_fname_2 = download_rename("beast_example_phat_seds.grid.hd5")

    # run create_physicsmodel
    # (the import is here because it'll break without datamodel being created first)
    from beast.tools.run import create_physicsmodel
    create_physicsmodel.create_physicsmodel(nsubs=1, nprocs=1)

    # check that files match
    # - isochrones
    table_cache = Table.read(
        iso_fname_cache, format="ascii.csv", comment="#", delimiter=",",
    )
    table_new = Table.read(
        "./beast_no_subgrid/beast_no_subgrid_iso.csv",
        format="ascii.csv", comment="#", delimiter=",",
    )
    compare_tables(table_cache, table_new)
    # - spectra with priors
    compare_hdf5(priors_fname_cache, "./beast_no_subgrid/beast_no_subgrid_spec_w_priors.grid.hd5")
    # - SEDs grid
    compare_hdf5(seds_fname_cache, "./beast_no_subgrid/beast_no_subgrid_seds.grid.hd5")
