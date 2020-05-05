from beast.tools.run import create_physicsmodel
from beast.tests.helpers import download_rename
from astropy.tests.helper import remote_data
import os

@remote_data
def test_create_physicsmodel_no_subgrid():
    """
    Test create_physicsmodel, assuming no subgrids
    """

    # download files
    # - datamodel settings
    datamodel_fname = download_rename("datamodel_no_subgrid.py")
    # - input files (referenced in datamodel)
    filter_fname = download_rename("filters.hd5")
    # - intermediate files (if they exist, create_physicsmodel won't
    #   regenerate them, which is good -- that part doesn't need testing, and
    #   there's a time limit for running things on Travis)
    #priors_fname = download_rename("beast_example_phat_spec_w_priors.grid.hd5")
    # - anticipated output files
    seds_fname_1 = download_rename("beast_example_phat_seds.grid.hd5")
    #seds_fname_2 = download_rename("beast_example_phat_seds.grid.hd5")

    # rename datamodel so it can be imported
    os.rename(datamodel_fname, 'datamodel.py')
