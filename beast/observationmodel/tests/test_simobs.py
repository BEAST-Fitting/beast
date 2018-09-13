import os.path

import numpy as np

from astropy.utils.data import download_file
from astropy.tests.helper import remote_data
from astropy.table import Table

from ..noisemodel import generic_noisemodel as noisemodel
from ...physicsmodel.grid import FileSEDGrid
from beast.observationmodel.observations import gen_SimObs_from_sedgrid


def _download_rename(filename):
    """
    Download a file and rename it to have the right extension

    Otherwise, downloaded file will not have an extension at all
    """
    url_loc = 'http://www.stsci.edu/~kgordon/beast/'
    fname_dld = download_file('%s%s' % (url_loc, filename))
    extension = filename.split('.')[-1]
    fname = '%s.%s' % (fname_dld, extension)
    os.rename(fname_dld, fname)
    return fname


@remote_data
def test_trim_grid():

    # download the needed files
    vega_fname = _download_rename('vega.hd5')
    seds_fname = _download_rename('beast_example_phat_seds.grid.hd5')
    noise_fname = _download_rename('beast_example_phat_noisemodel.hd5')

    # download cached version of noisemodel on the sed grid
    simobs_fname_cache = _download_rename('beast_example_phat_simobs.fits')

    ################

    # get the physics model grid - includes priors
    modelsedgrid = FileSEDGrid(seds_fname)

    # read in the noise model - includes bias, unc, and completeness
    noisegrid = noisemodel.get_noisemodelcat(noise_fname)

    table_new = gen_SimObs_from_sedgrid(modelsedgrid, noisegrid,
                                        nsim=100,
                                        compl_filter="f475w",
                                        ranseed=1234,
                                        vega_fname=vega_fname)

    # check that the simobs files are exactly the same
    table_cache = Table.read(simobs_fname_cache)

    assert len(table_new) == len(table_cache)

    for tcolname in table_new.colnames:
        np.testing.assert_equal(table_new[tcolname], table_cache[tcolname],
                                '%s columns not equal' % tcolname)
