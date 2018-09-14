from astropy.tests.helper import remote_data
from astropy.table import Table

from ..noisemodel import generic_noisemodel as noisemodel
from ...physicsmodel.grid import FileSEDGrid
from beast.observationmodel.observations import gen_SimObs_from_sedgrid
from beast.tests.helpers import (download_rename, compare_tables)


@remote_data
def test_trim_grid():

    # download the needed files
    vega_fname = download_rename('vega.hd5')
    seds_fname = download_rename('beast_example_phat_seds.grid.hd5')
    noise_fname = download_rename('beast_example_phat_noisemodel.hd5')

    # download cached version of noisemodel on the sed grid
    simobs_fname_cache = download_rename('beast_example_phat_simobs.fits')

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

    compare_tables(table_cache, table_new)
