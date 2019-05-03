import numpy as np

from astropy.tests.helper import remote_data

from ...physicsmodel.grid import FileSEDGrid
from ...observationmodel.observations import Observations
from ...observationmodel.vega import Vega
from ...observationmodel.noisemodel import generic_noisemodel as noisemodel
from ..trim_grid import trim_models
from beast.tests.helpers import (download_rename, compare_hdf5)


class GenFluxCatalog(Observations):
    """Generic n band filter photometry
    This class implements a direct access to the Generic HST measured fluxes.

    ..note::
        it does not implement uncertainties as in this model, the noise is
        given through artificial star tests
    """
    def __init__(self, inputFile, filters, obs_colnames,
                 vega_fname=None):
        """ Construct the interface """
        desc = 'GENERIC star: %s' % inputFile
        Observations.__init__(self, inputFile, desc=desc)
        self.setFilters(filters, vega_fname=vega_fname)
        # some bad values smaller than expected
        # in physical flux units
        self.setBadValue(6e-40)

        # rate column needed as this is the *flux* column
        for ik, k in enumerate(filters):
            self.data.set_alias(k, obs_colnames[ik])

    def getFlux(self, num, units=False):
        """returns the absolute flux of an observation

        Parameters
        ----------
        num: int
            index of the star in the catalog to get measurement from

        units: bool
            if set returns the fluxes with a unit capsule

        Returns
        -------
        flux: ndarray[dtype=float, ndim=1]
            Measured integrated flux values throughout the filters
            in erg/s/cm^2/A
        """

        # case for using '_flux' result
        d = self.data[num]

        flux = np.array([d[self.data.resolve_alias(ok)]
                         for ok in self.filters])*self.vega_flux

        if units is True:
            return flux*units.erg/(units.s*units.cm*units.cm*units.angstrom)
        else:
            return flux

    def setFilters(self, filters, vega_fname=None):
        """ set the filters and update the vega reference for the conversions

        Parameters
        ----------
        filters: sequence
            list of filters using the internally normalized namings
        """
        self.filters = filters

        # Data "rates" are normalized to Vega already, fits are not using vega

        # for optimization purpose: pre-compute
        #  getting vega mags, require to open and read the content of one file.
        #  since getObs, calls getFlux, for each star you need to do this
        #  expensive operation
        with Vega(source=vega_fname) as v:
            _, vega_flux, _ = v.getFlux(filters)

        self.vega_flux = vega_flux


def get_obscat(obsfile, filters, obs_colnames, vega_fname=None,
               *args, **kwargs):
    """ Function that generates a data catalog object with the correct
    arguments

    Parameters
    ----------
    obsfile: str, optional (default datamodel.obsfile)
        observation file

    filters: sequence(str), optional, datamodel.filters
        seaquence of filters of the data

    returns
    -------
    obs: GenFluxCatalog
        observation catalog
    """
    obs = GenFluxCatalog(obsfile, filters, obs_colnames, vega_fname=vega_fname)
    return obs


@remote_data
# @pytest.mark.skip(reason="temporarily disable")
def test_trim_grid():

    # download the needed files
    vega_fname = download_rename('vega.hd5')
    seds_fname = download_rename('beast_example_phat_seds.grid.hd5')
    noise_fname = download_rename('beast_example_phat_noisemodel.hd5')
    obs_fname = download_rename('b15_4band_det_27_A.fits')

    # download cached version of noisemodel on the sed grid
    noise_trim_fname_cache = download_rename(
                                'beast_example_phat_noisemodel_trim.grid.hd5')
    seds_trim_fname_cache = download_rename(
                                'beast_example_phat_seds_trim.grid.hd5')

    ################

    # read in the observed data
    filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
               'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
    basefilters = ['F275W', 'F336W', 'F475W',
                   'F814W', 'F110W', 'F160W']
    obs_colnames = [f.lower() + '_rate' for f in basefilters]

    obsdata = get_obscat(obs_fname,
                         filters,
                         obs_colnames,
                         vega_fname=vega_fname)

    # get the modesedgrid
    modelsedgrid = FileSEDGrid(seds_fname)

    # read in the noise model just created
    noisemodel_vals = noisemodel.get_noisemodelcat(noise_fname)

    # trim the model sedgrid
    seds_trim_fname = 'beast_example_phat_sed_trim.grid.hd5'
    noise_trim_fname = seds_trim_fname.replace('_sed', '_noisemodel')

    trim_models(modelsedgrid, noisemodel_vals, obsdata,
                seds_trim_fname, noise_trim_fname, sigma_fac=3.)

    # compare the new to the cached version
    compare_hdf5(seds_trim_fname_cache, seds_trim_fname, ctype='seds')
    compare_hdf5(noise_trim_fname_cache, noise_trim_fname, ctype='noise')
