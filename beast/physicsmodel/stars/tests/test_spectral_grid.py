from astropy.tests.helper import remote_data
from astropy import units
from astropy import constants as const

from ...stars import stellib
from ...stars.isochrone import ezIsoch
from ...dust import extinction
from ... import grid
from ...model_grid import (make_spectral_grid,
                           add_stellar_priors)
from beast.tests.helpers import (download_rename, compare_hdf5)


@remote_data
def test_make_kurucz_tlusty_spectral_grid():

    # download the needed files
    kurucz_fname = download_rename('kurucz2004.grid.fits')
    tlusty_fname = download_rename('tlusty.lowres.grid.fits')
    filter_fname = download_rename('filters.hd5')
    iso_fname = download_rename('beast_example_phat_iso.csv')

    # download cached version of spectral grid
    spec_fname_cache = download_rename('beast_example_phat_spec_grid.hd5')

    ################
    # generate the same spectral grid from the code

    # read in the cached isochrones
    oiso = ezIsoch(iso_fname)

    # define the distance
    distances = [24.47]
    distance_unit = units.mag

    velocity = -300 * units.km / units.s
    redshift = (velocity / const.c).decompose().value

    # define the spectral libraries to use
    osl = stellib.Tlusty(filename=tlusty_fname) \
        + stellib.Kurucz(filename=kurucz_fname)

    # define the extinction curve to use
    extLaw = extinction.Gordon16_RvFALaw()


    filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
               'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
    add_spectral_properties_kwargs = dict(filternames=filters)

    spec_fname = '/tmp/beast_example_phat_spec_grid.hd5'
    spec_fname, g = make_spectral_grid(
        'test',
        oiso,
        osl=osl,
        redshift=redshift,
        distance=distances,
        distance_unit=distance_unit,
        spec_fname=spec_fname,
        filterLib=filter_fname,
        extLaw=extLaw,
        add_spectral_properties_kwargs=add_spectral_properties_kwargs)

    # compare the new to the cached version
    compare_hdf5(spec_fname_cache, spec_fname)


@remote_data
def test_add_stellar_priors_to_spectral_grid():

    # download the needed files
    gspec_fname = download_rename('beast_example_phat_spec_grid.hd5')

    # download cached version of spectral grid with priors
    priors_fname_cache = download_rename(
        'beast_example_phat_spec_w_priors.grid.hd5')

    ###############
    # generate the spectral grid with stellar priors from the code

    gspec_fname = '/tmp/beast_example_phat_spec_grid.hd5'
    specgrid = grid.FileSpectralGrid(gspec_fname, backend='memory')

    priors_fname = '/tmp/beast_example_phat_spec_w_priors.grid.hd5'
    priors_fname, g = add_stellar_priors('test', specgrid,
                                         priors_fname=priors_fname)

    # compare the new to the cached version
    compare_hdf5(priors_fname_cache, priors_fname)
