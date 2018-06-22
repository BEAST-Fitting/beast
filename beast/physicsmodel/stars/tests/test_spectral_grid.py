import os.path

import numpy as np
import h5py

from astropy.utils.data import download_file
from astropy.tests.helper import remote_data
from astropy import units
from astropy import constants as const

from ...stars import stellib
from ...stars.isochrone import ezIsoch
from ... import grid
from ...model_grid import (make_spectral_grid,
                           add_stellar_priors)


@remote_data
# @pytest.mark.skip(reason="temporarily disable")
def test_make_kurucz_tlusty_spectral_grid():

    # download the needed files
    url_loc = 'http://www.stsci.edu/~kgordon/beast/'
    kurucz_fname_dld = download_file('%s%s'
                                     % (url_loc, 'kurucz2004.grid.fits'))
    tlusty_fname_dld = download_file('%s%s'
                                     % (url_loc, 'tlusty.lowres.grid.fits'))
    filter_fname_dld = download_file('%s%s'
                                     % (url_loc, 'filters.hd5'))
    iso_fname_dld = download_file('%s%s'
                                  % (url_loc, 'beast_example_phat_iso.csv'))

    # rename files to have the correct extensions
    kurucz_fname = '%s.fits' % (kurucz_fname_dld)
    os.rename(kurucz_fname_dld, kurucz_fname)
    tlusty_fname = '%s.fits' % (tlusty_fname_dld)
    os.rename(tlusty_fname_dld, tlusty_fname)
    filter_fname = '%s.hd5' % (filter_fname_dld)
    os.rename(filter_fname_dld, filter_fname)
    iso_fname = '%s.csv' % (iso_fname_dld)
    os.rename(iso_fname_dld, iso_fname)

    # download cached version of spectral grid
    filename = download_file('%s%s' % (url_loc,
                                     'beast_example_phat_spec_grid.hd5'))

    hdf_cache = h5py.File(filename, 'r')

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

    filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
               'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
    add_spectral_properties_kwargs = dict(filternames=filters)

    spec_fname = '/tmp/beast_example_phat_spec_grid.hd5'
    spec_fname, g = make_spectral_grid('test',
                                       oiso,
                                       osl=osl,
                                       redshift=redshift,
                                       distance=distances,
                                       distance_unit=distance_unit,
                                       spec_fname=spec_fname,
                                       filterLib=filter_fname,
         add_spectral_properties_kwargs=add_spectral_properties_kwargs)

    # open the hdf file with the specral grid
    hdf_new = h5py.File(spec_fname, 'r')

    # go through the file and check if it is exactly the same
    for sname in hdf_cache.keys():
        if isinstance(hdf_cache[sname], h5py.Dataset):
            cvalue = hdf_cache[sname]
            cvalue_new = hdf_new[sname]
            if cvalue.dtype.fields is None:
                np.testing.assert_equal(cvalue.value, cvalue_new.value,
                                        'testing %s' % (sname))
            else:
                for ckey in cvalue.dtype.fields.keys():
                    np.testing.assert_equal(cvalue.value[ckey],
                                            cvalue_new.value[ckey],
                                            'testing %s/%s' % (sname, ckey))


@remote_data
# @pytest.mark.skip(reason="temporarily disable")
def test_add_stellar_priors_to_spectral_grid():

    # download the needed files
    url_loc = 'http://www.stsci.edu/~kgordon/beast/'
    gspec_fname_dld = download_file('%s%s' % (url_loc,
                                        'beast_example_phat_spec_grid.hd5'))
    # rename files to have the correct extensions
    gspec_fname = '%s.hd5' % (gspec_fname_dld)
    os.rename(gspec_fname_dld, gspec_fname)

    filename = download_file('%s%s' % (url_loc,
                                'beast_example_phat_spec_w_priors.grid.hd5'))
    hdf_cache = h5py.File(filename, 'r')

    ###############
    # generate the spectral grid with stellar priors from the code

    gspec_fname = '/tmp/beast_example_phat_spec_grid.hd5'
    specgrid = grid.FileSpectralGrid(gspec_fname, backend='memory')

    priors_fname = '/tmp/beast_example_phat_spec_w_priors.grid.hd5'
    priors_fname, g = add_stellar_priors('test', specgrid,
                                         priors_fname=priors_fname)

    # open the hdf file with the specral grid with priors
    hdf_new = h5py.File(priors_fname, 'r')

    # go through the file and check if it is exactly the same
    for sname in hdf_cache.keys():
        if isinstance(hdf_cache[sname], h5py.Dataset):
            cvalue = hdf_cache[sname]
            cvalue_new = hdf_new[sname]
            if cvalue.dtype.fields is None:
                np.testing.assert_equal(cvalue.value, cvalue_new.value,
                                        'testing %s' % (sname))
            else:
                for ckey in cvalue.dtype.fields.keys():
                    np.testing.assert_equal(cvalue.value[ckey],
                                            cvalue_new.value[ckey],
                                            'testing %s/%s' % (sname, ckey))


if __name__ == '__main__':

    test_add_stellar_priors_to_spectral_grid()
