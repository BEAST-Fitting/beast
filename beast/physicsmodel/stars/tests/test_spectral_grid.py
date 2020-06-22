from astropy.tests.helper import remote_data
from astropy import units
from astropy import constants as const

from beast.physicsmodel.stars import stellib
from beast.physicsmodel.stars.isochrone import ezIsoch
from beast.physicsmodel.dust import extinction
from beast.physicsmodel.grid import SpectralGrid
from beast.physicsmodel.model_grid import make_spectral_grid, add_stellar_priors
from beast.tests.helpers import download_rename, compare_hdf5
from beast.tools import get_libfiles


@remote_data
class TestBeast:

    # download the BEAST library files
    get_libfiles.get_libfiles()

    # download the cached version for use and comparision
    iso_fname_cache = download_rename("beast_example_phat_iso.csv")
    spec_fname_cache = download_rename("beast_example_phat_spec_grid.hd5")
    priors_fname_cache = download_rename("beast_example_phat_spec_w_priors.grid.hd5")

    def test_make_kurucz_tlusty_spectral_grid(self):

        # read in the cached isochrones
        oiso = ezIsoch(self.iso_fname_cache)

        # define the distance
        distances = [24.47]
        distance_unit = units.mag

        velocity = -300 * units.km / units.s
        redshift = (velocity / const.c).decompose().value

        # define the spectral libraries to use
        osl = stellib.Tlusty() + stellib.Kurucz()

        # define the extinction curve to use
        extLaw = extinction.Gordon16_RvFALaw()

        filters = [
            "HST_WFC3_F275W",
            "HST_WFC3_F336W",
            "HST_ACS_WFC_F475W",
            "HST_ACS_WFC_F814W",
            "HST_WFC3_F110W",
            "HST_WFC3_F160W",
        ]
        add_spectral_properties_kwargs = dict(filternames=filters)

        spec_fname = "/tmp/beast_example_phat_spec_grid.hd5"
        spec_fname, g = make_spectral_grid(
            "test",
            oiso,
            osl=osl,
            redshift=redshift,
            distance=distances,
            distance_unit=distance_unit,
            spec_fname=spec_fname,
            # filterLib=filter_fname,
            extLaw=extLaw,
            add_spectral_properties_kwargs=add_spectral_properties_kwargs,
        )

        # compare the new to the cached version
        compare_hdf5(self.spec_fname_cache, spec_fname)

    def test_add_stellar_priors_to_spectral_grid(self):

        # gspec_fname = "/tmp/beast_example_phat_spec_grid.hd5"
        specgrid = SpectralGrid(self.spec_fname_cache, backend="memory")

        priors_fname = "/tmp/beast_example_phat_spec_w_priors.grid.hd5"
        priors_fname, g = add_stellar_priors("test", specgrid, priors_fname=priors_fname)

        # compare the new to the cached version
        compare_hdf5(self.priors_fname_cache, priors_fname)
