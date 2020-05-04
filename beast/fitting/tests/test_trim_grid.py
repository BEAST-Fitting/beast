from astropy.tests.helper import remote_data

from beast.physicsmodel.grid import FileSEDGrid
from beast.observationmodel.noisemodel import generic_noisemodel as noisemodel
from beast.fitting.trim_grid import trim_models
from beast.tests.helpers import download_rename, compare_hdf5
from beast.fitting.tests.helpers import get_obscat


@remote_data
# @pytest.mark.skip(reason="temporarily disable")
def test_trim_grid():

    # download the needed files
    vega_fname = download_rename("vega.hd5")
    seds_fname = download_rename("beast_example_phat_seds.grid.hd5")
    noise_fname = download_rename("beast_example_phat_noisemodel.grid.hd5")
    obs_fname = download_rename("b15_4band_det_27_A.fits")

    # download cached version of noisemodel on the sed grid
    noise_trim_fname_cache = download_rename(
        "beast_example_phat_noisemodel_trim.grid.hd5"
    )
    seds_trim_fname_cache = download_rename("beast_example_phat_seds_trim.grid.hd5")

    ################

    # read in the observed data
    filters = [
        "HST_WFC3_F275W",
        "HST_WFC3_F336W",
        "HST_ACS_WFC_F475W",
        "HST_ACS_WFC_F814W",
        "HST_WFC3_F110W",
        "HST_WFC3_F160W",
    ]
    basefilters = ["F275W", "F336W", "F475W", "F814W", "F110W", "F160W"]
    obs_colnames = [f.lower() + "_rate" for f in basefilters]

    obsdata = get_obscat(obs_fname, filters, obs_colnames, vega_fname=vega_fname)

    # get the modesedgrid
    modelsedgrid = FileSEDGrid(seds_fname)

    # read in the noise model just created
    noisemodel_vals = noisemodel.get_noisemodelcat(noise_fname)

    # trim the model sedgrid
    seds_trim_fname = "beast_example_phat_seds_trim.grid.hd5"
    noise_trim_fname = seds_trim_fname.replace("_seds", "_noisemodel")

    trim_models(
        modelsedgrid,
        noisemodel_vals,
        obsdata,
        seds_trim_fname,
        noise_trim_fname,
        sigma_fac=3.0,
    )

    # compare the new to the cached version
    compare_hdf5(seds_trim_fname_cache, seds_trim_fname, ctype="seds")
    compare_hdf5(noise_trim_fname_cache, noise_trim_fname, ctype="noise")
