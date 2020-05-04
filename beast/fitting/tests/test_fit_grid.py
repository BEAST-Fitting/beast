from astropy.table import Table
from astropy.tests.helper import remote_data

import beast.observationmodel.noisemodel.generic_noisemodel as noisemodel
from beast.fitting import fit
from beast.tests.helpers import download_rename, compare_tables, compare_fits
from beast.fitting.tests.helpers import get_obscat


@remote_data
def test_fit_grid():

    # download the needed files
    vega_fname = download_rename("vega.hd5")
    obs_fname = download_rename("b15_4band_det_27_A.fits")
    noise_trim_fname = download_rename("beast_example_phat_noisemodel_trim.grid.hd5")
    seds_trim_fname = download_rename("beast_example_phat_seds_trim.grid.hd5")

    # download cached version of fitting results
    stats_fname_cache = download_rename("beast_example_phat_stats.fits")
    pdf1d_fname_cache = download_rename("beast_example_phat_pdf1d.fits")

    ################

    # read in the the AST noise model
    noisemodel_vals = noisemodel.get_noisemodelcat(noise_trim_fname)

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
    # output files
    stats_fname = "/tmp/beast_example_phat_stats.fits"
    pdf1d_fname = "/tmp/beast_example_phat_pdf1d.fits"
    lnp_fname = "/tmp/beast_example_phat_lnp.hd5"

    fit.summary_table_memory(
        obsdata,
        noisemodel_vals,
        seds_trim_fname,
        threshold=-10.0,
        save_every_npts=100,
        lnp_npts=60,
        stats_outname=stats_fname,
        pdf1d_outname=pdf1d_fname,
        lnp_outname=lnp_fname,
    )

    # check that the stats files are exactly the same
    table_cache = Table.read(stats_fname_cache)
    table_new = Table.read(stats_fname)

    compare_tables(table_cache, table_new)

    # lnp files not checked as they are randomly sparsely sampled
    #   hence will be different every time the fitting is run

    # check that the pdf1d files are exactly the same
    compare_fits(pdf1d_fname_cache, pdf1d_fname)
