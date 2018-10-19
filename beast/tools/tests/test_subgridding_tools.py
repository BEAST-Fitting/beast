import os

import numpy as np
from astropy.tests.helper import remote_data
from astropy.table import Table
from astropy.io import fits
import tables

from beast.tools import subgridding_tools
from beast.tests.helpers import download_rename
from beast.physicsmodel.grid import FileSEDGrid
from beast.observationmodel.noisemodel.generic_noisemodel import get_noisemodelcat
from beast.fitting.tests.test_fit_grid import get_obscat
from beast.fitting import fit


def split_and_check(grid_fname, num_subgrids):
    complete_g = FileSEDGrid(grid_fname)
    sub_fnames = subgridding_tools.split_grid(grid_fname, num_subgrids)

    # count the number of grid cells
    sub_seds = []
    sub_grids = []

    for sub_fname in sub_fnames:
        sub_g = FileSEDGrid(sub_fname)

        sub_seds.append(sub_g.seds)
        sub_grids.append(sub_g.grid.data)

        np.testing.assert_equal(complete_g.lamb, sub_g.lamb)
        assert complete_g.grid.columns.items() == sub_g.grid.columns.items()

    sub_seds_reconstructed = np.concatenate(sub_seds)
    np.testing.assert_equal(sub_seds_reconstructed, complete_g.seds)

    sub_grids_reconstructed = np.concatenate(sub_grids)
    np.testing.assert_equal(sub_grids_reconstructed, complete_g.grid.data)

    # the split method skips anything that already exists, so if we
    # want to use this function multiple times for the same test
    # grid, we need to do this.
    for f in sub_fnames:
        os.remove(f)


@remote_data
def test_split_grid():
    seds_trim_fname = download_rename('beast_example_phat_seds_trim.grid.hd5')
    split_and_check(seds_trim_fname, 1)  # an edge case
    split_and_check(seds_trim_fname, 3)  # an odd numer
    split_and_check(seds_trim_fname, 4)  # an even number


@remote_data
def test_reduce_grid_info():
    seds_trim_fname = download_rename('beast_example_phat_seds_trim.grid.hd5')
    sub_fnames = subgridding_tools.split_grid(seds_trim_fname, 3)

    complete_g_info = subgridding_tools.subgrid_info(seds_trim_fname)
    cap_unique = 50
    sub_g_info = subgridding_tools.reduce_grid_info(
        sub_fnames, nprocs=3, cap_unique=cap_unique)

    for q in complete_g_info:
        assert q in sub_g_info
        assert complete_g_info[q]['min'] == sub_g_info[q]['min']
        assert complete_g_info[q]['max'] == sub_g_info[q]['max']
        num_unique = len(complete_g_info[q]['unique'])
        if num_unique > cap_unique:
            # Cpan still be larger if one of the sub results during the
            # reduction is larger. This is as intended.
            assert sub_g_info[q]['num_unique'] >= cap_unique
        else:
            assert sub_g_info[q]['num_unique'] == num_unique


@remote_data
def test_merge_pdf1d_stats():
    ########################################
    ## STEP 1: GET SOME DATA TO WORK WITH ##
    ########################################
    vega_fname = download_rename('vega.hd5')
    obs_fname = download_rename('b15_4band_det_27_A.fits')
    noise_trim_fname = download_rename(
        'beast_example_phat_noisemodel_trim.grid.hd5')
    seds_trim_fname = download_rename(
        'beast_example_phat_seds_trim.grid.hd5')

    # download cached version of fitting results
    # stats_fname_cache = download_rename('beast_example_phat_stats.fits')
    # pdf1d_fname_cache = download_rename('beast_example_phat_pdf1d.fits')

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

    ###########################################################################################
    ## STEP 2: SPLIT THE GRIDS AND GENERATE THE GRID INFO DICT AS IN THE SUBGRIDDING EXAMPLE ##
    ###########################################################################################
    num_subgrids = 3

    # Split SED grid
    sub_seds_trim_fnames = subgridding_tools.split_grid(
        seds_trim_fname, num_subgrids, overwrite=True)

    # Split noise grid (a standardized function does not exist)
    sub_noise_trim_fnames = []

    noisemodel_vals = get_noisemodelcat(noise_trim_fname)
    slices = subgridding_tools.uniform_slices(len(noisemodel_vals.root.bias),
                                              num_subgrids)
    for i, slc in enumerate(slices):
        outname = noise_trim_fname.replace('.hd5', 'sub{}.hd5'.format(i))
        with tables.open_file(outname, 'w') as outfile:
            outfile.create_array(outfile.root, 'bias',
                                 noisemodel_vals.root.bias[slc])
            outfile.create_array(outfile.root, 'error',
                                 noisemodel_vals.root.error[slc])
            outfile.create_array(outfile.root, 'completeness',
                                 noisemodel_vals.root.completeness[slc])
        sub_noise_trim_fnames.append(outname)

    # Collect information about the parameter rangers, to make the pdf1d bins
    # consistent between subgrids
    grid_info_dict = subgridding_tools.reduce_grid_info(sub_seds_trim_fnames,
                                                        sub_noise_trim_fnames,
                                                        nprocs=1,
                                                        cap_unique=100)

    ####################################################
    ## STEP 3: GENERATE FILENAMES AND RUN THE FITTING ##
    ####################################################
    def make_gridsub_fnames(base_fname, num_subgrids, extension='.fits'):
        return [base_fname.replace(extension, 'gridsub{}{}'.format(i, extension))
                for i in range(num_subgrids)]

    stats_fname = '/tmp/beast_example_phat_stats.fits'
    pdf1d_fname = '/tmp/beast_example_phat_pdf1d.fits'
    lnp_fname = '/tmp/beast_example_phat_lnp.hd5'

    subgrid_pdf1d_fnames = make_gridsub_fnames(pdf1d_fname, num_subgrids)
    subgrid_stats_fnames = make_gridsub_fnames(stats_fname, num_subgrids)
    subgrid_lnp_fnames = make_gridsub_fnames(lnp_fname, num_subgrids,
                                             extension='.hd5')

    for i in range(num_subgrids):
        sub_noisemodel_vals = get_noisemodelcat(sub_noise_trim_fnames[i])
        fit.summary_table_memory(obsdata, sub_noisemodel_vals,
                                 sub_seds_trim_fnames[i], threshold=-40.,
                                 save_every_npts=100, lnp_npts=60,
                                 stats_outname=subgrid_stats_fnames[i],
                                 pdf1d_outname=subgrid_pdf1d_fnames[i],
                                 lnp_outname=subgrid_lnp_fnames[i],
                                 grid_info_dict=grid_info_dict,
                                 do_not_normalize=True)
        # The do_not_normalize option is absolutely crucial!

    # Now merge the results
    merged_pdf1d_fname, merged_stats_fname = \
        subgridding_tools.merge_pdf1d_stats(subgrid_pdf1d_fnames,
                                            subgrid_stats_fnames)

    # Do a full fit also
    normal_stats = 'normal_stats.fits'
    normal_pdf1d = 'normal_pdf1d.fits'
    normal_lnp = 'normal_lnp.hd5'
    fit.summary_table_memory(obsdata, noisemodel_vals, seds_trim_fname,
                             threshold=-40., save_every_npts=100,
                             lnp_npts=60, stats_outname=normal_stats,
                             pdf1d_outname=normal_pdf1d,
                             lnp_outname=normal_lnp,
                             do_not_normalize=True)
    # Here, we also need to use do_not_normalize, otherwise Pmax will be
    # different by a factor

    # CHECKS
    tolerance = 1e-6
    print("comparing pdf1d")
    # fits_cache = fits.open(pdf1d_fname_cache)
    fits_normal = fits.open(normal_pdf1d)
    fits_new = fits.open(merged_pdf1d_fname)

    assert len(fits_new) == len(fits_normal)

    # A similar problem to the above will also occur here
    for k in range(1, len(fits_new)):
        qname = fits_new[k].header['EXTNAME']
        print(qname)
        np.testing.assert_allclose(fits_new[k].data, fits_normal[qname].data,
                                   rtol=tolerance, atol=tolerance)

    print("comparing stats")
    # table_cache = Table.read(stats_fname_cache)
    table_normal = Table.read(normal_stats)
    table_new = Table.read(merged_stats_fname)

    assert len(table_normal) == len(table_new)

    # These will normally fail, as the merging process can not be made
    # bit-correct due do floating point math (exacerbated by exponentials)
    for c in table_new.colnames:
        print(c)
        if c == 'Name' or c == 'RA' or c == 'DEC':
            np.testing.assert_equal(table_normal[c], table_new[c],
                                    err_msg='column {} is not equal'.format(c))
        else:
            np.testing.assert_allclose(table_normal[c], table_new[c], rtol=tolerance,
                                       equal_nan=True, err_msg='column {} is not close enough'.format(c))
