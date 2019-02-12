#!/usr/bin/env python3
"""Create a background map or a source density map, depending on the
subcommand given. The background map is based on an input catalog and a
fits image on which the backgrounds are measured for each source using
extended source photometry. Several files are are saved to disk, and the
resulting map in *hd5 format can later be loaded and reused using the
tools.DensityMap class in other parts of the BEAST.

"""

import argparse
import astropy
from astropy import wcs
from astropy.table import Table
from astropy import units
from astropy.io import fits
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import numpy as np
import photutils as pu
from density_map import DensityMap
import itertools as it
import os


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='map_type')
    subparsers.required = True

    # Common options for both types of map
    commonparser = argparse.ArgumentParser(add_help=False)
    commonparser.add_argument('-catfile', type=str, required=True,
                              help='catalog FITS file')

    npix_or_pixsize = commonparser.add_mutually_exclusive_group()
    npix_or_pixsize.add_argument('--npix', type=int, default=None,
                                 help='resolution')
    npix_or_pixsize.add_argument('--pixsize', type=float, default=None)

    background_parser = subparsers.add_parser('background', parents=[commonparser],
        help="""Create a background intensity map based on annulus
        measurements around the sources listed in the catalog""")
    sourceden_parser = subparsers.add_parser('sourceden', parents=[commonparser],
        help="""Create a source density map, counting the number of
        sources in each tile of the map""")

    background_parser.add_argument('-reference', type=str,  metavar='FITSIMAGE', required=True,
                                   help='reference image (FITS)')
    background_parser.add_argument('--plot_map_on_image', action='store_true',
                                   help='plot the map using transparant colored tiles on top of the reference image')
    background_parser.add_argument('--colorbar', type=str, default=None, metavar='LABEL',
                                   help='put a colorbar next to the plot and label it with the given string')

    args = parser.parse_args()

    # Common actions: load the catalog and set up the grid
    cat = Table.read(args.catfile)
    for name in cat.colnames:
        cat.rename_column(name, name.upper())

    ra = cat['RA']
    dec = cat['DEC']

    if args.npix is not None:
        n_x, n_y = args.npix, args.npix
        ra_grid = np.linspace(ra.min(), ra.max(), n_x + 1)
        dec_grid = np.linspace(dec.min(), dec.max(), n_y + 1)
    elif args.pixsize is not None:
        pixsize_degrees = args.pixsize / 3600
        n_x, n_y, ra_delt, dec_delt = calc_nx_ny_from_pixsize(cat, pixsize_degrees)
        # the ra spacing needs to be larger, as 1 degree of RA ==
        # cos(DEC) degrees on the great circle
        ra_grid = ra.min() + ra_delt * np.arange(0, n_x + 1, dtype=float)
        dec_grid = dec.min() + dec_delt * np.arange(0, n_y + 1, dtype=float)
    else:
        n_x, n_y = 10, 10
        ra_grid = np.linspace(ra.min(), ra.max(), n_x + 1)
        dec_grid = np.linspace(dec.min(), dec.max(), n_y + 1)

    output_base=os.path.basename(args.catfile).replace('.fits', '')

    if args.map_type == 'sourceden':
        map_values_array = make_source_dens_map(cat, ra_grid, dec_grid,
                                                output_base)

    if args.map_type == 'background':
        hdul = astropy.io.fits.open(args.reference)
        image = hdul[1]
        result = make_background_map(cat, ra_grid, dec_grid, ref_im=image,
                                     output_base=output_base)
        map_values_array = result['background_map']

        # Parts of the plotting which depend on command line arguments.
        # Maybe make this kind of plot possible for the sourceden map
        # too.
        n_map = result['nsources_map']
        ra_grid = result['ra_grid']
        dec_grid = result['dec_grid']
        mask = result['mask']

        # Plot the maps directly
        _, ax = plt.subplots(1, 2)
        figs = []
        figs.append(ax[0].imshow(map_values_array))
        ax[0].set_title('bg_density_estimate')
        figs.append(ax[1].imshow(n_map))
        ax[1].set_title('number of sources')
        for f, a in zip(figs, ax):
            plt.colorbar(f, ax=a)
            plt.tight_layout()

        plt.savefig('{}_plot_bg_and_npts.png'.format(output_base))

        # Overplot the map on the image used for the calculation
        if args.plot_map_on_image:
            image_fig, image_ax, patch_col = plot_on_image(
                image, map_values_array, ra_grid, dec_grid, mask=mask)

            if args.colorbar is not None:
                cb = image_fig.colorbar(patch_col)
                cb.set_alpha(1)
                cb.draw_all()

                if len(args.colorbar) == 0:
                    label = 'background (arbitrary units)'
                else:
                    label = args.colorbar

                cb.set_label(label)

            image_fig.savefig('{}_plot_overlay.pdf'.format(output_base))

    # Save a file describing the properties of the bins in a handy format
    bin_details = astropy.table.Table(
        names=['i_ra', 'i_dec', 'value',
               'min_ra', 'max_ra',
               'min_dec', 'max_dec'])
    for x, y in xyrange(n_x, n_y):
        bin_details.add_row([x, y, map_values_array[x, y],
                             ra_grid[x], ra_grid[x + 1],
                             dec_grid[y], dec_grid[y + 1]])

    # Add the ra and dec grids as metadata
    bin_details.meta['ra_grid'] = ra_grid
    bin_details.meta['dec_grid'] = dec_grid

    # This works for both the background density map or the source
    # density map
    dm = DensityMap(bin_details)
    dm.write('{}_{}_map.hd5'.format(output_base, args.map_type))


def make_background_map(cat, ra_grid, dec_grid, ref_im, output_base):
    """
    Divide the image into a number of bins, and calculate the median
    background for the stars that fall within each bin. Create a new
    catalog on disk, which contains the individually measured and binned
    background estimate for each source. Create a new fits table on
    disk, containing the x,y position, bin edges (ra and dec in
    degrees), and the median background density in this bin.

    Parameters
    ----------
    cat: astropy Table
        The photometry catalog. The positions of the sources will be
        used to measure the backgrounds, and mask the sources
        themselves.

    ra_grid: 1D array-like of float
        the edges of the RA bins

    dec_grid: 1D array-like of float
        the edges of the DEC bins

    ref_im: imageHDU
        image which will be used for the background measurements

    output_base: string
        base name (without extension) to be used for the output files

    Returns
    -------
    results: dict
        this dict contains:
        'background_map': 2d ndarray,
        'nsources_map': 2d ndarray,
        'ra_grid': list of ra bin edges,
        'dec_grid': list of dec bin edges,
        'mask': 2d ndarray the pixels in ref_im that were ignored
    """
    # A list of background values for each source of the catalog will be
    # built up. Mask used is also returned.
    individual_backgrounds, mask = measure_backgrounds(cat, ref_im)

    w = make_wcs_for_map(ra_grid, dec_grid)
    pix_x, pix_y = get_pix_coords(cat, w)

    n_x = len(ra_grid)
    n_y = len(dec_grid)

    background_map = np.zeros((n_x, n_y))
    nsources_map = np.zeros((n_x, n_y))
    median_backgrounds = np.zeros((len(cat),))
    for x, y in xyrange(n_x, n_y):
        idxs = indices_for_pixel(pix_x, pix_y, x, y)
        n = len(idxs)
        nsources_map[x, y] = n
        if n:
            background_map[x, y] = np.median(individual_backgrounds[idxs])
        if n == 1:
            print('Only 1 source in bin {},{}'.format(x, y))

        # store the median background for each source
        median_backgrounds[idxs] = background_map[x, y]

    background_map[nsources_map == 0] = 0

    # Save the catalog with extra density info
    extra_columns = {'indiv_bg': individual_backgrounds,
                     'bin_median_bg': median_backgrounds}
    for k in extra_columns:
        c = astropy.table.Column(extra_columns[k], name=k)
        cat.add_column(c)
    cat.write(output_base + '_with_bg.fits', format='fits', overwrite=True)

    # Save the map as a fits file here. The map in hd5 format is written
    # by some common code, in the main() function .
    save_map_fits(background_map, w, output_base + '_background.fits')
    save_map_fits(nsources_map, w, output_base + '_nsources.fits')
    # Return a bunch of stuff, to be used for plots
    return {'background_map': background_map,
            'nsources_map': nsources_map,
            'ra_grid': ra_grid,
            'dec_grid': dec_grid,
            'mask': mask}


def measure_backgrounds(cat_table, ref_im):
    """
    Measure the background for all the sources in cat_table, using
    ref_im.

    Parameters
    ----------
    cat_table: astropy Table
        the catalog in astropy table form, as loaded from a .gst.fits
        photometry catalog

    ref_im: imageHDU
        fits image which will be used to estimate the background

    Returns
    -------

    measurements: list of float
        a metric for the background intensity at the star's position
        (photometry / area)

    mask: 2d ndarray of bool
        an array containing a bool for each pixel of the image, True if
        the pixel was ignored for the background calculations
    """
    w = wcs.WCS(ref_im.header)
    shp = ref_im.data.shape

    inner_rad = 30 * units.pixel
    outer_rad = inner_rad + 20 * units.pixel
    mask_rad = inner_rad

    # More elaborate way using an actual image (do not care about the
    # filter for the moment, just take the image (which was given on the
    # command line) at face value)
    ra = cat_table['RA']
    dec = cat_table['DEC']
    c = astropy.coordinates.SkyCoord(ra * units.degree, dec * units.degree)

    # Annuli, of which the counts per surface area will be used as
    # background measurements
    annuli = pu.SkyCircularAnnulus(c, r_in=inner_rad, r_out=outer_rad)
    area = annuli.to_pixel(w).area()

    # A mask to make sure that no sources end up in the background
    # calculation
    circles = pu.SkyCircularAperture(c, mask_rad)
    source_masks = circles.to_pixel(w).to_mask()
    mask_union = np.zeros(shp)
    for i, ap_mask in enumerate(source_masks):
        # the masks have a bounding box which we need to take into
        # account. Here we will calculate the overlap between the mask
        # and the image (important if the mask is near the edge, so that
        # the box might go outside of it).
        data_slices = list(ap_mask.bbox.slices)
        # These slices go outside of the box sometimes! In that case we
        # will override the slices on both sides.

        # Default slices (go over the whole mask)
        mask_slices = [slice(None, None), slice(None, None)]

        # Adjust the slices in each dimension
        for j in range(2):
            # DATA: . . a b c d e f
            # index - - 0 1 2 3 4 5
            # ---------------
            # MASK: - + + + -
            # index 0 1 2 3 4
            # --> DATA_SLICE [0:stop]
            # --> MASK_SLICE [2:]
            data_start = data_slices[j].start
            if data_start < 0:
                # Move the start from negative n to zero
                data_slices[j] = slice(0, data_slices[j].stop)
                # Move the start from 0 to positive n
                mask_slices[j] = slice(-data_start, mask_slices[j].stop)
                # --> we slice over the part of the small mask that
                # falls on the positive side of the axis

            data_stop = data_slices[j].stop
            if data_stop > shp[j]:
                overflow = data_stop - shp[j]
                # Move the stop 'overflow' to shp[j]
                data_slices[j] = slice(data_slices[j].start, shp[j])
                # Move the stop to 'overflow' pixels from the end
                mask_slices[j] = slice(mask_slices[j].start, -overflow)
                # --> slice over the part that falls below the maximum

        mask_union[tuple(data_slices)] += ap_mask.data[tuple(mask_slices)]

    # Threshold
    mask_union = mask_union > 0
    phot = pu.aperture_photometry(ref_im.data, annuli, wcs=w, mask=mask_union)
    return phot['aperture_sum'] / area, mask_union


def make_source_dens_map(cat,
                         ra_grid, dec_grid,
                         output_base,
                         mag_name='F475W_VEGA',
                         mag_cut=[24.5, 27]):
    """
    Computes the source density map and store it in a pyfits HDU
    Also writes a text file storing the source density for each source

    INPUTS:
    -------
    cat: astropy Table
        the photometry catalog

    ra_grid, dec_grid: 1D array-like of float
        the edges of the bins in RA and DEC space

    output_base: string

    mag_name: string
         name of magnitude column in table

    mag_cut: 2-element list
         magnitude range on which the source density is computed

    OUTPUT:
    -------
    FITS files written to disk.
    """
    # force filter magnitude name to be upper case to match column names
    mag_name = mag_name.upper()

    # get the columns with fluxes
    rate_cols = [s for s in cat.colnames if s[-4:] == 'RATE']
    n_filters = len(rate_cols)

    # create the indexs where any of the rates are zero and non-zero
    #   zero = missing data, etc. -> bad for fitting
    #   non-zero = good data, etc. -> great for fitting
    initialize_zero = False
    band_zero_indxs = {}
    print('band, good, zero')
    for cur_rate in rate_cols:
        cur_good_indxs, = np.where(cat[cur_rate] != 0.0)
        cur_indxs, = np.where(cat[cur_rate] == 0.0)
        print(cur_rate, len(cur_good_indxs), len(cur_indxs))
        if not initialize_zero:
            initialize_zero = True
            zero_indxs = cur_indxs
            nonzero_indxs = cur_good_indxs
        else:
            zero_indxs = np.union1d(zero_indxs, cur_indxs)
            nonzero_indxs = np.intersect1d(nonzero_indxs, cur_good_indxs)

        # save the zero indexs for each band
        band_zero_indxs[cur_rate] = zero_indxs

    print('all bands', len(nonzero_indxs), len(zero_indxs))

    N_stars = len(cat)

    w = make_wcs_for_map(ra_grid, dec_grid)
    pix_x, pix_y = get_pix_coords(cat, w)

    n_x = len(ra_grid) - 1
    n_y = len(dec_grid) - 1
    npts_map = np.zeros([n_x, n_y], dtype=float)
    npts_zero_map = np.zeros([n_x, n_y], dtype=float)
    npts_band_zero_map = np.zeros([n_x, n_y, n_filters], dtype=float)
    source_dens = np.zeros(N_stars, dtype=float)

    for i, j in xyrange(n_x, n_y):
        indxs = indices_for_pixel(pix_x, pix_y, i, j)
        indxs_for_SD, = np.where(np.logical_and(cat[mag_name][indxs] >= mag_cut[0],
                                                   cat[mag_name][indxs] <= mag_cut[1]))
        n_indxs = len(indxs_for_SD)
        if n_indxs > 0:
            pix_area = (ra_grid[i + 1] - ra_grid[i]) * (dec_grid[j + 1] - dec_grid[j]) * 3600 ** 2
            npts_map[i, j] = n_indxs / pix_area

            # now make a map of the sources with zero fluxes in
            #   at least one band
            zindxs, = np.where((pix_x[zero_indxs] > i)
                               & (pix_x[zero_indxs] <= i+1)
                               & (pix_y[zero_indxs] > j)
                               & (pix_y[zero_indxs] <= j+1))
            if len(zindxs) > 0:
                npts_zero_map[i, j] = len(zindxs)

            # do the same for each band
            for k, cur_rate in enumerate(rate_cols):
                tindxs = band_zero_indxs[cur_rate]
                zindxs, = np.where((pix_x[tindxs] > i)
                                   & (pix_x[tindxs] <= i+1)
                                   & (pix_y[tindxs] > j)
                                   & (pix_y[tindxs] <= j+1))
                if len(zindxs) > 0:
                    npts_band_zero_map[i, j, k] = len(zindxs)

        # save the source density as an entry for each source
        source_dens[indxs_for_SD] = npts_map[i, j]

    save_map_fits(npts_map, w, output_base + '_source_den_image.fits')
    save_map_fits(npts_zero_map, w, output_base + '_npts_zero_fluxes_image.fits')
    for k, cur_rate in enumerate(rate_cols):
        save_map_fits(npts_band_zero_map[:, :, k], w, output_base + cur_rate + '_image.fits')

    # Save the source density for individual stars in a new catalog file
    cat['SourceDensity'] = source_dens
    cat.write(output_base + '_with_sourceden_inc_zerofluxes.fits',
              overwrite=True)

    # Save the source density for individual stars in a new catalog file
    #   only those that have non-zero fluxes in all bands
    cat[nonzero_indxs].write(output_base + '_with_sourceden.fits',
                             overwrite=True)

    return npts_map


def plot_on_image(image, background_map, ra_grid, dec_grid, mask=None):
    """
    Plot the density grid as a collection of colored rectangles layered
    over the given fits image.

    Parameters
    ----------
    image: imageHDU
        the fits image, which should include a WCS

    background_map: 2d ndarray
        the values for the bins (indexed on [ra_index,dec_index])

    ra_grid: list of float
        edges of the right ascension bins in degrees

    dec_grid: list of float
        edges of the declination bins in degrees

    mask: 2d ndarray of bool
        Mask that blacks out pixels. Needs to be of same dimensions as
        image.data (indexed on y,x).

    Returns
    -------
    image_fig: matplotlib Figure
        the figure on which the plot was made

    image_ax: matplotlib Axes
        the axes on which the imshow and the patches were applied

    patch_col: matplotlib PatchCollection
        the patch collection that was added to the axes to plot the
        rectangles representing the background map tiles
    """
    # plot the image
    image_wcs = wcs.WCS(image.header)
    image_fig = plt.figure()
    image_ax = image_fig.add_subplot(1, 1, 1, projection=image_wcs)
    # image_ax = plt.subplot(1, 1, 1, projection=wcs.WCS(image))
    imdata = image.data.astype(float)
    vmin = np.percentile(imdata, 16)
    vmax = np.percentile(imdata, 99)
    if mask is not None:
        imdata = np.where(mask, vmin, imdata)
    plt.imshow(imdata, cmap='gray_r', interpolation='none', vmin=vmin,
               vmax=vmax)

    # Make a rectangular patch for each grid point
    n_x = len(ra_grid) - 1
    n_y = len(dec_grid) - 1
    rectangles = []
    values = []
    for ix, iy in xyrange(n_x, n_y):
        # bottom left, bottom right, upper left and upper right in the RA, DEC grid
        l = ix
        r = ix + 1
        b = iy
        u = iy + 1
        ra_corners = [ra_grid[i] for i in [l, r, r, l]]
        dec_corners = [dec_grid[i] for i in [b, b, u, u]]
        ra_dec_corners = np.column_stack((ra_corners, dec_corners))
        pix_corners = image_wcs.wcs_world2pix(ra_dec_corners, 0)
        rec = Polygon(pix_corners, closed=True)
        rectangles.append(rec)

        values.append(background_map[ix, iy])

    patch_col = PatchCollection(rectangles, cmap='viridis')
    patch_col.set_alpha(0.3)

    # Associate the values with the patches. They will be used to
    # pick colors from the colorbar. By passing the patch collection
    # as an argument, the patch collection will be treated as a
    # 'mappable', which works because we did set_array.
    image_ax.add_collection(patch_col)
    patch_col.set_array(np.array(values))
    return image_fig, image_ax, patch_col


def calc_nx_ny_from_pixsize(cat, pixsize_degrees):
    min_ra = cat['RA'].min()
    max_ra = cat['RA'].max()
    min_dec = cat['DEC'].min()
    max_dec = cat['DEC'].max()

    # Compute the required width of the bins expressed in RA and DEC to
    # reach the requested physical pixel size on the sky
    dec_delt = pixsize_degrees
    cos_avg_dec = math.cos((max_dec + min_dec) / 2 * math.pi / 180)
    # ra_delt * cos(dec) = requested physical size
    # --> ra_delt \approx requested physical size / cos(avg dec)
    ra_delt = dec_delt / cos_avg_dec

    n_x = np.fix(np.round((max_ra - min_ra) / ra_delt))
    n_y = np.fix(np.round((max_dec - min_dec) / dec_delt))
    n_x = int(np.max([n_x, 1]))
    n_y = int(np.max([n_y, 1]))
    print('# of x & y pixels = ', n_x, n_y)
    return n_x, n_y, ra_delt, dec_delt


def xyrange(n_x, n_y):
    return it.product(range(n_x), range(n_y))


def indices_for_pixel(pix_x, pix_y, x, y):
    """
    Return the indices of the sources for which the coordinates lie in
    the x, y pixel
    """
    indxs, = np.where(np.logical_and.reduce([pix_x > x, pix_x <= x + 1,
                                             pix_y > y, pix_y <= y + 1]))
    return indxs


def make_wcs_for_map(ra_grid, dec_grid):
    """make wcs corresponding to a linear ra_grid and dec_grid"""
    n_x = len(ra_grid) - 1
    n_y = len(dec_grid) - 1
    center_ra = (ra_grid.min() + ra_grid.max()) / 2.
    center_dec = (dec_grid.min() + dec_grid.max()) / 2.
    ra_delt = ra_grid[1] - ra_grid[0]
    dec_delt = dec_grid[1] - dec_grid[0]

    w = wcs.WCS(naxis=2)
    w.wcs.crpix = np.asarray([n_x, n_y], dtype=float) / 2.
    w.wcs.crval = np.array([center_ra, center_dec])
    w.wcs.cdelt = [ra_delt * math.cos(center_dec), dec_delt]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return w


def get_pix_coords(cat, map_wcs):
    """get the pixel coordinates of all the sources in the catalog
       according to the given wcs"""
    world = np.column_stack((cat['RA'], cat['DEC']))
    print('working on converting ra, dec to pix x,y')
    pixcrd = map_wcs.wcs_world2pix(world, 1)
    pix_x = pixcrd[:, 0]
    pix_y = pixcrd[:, 1]
    return pix_x, pix_y


def save_map_fits(map_data_xy, map_wcs, file_name):
    header = map_wcs.to_header()
    hdu = fits.PrimaryHDU(map_data_xy.T, header=header)
    hdu.writeto(file_name, overwrite=True)


if __name__ == '__main__':
    main()
