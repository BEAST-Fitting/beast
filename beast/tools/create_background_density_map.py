#!/usr/bin/env python3
"""
Create a background density map based on an input catalog and a fits
image. The results are saved to disk, and can be used to make sure that
areas of different background densities are properly sampled when
sampling SEDs and positions for the AST input file.

"""

import argparse
import astropy
from astropy import wcs
from astropy.table import Table
from astropy import units
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import numpy as np
import photutils as pu
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('catfile', type=str, help='catalog FITS file')
    parser.add_argument('--npix', type=int, default=10, help='resolution')
    parser.add_argument('--reference', type=str, help='reference image (FITS)',
                        default=None)
    parser.add_argument('--nointeract', action='store_true')
    args = parser.parse_args()

    ref_base = os.path.basename(args.reference).replace('.fits', '')
    hdul = astropy.io.fits.open(args.reference)
    image = hdul[1]

    result = make_background_map(
        args.catfile, args.npix, ref_im=image, outfile_base=ref_base)

    bg_map = result['background_map']
    n_map = result['nsources_map']
    ra_grid = result['ra_grid']
    dec_grid = result['dec_grid']
    mask = result['mask']

    # Plot the maps directly
    _, ax = plt.subplots(1, 2)
    figs = []
    figs.append(ax[0].imshow(bg_map))
    ax[0].set_title('density_estimate')
    figs.append(ax[1].imshow(n_map))
    ax[1].set_title('number of sources')
    for f, a in zip(figs, ax):
        plt.colorbar(f, ax=a)
    plt.savefig('{}-maps.png'.format(ref_base))

    # Overplot the map on the image used for the calculation
    if image:
        plot_on_image(image, bg_map, ra_grid, dec_grid, mask=mask,
                      title=ref_base)

    if not args.nointeract:
        plt.show()


def make_background_map(catfile, npix, ref_im, outfile_base):
    """
    Divide the image into a number of bins, and calculate the median
    background for the stars that fall within each bin. Create a new
    catalog on disk, which contains the individually measured and binned
    background estimate for each source. Create a new fits table on
    disk, containing the x,y position, bin edges (ra and dec in
    degrees), and the median background density in this bin.

    Parameters
    ----------
    catfile: str
        file name of the photometry catalog. The positions of the
        sources will be used to measure the backgrounds, and mask the
        sources themselves.

    ref_im: imageHDU
        image which will be used for the background measurements

    outfile: str
        name for the output file

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
    cat = Table.read(catfile)

    # Dimensions of the grid
    nx, ny = npix, npix

    # Coordinates
    ra = cat['RA']
    dec = cat['DEC']

    # The bin edges
    ra_grid = np.linspace(ra.min(), ra.max(), nx + 1)
    dec_grid = np.linspace(dec.min(), dec.max(), ny + 1)

    # The left sides of the bins
    ra_limits = ra_grid[:-1]
    dec_limits = dec_grid[:-1]

    # A list of background values for each source of the catalog will be
    # built up. Mask used is also returned.
    individual_backgrounds, mask = measure_backgrounds(cat, ref_im)

    # Dictionary indexed on (x,y). Contains lists of indices and
    # measurements for each bin.
    sources_foreach_bin = {}
    for x in range(nx):
        for y in range(ny):
            sources_foreach_bin[x, y] = {'indices': [], 'measurements': []}

    # Indexed on source nr i. Will contain [x,y] for each source
    bin_foreach_source = [None] * len(cat)

    # Go over all the sources, and put them in the right bins
    for i in range(len(cat)):
        # Find the correct bin in the map. With side='right', the points
        # that are exactly on the min or max get assigned the next
        # insertion point. When we do minus one, we are sure that we
        # have the point to the left of the value.
        x = np.searchsorted(ra_limits, ra[i], side='right') - 1
        y = np.searchsorted(dec_limits, dec[i], side='right') - 1

        sources_foreach_bin[x, y]['indices'].append(i)
        sources_foreach_bin[x, y]['measurements'].append(
            individual_backgrounds[i])
        bin_foreach_source[i] = [x, y]

    background_map = np.zeros((nx, ny))
    nsources_map = np.zeros((nx, ny))
    for x in range(nx):
        for y in range(ny):
            # For plotting the number sources in each bin
            n = len(sources_foreach_bin[x, y]['indices'])
            nsources_map[x, y] = n
            # Get the median background of all the sources in each bin
            if n:
                background_map[x, y] = np.median(
                    sources_foreach_bin[x, y]['measurements'])
            if n == 1:
                print('Only 1 source in bin {},{}'.format(x, y))

    background_map[nsources_map == 0] = 0

    # Save the catalog with extra density info
    median_background_foreach_source = [
        background_map[xy_bin[0], xy_bin[1]] for xy_bin in bin_foreach_source]
    extra_columns = {'indiv_bg': individual_backgrounds,
                     'bin_median_bg': median_background_foreach_source}
    for k in extra_columns:
        c = astropy.table.Column(extra_columns[k], name=k)
        cat.add_column(c)
    mod_catfile = catfile.replace(
        '.fits', '_with_{}_bg.fits'.format(outfile_base))
    cat.write(mod_catfile, format='fits', overwrite=True)

    # Save a file describing the properties of the bins in a handy format
    bin_details = astropy.table.Table(
        names=['i_ra', 'i_dec', 'median_bg',
               'min_ra', 'max_ra',
               'min_dec', 'max_dec'])
    for x in range(nx):
        for y in range(ny):
            bin_details.add_row([x, y, background_map[x, y],
                                 ra_grid[x], ra_grid[x + 1],
                                 dec_grid[y], dec_grid[y + 1]])
    bin_details.write(outfile_base + '_background_map.fits',
                      format='fits', overwrite=True)

    # Return a bunch of stuff, to be used for plots
    return {'background_map': background_map,
            'nsources_map': nsources_map,
            'ra_grid': ra_grid,
            'dec_grid': dec_grid,
            'mask': mask}


def plot_on_image(image, background_map, ra_grid, dec_grid, mask=None, title=None):
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

    title: str
        title of the resulting image, also used as prefix for the image
        file
    """
    # plot the image
    image_fig, image_ax = plt.subplots()
    imdata = image.data.astype(float)
    f = np.sort(imdata.flatten())
    p = len(f) // 32
    vmin = np.median(f[:p])
    vmax = np.median(f[-p:])
    if mask is not None:
        imdata = np.where(mask, vmin, imdata)
    plt.imshow(imdata, cmap='gray', interpolation='none', vmin=vmin,
               vmax=vmax)
    plt.colorbar()

    # If we want to rotate the rectangles so they align with the WCS
    rotation = (90. - image.header['ORIENTAT']) * math.pi / 180.

    # Make a rectangular patch for each grid point
    rectangles = []
    values = []
    image_wcs = wcs.WCS(image.header)

    height, width = imdata.shape
    height /= len(dec_grid)
    width /= len(ra_grid)

    for ix in range(len(ra_grid) - 1):
        for iy in range(len(dec_grid) - 1):
            # Current, one to the right, and one upwards. Remember
            # that the x and xup can be different because of the
            # orientation of the WCS.
            [x], [y] = image_wcs.wcs_world2pix(
                [ra_grid[ix]], [dec_grid[iy]], 0)

            rec = Rectangle((x, y), width, height)
            rot = mpl.transforms.Affine2D().rotate_around(x, y, rotation)
            rec.set_transform(rot)
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
    cb = image_fig.colorbar(patch_col)
    cb.set_alpha(1)
    cb.draw_all()
    plt.title(title)
    image_fig.savefig('{}_overlay.png'.format(title))


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
    if not ref_im:
        # Return a dumb value
        return cat_table['F814W_CHI']

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
        data_slices = list(ap_mask.bbox.slices)
        # These slices go outside of the box sometimes! In that case we
        # will override the slices on both sides.

        # Default slices (go over the whole mask)
        mask_slices = [slice(None, None), slice(None, None)]

        # Adjust the slices for x and y if necessary
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

        mask_union[data_slices] += ap_mask.data[mask_slices]

    # Threshold
    mask_union = mask_union > 0

    phot = pu.aperture_photometry(ref_im.data, annuli, wcs=w, mask=mask_union)

    return phot['aperture_sum'] / area, mask_union


if __name__ == '__main__':
    main()
