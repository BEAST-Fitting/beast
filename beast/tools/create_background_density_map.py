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
from astropy import units as u
from astropy.io import fits
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import numpy as np
import photutils as pu
from beast.tools.density_map import DensityMap
from beast.tools import cut_catalogs
import itertools as it
import os
from shapely import geometry


def main():  # pragma: no cover
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand")
    subparsers.required = True

    # common options for both types of map
    commonparser = argparse.ArgumentParser(add_help=False)
    commonparser.add_argument(
        "-catfile", type=str, required=True, help="catalog FITS file"
    )
    commonparser.add_argument(
        "-erode_boundary",
        type=float,
        required=False,
        help="number of arcsec to erode the SD map boundary by",
    )
    npix_or_pixsize = commonparser.add_mutually_exclusive_group()
    npix_or_pixsize.add_argument("--npix", type=int, default=None, help="resolution")
    npix_or_pixsize.add_argument("--pixsize", type=float, default=None)

    # map making subcommands
    background_parser = subparsers.add_parser(
        "background",
        parents=[commonparser],
        help="""Create a background intensity map based on annulus
        measurements around the sources listed in the catalog""",
    )
    sourceden_parser = subparsers.add_parser(
        "sourceden",
        parents=[commonparser],
        help="""Create a source density map, counting the number of
        sources in each tile of the map""",
    )

    # plot subcommand
    plot_parser = subparsers.add_parser(
        "tileplot",
        help="""Plot the given map using transparent colored tiles on
        top of a given reference image""",
    )

    # arguments unique to background map
    background_parser.add_argument(
        "-reference",
        type=str,
        metavar="FITSIMAGE",
        required=True,
        help="reference image (FITS)",
    )
    background_parser.add_argument(
        "--mask_radius",
        type=float,
        metavar="RADIUS",
        default=30,
        help="radius (in pixels) of mask for catalog sources",
    )
    background_parser.add_argument(
        "--ann_width",
        type=float,
        metavar="ANNULUS_WIDTH",
        default=20,
        help="width of annulus (in pixels) for calculating bkgd around each catalog source",
    )
    background_parser.add_argument(
        "--cat_filter",
        type=str,
        nargs=2,
        metavar="FILTER MAG",
        default=None,
        help="catalog entries with FILTER_VEGA > MAG will not be masked",
    )

    # arguments unique to sourceden map
    sourceden_parser.add_argument(
        "--mag_cut",
        type=float,
        nargs=2,
        metavar="MIN MAX",
        default=[24.5, 27],
        help="magnitude range on which the source density is computed",
    )
    sourceden_parser.add_argument(
        "--mag_name",
        type=str,
        default="F475W_VEGA",
        metavar="FILTER",
        help="name of magnitude column in table",
    )
    sourceden_parser.add_argument(
        "--flag_name",
        type=str,
        default=None,
        metavar="FILTER_FLAG",
        help="if set, ignore sources with flag >= 99",
    )

    # options unique to plot command
    plot_parser.add_argument(
        "densitymap",
        metavar="MAP.HD5",
        help="the map to plot (e.g. the output of this script)",
    )
    plot_parser.add_argument(
        "-image",
        metavar="IM.FITS",
        type=str,
        required=True,
        help="image to overplot the tiles onto",
    )
    plot_parser.add_argument(
        "--colorbar",
        metavar="LABEL",
        type=str,
        help="use colorbar, and use the given string as label",
    )
    plot_parser.add_argument("--dpi", type=float, help="dpi for the resulting pdf")

    args = parser.parse_args()
    if args.subcommand in ["background", "sourceden"]:
        main_make_map(args)
    elif args.subcommand == "tileplot":
        main_plot(args)


def main_make_map(args):
    # Common actions: load the catalog and set up the grid
    cat = Table.read(args.catfile)
    for name in cat.colnames:
        cat.rename_column(name, name.upper())

    ra = cat["RA"]
    dec = cat["DEC"]

    # if erode_boundary is set, erode catalog boundary to compute
    # source densities on eroded region
    if args.erode_boundary:
        erode_deg = args.erode_boundary / 3600
        min_ra = ra.min() + erode_deg
        max_ra = ra.max() - erode_deg
        min_dec = dec.min() + erode_deg
        max_dec = dec.max() - erode_deg
    else:
        min_ra = ra.min()
        max_ra = ra.max()
        min_dec = dec.min()
        max_dec = dec.max()

    if args.npix is not None:
        n_x, n_y = args.npix, args.npix
        ra_grid = np.linspace(min_ra, max_ra, n_x + 1)
        dec_grid = np.linspace(min_dec, max_dec, n_y + 1)
    elif args.pixsize is not None:
        pixsize_arcsecs = args.pixsize * u.arcsec
        pixsize_degrees = pixsize_arcsecs.to(u.degree)
        n_x, n_y, ra_delt, dec_delt = calc_nx_ny_from_pixsize(cat, pixsize_degrees)
        # the ra spacing needs to be larger, as 1 degree of RA ==
        # cos(DEC) degrees on the great circle
        ra_grid = min_ra + ra_delt * np.arange(0, n_x + 1, dtype=float)
        dec_grid = min_dec + dec_delt * np.arange(0, n_y + 1, dtype=float)
    else:
        n_x, n_y = 10, 10
        ra_grid = np.linspace(min_ra, max_ra, n_x + 1)
        dec_grid = np.linspace(min_dec, max_dec, n_y + 1)

    output_base = args.catfile.replace(".fits", "")

    if args.subcommand == "sourceden":
        map_values_array = make_source_dens_map(
            cat,
            ra_grid,
            dec_grid,
            output_base,
            mag_name=args.mag_name,
            mag_cut=args.mag_cut,
            flag_name=args.flag_name,
            erode_boundary=args.erode_boundary,
        )

    if args.subcommand == "background":
        hdul = astropy.io.fits.open(args.reference)
        image = hdul["SCI"]
        map_values_array, n_map = make_background_map(
            cat,
            ra_grid,
            dec_grid,
            ref_im=image,
            mask_radius=args.mask_radius,
            ann_width=args.ann_width,
            cat_filter=args.cat_filter,
            output_base=output_base,
        )

    # Save a file describing the properties of the bins in a handy format
    bin_details = astropy.table.Table(
        names=["i_ra", "i_dec", "value", "min_ra", "max_ra", "min_dec", "max_dec"]
    )
    for x, y in xyrange(n_x, n_y):
        bin_details.add_row(
            [
                x,
                y,
                map_values_array[x, y],
                ra_grid[x],
                ra_grid[x + 1],
                dec_grid[y],
                dec_grid[y + 1],
            ]
        )

    # Add the ra and dec grids as metadata
    bin_details.meta["ra_grid"] = ra_grid
    bin_details.meta["dec_grid"] = dec_grid

    # This works for both the background density map or the source
    # density map
    dm = DensityMap(bin_details)
    dm.write("{}_{}_map.hd5".format(output_base, args.subcommand))


def main_plot(args):
    dm = DensityMap(args.densitymap)

    hdul = astropy.io.fits.open(args.image)
    image = hdul["SCI"]

    image_fig, image_ax, patch_col = plot_on_image(dm, image)

    if args.colorbar is not None:
        cb = image_fig.colorbar(patch_col)
        cb.set_alpha(1)
        cb.draw_all()

        if len(args.colorbar) == 0:
            label = "background (arbitrary units)"
        else:
            label = args.colorbar

        cb.set_label(label)

    output_base = os.path.basename(args.densitymap).replace(".hd5", "")
    image_fig.savefig("{}_plot_overlay.pdf".format(output_base), dpi=args.dpi)


def make_background_map(
    cat, ra_grid, dec_grid, ref_im, mask_radius, ann_width, cat_filter, output_base
):
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

    mask_radius : float
        radius (in pixels) of mask for catalog sources

    ann_width : float
        width of annulus (in pixels) for calculating background around each catalog source

    cat_filter : list or None
        If list: Two elements in which the first is a filter (e.g. 'F475W') and
        the second is a magnitude.  Catalog entries with [filter]_VEGA > mag
        will not be masked.
        If None: all catalog entries will be considered.

    output_base: string
        base name (without extension) to be used for the output files

    Returns
    -------
    results: background_map, nsources_map: 2d ndarray, 2d ndarray

    """
    # A list of background values for each source of the catalog will be
    # built up. Mask used is also returned.
    individual_backgrounds = measure_backgrounds(
        cat, ref_im, mask_radius, ann_width, cat_filter
    )

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
            print("Only 1 source in bin {},{}".format(x, y))

        # store the median background for each source
        median_backgrounds[idxs] = background_map[x, y]

    background_map[nsources_map == 0] = 0

    # Save the catalog with extra density info
    extra_columns = {
        "indiv_bg": individual_backgrounds,
        "bin_median_bg": median_backgrounds,
    }
    for k in extra_columns:
        c = astropy.table.Column(extra_columns[k], name=k)
        cat.add_column(c)
    cat.write(output_base + "_with_bg.fits", format="fits", overwrite=True)

    # Save the map as a fits file here. The map in hd5 format is written
    # by some common code, in the main() function .
    save_map_fits(background_map, w, output_base + "_background.fits")
    save_map_fits(nsources_map, w, output_base + "_nsources.fits")
    return background_map, nsources_map


def measure_backgrounds(cat_table, ref_im, mask_radius, ann_width, cat_filter):
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

    mask_radius : float
        radius (in pixels) of mask for catalog sources

    ann_width : float
        width of annulus (in pixels) for calculating background around each catalog source

    cat_filter : list or None
        If list: Two elements in which the first is a filter (e.g. 'F475W') and
        the second is a magnitude.  Catalog entries with [filter]_VEGA > mag
        will not be masked.
        If None: all catalog entries will be considered.


    Returns
    -------

    measurements: list of float
        a metric for the background intensity at the star's position
        (photometry / area)

    """
    w = wcs.WCS(ref_im.header)
    shp = ref_im.data.shape

    inner_rad = mask_radius * u.pixel
    outer_rad = inner_rad + ann_width * u.pixel
    mask_rad = inner_rad

    # More elaborate way using an actual image (do not care about the
    # filter for the moment, just take the image (which was given on the
    # command line) at face value)
    ra = cat_table["RA"]
    dec = cat_table["DEC"]
    c = astropy.coordinates.SkyCoord(ra * u.degree, dec * u.degree)

    # Annuli, of which the counts per surface area will be used as
    # background measurements
    annuli = pu.SkyCircularAnnulus(c, r_in=inner_rad, r_out=outer_rad)
    area = annuli.to_pixel(w).area

    # A mask to make sure that no sources end up in the background
    # calculation
    if cat_filter is None:
        circles = pu.SkyCircularAperture(c, mask_rad)
    else:
        circles = pu.SkyCircularAperture(
            c[cat_table[cat_filter[0] + "_VEGA"] < float(cat_filter[1])], mask_rad
        )
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

    # also mask NaNs
    mask_union[np.isnan(ref_im.data)] = True

    # Save the masked reference image
    hdu = fits.PrimaryHDU(np.where(mask_union, 0, ref_im.data), header=ref_im.header)
    hdu.writeto("masked_reference_image.fits", overwrite=True)

    # Do the measurements
    phot = pu.aperture_photometry(ref_im.data, annuli, wcs=w, mask=mask_union)
    return phot["aperture_sum"] / area


def make_source_dens_map(
    cat, ra_grid, dec_grid, output_base, mag_name, mag_cut, flag_name, erode_boundary
):
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

    flag_name : string or None
        if set, ignore sources with flag >= 99

    erode_boundary : None, or float (default=None)
        If provided, this number of arcseconds will be eroded from the region
        over which SD is estimated.  The purpose is to make sure that the SD
        bins are the same ones used for placing ASTs (i.e., so that SD bins will not
        fall outisde the region used for AST placement, thus ensuring all SD bins
        are well-populated).

    OUTPUT:
    -------
    FITS files written to disk.
    """

    # force filter magnitude name to be upper case to match column names
    mag_name = mag_name.upper()

    N_stars = len(cat)

    w = make_wcs_for_map(ra_grid, dec_grid)
    pix_x, pix_y = get_pix_coords(cat, w)

    # make a convex hull of the catalog
    catalog_boundary = geometry.Polygon(
        cut_catalogs.convexhull_path(pix_x, pix_y).vertices
    )

    n_x = len(ra_grid) - 1
    n_y = len(dec_grid) - 1
    npts_map = np.zeros([n_x, n_y], dtype=float)
    source_dens = np.zeros(N_stars, dtype=float)

    # area of one pixel in square degrees
    pix_area = w.wcs.cdelt[0] * w.wcs.cdelt[1] * 3600 ** 2

    for i, j in xyrange(n_x, n_y):
        indxs = indices_for_pixel(pix_x, pix_y, i, j)

        if flag_name is None:
            (indxs_for_SD,) = np.where(
                (cat[mag_name][indxs] >= mag_cut[0])
                & (cat[mag_name][indxs] <= mag_cut[1])
            )
        else:
            flag_name = flag_name.upper()
            (indxs_for_SD,) = np.where(
                (cat[mag_name][indxs] >= mag_cut[0])
                & (cat[mag_name][indxs] <= mag_cut[1])
                & (cat[flag_name][indxs] < 99)
            )

        n_indxs = len(indxs_for_SD)
        if n_indxs > 0:
            # make a box for the current SD map pixel
            pix_box = geometry.box(i, j, i + 1, j + 1)
            # find fractional overlap area
            frac_area = catalog_boundary.intersection(pix_box).area

            if frac_area > 0.5:
                # stars per unit area
                npts_map[i, j] = n_indxs / (pix_area * frac_area)

        # save the source density as an entry for each source
        source_dens[indxs] = npts_map[i, j]

    save_map_fits(npts_map, w, output_base + "_source_den_image.fits")

    # Save the source density for individual stars in a new catalog file
    cat["SourceDensity"] = source_dens
    cat.write(output_base + "_with_sourceden.fits", overwrite=True)

    return npts_map


def plot_on_image(densitymap, image):
    """
    Plot the density grid as a collection of colored rectangles layered
    over the given fits image.

    Parameters
    ----------
    densitymap: DensityMap
        the density map instance

    image: imageHDU
        the fits image, which should include a WCS

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
    plt.imshow(imdata, cmap="gray_r", interpolation="mitchell", vmin=vmin, vmax=vmax)

    # Make a rectangular patch for each tile of the map
    rectangles = []
    for row in densitymap.tile_data:
        ll = row["min_ra"]  # left
        r = row["max_ra"]  # right
        b = row["min_dec"]  # bottom
        t = row["max_dec"]  # top
        ra_dec_corners = np.array([[ll, b], [r, b], [r, t], [ll, t]])
        pix_corners = image_wcs.wcs_world2pix(ra_dec_corners, 0)
        rec = Polygon(pix_corners, closed=True)
        rectangles.append(rec)

    # This will give each tile a color according to the given values.
    # Some of the plt.scatter functionality works with PatchCollection
    # under the hood, apparently.
    patch_col = PatchCollection(rectangles, cmap="viridis")
    patch_col.set_alpha(0.3)
    patch_col.set_array(densitymap.tile_vals())
    image_ax.add_collection(patch_col)
    return image_fig, image_ax, patch_col


def calc_nx_ny_from_pixsize(cat, pixsize_degrees):
    min_ra = cat["RA"].min()
    max_ra = cat["RA"].max()
    min_dec = cat["DEC"].min()
    max_dec = cat["DEC"].max()

    # Compute the required width of the bins expressed in RA and DEC to
    # reach the requested physical pixel size on the sky
    if type(pixsize_degrees) == astropy.units.quantity.Quantity:
        dec_delt = pixsize_degrees.value
    else:
        dec_delt = pixsize_degrees
    cos_avg_dec = math.cos(math.radians((max_dec + min_dec) / 2))
    # ra_delt * cos(dec) = requested physical size
    # --> ra_delt \approx requested physical size / cos(avg dec)
    ra_delt = dec_delt / cos_avg_dec

    n_x = np.ceil((max_ra - min_ra) / ra_delt)
    n_y = np.ceil((max_dec - min_dec) / dec_delt)
    n_x = int(np.max([n_x, 1]))
    n_y = int(np.max([n_y, 1]))
    print("# of x & y pixels = ", n_x, n_y)
    return n_x, n_y, ra_delt, dec_delt


def xyrange(n_x, n_y):
    return it.product(range(n_x), range(n_y))


def indices_for_pixel(pix_x, pix_y, x, y):
    """
    Return the indices of the sources for which the coordinates lie in
    the x, y pixel
    """
    (indxs,) = np.where(
        np.logical_and.reduce([pix_x > x, pix_x <= x + 1, pix_y > y, pix_y <= y + 1])
    )
    return indxs


def make_wcs_for_map(ra_grid, dec_grid):
    """make wcs corresponding to a linear ra_grid and dec_grid"""
    n_x = len(ra_grid) - 1
    n_y = len(dec_grid) - 1
    center_ra = (ra_grid.min() + ra_grid.max()) / 2.0
    center_dec = (dec_grid.min() + dec_grid.max()) / 2.0
    ra_delt = ra_grid[1] - ra_grid[0]
    phys_ra_delt = ra_delt * math.cos(math.radians(center_dec))
    dec_delt = dec_grid[1] - dec_grid[0]

    w = wcs.WCS(naxis=2)
    w.wcs.crpix = np.asarray([n_x, n_y], dtype=float) / 2.0 + 0.5
    w.wcs.crval = np.array([center_ra, center_dec])
    w.wcs.cdelt = np.abs([-phys_ra_delt, dec_delt])
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    return w


def get_pix_coords(cat, map_wcs):
    """get the pixel coordinates of all the sources in the catalog
    according to the given wcs"""
    world = np.column_stack((cat["RA"], cat["DEC"]))
    print("working on converting ra, dec to pix x,y")
    pixcrd = map_wcs.wcs_world2pix(world, 1) - 0.5
    pix_x = pixcrd[:, 0]
    pix_y = pixcrd[:, 1]
    return pix_x, pix_y


def save_map_fits(map_data_xy, map_wcs, file_name):
    header = map_wcs.to_header()
    hdu = fits.PrimaryHDU(map_data_xy.T, header=header)
    hdu.writeto(file_name, overwrite=True)


if __name__ == "__main__":  # pragma: no cover
    main()
