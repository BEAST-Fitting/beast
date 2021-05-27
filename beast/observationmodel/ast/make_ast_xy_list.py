import numpy as np
from matplotlib.path import Path
from tqdm import tqdm
import warnings

from astropy.io import ascii, fits
from astropy.table import Column, Table
from astropy import wcs

from shapely.geometry import box, Polygon

from beast.tools import density_map, cut_catalogs


def pick_positions_from_map(
    catalog,
    chosen_seds,
    input_map,
    bin_mode,
    N_bins,
    bin_width,
    custom_bins,
    Npermodel,
    outfile=None,
    refimage=None,
    refimage_hdu=1,
    wcs_origin=1,
    Nrealize=1,
    set_coord_boundary=None,
    region_from_filters=None,
    erode_boundary=None,
):
    """
    Spreads a set of fake stars across regions of similar values,
    given a map file generated by 'create background density map' or
    'create stellar density map' in the tools directory.

    The tiles of the given map are divided across a given
    number of bins. Each bin will then have its own set of tiles,
    which constitute a region on the image.

    Then, for each bin, the given set of fake stars is duplicated,
    and the stars are assigned random positions within this region.

    This way, it can be ensured that enough ASTs are performed for each
    regime of the map, making it possible to have a separate noise model
    for each of these regions.

    Parameters
    ----------

    catalog: Observations object
        Provides the observations

    chosen_seds: astropy Table
        Table containing fake stars to be duplicated and assigned positions

    input_map: str
        Path to a hd5 file containing the file written by a DensityMap

    bin_mode: str
        The convention for generating bins of source density. The options
        are "linear" (for linear binning) and "log" (for log binning). If "log",
        the number of bins (N_bins) must be set. If "linear", either N_bins
        or the bin width (bin_width), or neither (resulting in
        default integer binning by sources/arcsec^2) can be set.
        Default: "linear"

    N_bins: int
        The number of bins for the range of background density values.
        The bins will be picked on a linear grid or log grid (according to bin_mode)
        ranging from the minimum to the maximum value of the map. Then, each tile will be
        put in a bin, so that a set of tiles of the map is obtained for
        each range of source density/background values.

    bin_width: int
        The bin width for the range of  background density values, in units
        of number of sources per square arcsecond.
        The bins will be picked on a linear grid, ranging from the
        minimum to the maximum value of the map. Then, each tile will be
        put in a bin, so that a set of tiles of the map is obtained for
        each range of source density/background values.

    custom_bins: list (default=None)
        Custom values of bin edges for source or background density values.
        Each tile will be put into a bin, so that a set of tiles of the
        map is obtained for each range of source density/background values.

    refimage: str
        Path to fits image that is used for the positions. If none is
        given, the ra and dec will be put in the x and y output columns
        instead.

    refimage_hdu: int (default=1)
        index of the HDU from which to get the header, which will be used
        to extract WCS information

    wcs_origin : 0 or 1 (default=1)
        As described in the WCS documentation: "the coordinate in the upper
        left corner of the image. In FITS and Fortran standards, this is 1.
        In Numpy and C standards this is 0."

    Nrealize: integer
        The number of times each model should be repeated for each
        background regime. This is to sample the variance due to
        variations within each region, for each individual model.

    set_coord_boundary : None, or list of 2 numpy arrays
        If provided, these RA/Dec coordinates will be used to limit the
        region over which ASTs are generated.  Input should be list of two
        arrays, the first RA and the second Dec, ordered sequentially
        around the region (either CW or CCW).  If the input catalog only has x/y
        (no RA/Dec), a refimage is required.

    region_from_filters : None, list of filter name(s), or 'all'
        If provided, ASTs will only be placed in regions with this particular
        combination of filter(s).  Or, if 'all' is chosen, ASTs will only be
        placed where there is overlap with all filters.  In practice, this
        means creating a convex hull around the catalog RA/DEC of sources with
        valid values in these filters.  Note that if the region in question is
        a donut, this will put ASTs in the hole.  This will also only work
        properly if the region is a convex polygon.  A solution to these needs
        to be figured out at some point.

    erode_boundary : None, or float (default=None)
        If provided, this number of arcseconds will be eroded from the region
        over which ASTs are generated.  The purpose is to avoid placing ASTs
        near the image edge.  Erosion is applied to both the catalog boundary
        and the values from set_coord_boundary.  If the input catalog only has
        x/y (no RA/Dec), a refimage is required.

    Returns
    -------
    astropy Table: List of fake stars, with magnitudes and positions
    - optionally -
    ascii file of this table, written to outfile

    """

    # if refimage exists, extract WCS info
    if refimage is None:
        ref_wcs = None
    else:
        with fits.open(refimage) as hdu:
            imagehdu = hdu[refimage_hdu]
            ref_wcs = wcs.WCS(imagehdu.header)

    # if appropriate information is given, extract the x/y positions so that
    # there are no ASTs generated outside of the catalog footprint
    colnames = catalog.data.columns
    xy_pos = False
    radec_pos = False

    # if x/y in catalog, save them
    if ("X" in colnames) or ("x" in colnames):
        xy_pos = True
        if "X" in colnames:
            x_positions = catalog.data["X"][:]
            y_positions = catalog.data["Y"][:]
        if "x" in colnames:
            x_positions = catalog.data["x"][:]
            y_positions = catalog.data["y"][:]

    # if RA/Dec in catalog, save them
    if ("RA" in colnames) or ("ra" in colnames):
        radec_pos = True
        if "RA" in colnames:
            ra_positions = catalog.data["RA"][:]
            dec_positions = catalog.data["DEC"][:]
        if "ra" in colnames:
            ra_positions = catalog.data["ra"][:]
            dec_positions = catalog.data["dec"][:]

    # # if only one of those exists and there's a refimage, convert to the other
    # if xy_pos and not radec_pos and refimage:
    #     radec_pos = True
    #     x_positions, y_positions = ref_wcs.all_world2pix(
    #         ra_positions, dec_positions, wcs_origin
    #     )
    # if radec_pos and not xy_pos and refimage:
    #     xy_pos = True
    #     ra_positions, dec_positions = ref_wcs.all_pix2world(
    #         x_positions, y_positions, wcs_origin
    #     )

    # if only one of those exists and there's a refimage, convert to the other
    if xy_pos and not radec_pos and refimage:
        radec_pos = True
        ra_positions, dec_positions = ref_wcs.all_pix2world(
            x_positions, y_positions, wcs_origin

        )
    if radec_pos and not xy_pos and refimage:
        xy_pos = True
        x_positions, y_positions = ref_wcs.all_world2pix(
            ra_positions, dec_positions, wcs_origin
        )

    # if no x/y or ra/dec in the catalog, raise error
    if not xy_pos and not radec_pos:
        raise RuntimeError(
            "Your catalog does not supply X/Y or RA/DEC information to ensure ASTs are within catalog boundary"
        )

    # if erode_boundary is set, try to make a pixel version to go with xy positions
    erode_deg = None
    erode_pix = None
    if erode_boundary:
        erode_deg = erode_boundary / 3600
        if xy_pos and refimage:
            deg_per_pix = wcs.utils.proj_plane_pixel_scales(ref_wcs)[0]
            erode_pix = erode_deg / deg_per_pix

    # create path containing the positions (eroded if chosen)
    catalog_boundary_xy = None
    catalog_boundary_radec = None
    if xy_pos:
        catalog_boundary_xy = erode_path(
            cut_catalogs.convexhull_path(x_positions, y_positions), erode_pix
        )
    if radec_pos:
        catalog_boundary_radec = erode_path(
            cut_catalogs.convexhull_path(ra_positions, dec_positions), erode_deg
        )

    # if coord_boundary set, define an additional boundary for ASTs (eroded if chosen)
    if set_coord_boundary is not None:
        # initialize variables
        coord_boundary_xy = None
        coord_boundary_radec = None
        # evaluate one or both
        if xy_pos and refimage:
            bounds_x, bounds_y = ref_wcs.all_world2pix(
                set_coord_boundary[0], set_coord_boundary[1], wcs_origin
            )
            coord_boundary_xy = erode_path(
                Path(np.array([bounds_x, bounds_y]).T), erode_pix
            )
        if radec_pos:
            coord_boundary_radec = erode_path(
                Path(np.array([set_coord_boundary[0], set_coord_boundary[1]]).T),
                erode_deg,
            )

    # if region_from_filters is set, define an additional boundary for ASTs
    if region_from_filters is not None:

        # 1. find the sub-list of sources
        if isinstance(region_from_filters, list):
            # good stars with user-defined partial overlap
            _, good_stars = cut_catalogs.cut_catalogs(
                catalog.inputFile,
                "N/A",
                flagged=True,
                flag_filter=region_from_filters,
                no_write=True,
            )
        elif region_from_filters == "all":
            # good stars only with fully overlapping region
            _, good_stars = cut_catalogs.cut_catalogs(
                catalog.inputFile, "N/A", partial_overlap=True, no_write=True
            )
        else:
            raise RuntimeError("Invalid argument for region_from_filters")

        # 2. define the Path object for the convex hull
        # initialize variables
        filt_reg_boundary_xy = None
        filt_reg_boundary_radec = None
        # evaluate one or both
        if xy_pos:
            filt_reg_boundary_xy = cut_catalogs.convexhull_path(
                x_positions[good_stars == 1], y_positions[good_stars == 1]
            )
        if radec_pos:
            filt_reg_boundary_radec = cut_catalogs.convexhull_path(
                ra_positions[good_stars == 1], dec_positions[good_stars == 1]
            )

    # Load the background map
    print(Npermodel, " repeats of each model in each map bin")

    bdm = density_map.BinnedDensityMap.create(
        input_map,
        bin_mode=bin_mode,
        N_bins=N_bins,
        bin_width=bin_width,
        custom_bins=custom_bins,
    )

    tile_vals = bdm.tile_vals()
    max_val = np.amax(tile_vals)
    min_val = np.amin(tile_vals)
    tiles_foreach_bin = bdm.tiles_foreach_bin()

    # Remove any of the tiles that aren't contained within user-imposed
    # constraints (if any)
    if (set_coord_boundary is not None) or (region_from_filters is not None):

        tile_ra_min, tile_dec_min = bdm.min_ras_decs()
        tile_ra_delta, tile_dec_delta = bdm.delta_ras_decs()

        for i, tile_set in enumerate(tiles_foreach_bin):

            # keep track of which indices to discard
            keep_tile = np.ones(len(tile_set), dtype=bool)

            for j, tile in enumerate(tile_set):

                # corners of the tile
                ra_min = tile_ra_min[tile]
                ra_max = tile_ra_min[tile] + tile_ra_delta[tile]
                dec_min = tile_dec_min[tile]
                dec_max = tile_dec_min[tile] + tile_dec_delta[tile]

                # make a box object for the tile
                tile_box_radec = box(ra_min, dec_min, ra_max, dec_max)
                tile_box_xy = None
                if refimage:
                    bounds_x, bounds_y = ref_wcs.all_world2pix(
                        np.array([ra_min, ra_max]),
                        np.array([dec_min, dec_max]),
                        wcs_origin,
                    )

                    tile_box_xy = box(
                        np.min(bounds_x),
                        np.min(bounds_y),
                        np.max(bounds_x),
                        np.max(bounds_y),
                    )

                # discard tile if there's no overlap with user-imposed regions

                # - erode_boundary
                # if you only want to erode the boundary and not impose other
                # coordinate boundary constraints, still discard SD tiles that don't overlap
                if (set_coord_boundary is None) and (erode_boundary is not None):
                    if catalog_boundary_radec and tile_box_radec:
                        if (
                            Polygon(catalog_boundary_radec.vertices)
                            .intersection(tile_box_radec)
                            .area
                            == 0
                        ):
                            keep_tile[j] = False
                    elif catalog_boundary_xy and tile_box_xy:
                        if (
                            Polygon(catalog_boundary_xy.vertices)
                            .intersection(tile_box_xy)
                            .area
                            == 0
                        ):
                            keep_tile[j] = False
                # - set_coord_boundary
                if set_coord_boundary is not None:
                    # coord boundary is input in RA/Dec, and tiles are RA/Dec,
                    # so there's no need to check the x/y version of either
                    if (
                        Polygon(coord_boundary_radec.vertices)
                        .intersection(tile_box_radec)
                        .area
                        == 0
                    ):
                        keep_tile[j] = False

                # - region_from_filters
                if region_from_filters is not None:
                    if filt_reg_boundary_radec and tile_box_radec:
                        if (
                            Polygon(filt_reg_boundary_radec.vertices)
                            .intersection(tile_box_radec)
                            .area
                            == 0
                        ):
                            keep_tile[j] = False
                    elif filt_reg_boundary_xy and tile_box_xy:
                        if (
                            Polygon(filt_reg_boundary_xy.vertices)
                            .intersection(tile_box_xy)
                            .area
                            == 0
                        ):
                            keep_tile[j] = False
                    else:
                        warnings.warn(
                            "Unable to use regions_from_filters to remove SD/bg tiles"
                        )

            # remove anything that needs to be discarded
            tiles_foreach_bin[i] = tile_set[keep_tile]

    # Remove empty bins
    tile_sets = [tile_set for tile_set in tiles_foreach_bin if len(tile_set)]
    print(
        "{0} non-empty map bins (out of {1}) found between {2} and {3}".format(
            len(tile_sets), N_bins, min_val, max_val
        )
    )

    # Repeat the seds Nrealize times (sample each on at Nrealize
    # different positions, in each region)
    repeated_seds = np.repeat(chosen_seds, Nrealize)
    Nseds_per_region = len(repeated_seds)
    # For each set of tiles, repeat the seds and spread them evenly over
    # the tiles
    repeated_seds = np.repeat(repeated_seds, len(tile_sets))

    out_table = Table(repeated_seds, names=chosen_seds.colnames)
    ast_x_list = np.zeros(len(out_table))
    ast_y_list = np.zeros(len(out_table))
    bin_indices = np.zeros(len(out_table))

    tile_ra_min, tile_dec_min = bdm.min_ras_decs()
    tile_ra_delta, tile_dec_delta = bdm.delta_ras_decs()

    for bin_index, tile_set in enumerate(
        tqdm(
            tile_sets,
            desc="{:.2f} models per map bin".format(Nseds_per_region / Npermodel),
        )
    ):
        start = bin_index * Nseds_per_region
        stop = start + Nseds_per_region
        bin_indices[start:stop] = bin_index
        for i in range(Nseds_per_region):

            # keep track of whether we're still looking for valid coordinates
            x = None
            y = None

            while (x is None) or (y is None):
                # Pick a random tile in this tile set
                tile = np.random.choice(tile_set)
                # Within this tile, pick a random ra and dec
                ra = tile_ra_min[tile] + np.random.random_sample() * tile_ra_delta[tile]
                dec = (
                    tile_dec_min[tile]
                    + np.random.random_sample() * tile_dec_delta[tile]
                )

                # if we can't convert this to x/y, do everything in RA/Dec
                if ref_wcs is None:
                    x, y = ra, dec

                    # check that this x/y is within the catalog footprint
                    if catalog_boundary_radec:
                        # N,2 array of AST X and Y positions
                        inbounds = catalog_boundary_radec.contains_points([[x, y]])[0]

                        if not inbounds:
                            x = None

                    # check that this x/y is with any input boundary
                    if set_coord_boundary is not None:
                        if coord_boundary_radec:
                            inbounds = coord_boundary_radec.contains_points([[x, y]])[0]
                            if not inbounds:
                                x = None
                    if region_from_filters is not None:
                        if filt_reg_boundary_radec:
                            # fmt: off
                            inbounds = filt_reg_boundary_radec.contains_points([[x, y]])[0]
                            # fmt: on
                            if not inbounds:
                                x = None

                # if we can convert to x/y, do everything in x/y
                else:
                    [x], [y] = ref_wcs.all_world2pix(
                        np.array([ra]), np.array([dec]), wcs_origin
                    )

                    # check that this x/y is within the catalog footprint
                    # N,2 array of AST X and Y positions
                    inbounds = catalog_boundary_xy.contains_points([[x, y]])[0]
                    if not inbounds:
                        x = None

                    # check that this x/y is with any input boundary
                    if set_coord_boundary is not None:
                        if coord_boundary_xy:
                            inbounds = coord_boundary_xy.contains_points([[x, y]])[0]
                            if not inbounds:
                                x = None
                    if region_from_filters is not None:
                        if filt_reg_boundary_xy:
                            inbounds = filt_reg_boundary_xy.contains_points([[x, y]])[0]
                            if not inbounds:
                                x = None

            j = bin_index + i * len(tile_sets)
            ast_x_list[j] = x
            ast_y_list[j] = y

    # I'm just mimicking the format that is produced by the examples
    cs = []
    cs.append(Column(np.zeros(len(out_table), dtype=int), name="zeros"))
    cs.append(Column(np.ones(len(out_table), dtype=int), name="ones"))

    # positions were found using RA/Dec
    if ref_wcs is None:
        cs.append(Column(ast_x_list, name="RA"))
        cs.append(Column(ast_y_list, name="DEC"))
    # positions were found using x/y
    else:
        cs.append(Column(ast_x_list, name="X"))
        cs.append(Column(ast_y_list, name="Y"))

    for i, c in enumerate(cs):
        out_table.add_column(c, index=i)  # insert these columns from the left

    # Write out the table in ascii
    if outfile:
        formats = {k: "%.5f" for k in out_table.colnames[2:]}
        ascii.write(out_table, outfile, overwrite=True, formats=formats)

    return out_table


def erode_path(path_object, erode_amount):
    """
    Returns the original Path object, but eroded by the defined amount.

    Parameters
    ----------
    path_object : Path object
        the Path object to be eroded

    erode_amount : float or None
        If float, amount to erode. Units (pixels or degrees) should match units
        of path_object. Absolute value will be used.
        If None, don't do any eroding - return original path_object.

    Returns
    -------
    Path object : the original path object (erode_amount=None) or the eroded
    path object (erode_amount=float)
    """

    if erode_amount is None:
        return path_object
    else:
        eroded_polygon = Polygon(path_object.vertices).buffer(-np.abs(erode_amount))
        bounds_x = [float(i) for i in eroded_polygon.exterior.coords.xy[0]]
        bounds_y = [float(i) for i in eroded_polygon.exterior.coords.xy[1]]
        return Path(np.array([bounds_x, bounds_y]).T)


def pick_positions(
    catalog, chosen_sed_file, outfile, separation, refimage=None, wcs_origin=1
):
    """
    Assigns positions to fake star list generated by pick_models

    Parameters
    ----------

    catalog : Observations object
        Provides the observations

    chosen_sed_file : string
        file containing fake stars to be assigned positions

    outfile : string
        Name of the file to save the fake star list

    separation : float
        Minimum pixel separation between AST and star in photometry
        catalog provided in the beast settings.

    refimage : string
        Name of the reference image.  If supplied, the method will use the
        reference image header to convert from RA and DEC to X and Y.

    wcs_origin : 0 or 1 (default=1)
        As described in the WCS documentation: "the coordinate in the upper
        left corner of the image. In FITS and Fortran standards, this is 1.
        In Numpy and C standards this is 0."

    Returns
    -------

    Ascii table that replaces [chosen_sed_file] with a new version, [outfile],
    that contains the necessary position columns for running the ASTs though
    DOLPHOT
    """

    noise = 3.0  # Spreads the ASTs in a circular annulus of 3 pixel width instead of all being
    # precisely [separation] from an observed star.

    colnames = catalog.data.columns

    if "X" or "x" in colnames:
        if "X" in colnames:
            x_positions = catalog.data["X"][:]
            y_positions = catalog.data["Y"][:]
        if "x" in colnames:
            x_positions = catalog.data["x"][:]
            y_positions = catalog.data["y"][:]
    else:
        if refimage:
            if ("RA" in colnames) or ("ra" in colnames):
                if "RA" in colnames:
                    ra_positions = catalog.data["RA"][:]
                    dec_positions = catalog.data["DEC"][:]
                if "ra" in colnames:
                    ra_positions = catalog.data["ra"][:]
                    dec_positions = catalog.data["dec"][:]
            else:
                raise RuntimeError(
                    "Your catalog does not supply X, Y or RA, DEC information for spatial AST distribution"
                )

        else:
            raise RuntimeError(
                "You must supply a Reference Image to determine spatial AST distribution."
            )
        ref_wcs = wcs.WCS(refimage)

        x_positions, y_positions = ref_wcs.all_world2pix(
            ra_positions, dec_positions, wcs_origin
        )

    astmags = ascii.read(chosen_sed_file)

    n_asts = len(astmags)

    # keep is defined to ensure that no fake stars are put outside of the image boundaries

    keep = (
        (x_positions > np.min(x_positions) + separation + noise)
        & (x_positions < np.max(x_positions) - separation - noise)
        & (y_positions > np.min(y_positions) + separation + noise)
        & (y_positions < np.max(y_positions) - separation - noise)
    )

    x_positions = x_positions[keep]
    y_positions = y_positions[keep]

    ncat = len(x_positions)
    ind = np.random.random(n_asts) * ncat
    ind = ind.astype("int")

    # Here we generate the circular distribution of ASTs surrounding random observed stars

    separation = np.random.random(n_asts) * noise + separation
    theta = np.random.random(n_asts) * 2.0 * np.pi
    xvar = separation * np.cos(theta)
    yvar = separation * np.sin(theta)

    new_x = x_positions[ind] + xvar
    new_y = y_positions[ind] + yvar
    column1 = 0 * new_x
    column2 = column1 + 1
    column1 = Column(name="zeros", data=column1.astype("int"))
    column2 = Column(name="ones", data=column2.astype("int"))
    column3 = Column(name="X", data=new_x, format="%.2f")
    column4 = Column(name="Y", data=new_y, format="%.2f")
    astmags.add_column(column1, 0)
    astmags.add_column(column2, 1)
    astmags.add_column(column3, 2)
    astmags.add_column(column4, 3)

    ascii.write(astmags, outfile, overwrite=True)
