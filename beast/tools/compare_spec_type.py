import numpy as np
from collections import defaultdict
from scipy.interpolate import RegularGridInterpolator
import pkg_resources

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u

import beast.plotting.plot_compare_spec_type as plot_match


def compare_spec_type(
    phot_cat_file,
    beast_stats_file,
    spec_ra,
    spec_dec,
    spec_type,
    spec_subtype,
    lumin_class,
    match_radius=1.0,
    bright_filter="F475W",
    output_filebase=None,
    plot=False
):
    """
    BEAST verification: compare spectrally-typed stars to BEAST fits

    Parameters
    ----------
    phot_cat_file : string
        name of the file with the photometry

    beast_stats_file : string
        name of the file with the BEAST stats

    spec_ra, spec_dec : array of floats
        coordinates for the spectrally typed stars (in degrees)

    spec_type : list of strings
        spectral type for each star (O/B/A/F/G/K/M)

    spec_subtype : list of floats
        spectral subtype for each star (e.g., 7 for a B7 star, or 3.5 for a
        G3.5 star)

    lumin_class : list of strings
        luminosity class for each star (I/II/III/IV/V)

    match_radius : float (default=1)
        radius (arcsec) for searching for the BEAST star matching the spectrally
        typed star

    bright_filter : string (default='F475W')
        the star within match_radius that is brightest in this filter will be
        chosen as the match

    output_filebase : string (default=None)
        Prefix for saving the file of match info ("_spectype_match.fits" will be
        appended).  If None, will return the table rather than saving it.

    plot : boolean (default=False)
        Plot stars found to be a match.

    Returns
    -------
    spec_match : dict
        if output_filebase is None, a dictionary with match info is returned.

    """

    # read in the photometry catalog
    beast_phot = Table.read(phot_cat_file)
    ra_col = [x for x in beast_phot.colnames if "RA" in x.upper()][0]
    dec_col = [x for x in beast_phot.colnames if "DEC" in x.upper()][0]
    beast_phot_catalog = SkyCoord(
        ra=beast_phot[ra_col] * u.degree, dec=beast_phot[dec_col] * u.degree
    )
    # read in BEAST results
    beast_stats = Table.read(beast_stats_file, hdu=1)
    beast_stats_catalog = SkyCoord(
        ra=beast_stats["RA"] * u.degree, dec=beast_stats["DEC"] * u.degree
    )

    # setup the effective temperature table
    teff_function = setup_teff_table()
    # setup the surface gravity table
    logg_function = setup_logg_table()

    # dictionary to hold the matching results
    spec_match = defaultdict(list)

    for i in range(len(spec_ra)):

        # save star info
        spec_match["spec_ra"].append(spec_ra[i])
        spec_match["spec_dec"].append(spec_dec[i])
        spec_match["spec_type"].append(
            "{0} {1} {2}".format(spec_type[i], spec_subtype[i], lumin_class[i])
        )

        # Finding distances to every star in the catalog will take a while.
        # Therefore, first check to see if the spectrally typed star is even in
        # the field, by finding the distance for the closest star.  If it's a
        # small enough number, then proceed with calculating all separations.

        c = SkyCoord(spec_ra[i], spec_dec[i], unit="deg")
        _, sep_test, _ = c.match_to_catalog_sky(beast_stats_catalog)

        # print(sep_test[0].arcmin)

        # if the star is close enough, continue matching

        if sep_test[0].arcsec <= match_radius:

            # calculate separations for full catalog
            sep = c.separation(beast_stats_catalog)

            # find sources within user-defined radius
            small_sep_ind = np.where(sep.arcsec < match_radius)[0]

            # get photometry for those sources
            phot_list = np.zeros(len(small_sep_ind))
            phot_ind_list = np.zeros(len(small_sep_ind)).astype(int)

            for j in range(len(small_sep_ind)):
                phot_ind, _, _ = beast_stats_catalog[
                    small_sep_ind[j]
                ].match_to_catalog_sky(beast_phot_catalog)
                phot_ind_list[j] = phot_ind

                phot_ref_col = bright_filter + "_RATE"
                if phot_ref_col.lower() in beast_phot.colnames:
                    phot_list[j] = beast_phot[phot_ref_col.lower()][phot_ind]
                elif phot_ref_col.upper() in beast_phot.colnames:
                    phot_list[j] = beast_phot[phot_ref_col.upper()][phot_ind]
                else:
                    raise ValueError("{} not in catalog file".format(bright_filter))

            # find brightest match
            best_ind = small_sep_ind[phot_list == np.max(phot_list)][0]
            best_phot_ind = phot_ind_list[phot_list == np.max(phot_list)][0]

            # grab physical parameters for that source
            # - temperature
            teff_p16 = 10 ** beast_stats["logT_p16"][best_ind]
            teff_p50 = 10 ** beast_stats["logT_p50"][best_ind]
            teff_p84 = 10 ** beast_stats["logT_p84"][best_ind]
            # - log(g)
            logg_p16 = beast_stats["logg_p16"][best_ind]
            logg_p50 = beast_stats["logg_p50"][best_ind]
            logg_p84 = beast_stats["logg_p84"][best_ind]

            # grab physical parameters for the spectral type
            star_teff = lookup_param(
                teff_function, spec_type[i], spec_subtype[i], lumin_class[i]
            )
            star_logg = lookup_param(
                logg_function, spec_type[i], spec_subtype[i], lumin_class[i]
            )

            # calculate the number of sigmas away from the "true" value
            # - temperature
            if star_teff > teff_p50:
                teff_sigma = (star_teff - teff_p50) / (teff_p84 - teff_p50)
            else:
                teff_sigma = -(teff_p50 - star_teff) / (teff_p50 - teff_p16)
            # - log(g)
            if star_logg > logg_p50:
                logg_sigma = (star_logg - logg_p50) / (logg_p84 - logg_p50)
            else:
                logg_sigma = -(logg_p50 - star_logg) / (logg_p50 - logg_p16)

            # save the results
            spec_match["spec_teff"].append(star_teff)
            spec_match["spec_logg"].append(star_logg)
            spec_match["phot_cat_ind"].append(best_phot_ind)
            spec_match["stats_cat_ind"].append(best_ind)
            spec_match["beast_teff_p50"].append(teff_p50)
            spec_match["beast_teff_p16"].append(teff_p16)
            spec_match["beast_teff_p84"].append(teff_p84)
            spec_match["beast_logg_p50"].append(logg_p50)
            spec_match["beast_logg_p16"].append(logg_p16)
            spec_match["beast_logg_p84"].append(logg_p84)
            spec_match["teff_sigma"].append(teff_sigma)
            spec_match["logg_sigma"].append(logg_sigma)

        # if there's not a close enough match, put in placeholder values
        else:
            spec_match["spec_teff"].append(np.nan)
            spec_match["spec_logg"].append(np.nan)
            spec_match["phot_cat_ind"].append(np.nan)
            spec_match["stats_cat_ind"].append(np.nan)
            spec_match["beast_teff_p50"].append(np.nan)
            spec_match["beast_teff_p16"].append(np.nan)
            spec_match["beast_teff_p84"].append(np.nan)
            spec_match["beast_logg_p50"].append(np.nan)
            spec_match["beast_logg_p16"].append(np.nan)
            spec_match["beast_logg_p84"].append(np.nan)
            spec_match["teff_sigma"].append(np.nan)
            spec_match["logg_sigma"].append(np.nan)

    # write out the table
    if output_filebase is not None:
        spec_match_table = output_filebase + "_spectype_match.fits"
        Table(spec_match).write(spec_match_table, overwrite=True)

        if plot:
            plot_match.plot_compare_spec_type(spec_match_table, savefig='png')
    else:
        return dict(spec_match)


def lookup_param(input_function, spec_type, spec_subtype, lumin_class):
    """
    For given spectral classification, evaluate a function to find the
    physical parameter (e.g., T_eff, log(g))

    Note: the `.item()` at the end is because it would otherwise return, e.g.,
    `array(99.9)` instead of `99.9`.

    Parameters
    ----------
    input_function : RegularGridInterpolator object
        function to interpolate across a table

    spec_type : string
        spectral type for the star (O/B/A/F/G/K/M)

    spec_subtype : float
        spectral subtype for the star (e.g., 7 for a B7 star, or 3.5 for a
        G3.5 star)

    lumin_class : string
        luminosity class for the star (I/II/III/IV/V)

    Returns
    -------
    function_value : float
        the physical parameter for the star

    """

    return input_function(
        (_lumin_class_to_number(lumin_class), _spectype_to_id(spec_type, spec_subtype))
    ).item()


def setup_teff_table():
    """
    Read in and interpolate the T_eff table

    Ref: Table B.3 and B.4 in Stellar Spectral Classification

    Returns
    -------
    teff_interp_function : RegularGridInterpolator object
        function to interpolate across the T_eff table
    """

    # read in the table
    data_path = pkg_resources.resource_filename("beast", "tools/data/")
    teff_table = Table.read(data_path + "effective_temperature.txt", format="ascii")
    # make each row a number rather than a letter+number
    teff_table["row_id"] = np.zeros(len(teff_table))
    for i in range(len(teff_table)):
        teff_table["row_id"][i] = _spectype_to_id(
            teff_table["spec_type"][i], teff_table["spec_subtype"][i]
        )
    # interpolate the columns with missing data
    nans = np.isnan(teff_table["I"])
    teff_table["I"][nans] = np.interp(
        teff_table["row_id"][nans], teff_table["row_id"][~nans], teff_table["I"][~nans]
    )
    nans = np.isnan(teff_table["III"])
    teff_table["III"][nans] = np.interp(
        teff_table["row_id"][nans],
        teff_table["row_id"][~nans],
        teff_table["III"][~nans],
    )
    # now make a general interpolator
    teff_interp_function = RegularGridInterpolator(
        (np.array([1, 3, 5]), teff_table["row_id"]),
        np.array([teff_table["I"], teff_table["III"], teff_table["V"]]),
    )

    return teff_interp_function


def setup_logg_table():
    """
    Read in and interpolate the log(g) table

    Ref: Table 15.8 in Allen's Astrophysical Quantities

    Returns
    -------
    logg_interp_function : RegularGridInterpolator object
        function to interpolate across the log(g) table
    """

    # read in the table
    data_path = pkg_resources.resource_filename("beast", "tools/data/")
    logg_table = Table.read(data_path + "surface_gravity.txt", format="ascii")
    # make each row a number rather than a letter+number
    logg_table["row_id"] = np.zeros(len(logg_table))
    for i in range(len(logg_table)):
        logg_table["row_id"][i] = _spectype_to_id(
            logg_table["spec_type"][i], logg_table["spec_subtype"][i]
        )
    # interpolate the columns with missing data
    nans = np.isnan(logg_table["I"])
    logg_table["I"][nans] = np.interp(
        logg_table["row_id"][nans], logg_table["row_id"][~nans], logg_table["I"][~nans]
    )
    nans = np.isnan(logg_table["III"])
    logg_table["III"][nans] = np.interp(
        logg_table["row_id"][nans],
        logg_table["row_id"][~nans],
        logg_table["III"][~nans],
    )

    # convert from g/g_sun units to cm/s^2
    logg_sun = 27444  # cm/s^2
    for col in ["I", "III", "V"]:
        logg_table[col] = np.log10(10 ** logg_table[col] * logg_sun)

    # now make a general interpolator
    logg_interp_function = RegularGridInterpolator(
        (np.array([1, 3, 5]), logg_table["row_id"]),
        np.array([logg_table["I"], logg_table["III"], logg_table["V"]]),
    )

    return logg_interp_function


def _spectype_to_id(spec_type, spec_subtype):
    """
    Convert a spectral type to a numerical ID, so they can be interpolated
    """

    if (spec_subtype < 0) or (spec_subtype >= 10):
        raise ValueError("spec_subtype must be >=0 and <10")

    return ["O", "B", "A", "F", "G", "K", "M"].index(spec_type) * 10 + spec_subtype


def _lumin_class_to_number(lumin_class):
    """
    Convert the luminosity class from roman numeral to integer
    """

    return ["I", "II", "III", "IV", "V"].index(lumin_class) + 1
