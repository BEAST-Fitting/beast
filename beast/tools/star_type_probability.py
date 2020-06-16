import numpy as np
from collections import defaultdict

from astropy.table import Table
from astropy.io import fits


def star_type_probability(
    pdf1d_files,
    pdf2d_files,
    output_filebase=None,
    ext_O_star_params=None,
    dusty_agb_params=None,
):
    """
    Calculate the probabilities of a set of star types by integrating either the
    1D or 2D PDF across the relevant range of parameter values.

    Currently does probabilities for these types.  See the docstrings of their
    respective functions for more details.  If the required parameters are not
    present in the PDFs, the function returns None.
    * extinguished O star (M_ini, Av)
    * dusty AGB star (Av, logT)

    Note that if more functions for different stellar types are added (please
    add more!), `params_to_save` needs to be updated.  This variable ensures
    that only the necessary parameters are stored in memory.

    Parameters
    ----------
    pdf1d_files : string or list of strings
        Name of the file(s) with the 1D PDFs.  If a list, it's each part of the
        subgrid.

    pdf2d_files : string or list of strings
        Name of the file(s) with the 2D PDFs.  If a list, it's in the same
        order as the subgrids above.

    output_filebase : string (default=None)
        Prefix for saving the file of probabilities ("_startype.fits" will be
        appended).  If None, will return the table rather than saving it.

    ext_O_star_params : dict or None
        Set to a dictionary to override the default cuts for extinguished early
        type stars.  Allowed keywords are 'min_M_ini' (float), 'min_Av' (float),
        and 'max_Av' (float).

    dusty_agb_params : dict or None
        Set to a dictionary to override the default cuts for dusty AGB stars
        (the high Av failure mode).  Allowed keywords are 'min_Av' (float),
        'min_logT' (float), and 'max_logT' (float).

    Returns
    -------
    star_prob : dict
        if output_filebase is None, a dictionary of probabilities is returned.

    """

    # read in the data

    # - set up dictionaries to hold PDFs and bins
    pdf1d_data = defaultdict(list)
    pdf1d_bins = defaultdict(list)
    pdf2d_data = defaultdict(list)
    pdf2d_bins = defaultdict(list)

    # - parameters to save
    params_to_save = ["Av", "M_ini", "logT"]

    # - go through each pair of files
    for (pdf1d_file, pdf2d_file) in zip(
        np.atleast_1d(pdf1d_files), np.atleast_1d(pdf2d_files)
    ):

        # 1D PDF data
        with fits.open(str(pdf1d_file)) as hdu:
            for ext in hdu:
                # only save the data if the parameter is in params_to_save
                if ext.name in params_to_save:
                    pdf1d_data[ext.name].append(ext.data[:-1, :])
                    pdf1d_bins[ext.name].append(ext.data[-1, :])

        # 2D PDF data
        with fits.open(str(pdf2d_file)) as hdu:
            for ext in hdu:
                # skip extensions without '+'
                if "+" not in ext.name:
                    continue

                # break up the name into the two parameters
                p1, p2 = ext.name.split("+")
                # only save the data if both parameters are in params_to_save
                if (p1 in params_to_save) and (p2 in params_to_save):
                    pdf2d_data[ext.name].append(ext.data[:-2, :, :])
                    pdf2d_bins[ext.name].append(ext.data[-2:, :, :])

    # combine arrays from each file

    for key in pdf1d_data:
        # check that the bins are the same for all
        bin_list = pdf1d_bins[key]
        bin_check = [
            not np.array_equal(bin_list[i], bin_list[i + 1])
            for i in range(len(bin_list) - 1)
        ]
        if np.sum(bin_check) > 0:
            raise ValueError("1D PDF bins not the same for each input file")
        # if so, just save the first one
        pdf1d_bins[key] = pdf1d_bins[key][0]
        # concatenate the PDFs
        pdf1d_data[key] = np.concatenate(pdf1d_data[key])

    for key in pdf2d_data:
        # check that the bins are the same for all
        bin_list = pdf2d_bins[key]
        bin_check = [
            not np.array_equal(bin_list[i], bin_list[i + 1])
            for i in range(len(bin_list) - 1)
        ]
        if np.sum(bin_check) > 0:
            raise ValueError("2D PDF bins not the same for each input file")
        # if so, just save the first one
        pdf2d_bins[key] = pdf2d_bins[key][0]
        # concatenate the PDFs
        pdf2d_data[key] = np.concatenate(pdf2d_data[key])

    # evaluate probabilities of things
    star_prob = {}

    # - extinguished O star
    if ext_O_star_params is None:
        star_prob["ext_O_star"] = ext_O_star(pdf2d_data, pdf2d_bins)
    else:
        star_prob["ext_O_star"] = ext_O_star(
            pdf2d_data, pdf2d_bins, **ext_O_star_params
        )

    # - dusty AGB star (high Av failure mode)
    if dusty_agb_params is None:
        star_prob["dusty_agb"] = dusty_agb(pdf2d_data, pdf2d_bins)
    else:
        star_prob["dusty_agb"] = dusty_agb(pdf2d_data, pdf2d_bins, **dusty_agb_params)

    # - other things

    # write out the table
    if output_filebase is not None:
        Table(star_prob).write(output_filebase + "_startype.fits", overwrite=True)
    else:
        return star_prob


def ext_O_star(pdf2d_data, pdf2d_bins, min_M_ini=10, min_Av=0.5, max_Av=99):
    """
    Calculate the probability that each star is an extinguished O star:
    * initial mass >= 10 Msun
    * A_V >= 0.5 mag
    There's a max A_V option to avoid possible high-Av artifacts.

    Some useful references for O/B stars
    https://ui.adsabs.harvard.edu/abs/2019A%26A...625A.104R/abstract
    https://ui.adsabs.harvard.edu/abs/2018A%26A...615A..40R/abstract
    https://ui.adsabs.harvard.edu/abs/2018A%26A...609A...7R/abstract

    Parameters
    ----------
    pdf2d_data : dict
        2D PDF data, each key has an array with shape (n_stars, nbin1, nbin2)

    pdf2d_bins : dict
        dictionary with corresponding bin values

    min_M_ini : float (default=10)
        minimum mass (in solar masses)

    min_Av : float (default=0.5)
        minimum Av (magnitudes)

    max_Av : float (default=99)
        maximum Av (magnitudes)

    Returns
    -------
    star_prob : array
        probability for each star

    """

    if "Av+M_ini" in pdf2d_data.keys():
        prob_data = pdf2d_data["Av+M_ini"]
        av_bins = pdf2d_bins["Av+M_ini"][0, :, :]
        mass_bins = pdf2d_bins["Av+M_ini"][1, :, :]
    elif "M_ini+Av" in pdf2d_data.keys():
        prob_data = pdf2d_data["M_ini+Av"]
        av_bins = pdf2d_bins["M_ini+Av"][1, :, :]
        mass_bins = pdf2d_bins["M_ini+Av"][0, :, :]
    else:
        print("2D PDFs don't contain M_ini and Av data")
        tot_stars = pdf2d_data[list(pdf2d_data)[0]].shape[0]
        return [np.nan] * tot_stars

    # reshape the arrays
    prob_data = prob_data.reshape(prob_data.shape[0], -1)
    av_bins = av_bins.reshape(-1)
    mass_bins = mass_bins.reshape(-1)

    keep = np.where(
        (mass_bins >= min_M_ini) & (av_bins >= min_Av) & (av_bins <= max_Av)
    )[0]

    return np.sum(prob_data[:, keep], axis=1)


def dusty_agb(pdf2d_data, pdf2d_bins, min_Av=7, min_logT=3.7, max_logT=4.2):
    """
    Calculate the probability that each star is a dusty AGB star, using the high
    Av failure mode:
    * A_V >= 7 mag
    * Log T_eff from 3.7 to 4.2

    Parameters
    ----------
    pdf2d_data : dict
        2D PDF data, each key has an array with shape (n_stars, nbin1, nbin2)

    pdf2d_bins : dict
        dictionary with corresponding bin values

    min_Av : float (default=0.5)
        minimum Av (magnitudes)

    min_logT, max_logT : float (default=3.7, 4.2)
        minimum and maximum logT


    Returns
    -------
    star_prob : array
        probability for each star

    """

    if "Av+logT" in pdf2d_data.keys():
        prob_data = pdf2d_data["Av+logT"]
        av_bins = pdf2d_bins["Av+logT"][0, :, :]
        logT_bins = pdf2d_bins["Av+logT"][1, :, :]
    elif "logT+Av" in pdf2d_data.keys():
        prob_data = pdf2d_data["logT+Av"]
        av_bins = pdf2d_bins["logT+Av"][1, :, :]
        logT_bins = pdf2d_bins["logT+Av"][0, :, :]
    else:
        print("2D PDFs don't contain Av and logT (T_eff) data")
        tot_stars = pdf2d_data[list(pdf2d_data)[0]].shape[0]
        return [np.nan] * tot_stars

    # reshape the arrays
    prob_data = prob_data.reshape(prob_data.shape[0], -1)
    av_bins = av_bins.reshape(-1)
    logT_bins = logT_bins.reshape(-1)

    keep = np.where(
        (av_bins >= min_Av) & (logT_bins >= min_logT) & (logT_bins <= max_logT)
    )[0]

    return np.sum(prob_data[:, keep], axis=1)
