import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt

from beast.tools import read_beast_data

from astropy.table import Table
from astropy.io import fits


def star_type_probability(
    pdf1d_files,
    pdf2d_files,
    output_filebase,
):
    """
    Parameters
    ----------
    pdf1d_files : string or list of strings
        Name of the file(s) with the 1D PDFs.  If a list, it's each part of the
        subgrid.

    pdf2d_files : string or list of strings
        Name of the file(s) with the 2D PDFs.  If a list, it's in the same
        order as the subgrids above.

    output_filebase : string
        prefix for saving the file of probabilities ("_startype.fits" will be
        appended)

    """

    # read in the data

    # - set up dictionaries to hold PDFs and bins
    pdf1d_data = defaultdict(list)
    pdf1d_bins = defaultdict(list)
    pdf2d_data = defaultdict(list)
    pdf2d_bins = defaultdict(list)

    # - parameters to save
    param_list = ['Av', 'M_ini', 'logA']

    # - go through each pair of files
    for (pdf1d_file, pdf2d_file) in zip(np.atleast_1d(pdf1d_files), np.atleast_1d(pdf2d_files)):

        # 1D PDF data
        with fits.open(str(pdf1d_file)) as hdu:
            for ext in hdu:
                # only save the data if the parameter is in param_list
                if ext.name in param_list:
                    pdf1d_data[ext.name].append(ext.data[:-1,:])
                    pdf1d_bins[ext.name].append(ext.data[-1,:])

        # 2D PDF data
        with fits.open(str(pdf2d_file)) as hdu:
            for ext in hdu:
                # skip extensions without '+'
                if '+' not in ext.name:
                    continue

                # break up the name into the two parameters
                p1, p2 = ext.name.split('+')
                # only save the data if both parameters are in param_list
                if (p1 in param_list) and (p2 in param_list):
                    pdf2d_data[ext.name].append(ext.data[:-2,:,:])
                    pdf2d_bins[ext.name].append(ext.data[-2:,:,:])


    # combine arrays from each file

    for key in pdf1d_data:
        # check that the bins are the same for all
        bin_list = pdf1d_bins[key]
        bin_check = [
            not np.array_equal(bin_list[i], bin_list[i+1])
            for i in range(len(bin_list)-1)
        ]
        if np.sum(bin_check) > 0:
            raise ValueError('1D PDF bins not the same for each input file')
        # if so, just save the first one
        pdf1d_bins[key] = pdf1d_bins[key][0]
        # concatenate the PDFs
        pdf1d_data[key] = np.concatenate(pdf1d_data[key])

    for key in pdf2d_data:
        # check that the bins are the same for all
        bin_list = pdf2d_bins[key]
        bin_check = [
            not np.array_equal(bin_list[i], bin_list[i+1])
            for i in range(len(bin_list)-1)
        ]
        if np.sum(bin_check) > 0:
            raise ValueError('2D PDF bins not the same for each input file')
        # if so, just save the first one
        pdf2d_bins[key] = pdf2d_bins[key][0]
        # concatenate the PDFs
        pdf2d_data[key] = np.concatenate(pdf2d_data[key])



    # evaluate probabilities of things
    star_prob = {}

    # - extinguished O star
    star_prob['ext_O_star'] = ext_O_star(pdf2d_data, pdf2d_bins)

    # - dusty AGB star (high Av failure mode)
    star_prob['dusty_agb'] = dusty_agb(pdf2d_data, pdf2d_bins)

    # - other things


    # write out the table
    Table(star_prob).write(output_filebase+"_startype.fits", overwrite=True)



def ext_O_star(pdf2d_data, pdf2d_bins):
    """
    Calculate the probability that each star is an extinguished O star:
    * initial mass >= 10 Msun
    * A_V >= 0.5 mag

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

    Returns
    -------
    star_prob : array
        probability for each star

    """

    if 'Av+M_ini' in pdf2d_data.keys():
        prob_data = pdf2d_data['Av+M_ini']
        av_bins = pdf2d_bins['Av+M_ini'][0,:,:]
        mass_bins = pdf2d_bins['Av+M_ini'][1,:,:]
    elif 'M_ini+Av' in pdf2d_data.keys():
        prob_data = pdf2d_data['M_ini+Av']
        av_bins = pdf2d_bins['M_ini+Av'][1,:,:]
        mass_bins = pdf2d_bins['M_ini+Av'][0,:,:]
    else:
        raise ValueError("2D PDFs don't contain M_ini and Av data")

    # reshape the arrays
    prob_data = prob_data.reshape(prob_data.shape[0], -1)
    av_bins = av_bins.reshape(-1)
    mass_bins = mass_bins.reshape(-1)

    keep = np.where((mass_bins > 10) & (av_bins > 0.5))[0]

    return np.sum(prob_data[:,keep], axis=1)


def dusty_agb(pdf2d_data, pdf2d_bins):
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

    Returns
    -------
    star_prob : array
        probability for each star

    """

    if 'Av+logT' in pdf2d_data.keys():
        prob_data = pdf2d_data['Av+logT']
        av_bins = pdf2d_bins['Av+logT'][0,:,:]
        logT_bins = pdf2d_bins['Av+logT'][1,:,:]
    elif 'logT+Av' in pdf2d_data.keys():
        prob_data = pdf2d_data['logT+Av']
        av_bins = pdf2d_bins['logT+Av'][1,:,:]
        logT_bins = pdf2d_bins['logT+Av'][0,:,:]
    else:
        raise ValueError("2D PDFs don't contain Av and logT (T_eff) data")

    # reshape the arrays
    prob_data = prob_data.reshape(prob_data.shape[0], -1)
    av_bins = av_bins.reshape(-1)
    logT_bins = logT_bins.reshape(-1)

    keep = np.where((av_bins >= 7) & (logT_bins >= 3.7) & (logT_bins <= 4.2))[0]

    return np.sum(prob_data[:,keep], axis=1)
