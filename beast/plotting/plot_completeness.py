import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import copy
from scipy.stats import binned_statistic, binned_statistic_2d
from astropy.table import Table, vstack

from beast.physicsmodel.grid import FileSEDGrid
import beast.observationmodel.noisemodel.generic_noisemodel as noisemodel


def plot_completeness(
    physgrid_list,
    noise_model_list,
    output_plot_filename,
    param_list=['Av', 'Rv', 'logA', 'f_A', 'M_ini', 'Z', 'distance'],
    compl_filter='F475W',
):
    """
    Make visualization of the completeness
    Parameters
    ----------
    physgrid_list : string or list of strings
        Name of the physics model file.  If there are multiple physics model
        grids (i.e., if there are subgrids), list them all here.

    noise_model_list : string or list of strings
        Name of the noise model file.  If there are multiple files for
        physgrid_list (because of subgrids), list the noise model file
        associated with each physics model file.

    param_list : list of strings
        names of the parameters to plot

    compl_filter : str
        filter to use for completeness (required for toothpick model)

    output_plot_filename : string
        name of the file in which to save the output plot

    """

    n_params = len(param_list)

    # If there are subgrids, we can't read them all into memory.  Therefore,
    # we'll go through each one and just grab the relevant parts.
    compl_table_list = []

    # make a table for each physics model + noise model
    for physgrid, noise_model in zip(
        np.atleast_1d(physgrid_list), np.atleast_1d(noise_model_list)
    ):

        # get the physics model grid - includes priors
        modelsedgrid = FileSEDGrid(str(physgrid))
        # get list of filters
        short_filters = [filter.split(sep="_")[-1].upper() for filter in modelsedgrid.filters]
        if compl_filter.upper() not in short_filters:
            raise ValueError("requested completeness filter not present")
        filter_k = short_filters.index(compl_filter.upper())
        print("Completeness from {0}".format(modelsedgrid.filters[filter_k]))


        # read in the noise model
        noisegrid = noisemodel.get_noisemodelcat(str(noise_model))
        # get the completeness
        model_compl = noisegrid["completeness"]
        # close the file to save memory
        noisegrid.close()

        # put it all into a table
        table_dict = {x:modelsedgrid[x] for x in param_list}
        table_dict['compl'] = model_compl[:, filter_k]

        # append to the list
        compl_table_list.append(Table(table_dict))

    # stack all the tables into one
    compl_table = vstack(compl_table_list)

    #import pdb; pdb.set_trace()

    # figure
    fig = plt.figure(figsize=(4*n_params, 4*n_params))

    # label font sizes
    label_font = 25
    tick_font = 22

    # load in color map
    cmap = matplotlib.cm.get_cmap('magma')

    # iterate through the panels
    for i,pi in enumerate(param_list):
        for j,pj in enumerate(param_list[i:], i):

            print('plotting {0} and {1}'.format(pi, pj))

            # not along diagonal
            if i != j:

                # set up subplot
                plt.subplot(n_params, n_params, i + j*(n_params) + 1)
                ax = plt.gca()

                # create image and labels
                x_col, x_bins, x_label = setup_axis(compl_table, pi)
                y_col, y_bins, y_label = setup_axis(compl_table, pj)
                compl_image, _, _, _ = binned_statistic_2d(
                    x_col,
                    y_col,
                    compl_table['compl'],
                    statistic='mean',
                    bins=(x_bins, y_bins),
                )

                # plot points
                im = plt.imshow(
                    compl_image.T,
                    #np.random.random((4,4)),
                    extent=(
                        np.min(x_bins), np.max(x_bins),
                        np.min(y_bins), np.max(y_bins)
                    ),
                    cmap='magma',
                    vmin=0,
                    vmax=1,
                    aspect='auto',
                    origin='lower',
                )

                ax.tick_params(axis='both', which='both', direction='in', labelsize=tick_font,
                                   bottom=True, top=True, left=True, right=True)

                # axis labels and ticks
                if i == 0:
                    ax.set_ylabel(y_label, fontsize=label_font)
                    #ax.get_yaxis().set_label_coords(-0.35,0.5)
                else:
                    ax.set_yticklabels([])
                if j == n_params-1:
                    ax.set_xlabel(x_label, fontsize=label_font)
                    plt.xticks(rotation=-45)
                else:
                    ax.set_xticklabels([])


            # along diagonal
            if i == j:

                # set up subplot
                plt.subplot(n_params, n_params, i + j*(n_params) + 1)
                ax = plt.gca()

                # create histogram and labels
                x_col, x_bins, x_label = setup_axis(compl_table, pi)
                compl_hist, _, _ = binned_statistic(
                    x_col,
                    compl_table['compl'],
                    statistic='mean',
                    bins=x_bins,
                )
                # make histogram
                _, _, patches = plt.hist(x_bins[:-1], x_bins, weights=compl_hist)
                # color each bar by its completeness
                for c,comp in enumerate(compl_hist):
                    patches[c].set_color(cmap(comp))
                    patches[c].set_linewidth=0.1
                # make a black outline so it stands out as a histogram
                plt.hist(x_bins[:-1], x_bins, weights=compl_hist, histtype='step', color='k')
                # axis ranges
                plt.xlim(np.min(x_bins), np.max(x_bins))
                plt.ylim(0, 1.05)

                ax.tick_params(axis='y',which='both',length=0, labelsize=tick_font)
                ax.tick_params(axis='x',which='both',direction='in', labelsize=tick_font)

                # axis labels and ticks
                ax.set_yticklabels([])
                if i < n_params-1:
                    ax.set_xticklabels([])
                if i == n_params-1:
                    ax.set_xlabel(x_label, fontsize=label_font)
                    plt.xticks(rotation=-45)

    #plt.subplots_adjust(wspace=0.05, hspace=0.05)
    plt.tight_layout()

    # add a colorbar
    gs = GridSpec(nrows=20, ncols=n_params)
    cax = fig.add_subplot(gs[0, 2:])
    cbar = plt.colorbar(im, cax=cax, orientation='horizontal')
    cbar.set_label('Completeness', fontsize=label_font)
    cbar.ax.tick_params(labelsize=tick_font)
    gs.tight_layout(fig)

    fig.savefig(output_plot_filename)
    plt.close(fig)


def setup_axis(compl_table, param):
    """
    Set up the bins and labels for a parameter

    Parameters
    ----------
    compl_table : astropy table
        table with each set of physical parameters and their completeness

    param : string
        name of the parameter we're binning/labeling

    Returns
    -------
    col : numpy array
        column to plot

    bins : numpy array
        bin edges

    label : string
        the axis label to use

    """

    # mass isn't reguarly spaced, so take log and manually define bins
    if 'M_' in param:
        col = np.log10(compl_table[param])
        bins = np.linspace(np.min(col), np.max(col), 20)
        label = 'log '+param
    # metallicity just needs to be log
    elif param == 'Z':
        col = np.log10(compl_table[param])
        bins = np.linspace(np.min(col), np.max(col), len(np.unique(col))+1)
        label = 'log '+param
    # for all others, standard linear spacing is ok
    else:
        col = copy.copy(compl_table[param])
        bins = np.linspace(np.min(col), np.max(col), len(np.unique(col))+1)
        label = copy.copy(param)

    return col, bins, label
