import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import copy
from scipy.stats import binned_statistic_2d
from astropy.io import fits
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
        model_compl = noisegrid.root.completeness[:]
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
                compl_data = make_image(compl_table, pi, pj)

                # plot points
                im = plt.imshow(
                    compl_data['compl_image'].T,
                    #np.random.random((4,4)),
                    extent=(
                        np.min(compl_data['x_bins']), np.max(compl_data['x_bins']),
                        np.min(compl_data['y_bins']), np.max(compl_data['y_bins'])
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
                    ax.set_ylabel(compl_data['y_label'], fontsize=label_font)
                    #ax.get_yaxis().set_label_coords(-0.35,0.5)
                else:
                    ax.set_yticklabels([])
                if j == n_params-1:
                    ax.set_xlabel(compl_data['x_label'], fontsize=label_font)
                    plt.xticks(rotation=-45)
                else:
                    ax.set_xticklabels([])


            # along diagonal
            if i == j:

                # set up subplot
                plt.subplot(n_params, n_params, i + j*(n_params) + 1)
                ax = plt.gca()

                # make histogram
                plt.hist(np.random.randint(0,10,100), bins=10,
                         facecolor='grey', linewidth=0.25, edgecolor='grey')

                ax.tick_params(axis='y',which='both',length=0, labelsize=tick_font)
                ax.tick_params(axis='x',which='both',direction='in', labelsize=tick_font)

                # axis labels and ticks
                ax.set_yticklabels([])
                if i < n_params-1:
                    ax.set_xticklabels([])
                if i == n_params-1:
                    ax.set_xlabel(compl_data['x_label'], fontsize=label_font)
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



def make_image(compl_table, x_param, y_param):
    """
    Make the completeness image

    compl_table : astropy table
        table with each set of physical parameters and their completeness

    x_param, y_param : strings
        parameters for which to create a completeness image
    """

    # ---- set up things for x
    # mass isn't reguarly spaced, so take log and manually define bins
    if 'M_' in x_param:
        x_col = np.log(compl_table[x_param])
        x_bins = np.linspace(np.min(x_col), np.max(x_col), 20)
        x_label = 'log '+x_param
    # metallicity just needs to be log
    elif x_param == 'Z':
        x_col = np.log(compl_table[x_param])
        x_bins = np.linspace(np.min(x_col), np.max(x_col), len(np.unique(x_col))+1)
        x_label = 'log '+x_param
    # for all others, standard linear spacing is ok
    else:
        x_col = copy.copy(compl_table[x_param])
        x_bins = np.linspace(np.min(x_col), np.max(x_col), len(np.unique(x_col))+1)
        x_label = copy.copy(x_param)

    # ---- set up things for y
    # mass isn't reguarly spaced, so take log and manually define bins
    if 'M_' in y_param:
        y_col = np.log(compl_table[y_param])
        y_bins = np.linspace(np.min(y_col), np.max(y_col), 20)
        y_label = 'log '+y_param
    # metallicity just needs to be log
    elif y_param == 'Z':
        y_col = np.log(compl_table[y_param])
        y_bins = np.linspace(np.min(y_col), np.max(y_col), len(np.unique(y_col))+1)
        y_label = 'log '+y_param
    # for all others, standard linear spacing is ok
    else:
        y_col = copy.copy(compl_table[y_param])
        y_bins = np.linspace(np.min(y_col), np.max(y_col), len(np.unique(y_col))+1)
        y_label = copy.copy(y_param)


    # ---- make the image

    # compl_image = np.zeros((len(x_bins)-1, len(y_bins)-1))
    #
    # for i in range(len(x_bins)-1):
    #     for j in range(len(y_bins)-1):
    #
    #         # find the indices in this bin
    #         match = np.where(
    #             (x_col >= x_bins[i]) & (x_col <= x_bins[i+1]) &
    #             (y_col >= y_bins[j]) & (y_col <= y_bins[j+1])
    #         )
    #         # calculate average completeness for those indices
    #         compl_image[i,j] = np.average(compl_table['compl'][match])

    compl_image, _, _, _ = binned_statistic_2d(
        x_col,
        y_col,
        compl_table['compl'],
        statistic='mean',
        bins=(x_bins, y_bins),
    )

    return {
        'compl_image':compl_image,
        'x_bins':x_bins,
        'x_label':x_label,
        'y_bins':y_bins,
        'y_label':y_label,
    }
