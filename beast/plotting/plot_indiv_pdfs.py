import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from astropy.io import fits
import copy

__all__ = ["plot_indiv_pdfs"]


def plot_indiv_pdfs(pdf1d_file, pdf2d_file, starnum):
    """
    Make a triangle/corner plot with the 1D and 2D PDFs of a star

    Parameters
    ----------
    pdf1d_file : str
        path+file for the BEAST 1D PDFs

    pdf2d_file : list of strings
        path+file for the BEAST 2D PDFs

    starnum : int
        the index of the star to plot

    """

    with fits.open(pdf2d_file) as hdu_2d, fits.open(pdf1d_file) as hdu_1d:

        # start with the 2D PDF file to figure out which parameters to plot

        # list of extension names
        # - skip the extension named "PRIMARY"
        # - skip any extension with a dimension of 1 (means that param wasn't fit)
        ext_list = [
            hdu_2d[i].name
            for i in range(len(hdu_2d))
            if ("PRIMARY" not in hdu_2d[i].name) and (1 not in hdu_2d[i].data.shape)
        ]

        # grab the parameter names
        param_list_temp = [i for x in ext_list for i in x.split("+")]

        # put them in the order we want
        _plot_order = ["M_ini", "logA", "distance", "Z", "Av", "Rv", "f_A"]
        param_list = [p for p in _plot_order if p in param_list_temp]
        # if there are parameters not in the predetermined plot order, append them
        for p in param_list_temp:
            if p not in param_list:
                param_list.append(p)
        # total number of parameters
        n_params = len(param_list)

        # figure
        fig = plt.figure(figsize=(4 * n_params, 4 * n_params))

        # label font sizes
        label_font = 25
        tick_font = 22

        # iterate through the panels
        for i, pi in enumerate(param_list):
            for j, pj in enumerate(param_list[i:], i):

                # print('plotting {0} and {1}'.format(pi, pj))

                # not along diagonal
                if i != j:

                    # set up subplot
                    plt.subplot(n_params, n_params, i + j * (n_params) + 1)
                    ax = plt.gca()

                    # find the 2D PDF and make sure it's properly rotated
                    try:
                        image = hdu_2d[pi + "+" + pj].data[starnum, :, :].T
                    except KeyError:
                        image = hdu_2d[pj + "+" + pi].data[starnum, :, :]
                    except Exception:
                        raise

                    # create axis/labels
                    x_bins, x_label = setup_axis(pi, hdu_1d)
                    y_bins, y_label = setup_axis(pj, hdu_1d)

                    # plot 2D PDF image
                    im = plt.imshow(
                        np.log(image / np.max(image)),
                        extent=(
                            np.min(x_bins),
                            np.max(x_bins),
                            np.min(y_bins),
                            np.max(y_bins),
                        ),
                        cmap="magma",
                        vmin=-10,
                        vmax=0,
                        aspect="auto",
                        origin="lower",
                    )

                    # attempt to plot 1/2/3 sigma contours
                    # (doesn't work if the probability is super concentrated,
                    # which is generally due to the grid being really coarse)
                    try:
                        im_sort = np.sort(image, axis=None)[::-1]
                        cumsum = np.cumsum(im_sort)
                        cumsum /= np.max(cumsum)
                        clevels = [
                            np.log(im_sort[np.where(cumsum <= p)[0][-1]])
                            for p in [0.68, 0.95, 0.997]
                        ]
                        plt.contour(
                            x_bins,
                            y_bins,
                            image,
                            levels=clevels,
                            colors="k",
                            linestyles="-",
                        )
                    except IndexError:
                        # print("  can't make contours for this")
                        pass
                    except ValueError:
                        pass
                    except Exception:
                        raise

                    ax.tick_params(
                        axis="both",
                        which="both",
                        direction="in",
                        labelsize=tick_font,
                        bottom=True,
                        top=True,
                        left=True,
                        right=True,
                    )

                    # axis labels and ticks
                    if i == 0:
                        ax.set_ylabel(y_label, fontsize=label_font)
                        # ax.get_yaxis().set_label_coords(-0.35,0.5)
                    else:
                        ax.set_yticklabels([])
                    if j == n_params - 1:
                        ax.set_xlabel(x_label, fontsize=label_font)
                        plt.xticks(rotation=-45)
                    else:
                        ax.set_xticklabels([])

                # along diagonal
                if i == j:

                    # set up subplot
                    plt.subplot(n_params, n_params, i + j * (n_params) + 1)
                    ax = plt.gca()

                    # create axis/labels
                    x_bins, x_label = setup_axis(pi, hdu_1d)

                    # make histogram
                    _pdf = hdu_1d[pi].data[starnum, :]
                    plt.plot(
                        x_bins,
                        _pdf / np.max(_pdf),
                        marker="o",
                        mew=0,
                        color="black",
                        markersize=2,
                        linestyle="-",
                    )

                    # axis ranges
                    plt.xlim(np.min(x_bins), np.max(x_bins))
                    plt.ylim(0, 1.05)

                    ax.tick_params(
                        axis="y", which="both", length=0, labelsize=tick_font
                    )
                    ax.tick_params(
                        axis="x", which="both", direction="in", labelsize=tick_font
                    )

                    # axis labels and ticks
                    ax.set_yticklabels([])
                    if i < n_params - 1:
                        ax.set_xticklabels([])
                    if i == n_params - 1:
                        ax.set_xlabel(x_label, fontsize=label_font)
                        plt.xticks(rotation=-45)

    # plt.subplots_adjust(wspace=0.05, hspace=0.05)
    plt.tight_layout()

    # add a colorbar
    gs = GridSpec(nrows=20, ncols=n_params)
    cax = fig.add_subplot(gs[0, 2:])
    cbar = plt.colorbar(im, cax=cax, orientation="horizontal")
    cbar.set_label("Log Likelihood", fontsize=label_font)
    cbar.ax.tick_params(labelsize=tick_font)
    gs.tight_layout(fig)

    fig.savefig(
        pdf1d_file.replace("pdf1d.fits", "pdfs_starnum_{0}.pdf".format(starnum))
    )
    plt.close(fig)


def setup_axis(param, pdf1d_hdu):
    """
    Set up the bins and labels for a parameter

    Parameters
    ----------
    param : string
        name of the parameter we're binning/labeling

    pdf1d_hdu : HDU object
        the HDU for param (to get the bins)

    Returns
    -------
    bins : numpy array
        bin edges

    label : string
        the axis label to use

    """

    # mass and metallicity get log spaced
    if ("M_" in param) or (param == "Z"):
        bins = np.log10(pdf1d_hdu[param].data[-1, :])
        label = "log " + param
    # for all others, standard linear spacing is ok
    else:
        bins = pdf1d_hdu[param].data[-1, :]
        label = copy.copy(param)

    return bins, label
