from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


def plot_mag_hist(data_file, stars_per_bin=100, max_bins=75):
    """
    Make histograms of magnitudes.  This only uses the [filter]_VEGA, so
    sources removed with quality cuts are not included.

    Parameters
    ----------
    data_file : str
        path+file for the stellar photometry

    stars_per_bin : float (default=100)
        This is the average number of stars per histogram bin. Calculate
        the number of bins to use for each histogram by dividing the total
        number of stars by this value (this ensures each histogram is
        reasonably smooth).  The total number of bins is capped at max_bins.

    max_bins : int (default=75)
        maximum number of bins for each histogram.

    Returns
    -------
    peak_mags : dict
        dictionary with the peak magnitudes, for possible later use

    """

    # read in data
    with fits.open(data_file) as hdu:
        data_table = hdu[1].data
    filter_list = [col[:-5] for col in data_table.columns.names if 'VEGA' in col]
    n_filter = len(filter_list)

    # save the peak mags
    peak_mags = {}
    
    # figure
    fig = plt.figure(figsize=(5,4*n_filter))

    
    # make histograms
    for f,filt in enumerate(filter_list):

        # subplot region
        ax = plt.subplot(n_filter, 1, f+1)

        # histogram
        plot_this = data_table[filt+'_VEGA'][np.where(data_table[filt+'_VEGA'] < 90)]
        n_bins = np.min([ int( len(plot_this) / stars_per_bin ), max_bins ])
        
        hist = plt.hist(plot_this, bins=n_bins,
                        facecolor='grey', linewidth=0.25, edgecolor='grey')

        # peak magnitude
        peak_mags[filt] = hist[1][np.where(hist[0] == np.max(hist[0]))][0] + (hist[1][1]-hist[1][0])/2
        hist_ylim = ax.get_ylim()
        
        plt.plot([peak_mags[filt], peak_mags[filt]], [-100,1.2*hist_ylim[1]],
                    linestyle='--', linewidth=2, color='black', alpha=0.75)
        ax.set_ylim(hist_ylim)

        # label peak mag and total number of stars
        plt.text(0.65, 0.93, r'N$_{\mathrm{tot}}$: '+'{}'.format(len(plot_this)),
                     ha='left', va='center', transform=ax.transAxes, fontsize=12)
        plt.text(0.65, 0.85, r'M$_{\mathrm{peak}}$: '+'{:.2f}'.format(peak_mags[filt]),
                     ha='left', va='center', transform=ax.transAxes, fontsize=12)


        
        #pdb.set_trace()

        # axis labels
        ax.tick_params(axis='both', which='major', labelsize=13)
        ax.set_xlim(ax.get_xlim()[::-1])
        plt.xlabel(filt+' (Vega mag)', fontsize=14)
        plt.ylabel('N', fontsize=14)


    plt.tight_layout()

    fig.savefig(data_file.replace('.fits', '_maghist.pdf'))
    plt.close(fig)


    # return peak magnitudes
    return peak_mags
