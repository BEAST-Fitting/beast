from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


def plot_mag_hist(data_file, n_bins=50):
    """
    Make histograms of magnitudes.  This only uses the [filter]_VEGA, so
    sources removed with quality cuts are not included.

    Parameters
    ----------
    data_file : str
        path+file for the stellar photometry

    n_bins : int (default=50)
        number of bins to use in the histogram

    """

    # read in data
    with fits.open(data_file) as hdu:
        data_table = hdu[1].data
    filter_list = [col[:-5] for col in data_table.columns.names if 'VEGA' in col]
    n_filter = len(filter_list)

    
    # figure
    fig = plt.figure(figsize=(5,4*n_filter))

    
    # make histograms
    for f,filt in enumerate(filter_list):

        # subplot region
        ax = plt.subplot(n_filter, 1, f+1)

        # histogram
        plot_this = data_table[filt+'_VEGA'][np.where(data_table[filt+'_VEGA'] < 90)]
        
        hist = plt.hist(plot_this, bins=n_bins,
                        facecolor='grey', linewidth=0.25, edgecolor='grey')

        # peak magnitude
        mag_peak = hist[1][np.where(hist[0] == np.max(hist[0]))][0] + (hist[1][1]-hist[1][0])/2
        hist_ylim = ax.get_ylim()
        
        plt.plot([mag_peak, mag_peak], [-100,1.2*hist_ylim[1]],
                    linestyle='--', linewidth=2, color='black', alpha=0.75)
        ax.set_ylim(hist_ylim)

        # label peak mag and total number of stars
        plt.text(0.65, 0.93, r'N$_{\mathrm{tot}}$: '+'{}'.format(len(plot_this)),
                     ha='left', va='center', transform=ax.transAxes, fontsize=12)
        plt.text(0.65, 0.85, r'M$_{\mathrm{peak}}$: '+'{:.2f}'.format(mag_peak),
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

