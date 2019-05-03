from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


def plot(beast_stats_file, n_bins=50):
    """
    Make a histogram of the chi2 values
    
    Parameters
    ----------
    beast_stats_file : str
        path+file for the BEAST fitting results

    n_bins : int (default=50)
        number of bins to use in the histogram

    """

    # read in stats file
    with fits.open(beast_stats_file) as beast_hdu:
        beast_table = beast_hdu[1].data

    # grab some info
    # - chi2
    chi2 = beast_table['chi2min']
    # - filters
    #filter_list = list(set( [re.search('[F]...[W]',col).group(0) for col in beast_table.columns.names if
    #                             re.search('[F]...[W]',col) is not None] ))
    #n_filter = len(filter_list)
    # - fit parameters
    #param_list = [col[:-4] for col in beast_table.columns.names if 'p50' in col and '_F' not in col]
    #n_param = 6

    # figure
    fig = plt.figure(figsize=(5,4))
    ax = plt.gca()

    # make histogram
    plt.hist(np.log10(chi2), bins=n_bins,
                 facecolor='grey', linewidth=0.25, edgecolor='grey')
    hist_ylim = ax.get_ylim()
    
    # reduced chi2
    # ... well, that's not so straightforward

    # mark chi2=1
    plt.plot([0,0], [-100,1.2*hist_ylim[1]],
                 linestyle=':', linewidth=2, color='black', alpha=0.75)

    # mark median chi2
    plt.plot([np.log10(np.median(chi2)), np.log10(np.median(chi2))], [-100,1.2*hist_ylim[1]],
                 linestyle='--', linewidth=2, color='black', alpha=0.75)
    plt.text(np.log10(np.median(chi2)) + 0.2, 0.93*hist_ylim[1],
                          r'$\chi^2_{\mathrm{p50}}$ = ' + '{:.2f}'.format(np.median(chi2)),
                          ha='left', va='center')
    
    

    ax.set_xlabel(r'Log $\chi^2$')
    ax.set_ylim(hist_ylim)

    plt.tight_layout()

    # save figure
    fig.savefig(beast_stats_file.replace('.fits', '_chi2hist.pdf'))

    plt.close('all')
