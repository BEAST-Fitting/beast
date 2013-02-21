import numpy as np
import pylab as plt
from sedfitter import eztables


def plot_keys(keys, outdir, outnames, show=True):
    """ Plot average offset per unique values in keys
    INPUTS:
        keys        list    columns to process
        outdir      str     directory where to find the summary fits tables
        outnames    list    list of names for tables

    TODO: replace inputs by a list of files and list of keys
    """
    for i in range(len(outnames)):
        summary_table = eztables.Table(outdir + 'summary_' + outnames[i] + '.fits')
        for j, key in enumerate(keys):
            #rec_vals = summary_table.data[key+'_recovered']  # recovered values
            true_vals = summary_table.data[key]              # true values
            #rec_vals = rec_vals - true_vals                  # offset
            rec_vals = summary_table.evalexpr('{}_recovered - {}'.format(key, key))  # offset
            uniq_vals = np.unique(true_vals)  # unique true values
            avg_offs = np.zeros(uniq_vals.size)  # Mean of recovered params for given input param
            for k in range(avg_offs.size):
                sel = np.where(true_vals == uniq_vals[k])
                avg_offs[k] = rec_vals[sel].mean()
            #figure
            plt.figure(j)
            ax = plt.subplot(111)
            ax.plot(uniq_vals, avg_offs, label=(outnames[i]).replace('_', ' '), marker='+', markersize=10)

    for j, key in enumerate(keys):
        plt.figure(j)
        ax = plt.subplot(111)
        ax.set_xlabel(key.replace('_', ''))
        ax.set_ylabel(key.replace('_', '') + ' (out - in)')
        ax.legend(loc='best')
        ax.plot(ax.get_xlim(), [0, 0], linestyle='dotted', color='black')
    if show:
        plt.show()
