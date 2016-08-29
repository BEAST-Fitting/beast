""" Script to plot the results of the single band noise model

.. history::
    Written 28Oct14 to allow the single band noise model generation
       from ASTs to be visually inspected.

"""
import numpy as np
import tables

from beast.core.noisemodel import toothpick
from beast.external.ezpipe.helpers import RequiredFile, task_decorator

class NGC4214_ToothPick_Noisemodel(toothpick.MultiFilterASTs):

    def set_data_mappings(self):
        """ hard code mapping directly with the interface to ASTs format

        .. note::

            As noted in the class documentation, it is trivial to adapt the
            class to specific formats
        """
        for k in self.filters:
            try:
                #self.data.set_alias(k, k.split('_')[-1].lower() + '_rate')
                self.data.set_alias(k + '_out', k.split('_')[-1].upper() + '_VEGA')
                self.data.set_alias(k + '_in', k.split('_')[-1].upper() + '_IN')
            except Exception as e:
                print(e)
                print('Warning: Mapping failed. This could lead to wrong results')

if __name__ == '__main__':

    astfile = 'data/11360_NGC-4214.gst.fake_new.fits'
    filters = ['HST_WFC3_F225W', 'HST_WFC3_F336W', 'HST_WFC3_F438W', 'HST_WFC3_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
    plot_filter_num = 2
    
    model = NGC4214_ToothPick_Noisemodel(astfile, filters)

    mag_in = model.data[filters[plot_filter_num] + '_in']
    mag_out = model.data[filters[plot_filter_num] + '_out']

    completeness_mag_cut = 80

    d = model._compute_stddev_mag(mag_in, mag_out, k=10, eps=0, completeness_mag_cut=completeness_mag_cut)

    d2 = model._compute_stddev_wrong(mag_in, mag_out, k=10, eps=0, completeness_mag_cut=completeness_mag_cut)

    #good_indxs,= np.where(mag_out < completeness_mag_cut)
    #flux_in = 10 ** (-0.4*mag_in[good_indxs])
    #flux_out = 10 ** (-0.4*mag_out[good_indxs])

    d3 = model._compute_stddev(mag_in, mag_out, k=10, eps=0, completeness_mag_cut=completeness_mag_cut)

    flux_out = model.data[filters[plot_filter_num].split('_')[-1].upper() + '_RATE']

    d4 = model._compute_stddev(mag_in, flux_out, k=10, eps=0, completeness_mag_cut=-1)

    import matplotlib.pyplot as pyplot

    fig, ax = pyplot.subplots(3, sharex=True, figsize=(10,13))

    ax[0].plot(d2['FLUX_IN'],d2['COMPLETENESS'],'ro',markerfacecolor='none',markeredgecolor='r')
    ax[0].plot(d['FLUX_IN'],d['COMPLETENESS'],'bo')
    ax[0].plot(d3['FLUX_IN'],d3['COMPLETENESS'],'go')
    ax[0].plot(d4['FLUX_IN'],d4['COMPLETENESS'],'co',markerfacecolor='none',markeredgecolor='c')
    ax[0].set_xscale('log')
    ax[0].set_yscale('linear')
    ax[0].set_ylabel('Completeness')
    ax[0].set_ylim(-0.1,1.1)

    ax[1].plot(d2['FLUX_IN'],d2['FLUX_BIAS']/d2['FLUX_IN'],'ro',markerfacecolor='none',markeredgecolor='r',label='old(mag)')
    ax[1].plot(d['FLUX_IN'],d['FLUX_BIAS']/d['FLUX_IN'],'bo',label='new(mag)')
    ax[1].plot(d3['FLUX_IN'],d3['FLUX_BIAS']/d3['FLUX_IN'],'go',label='newnew(flux)')
    ax[1].plot(d4['FLUX_IN'],d4['FLUX_BIAS']/d4['FLUX_IN'],'co',markerfacecolor='none',markeredgecolor='c',label='newrate(flux)')
    ax[1].plot([1e-13,1e-6],[0.0,0.0],'g--')
    ax[1].set_xscale('log')
    ax[1].set_yscale('linear')
    ax[1].set_ylabel('Bias/Flux(Input)')
    ax[1].set_ylim(-1.1,1.1)

    ax[2].plot(d2['FLUX_IN'],d2['FLUX_STD']/d2['FLUX_IN'],'ro',markerfacecolor='none',markeredgecolor='r')
    ax[2].plot(d['FLUX_IN'],d['FLUX_STD']/d['FLUX_IN'],'bo')
    ax[2].plot(d3['FLUX_IN'],d3['FLUX_STD']/d3['FLUX_IN'],'go')
    ax[2].plot(d4['FLUX_IN'],d4['FLUX_STD']/d4['FLUX_IN'],'co',markerfacecolor='none',markeredgecolor='c')
    ax[2].set_xscale('log')
    ax[2].set_yscale('log')
    ax[2].set_xlabel('Flux(Input) [Normalized to Vega]')
    ax[2].set_ylabel('Sigma/Flux(Input)')

    ax[1].legend()

    fig.subplots_adjust(hspace=0.05)

    pyplot.show()
