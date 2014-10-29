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

    d = model._compute_stddev(mag_in, mag_out, k=10, eps=0, completeness_mag_cut=80)

    d2 = model._compute_stddev_wrong(mag_in, mag_out, k=10, eps=0, completeness_mag_cut=80)

    import matplotlib.pyplot as pyplot

    fig, ax = pyplot.subplots(3, sharex=True, figsize=(10,13))

    ax[0].plot(d2['FLUX_IN'],d2['COMPLETENESS'],'ro',markerfacecolor='none',markeredgecolor='r')
    ax[0].plot(d['FLUX_IN'],d['COMPLETENESS'],'bo')
    ax[0].set_xscale('log')
    ax[0].set_yscale('linear')
    ax[0].set_ylabel('Completeness')
    ax[0].set_ylim(-0.1,1.1)

    ax[1].plot(d2['FLUX_IN'],d2['FLUX_BIAS']/d2['FLUX_IN'],'ro',markerfacecolor='none',markeredgecolor='r',label='old')
    ax[1].plot(d['FLUX_IN'],d['FLUX_BIAS']/d['FLUX_IN'],'bo',label='new')
    ax[1].plot([1e-13,1e-6],[0.0,0.0],'g--')
    ax[1].set_xscale('log')
    ax[1].set_yscale('linear')
    ax[1].set_ylabel('Bias/Flux(Input)')
    ax[1].set_ylim(-1.1,1.1)
    ax[1].legend()

    ax[2].plot(d2['FLUX_IN'],d2['FLUX_STD']/d2['FLUX_IN'],'ro',markerfacecolor='none',markeredgecolor='r')
    ax[2].plot(d['FLUX_IN'],d['FLUX_STD']/d['FLUX_IN'],'bo')
    ax[2].set_xscale('log')
    ax[2].set_yscale('log')
    ax[2].set_xlabel('Flux(Input)')
    ax[2].set_ylabel('Sigma/Flux(Input)')

    pyplot.show()

