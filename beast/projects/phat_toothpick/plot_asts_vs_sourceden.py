""" Script to plot the variation in the noise model versus source density

.. history::
    Written 22dec14 for the BEAST paper

"""
import numpy as np
import tables
import matplotlib.pyplot as pyplot

from beast.core.noisemodel import toothpick
from beast.external.ezpipe.helpers import RequiredFile, task_decorator
from beast.external.eztables import Table
from merge_phat_asts import merge_phat_asts

class PHAT_ToothPick_Noisemodel(toothpick.MultiFilterASTs):

    def set_data_mappings(self):
        """ hard code mapping directly with the interface to ASTs format

        .. note::

            As noted in the class documentation, it is trivial to adapt the
            class to specific formats
        """
        for k in self.filters:
            try:
                self.data.set_alias(k + '_out', k.split('_')[-1].upper() + '_VEGA')
                self.data.set_alias(k + '_in', k.split('_')[-1].upper() + '_IN')
            except Exception as e:
                print(e)
                print('Warning: Mapping failed. This could lead to wrong results')

        # hard coded for testing single camera ASTs
        # later need to write code to create to appropriate FITS file from
        #   the single camera AST files
        #self.data.set_alias(filters[0] + '_out', 'MAG1OUT')
        #self.data.set_alias(filters[0] + '_in', 'MAG1IN')
        #self.data.set_alias(filters[1] + '_out', 'MAG2OUT')
        #self.data.set_alias(filters[1] + '_in', 'MAG2IN')

if __name__ == '__main__':

    vega_flux = [3.73563424e-09,3.25401916e-09,5.31980011e-09,1.13432076e-09,4.08364770e-10,1.44227092e-10]
    filters = ['HST_F275W','HST_F336W','HST_F475W','HST_F814W','HST_F110W','HST_F160W']
    filters_tag = ['F275W','F336W','F475W','F814W','F110W','F160W']

    fig, ax = pyplot.subplots(6, sharex=True, figsize=(10,13))
    
    blabels = ['1_2','10_11','20_21','29_30']
    #blabels = ['b15_0','b15_10','b15_20','b15_27']
    #blabels = ['b15_0','b15_20']
    #blabels = ['b15_0']

    colsym = ['ro','bo','go','yo']
    col = ['r','b','g','y']

    for k in range(len(blabels)):
        blabel = blabels[k]
        basename = 'BEAST_production/merged_asts/PHAT_fake_stars_SD_' + blabel
	astfile = basename + '.fits'
        #astfile = basename+'_all.hd5'
        #merge_phat_asts(basename+'_uv.fits',basename+'_opt.fits',basename+'_ir.fits',astfile)
        model = PHAT_ToothPick_Noisemodel(astfile, filters)
        
        for i in range(6):
            mag_in = model.data[filters[i] + '_in']
            mag_out = model.data[filters[i] + '_out']
            
            completeness_mag_cut = 80
            
            d = model._compute_stddev_bins(mag_in, mag_out, nbins=30, completeness_mag_cut=completeness_mag_cut)
            
	    if i == 0:
                ax[i].plot(vega_flux[i]*d['FLUX_IN'],d['FLUX_BIAS']/d['FLUX_IN'],colsym[k], label=blabel + ' $\mu$')
	    else:
                ax[i].plot(vega_flux[i]*d['FLUX_IN'],d['FLUX_BIAS']/d['FLUX_IN'],colsym[k])

	    if i == 1:
                ax[i].plot(vega_flux[i]*d['FLUX_IN'],(d['FLUX_STD']/d['FLUX_IN']),colsym[k], markerfacecolor='none',markeredgecolor=col[k], label=blabel + ' $\sigma$')
            else:
                ax[i].plot(vega_flux[i]*d['FLUX_IN'],(d['FLUX_STD']/d['FLUX_IN']),colsym[k], markerfacecolor='none',markeredgecolor=col[k])

            ax[i].text(1e-17,0.6,filters_tag[i])

    for i in range(6):
        ax[i].set_yscale('linear')
        #ax[i].set_ylabel('($\mu$ or $\sigma$)/F(band)') 
        ax[i].set_ylabel('$\mu$/F(band)') 
        #ax[i].set_ylabel('Sigma/Flux(Input)')
        #ax[i].set_ylim(-0.5,2.0)
        ax[i].set_ylim(-0.1,0.75)
        ax[i].set_xscale('log')
        ax[i].set_xlim(1e-20,2e-14)

    ax[5].set_xlabel('F(band) [ergs cm$^{-2}$ s$^{-1}$ $\AA^{-1}$]')
    ax[0].legend()
    ax[1].legend()


    fig.subplots_adjust(hspace=0.05)

    pyplot.show()

