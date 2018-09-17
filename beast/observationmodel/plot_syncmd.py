# plot_syncmd.py
# Plots a generic CMD from real or simulated BEAST fitting data
# PYMJ
# Created 9/13/18
# Updated 9/16/18

from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from functools import reduce

def main(fitsfile, mag1_filter='F475W', mag2_filter='F814W', 
         mag3_filter='F475W', showplot=False, saveplot=True, 
         xlim=None, ylim=None):
    """ 
    Read in flux from real or simulated data in fitsfile and plot a 
    color-magnitude diagram based on specified filters.
    
    fitsfile:           str
        input fitsfile (includes full path to file); format = .fits
    mag1_filter:        str
        1st color filter; default = 'F475W'
    mag2_filter:        str
        2nd color filter; default = 'F814W'
    mag3_filter:        str
        magnitude; default = 'F475W'
    showplot:           boolean
        Keep plot after saving it; default = False
    saveplot:           boolean
        Save plot; default = True
    xlim:               tuple
        color limit; default = None
    ylim:               tuple
        mag limit; default = None
    """
    
    fits_data = fits.open(fitsfile)
    table = fits_data[1].data
    
    mag1_flux = table['%s' % (mag1_filter + '_rate')]
    mag2_flux = table['%s' % (mag2_filter + '_rate')]
    mag_flux = table['%s' % (mag3_filter + '_rate')]


    # Exclude negative or 0 fluxes (also to avoid division by zero):    
    m1_pos_inds = np.where(mag1_flux > 0.0)
    m2_pos_inds = np.where(mag2_flux > 0.0)
    m_pos_inds = np.where(mag_flux > 0.0)
    pos_inds = reduce(np.intersect1d, (m1_pos_inds, m2_pos_inds, m_pos_inds))    
    mag1_flux_pos = mag1_flux[pos_inds]
    mag2_flux_pos = mag2_flux[pos_inds]
    mag_flux_pos = mag_flux[pos_inds]

# ===================== Convert from flux to mags ============================    
    mag1 = ((-2.5)*np.log10(mag1_flux_pos))
    mag2 = ((-2.5)*np.log10(mag2_flux_pos))
    mag = ((-2.5)*np.log10(mag_flux_pos))
    
    col = mag1 - mag2
    
    # Make cuts on col/mag if needed
    if xlim and ylim:
	lim_inds = np.where((col >= xlim[0]) & (col <= xlim[1]) & 
        	        (mag >= ylim[0]) & (mag <= ylim[1]))
	# May produce a warning if NaNs exist
        col = col[lim_inds]
        mag = mag[lim_inds]

    plt.figure(figsize=(9,9))
    plt.plot(col, mag, '.')

    if xlim: plt.xlim(xlim); 
    if ylim: plt.ylim(ylim)
    plt.gca().invert_yaxis()

    plt.xlabel('%s - %s' % (mag1_filter, mag2_filter))
    plt.ylabel(mag3_filter)

    if '/' in fitsfile:
        plottitle = fitsfile.rpartition('/')[-1].replace('fits','png')
    else:
        plottitle = fitsfile.replace('fits','png')
    plt.title(plottitle)
                
    if saveplot:
        figname = fitsfile.replace('.fits','.png')
        plt.savefig('%s' % figname)
    
    if not showplot: plt.close()
