"""
Create the catalog from the nD likelihood functions

Jan 2013: Written Karl G.

"""

__version__ = '0.1dev'

import output
import glob

if __name__ == '__main__':

    grid_filename='libs/stellib_kurucz2004_padovaiso.spectralgrid_sed_extinguished.grid.fits'
	
    filters  = 'hst_wfc3_f225w hst_wfc3_f336w hst_acs_hrc_f475w hst_acs_hrc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split()

    files = glob.glob('Tests/fake_many_0/fake_star*.fits')

    output.summaryTable(files, grid_filename, filters, 'Tests/fake_many_0.fits')
