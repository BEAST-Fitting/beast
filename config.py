# Set to use some C code instead of pure python to speed up the computations.
# If False, only numpy and python code are used.
#__WITH_C_LIBS__ = True
__WITH_C_LIBS__ = False

#numexpr -- optimized multi-threaded numpy evaluations
__USE_NUMEXPR__ = True

# Default number of threads to use
__NTHREADS__ = 25

# Online libraries
# will be replaced by a more flexible support (JSON is easy!)
libs_server = 'http://chex.astro.washington.edu:7777/beastlibs/'
libs = dict(
    vega     = 'vega.hd5',
    filters  = 'filters.hd5',
    #basel22p = 'BaSeL_v2.2.pegase.grid.fits',
    #elodie31 = 'Elodie_v3.1.grid.fits',
    kurucz04 = 'kurucz2004.grid.fits',
    tlusty09 = 'tlusty.grid.fits'
    #specgrid = 'basel_padova2010.spectralgrid.fits',
    #sedgrid  = 'PHATSEDs_basel_padovaiso.fits'
    #kuruczgrid = 'kurucz2004.grid.fits',
    #padovaiso  = 'padova2010.iso.fits',
    #kuruczisog  = 'stellib_kurucz2004_padovaiso.spectralgrid.fits'
)

#Make sure the configuration is coherent for the python installation
try:
    import numexpr
    numexpr.set_num_threads(__NTHREADS__)
except ImportError:
    __USE_NUMEXPR__ = False

import os
import inspect


__ROOT__ = '/'.join(os.path.abspath(inspect.getfile(inspect.currentframe())).split('/')[:-1])


def printConfig():
	print """ ============ pyPEGASE defaut configuration ===========
	* Including C-code during computations: %s
	* Parallel processing using %d threads
	""" % (__WITH_C_LIBS__, __NTHREADS__)
