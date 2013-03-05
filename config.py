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
libs_server = 'http://chex.astro.washington.edu:8899/sedfitter/'
libs = dict(
    vega     = 'vega.hd5',
    filters  = 'filters.hd5',
    #specgrid = 'basel_padova2010.spectralgrid.fits',
    #sedgrid  = 'PHATSEDs_basel_padovaiso.fits'
    kuruczgrid = 'kurucz2004.grid.fits',
    padovaiso  = 'padova2010.iso.fits',
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
