# Set to use some C code instead of pure python to speed up the computations.
# If False, only numpy and python code are used.
__WITH_C_LIBS__ = True
#__WITH_C_LIBS__ = False

# Default number of threads to use
__NTHREADS__ = 3

# Online libraries
# will be replaced by a more flexible support (JSON is easy!)
libs_server = 'http://chex.astro.washington.edu:8899/sedfitter/'
libs = dict(
    vega     = 'vega.hd5',
    filters  = 'filters.hd5',
    #specgrid = 'basel_padova2010.spectralgrid.fits',
    #sedgrid  = 'PHATSEDs_basel_padovaiso.fits'
    kuruczgrid  = 'kurucz2004.grid.fits'
)

def printConfig():
	print """ ============ pyPEGASE defaut configuration ===========
	* Including C-code during computations: %s
	* Parallel processing using %d threads
	""" % (__WITH_C_LIBS__, __NTHREADS__)
