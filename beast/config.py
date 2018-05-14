from __future__ import print_function
import os
from os.path import expanduser
import inspect

# Set to use some C code instead of pure python to speed up the computations.
# If False, only numpy and python code are used.
# __WITH_C_LIBS__ = True
__WITH_C_LIBS__ = False

# numexpr -- optimized multi-threaded numpy evaluations
__USE_NUMEXPR__ = True

# Default number of threads to use when parallel computing
# __NTHREADS__ = 6
__NTHREADS__ = 1

# library directory
beast_envvar = "BEAST_LIBS"
userhome = expanduser("~")
ploc = userhome + "/.beast/"
if beast_envvar in os.environ:
    __ROOT__ = os.environ[beast_envvar]
elif os.path.isdir(ploc):
    __ROOT__ = ploc
else:
    __ROOT__ = '/'.join(os.path.abspath(inspect.getfile(inspect.currentframe())).split('/')[:-1])
    __ROOT__ += '/libs/'

# Online libraries
# will be replaced by a more flexible support (JSON is easy!)
libs_server = 'http://www.stsci.edu/~kgordon/beast/'
libs = dict(
    vega='vega.hd5',
    filters='filters.hd5',
    kurucz04='kurucz2004.grid.fits',
    tlusty09='tlusty.lowres.grid.fits',
    hstcovar='hst_whitedwarf_frac_covar.fits'
    # basel22p = 'BaSeL_v2.2.pegase.grid.fits',
    # elodie31 = 'Elodie_v3.1.grid.fits'
)

# Make sure the configuration is coherent for the python installation
try:
    import numexpr
    if not __USE_NUMEXPR__:
        numexpr.set_num_threads(1)
        numexpr.set_vml_num_threads(1)
    else:
        numexpr.set_num_threads(__NTHREADS__)
        numexpr.set_vml_num_threads(__NTHREADS__)
except ImportError:
    __USE_NUMEXPR__ = False

try:
    import tables
    tables.parameters.MAX_NUMEXPR_THREADS = __NTHREADS__
    tables.parameters.MAX_BLOSC_THREADS = __NTHREADS__
    tables.set_blosc_max_threads(__NTHREADS__)
except ImportError:
    pass


def printConfig():
    print(""" ============ BEAST defaut configuration ===========
    * Including C-code during computations: %s
    * Parallel processing using %d threads
    """ % (__WITH_C_LIBS__, __NTHREADS__))
