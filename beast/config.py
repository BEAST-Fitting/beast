import os
from os.path import expanduser

# Set to use some C code instead of pure python to speed up the computations.
# If False, only numpy and python code are used.
# __WITH_C_LIBS__ = True
__WITH_C_LIBS__ = False

# numexpr -- optimized multi-threaded numpy evaluations
__USE_NUMEXPR__ = True

# Default number of threads to use when parallel computing
# __NTHREADS__ = 6
__NTHREADS__ = 1

# useful constant(s)
solar_metalicity = 0.0153

# library directory
beast_envvar = "BEAST_LIBS"
userhome = expanduser("~")
ploc = userhome + "/.beast/"
if beast_envvar in os.environ:
    __ROOT__ = os.environ[beast_envvar]
else:
    if not os.path.exists(ploc):
        os.mkdir(ploc)
    __ROOT__ = ploc

# Online libraries
# will be replaced by a more flexible support (JSON is easy!)
libs_server = "http://www.stsci.edu/~kgordon/beast/"
libs = dict(
    vega="vega.hd5",
    filters="filters.hd5",
    kurucz04="kurucz2004.grid.fits",
    tlusty09="tlusty.lowres.grid.fits",
    hstcovar="hst_whitedwarf_frac_covar.fits",
    mist1="MIST_FeH-4.00_vvcrit0.4.fits",
    mist2="MIST_FeH-3.50_vvcrit0.4.fits",
    mist3="MIST_FeH-3.00_vvcrit0.4.fits",
    mist4="MIST_FeH-2.50_vvcrit0.4.fits",
    mist5="MIST_FeH-2.00_vvcrit0.4.fits",
    mist6="MIST_FeH-1.75_vvcrit0.4.fits",
    mist7="MIST_FeH-1.50_vvcrit0.4.fits",
    mist8="MIST_FeH-1.25_vvcrit0.4.fits",
    mist9="MIST_FeH-1.00_vvcrit0.4.fits",
    mist10="MIST_FeH-0.75_vvcrit0.4.fits",
    mist11="MIST_FeH-0.50_vvcrit0.4.fits",
    mist12="MIST_FeH-0.25_vvcrit0.4.fits",
    mist13="MIST_FeH0.00_vvcrit0.4.fits",
    mist14="MIST_FeH0.25_vvcrit0.4.fits",
    mist15="MIST_FeH0.50_vvcrit0.4.fits",
    # basel22="stellib_BaSeL_v2.2.grid.fits",
    # munari="atlas9-munari.hires.grid.fits",
    # btsettl="bt-settl.lowres.grid.fits",
    # elodie31 = 'Elodie_v3.1.grid.fits'
    # aringer="Aringer.AGB.grid.fits"
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
    print(
        """ ============ BEAST defaut configuration ===========
    * Including C-code during computations: %s
    * Parallel processing using %d threads
    """
        % (__WITH_C_LIBS__, __NTHREADS__)
    )
