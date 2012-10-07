# Set to use some C code instead of pure python to speed up the computations.
# If False, only numpy and python code are used.
__WITH_C_LIBS__ = True
#__WITH_C_LIBS__ = False

# Default number of threads to use
__NTHREADS__ = 3


def printConfig():
	print """ ============ pyPEGASE defaut configuration ===========
	* Including C-code during computations: %s
	* Parallel processing using %d threads
	""" % (__WITH_C_LIBS__, __NTHREADS__)
