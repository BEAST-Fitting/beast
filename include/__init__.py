""" C extensions -- replace some python code when C interface is available """

try:
	from .interp import *
except ImportError:
	pass
	#import setup
	#setup.main()
	#from .ctools import *	


