# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This is an Astropy affiliated package.
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

# For egg_info test builds to pass, put package imports here.
#if not _ASTROPY_SETUP_:
#    from .example_mod import *

# BEAST specific
#   this is all old.  Need to check if it is needed.  Currently commented out during reorg.

#__all__ = ['external', 'tools', 'config', 'anased', 'creategrid', 'extinction', 'grid',
#           'isochrone', 'observations', 'phot', 'photometry', 'stellib', 'vega' ]

#from . import external
#from . import tools
#from . import config
#from .core import anased
#from .core import creategrid
#from .core import extinction
#from .core import grid
#from .core import isochrone
#from .core import observations
#from .core import phot
#from .core import photometry
#from .core import stellib
#from .core import vega
#from .core import noisemodel
