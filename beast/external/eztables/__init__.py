from .table import *
from .astro import AstroTable

__version__ = "1.0"
__author__ = "MF"
import os
import inspect

localpath = "/".join(
    os.path.abspath(inspect.getfile(inspect.currentframe())).split("/")[:-1]
)
# __doc__ = open(localpath+'/README.md').read()
