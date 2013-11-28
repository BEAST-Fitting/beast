from . import stellib
from matplotlib.path import Path


class CompositeStellib(object):
    def __init__(self, osllist, *args, **kwargs):
        self._olist = osllist

    def which_osl(self):
        pass

    def genQ(self, qname, r, **kwargs):
        pass

    def genSpectrum(self, T0, g0=None, Z0=None, weights=None, **kwargs):
        pass

    def interpMany(self, T0, g0, Z0, L0, dT_max=0.1, eps=1e-6, weights=None, **kwargs):
        pass

    def get_boundaries(self, dlogT=0.1, dlogg=0.3, closed=True):
        pass

    def __repr__(self):
        return "CompositeStellib, {0}\n{1}".format(object.__repr__(self), [k.name for k in self._olist])

