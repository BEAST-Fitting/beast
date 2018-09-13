""" Defines a generic interface to observation catalog
    This enables to handle non detections, (upper limits one day?), flux and
    magnitude conversions to avoid painful preparation of the dataset

    Data model v2 with limited quantity units handling
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from scipy.interpolate import interp1d

__all__ = ['Observations', 'FakeObs', 'PhotCharact']


class Observations(object):
    """ A generic class that interfaces observation catalog in a standardized way

    Attributes
    ----------
    inputFile: str
        catalog source file

    filters: sequence(str)
        list of filter names (internal standards)

    desc: str, optional
        description of the observations

    badvalue: float, optional
        value that tags a bad measurement that should not be used in the
        fitting.

    nObs: int
        number of observations in the catalog
    """
    def __init__(self, inputFile, desc=None):
        """ Generate a data interface object """
        self.inputFile = inputFile
        self.filters = None
        self.desc = desc
        self.readData()
        self.badvalue = None

    @property
    def nObs(self):
        return self.data.nrows

    def __len__(self):
        return self.nObs

    def __call__(self):
        """ Calling the object will show info """
        self.info()

    def info(self):
        """ Prints some information about the catalog """
        txt = "Data read from {s.inputFile:s}\n"
        if self.desc is not None:
            txt += "Description: {s.desc:s}\n"
        txt += "Number of records: {s.nObs:d}\n\n"
        txt += "Dataset contains:"

        print("Data read from %s " % self.inputFile)
        if self.desc is not None:
            print("Description: %s" % self.desc)
            print("Number of records: %d" % self.nObs)
            print("")
            print("Dataset contains:")

        for k in list(self.data.keys()):
            txt += "\t {0:s}\n".format(k)

        if self.filters is None:
            txt += '\n No filters given yet!'
        else:
            txt += '\n Using Filters: {s.filters}\n'

        print(txt.format(s=self))

    def __getitem__(self, *args, **kwargs):
        """ get item will generate a subsample """
        return self.data.__getitem__(*args, **kwargs)

    def keys(self):
        """ Returns dataset content names """
        return self.data.keys()

    def setDescription(self, txt):
        self.desc = txt

    def setBadValue(self, val):
        self.badvalue = val

    def getFilters(self):
        return self.filters

    def setFilters(self, filters):
        self.filters = filters

    def getMags(self, num, filters):
        raise Exception('Do not use as magnitudes')
        return np.array([self.data[tt][num] for tt in filters])

    def getErrors(self, num, filters):
        raise Exception('Do not use as magnitudes')
        return np.array([self.data[tt + 'err'][num] for tt in filters])

    def getFlux(self, num):
        """returns the flux of an observation from the number of counts"""

        flux = np.empty(len(self.filters), dtype=float)
        for ek, ok in enumerate(self.filters):
            flux[ek] = self.data[ok][num]

        return flux

    def getFluxerr(self, num):
        """returns the error on the flux of an observation from the number of
        counts (not used in the analysis)"""

        fluxerr = np.empty(len(self.filters), dtype=float)

        for ek, ok in enumerate(self.filters):
            fluxerr[ek] = self.data[ok + '_err'][num]

        return fluxerr

    def getObs(self, num=0):
        """ returns the flux"""
        if self.filters is None:
            raise AttributeError('No filter set provided.')

        flux = self.getFlux(num, self.filters)

        return flux

    def readData(self):
        """ read the dataset from the original source file """
        from ..external.eztables import AstroTable
        if type(self.inputFile) == str:
            self.data = AstroTable(self.inputFile)
        else:
            self.data = self.inputFile

    def iterobs(self):
        """ yield getObs """
        for k in range(self.nObs):
            yield self.getObs(k)

    def enumobs(self):
        for k in range(self.nObs):
            yield k, self.getObs(k)


def gen_SimObs_from_sedgrid(sedgrid, noisemodel):
    """
    Generate simulated observations using the physics and observation grids.
    The priors are sampled as they give the ensemble model for the stellar
    and dust distributions (IMF, Av distribution etc.).
    The physics model gives the SEDs based on the priors.
    The observation model gives the noise, bias, and completeness all of
    which are used in simulating the observations.

    Currently written to only work for the toothpick noisemodel.

    Keywords
    ---------
    sedgrid: grid.SEDgrid instance
        model grid

    sedgrid_noisemodel: beast noisemodel instance
        noise model data
    """


def gen_FakeObs_from_sedgrid(sedgrid, nrows, err=0.05,
                             filters=None, save=False):
    from ..external.eztables import Table
    from . import grid
    if type(sedgrid) == str:
        sedgrid = grid.FileSEDGrid(sedgrid)

    inds = np.random.randint(0, high=sedgrid.grid.nrows, size=nrows)
    obsTab = Table()
    if filters is None:
        filters = sedgrid.grid.header.FILTERS.split()
    for e, filt in enumerate(filters):
        errs = np.random.normal(loc=0., scale=err, size=nrows)
        obsTab.addCol(filt, (1. + errs) * sedgrid.seds[inds, e])
        obsTab.addCol(filt + 'err', err * sedgrid.seds[inds, e])
    for key in list(sedgrid.grid.keys()):
        obsTab.addCol(key, sedgrid.grid[key][inds])

    if save is True:
        return obsTab
    else:
        obsTab.write(save, clobber=True, append=False)
