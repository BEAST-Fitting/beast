"""
Faking data
"""
from beast.core import grid
from beast.core.observations import Observations
from eztables import Table
from eztables.core.odict import odict
import numpy as np
from ezpipe.helpers import RequiredFile, task_decorator
import sys


__all__ = ['make_fake_catalog', 'FakeData', 't_fakedata' ]


def make_fake_catalog(sed_grid, nstars=10, ferr=0.05, nsamp=1, outname=False):
    """ generate stars from a sed grid """

    if type(sed_grid) == str:
        g0 = grid.FileSEDGrid(sed_grid)
    else:
        g0 = sed_grid

    ind = np.random.randint(0, g0.grid.nrows, nstars)

    d = odict()

    for key in g0.grid.keys():
        d[key] = np.zeros(nstars * nsamp, dtype=float)

    for key in g0.filters:
        d[key] = np.zeros(nstars * nsamp, dtype=float)
        d[key + 'err'] = np.zeros(nstars * nsamp, dtype=float)

    # oh that's ugly...
    for ek, k in enumerate(ind):
        datak = g0.grid[k]
        sedk = g0.seds[k]

        for sampk in range(nsamp):
            rowk = ek * nsamp + sampk
            for key in g0.grid.keys():
                d[key][rowk] = datak[key]

            for ef, f in enumerate(g0.filters):
                _err = np.random.normal(0, ferr, 1)
                d[f][rowk] = sedk[ef] * (1. + _err)
                d[f + 'err'][rowk] = abs(_err * sedk[ef])

    t = Table(d, name='Fake_data')
    if outname is not False:
        t.write(outname)
    return t


#---------------------------------------------------------
# Pipeline interface                        [sec:pipeline]
#---------------------------------------------------------

@task_decorator(logger=sys.stdout)
def t_fakedata(project, sed_grid, nstars=10, ferr=0.05, nsamp=1):
    outname = '{}_fakedata.fits'.format(project)
    fake_source = RequiredFile(outname, make_fake_catalog, sed_grid, nstars=nstars, ferr=ferr, nsamp=nsamp, outname=outname)
    return project, FakeData(fake_source(), sed_grid.filters)


#---------------------------------------------------------
# Data interface                                [sec:data]
#---------------------------------------------------------

# derive the global class and update what's needed
class FakeData(Observations):
    """ PHAT catalog for clusters in M31 """
    def __init__(self, inputFile, filters, distanceModulus=0.):
        desc = 'PHAT star: %s' % inputFile
        Observations.__init__(self, inputFile, distanceModulus, desc=desc)
        self.setFilters( filters )
        self.setBadValue(50.0)  # some bad values smaller than expected
        self.minError = 0.000
        self.floorError = 0.00  # constant error term

    def getObs(self, num):
        """ Using the decorator @from_Vegamag_to_Flux
            Hence, results are in flux (not in flux/flux_vega)
            Returns the fluxes, errors and mask of an observation.
        """
        return Observations.getObs(self, num)

    def getObsinMag(self, num):
        """ Returns the original catalog magnitudes """
        pass

    def getErrors(self, num, filters):
        """ Redifined to impose a minimal error """
        err = np.array([ self.data[tt + 'err'][num] for tt in filters])
        if self.floorError > 0.:
            err = np.sqrt(err ** 2 + self.floorError ** 2)
        if self.minError > min(err):
            err[ err < self.minError ] = self.minError
        return err
