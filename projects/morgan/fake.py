"""
Faking data

Generate a fake catalog of stars
"""
import numpy as np
import sys
from beast.core import grid
from beast.core.observations import Observations
from beast.external.eztables import Table
from eztables.core.odict import odict
from beast.external.ezpipe.helpers import RequiredFile, task_decorator


__all__ = ['make_fake_catalog', 'FakeData', 't_fakedata' ]


def make_fake_catalog(sed_grid, nstars=10, ferr=0.05, nsamp=1, outname=False, gridbackend='cache'):
    """make_fake_catalog -- generate a list of stars from a sed grid

    keywords
    --------

    sed_grid:

    nstars: int
        number of stars to generate

    ferr: float
        fractional error in flux (0.05=5%)

    nsamp:
        number of random sampling per star. Each draw will include white noise
        on the flux measurements according to ferr

    outname: False or str
        if define, the catalog table will be saved into a file (see eztables.Table formats)

    gridbackend: str or grid.GridBackend
        backend to use to load the grid if necessary (memory, cache, hdf)
        (see beast.core.grid)

    returns
    -------
    t: eztables.Table
        Table of the fake catalog
    """
    if type(sed_grid) == str:
        g0 = grid.FileSEDGrid(sed_grid, backend=gridbackend)
    else:
        g0 = sed_grid

    ind = np.random.randint(0, len(g0.grid), nstars)

    d = odict()

    for key in g0.keys():
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
            for key in g0.keys():
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
    """t_fakedata  -- task makeing  fake data

    keywords
    --------

    project:
    sed_grid:
    nstars:
    ferr:
    nsamp:

    """

    outname = '{}_fakedata.fits'.format(project)
    fake_source = RequiredFile(outname, make_fake_catalog, sed_grid, nstars=nstars, ferr=ferr, nsamp=nsamp, outname=outname)
    return project, FakeData(fake_source(), sed_grid.filters)


#---------------------------------------------------------
# Data interface                                [sec:data]
#---------------------------------------------------------

# derive the global class and update what's needed
class FakeData(Observations):
    """ class to make a catalog of fake stars from the common class """
    def __init__(self, inputFile, filters, distanceModulus=0.):
        desc = 'PHAT star: %s' % inputFile
        Observations.__init__(self, inputFile, distanceModulus, desc=desc)
        self.setFilters( filters )
        self.setBadValue(50.0)  # some bad values smaller than expected
        self.minError = 0.000
        self.floorError = 0.00  # constant error term

    def getObs(self, num=0):
        """ returns the dictionnary used during the analysis """
        assert ( not self.filters is None), "No filter set."
        mags = self.getMags(num, self.filters)
        errs = self.getErrors(num, self.filters)
        if self.badvalue is not None:
            mask = (mags >= self.badvalue)
        else:
            mask = np.zeros(len(mags), dtype=bool)

        #faking non symmetric errors
        return mags, errs, errs, mask

    def getObsinMag(self, num):
        """ Using the decorator @from_Vegamag_to_Flux
            Hence, results are in flux (not in flux/flux_vega)
            Returns the fluxes, errors and mask of an observation.
        """
        pass

    def getErrors(self, num, filters):
        """ Redifined to impose a minimal error """
        err = np.array([ self.data[tt + 'err'][num] for tt in filters])
        if self.floorError > 0.:
            err = np.sqrt(err ** 2 + self.floorError ** 2)
        if self.minError > min(err):
            err[ err < self.minError ] = self.minError
        return err
