import numpy as np
from scipy import interpolate
from numpy import interp
from sedfitter import isochrone
from sedfitter.ezunits import unit, hasUnit
from sedfitter.eztables import Table


class ezIsoch(isochrone.Isochrone):
    """ Trying to make something that is easy to manipulate
    This class is basically a proxy to a table (whatever format works best)
    and tries to keep things coherent.
    """
    def __init__(self, source, interp=False):
        self.name = '<auto>'
        self.source = source
        self._load_table_(self.source)
        self.logages = np.unique(np.round(self.data['logA'], 6))  # round because of precision noise
        self.ages = np.round(10 ** self.logages)
        self.Z    = np.unique(np.round(self.data['Z'], 6))
        self.interpolation(interp)

    def selectWhere(self, *args, **kwargs):
        return self.data.selectWhere(*args, **kwargs)

    def interpolation(self, b=None):
        if b is not None:
            if hasattr(self, 'interp'):
                print "Do not use interpolation yet, at your own risks!!"
            self.interp = bool(b)
        else:
            return self.interp

    def _load_table_(self, source):
        self.data = Table(self.source).selectWhere('*', 'isfinite(logA)')

    def __getitem__(self, key):
        return self.data[key]

    def _get_isochrone(self, age, metal=None, FeH=None, inputUnit=unit['yr'], masses=None, *args, **kwargs):
        """ Retrieve isochrone from the original source
            internal use to adapt any library
        """
        if hasUnit(age):
            _age = int(age.to('yr').magnitude)
        else:
            _age = int(age * inputUnit.to('yr').magnitude)

        _logA = np.log10(_age)

        assert ((metal is not None) | (FeH is not None)), "Need a chemical par. value."

        if (metal is not None) & (FeH is not None):
            print "Warning: both Z & [Fe/H] provided, ignoring [Fe/H]."

        if metal is None:
            metal = self.FeHtometal(FeH)

        if self.interpolation():
            #Do the actual nd interpolation

            #Maybe already exists?
            if (metal in self.Z) & (_age in self.ages):
                t = self.selectWhere('*', '(round(Z, 6) == {}) & (round(logA, 6) == {})'.format(metal, _logA))
                if t.nrows > 0:
                    return t
            #apparently not
            #find 2 closest metal values
            ca1 = (self.ages <= _age)
            ca2 = (self.ages > _age)
            cz1 = (self.Z <= metal)
            cz2 = (self.Z > metal)
            if (metal in self.Z):
                #perfect match in metal, need to find ages
                if (_age in self.ages):
                    return self.selectWhere('*', '(round(Z, 6) == {}) & (round(logA, 6) == {})'.format(metal, _logA))
                elif ( True in ca1) & ( True in ca2 ):
                    # bracket on _age: closest values
                    a1, a2 = np.log10(max(self.ages[ca1])), np.log10(min(self.ages[ca2]))
                    iso = self.selectWhere('*', '(Z == 0.02) & ( (abs(logA - {}) < 1e-4) | (abs(logA - {}) < 1e-4 )  )'.format(a1, a2) )
                    if masses is None:
                        _logM = np.unique(iso['logM'])
                    else:
                        _logM = masses

                    #define interpolator
                    points = np.array([self[k] for k in 'logA logM Z'.split()]).T
                    values = np.array([ self[k] for k in self.data.keys() ]).T
                    _ifunc = interpolate.LinearNDInterpolator(points, values)

                    pts = np.array([ (_logA, logMk, metal) for logMk in _logM ])
                    r = _ifunc(pts)
                    return Table(r)
                else:
                    raise Exception('Age not covered by the isochrones')
            elif ( True in cz1 ) & ( True in cz2 ):
                #need to find closest Z
                pass
            return
        else:
            # find the closest match
            _Z = self.Z[((metal - self.Z) ** 2).argmin()]
            #_logA = np.log10(self.ages[((_age - self.ages) ** 2).argmin()])
            _logA = self.logages[ ((np.log10(_age) - self.logages) ** 2).argmin() ]
            tab = self.data.selectWhere('*', "(round(Z, 6) == {}) & (round(logA,6) == {})".format(_Z, _logA))
            #mass selection
            if masses is not None:
                #masses are expected in logM for interpolation
                #if masses.max() > 2.3:
                #    _m = np.log10(masses)
                #else:
                _m = masses
                data_logM = tab['logM'][:]
                # refuse extrapolation!
                #ind = np.where(_m <= max(data_logM))
                data = {}
                for kn in tab.keys():
                    data[kn] = interp(_m, data_logM, tab[kn], left=np.nan, right=np.nan)
                return Table(data)
