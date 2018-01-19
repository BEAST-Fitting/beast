from . import astrohelpers
from ..table import Table, __indent__
import numpy as np
from copy import deepcopy


class AstroTable(Table):
    """ Derived from the Table, this class add implementations of common astro tools especially conesearch """
    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        self.ra_name, self.dec_name = self.__autoRADEC__()
        if (len(args) > 0):
            if isinstance(args[0], AstroTable):
                self.ra_name = args[0].ra_name
                self.dec_name = args[0].dec_name
        self.ra_name = kwargs.get('ra_name', self.ra_name)
        self.dec_name = kwargs.get('dec_name', self.dec_name)

    def __autoRADEC__(self):
        """ Tries to identify the columns containing RA and DEC coordinates """
        if 'ra' in self:
            ra_name = 'ra'
        elif 'RA' in self:
            ra_name = 'RA'
        else:
            ra_name = None
        if 'dec' in self:
            dec_name = 'dec'
        elif 'DEC' in self:
            dec_name = 'DEC'
        else:
            dec_name = None
        return ra_name, dec_name

    def setRA(self, val):
        """ Set the column that defines RA coordinates """
        assert(val in self), 'column name {} not found in the table'.format(val)
        self.ra_name = val

    def setDEC(self, val):
        """ Set the column that defines DEC coordinates """
        assert(val in self), 'column name {} not found in the table'.format(val)
        self.dec_name = val

    def getRA(self, degree=True):
        """ Returns RA, converted from hexa/sexa into degrees """
        if self.ra_name is None:
            return None
        if (not degree) or (self.dtype[self.ra_name].kind != 'S'):
            return self[self.ra_name]
        else:
            if (len(str(self[0][self.ra_name]).split(':')) == 3):
                return np.asarray(astrohelpers.hms2deg(self[self.ra_name], delim=':'))
            elif (len(str(self[0][self.ra_name]).split(' ')) == 3):
                return np.asarray(astrohelpers.hms2deg(self[self.ra_name], delim=' '))
            else:
                raise Exception('RA Format not understood')

    def getDEC(self, degree=True):
        """ Returns RA, converted from hexa/sexa into degrees """
        if self.dec_name is None:
            return None
        if (not degree) or (self.dtype[self.dec_name].kind != 'S'):
            return self[self.dec_name]
        else:
            if (len(str(self[0][self.dec_name]).split(':')) == 3):
                return np.asarray(astrohelpers.dms2deg(self[self.dec_name], delim=':'))
            elif (len(str(self[0][self.dec_name]).split(' ')) == 3):
                return np.asarray(astrohelpers.dms2deg(self[self.dec_name], delim=' '))
            else:
                raise Exception('RA Format not understood')

    def info(self):
        print(self.header)
        print("Table contains: %i row(s) in %i column(s)\n" % (self.nrows, self.ncols))
        if (self.ra_name is not None) & (self.dec_name is not None):
            print("Position coordinate columns: {}, {}\n".format(self.ra_name, self.dec_name))
        if self._aliases is not None:
                if len(self._aliases) > 0:
                        print("Table contains alias(es):")
                        for k, v in self._aliases.items():
                                print('\t %s --> %s' % (k, v))
                        print('')
        fields = 'columns unit format description'.split()
        row    = [ (k, self.columns[k].unit, self.columns[k].format, self.columns[k].description) for k in list(self.keys()) ]
        out    = __indent__([fields] + row, hasHeader=True, hasUnits=False, delim=' ')
        print(out)

    def coneSearch(self, ra, dec, r, outtype=0):
        """ Perform a cone search on a table
        INPUTS:
            ra0 	ndarray[ndim=1, dtype=float]	column name to use as RA source in degrees
            dec0	ndarray[ndim=1, dtype=float]	column name to use as DEC source in degrees
            ra		float                       	ra to look for (in degree)
            dec		float	                        ra to look for (in degree)
            r		float		                    distance in degrees
        KEYWORDS:
            outtype int                             0 -- minimal, indices of matching coordinates
                                                    1 -- indices and distances of matching coordinates
                                                    2 -- full, boolean filter and distances
        """

        assert( (self.ra_name is not None) & (self.dec_name is not None) ), 'Coordinate columns not set.'

        ra0  = self.getRA()
        dec0 = self.getDEC()
        return astrohelpers.conesearch(ra0, dec0, ra, dec, r, outtype=outtype)

    def zoneSearch(self, ramin, ramax, decmin, decmax, outtype=0):
        """ Perform a zone search on a table, i.e., a rectangular selection
        INPUTS:
            ra0 	ndarray[ndim=1, dtype=float]	column name to use as RA source in degrees
            dec0	ndarray[ndim=1, dtype=float]	column name to use as DEC source in degrees
            ra		float                       	ra to look for (in degree)
            dec		float	                        ra to look for (in degree)
            r		float		                    distance in degrees
        KEYWORDS:
            outtype int                             0 -- minimal, indices of matching coordinates
                                                    1 -- indices and distances of matching coordinates
                                                    2 -- full, boolean filter and distances
        """

        assert( (self.ra_name is not None) & (self.dec_name is not None) ), 'Coordinate columns not set.'

        ra0  = self.getRA()
        dec0 = self.getDEC()
        ind = (ra0 >= ramin) & (ra0 <= ramax) & (dec0 >= decmin) & (dec0 <= decmax)
        if outtype <= 2:
            return ind
        else:
            return np.where(ind)

    def selectWhere(self, fields, condition=None, condvars=None, cone=None, zone=None, **kwargs):
        """ Read table data fulfilling the given `condition`.
            Only the rows fulfilling the `condition` are included in the result.
            conesearch is also possible through the keyword cone formatted as (ra, dec, r)
            zonesearch is also possible through the keyword zone formatted as (ramin, ramax, decmin, decmax)

            Combination of multiple selections is also available.
        """
        if cone is not None:
            if len(cone) != 3:
                raise ValueError('Expecting cone keywords as a triplet (ra, dec, r)')
        if zone is not None:
            if len(zone) != 4:
                raise ValueError('Expecting zone keywords as a tuple of 4 elements (ramin, ramax, decmin, decmax)')

        # only a cone search
        # make a copy without the data itself (memory gentle)
        tab = self.__class__()
        for k in list(self.__dict__.keys()):
            if k != 'data':
                setattr(tab, k, deepcopy(self.__dict__[k]))
        if (condition is None) & (cone is None):
            tab.data = deepcopy(self.data)

        if fields.count(',') > 0:
            _fields = fields.split(',')
        elif fields.count(' ') > 0:
            _fields = fields.split()
        else:
            _fields = fields

        if (condition in [True, 'True', None]):
            if (cone is None) & (zone is None):
                ind = None
            elif (cone is not None) & (zone is None):
                ra, dec, r = cone
                ind, d = self.coneSearch(ra, dec, r, outtype=1)
            elif (cone is None) & (zone is not None):
                ind = self.zoneSearch(zone[0], zone[1], zone[2], zone[3], outtype=1)
            else:  # cone + zone
                ra, dec, r = cone
                ind, d = self.coneSearch(ra, dec, r, outtype=2)
                ind = ind & self.zoneSearch(zone[0], zone[1], zone[2], zone[3], outtype=2)
                d = d[ind]
                ind = np.where(ind)
        else:
            if condvars is None:
                condvars = {}
            condi = '({0:s})'.format(condition)
            if cone is not None:
                ra, dec, r = cone
                cind, d = self.coneSearch(ra, dec, r, outtype=2)
                d = d[cind]
                condvars['_cone_'] = cind
                condi += ' & (_cone_)'
            if zone is not None:
                zind = self.zoneSearch(zone[0], zone[1], zone[2], zone[3], outtype=2)
                condvars['_zone_'] = zind
                condi += ' & (_zone_)'
            ind = self.where(condi, condvars=condvars, **kwargs)

        if _fields == '*':
            if ind is not None:
                tab.data = self.data[ind]
            else:
                tab.data = deepcopy(self.data)
        else:
            if ind is not None:
                tab.data = self.data[tab.resolve_alias(_fields)][ind]
            else:
                tab.data = self.data[tab.resolve_alias(_fields)]
            names = tab.data.dtype.names
            #cleanup aliases and columns
            for k in list(self.keys()):
                if k not in names:
                    al = self.reverse_alias(k)
                    for alk in al:
                        tab.delCol(alk)
                    if k in list(tab.keys()):
                        tab.delCol(k)
        if cone is not None:
            tab.add_column('separation', np.asarray(d), unit='degree')
            tab.header['COMMENT'] = 'SELECT %s FROM %s WHERE %s' % (fields, self.header['NAME'], 'distance from (%0.3f, %0.3f) <= %0.3f' % (ra, dec, r) )
        if zone is not None:
            tab.header['COMMENT'] = 'SELECT %s FROM %s WHERE %s' % (fields, self.header['NAME'], 'coordinates inside box(%0.3f, %0.3f, %0.3f, %0.3f)' % (zone[0], zone[1], zone[2], zone[3]) )
        if condition not in [True, 'True', None]:
            tab.header['COMMENT'] = 'SELECT %s FROM %s WHERE %s' % (fields, self.header['NAME'], str(condition) )

        tab.setRA(self.ra_name)
        tab.setDEC(self.dec_name)
        return tab
