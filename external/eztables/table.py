"""
This module contains the base of any table regardless of their storage
"""

from __future__ import (absolute_import, division, print_function)

import numpy as np
from .core.odict import odict
from .core.helpers import *
from .core.tableheader import TableHeader
from .core.columnheader import ColumnHeader
import operator
import io
from numpy.lib import recfunctions
from .registered_backends import *
from copy import deepcopy
from .core.decorators import warning
import sys
from functools import reduce

try:
    unicode = unicode
except NameError:
    # 'unicode' is undefined, must be Python 3
    str = str
    unicode = str
    bytes = bytes
    basestring = (str,bytes)
else:
    # 'unicode' exists, must be Python 2
    str = str
    unicode = unicode
    bytes = str
    basestring = basestring

"""
Notes

a table is
   a group of *Ordered* columns
   a header for description of the table (incl. title, name)
   a description of each column

"""

__all__ = ['Table', 'ColumnHeader', 'TableHeader']


def pretty_size_print(num_bytes):
    """
    Output number of bytes in a human readable format
    """
    if num_bytes is None:
        return

    KiB = 1024
    MiB = KiB * KiB
    GiB = KiB * MiB
    TiB = KiB * GiB
    PiB = KiB * TiB
    EiB = KiB * PiB
    ZiB = KiB * EiB
    YiB = KiB * ZiB

    if num_bytes > YiB:
        output = '%.3g YB' % (num_bytes / YiB)
    elif num_bytes > ZiB:
        output = '%.3g ZB' % (num_bytes / ZiB)
    elif num_bytes > EiB:
        output = '%.3g EB' % (num_bytes / EiB)
    elif num_bytes > PiB:
        output = '%.3g PB' % (num_bytes / PiB)
    elif num_bytes > TiB:
        output = '%.3g TB' % (num_bytes / TiB)
    elif num_bytes > GiB:
        output = '%.3g GB' % (num_bytes / GiB)
    elif num_bytes > MiB:
        output = '%.3g MB' % (num_bytes / MiB)
    elif num_bytes > KiB:
        output = '%.3g KB' % (num_bytes / KiB)
    else:
        output = '%.3g Bytes' % (num_bytes)

    return output


#================================================================================
class Table(object):
    """ This class implements a Table object which aims at being able to
    manage dataset table independently of its storage format.

    Tables allows row and column access using the __getitem__ intrinsic
    method. Columns can be added to or deleted from the table.
    Any columns can be described by its own header which includes units.

    Besides Table object can be considered as dict objects and so are
    iterable over the columns,  Table's data are stored under np.recarray
    format for memory footprint and flexibility.


    reading and writing methods are not directly defined in this class but
    apart using TableManager objects that are then registered with the
    "register_extension" function.
    """
#================================================================================
    def __init__(self, *args, **kwargs):
        """
        Create a table instance

        *args: Optional Arguments:

            If no arguments are given, and empty table is created

            If one or more arguments are given they are passed to the
            Table.read() method.

        **kwargs: Optional Keyword Arguments (independent of table type):

            name: [ string ]
                The table name
        """
        self.__set_defaults__()

        if len(args) > 0:
            if isinstance(args[0], self.__class__):
                self.__set_from_table__(from_Table(*args, **kwargs))
            elif (type(args[0]) in [ np.core.records.recarray, np.ndarray]):
                self.__set_from_table__(from_ndArray(*args, **kwargs))
            elif hasattr(args[0], 'iteritems') or isinstance(args[0], dict):
                self.__set_from_table__(from_dict(*args, **kwargs))
            else:
                self.read(*args, **kwargs)
        self.set_name( kwargs.get('name', None) )
        self.caseless = kwargs.get('caseless', False)

    def read(self, filename, type=None, manager=None, silent=False, **kwargs):
        """ This function is a general function aiming at reading files.
        it uses the registered extensions to use the appropriate reading
        function.
        inputs:
            filename -- [ string | buffer ]
                        adress of the source file to be used
        keywords:
            type     -- [ string ]
                    if specified, this will force the function
                    to use considered the file to be of this
                    type.
            manager  -- [ Backend ]
                    If specified, it will use this format
                    manager (even if not registered)

            **kwargs are sent to the Backend.read function
        """
        if type is None:
            manager = determine_type(filename, verbose=not silent)
        else:
            manager = determine_type(type, verbose=False)

        t = manager().read(filename, **kwargs)
        self.data     = t.data
        self.header   = t.header
        self._aliases = t._aliases
        self.columns  = t.columns

    def write(self, filename, type=None, manager=None, silent=False, **kwargs):
        """ This function is a general function aiming at reading files.
        it uses the registered extensions to use the appropriate reading
        function.
        inputs:
            filename -- [ string | buffer ]
                        adress of the destination file to be used
        keywords:
            type     -- [ string ]
                    if specified, this will force the function
                    to use considered the file to be of this
                    type.
            manager  -- [ Backend ]
                    If specified, it will use this format
                    manager (even if not registered)

            **kwargs are sent to the Backend.read function
        """
        if type is None:
            manager = determine_type(filename, verbose=not silent)
        else:
            manager = determine_type(type, verbose=False)

        return manager().write(self, filename, **kwargs)

    def __set_from_table__(self, t):
        """
        Set defaults from another Table obj
        """
        for k, v in t.__dict__.items():
            self.__setattr__(k, v)

    def __set_defaults__(self):
        """
        Empty the table
        """
        self.header       = TableHeader()
        self.columns      = odict()
        self._aliases     = dict()
        self.data         = None
        self._primary_key = None
        self.caseless     = False

    def __call__(self, args=None):
        if args is None:
            return self.info()
        else:
            return self[args]

    def keys(self):
        return self.colnames

    @property
    def colnames(self):
        return list(self.columns.keys())

    @property
    def ncols(self):
        return len(self.colnames)

    @property
    def nrows(self):
        if self.ncols > 0:
            return self.data.shape[0]
        else:
            return 0

    @property
    def shape(self):
        return (self.nrows, self.ncols)

    @property
    def nbytes(self):
        """ return the number of bytes of the object """
        n = sum(k.nbytes if hasattr(k, 'nbytes') else sys.getsizeof(k) for k in list(self.__dict__.values()))
        return n

    def __len__(self):
        return self.nrows

    def set_name(self, name):
        """
        Set table's name
        name: [ string ] The table name

        """
        self.header['NAME'] = name

    def ravel(self, order='C'):
        """Return a flattened array.
            see np.ravel
        """
        return self.data.ravel(order=order)

    @property
    def dtype(self):
        return self.data.dtype

    def match(self, r2, key):
        """ Returns the indices at which the tables match
        matching uses 2 columns that are compared in values
        INPUTS:
            r2  [Table] second table to use
            key [str]   fields used for comparison.

        OUTPUS:
            tuple of both indices list where the two columns match.
        """
        return np.where( np.equal.outer( self[key], r2[key] ) )

    def join_by(self, r2, key, jointype='inner', r1postfix='1', r2postfix='2',
                defaults=None, asrecarray=False, asTable=True):
        """
        Join arrays `r1` and `r2` on key `key`.

        The key should be either a string or a sequence of string corresponding
        to the fields used to join the array.
        An exception is raised if the `key` field cannot be found in the two input
        arrays.
        Neither `r1` nor `r2` should have any duplicates along `key`: the presence
        of duplicates will make the output quite unreliable. Note that duplicates
        are not looked for by the algorithm.

        INPUTS:
            key        {str, seq}   A string or a sequence of strings
                        corresponding to the fields used for comparison.
            r2     [Table]  Table to join with

        KEYWORDS:
            jointype   [str]    {'inner', 'outer', 'leftouter'}
                'inner'     : returns the elements common to both r1 and r2.
                'outer'     : returns the common elements as well as the elements
                      of r1 not in r2 and the elements of not in r2.
                'leftouter' : returns the common elements and the elements of r1 not in r2.

            r1postfix  [str]    String appended to the names of the fields of r1 that are present in r2
            r2postfix  [str]    String appended to the names of the fields of r2 that are present in r1
            defaults   [dict]   Dictionary mapping field names to the corresponding default values.
            asrecarray [bool]   Whether to return a recarray or just a flexible-type ndarray.
            asTable    [bool]   Whether to return a Table (default).

        *Notes*:
        - The output is sorted along the key.
        - A temporary array is formed by dropping the fields not in the key for the
          two arrays and concatenating the result. This array is then sorted, and
          the common entries selected. The output is constructed by filling the fields
          with the selected entries. Matching is not preserved if there are some
          duplicates...

        """
        #TODO: return a Table by default
        if asTable:
            asrecarray = True
        arr = recfunctions.join_by(key, self, r2, jointype=jointype,
                                   r1postfix=r1postfix, r2postfix=r2postfix,
                                   defaults=defaults, usemask=False,
                                   asrecarray=asrecarray)

        return arr

    def dtype_for(self, cols):
        """
        Return the dtype of a row of a subset of columns.

        INPUTS:
            cols    [ iterable ] The subset of column names for which to return the dtype.

        OUTPUTS:
            dtype   [ np.dtype ]
        """
        return np.dtype([ (name, self.columns[name].dtype) for name in cols ])

    def tolist(self):
        return self.data.tolist()

    def set_alias(self, alias, colname):
        """
        Define an alias to a column

        INPUTS:
        alias   [ string ]   The new alias of the column
        colname [ string ]   The column being aliased
        """
        assert (colname in self.colnames), "Column %s does not exist" % colname
        self._aliases[alias] = colname

    def reverse_alias(self, colname):
        """
        Return aliases of a given column.

        Given a colname, return a sequence of aliases associated to this column
        Aliases are defined by using .define_alias()
        """
        _colname = self.resolve_alias(colname)
        # User aliases
        assert _colname in self.colnames

        return tuple([ k for k, v in self._aliases.items() if (v == _colname) ])

    def resolve_alias(self, colname):
        """
        Return the name of an aliased column.

        Given an alias, return the column name it aliases. This
        function is a no-op if the alias is a column name itself.

        Aliases are defined by using .define_alias()
        """
        # User aliases
        if hasattr(colname, '__iter__') and not isinstance(colname, basestring):
            return [ self.resolve_alias(k) for k in colname ]
        else:
            if self.caseless is True:
                maps = dict( [ (k.lower(), v) for k, v in list(self._aliases.items()) ] )
                maps.update( (k.lower(), k) for k in list(self.keys()) )
                return maps.get(colname.lower(), colname)
            else:
                return self._aliases.get(colname, colname)

    def add_empty_column(self, name, dtype, unit='', null='', description='', format=None,
                         shape=None, before=None, after=None, position=None, col_hdr=None):
        '''
        Add an empty column to the table. This only works if there
        are already existing columns in the table.

        INPUTS:
            name    [ string ]  The name of the column to add
            dtype   [ np.dtype ]    Numpy type of the column.

        KEYWORDS:
            unit    [ string ]  The unit of the values in the column
            null    [ datatype ]    The values corresponding to 'null', if not NaN
            description [ string ]  A description of the content of the column
            format  [ string ]  ASCII printing format
            shape   [ tuple ]   like shape of ndarray (nrow x ncols)
            before  [ string ]  Column before which the new column should be inserted
            after   [ string ]  Column after which the new column should be inserted
            position    [ integer ] Position at which the new column should be inserted (0 = first)

            col_hdr     [ ColumnHeader ] The metadata from an existing column to copy over.
            Metadata includes the unit, null value, description,
            format, and dtype.  For example, to create a column 'b'
            with identical metadata to column 'a' in table 't', use::
        
                t.add_column('b', column_header=t.columns[a])

        '''
        if shape:
            data = np.empty(shape, dtype=name2dtype(dtype))
        elif self.__len__() > 0:
            data = np.empty(self.__len__(), dtype=name2dtype(dtype))
        else:
            raise Exception("Table is empty, shape keyword is required for the first column")

        self.add_column(name, data, unit=unit, null=null,
                        description=description, format=format,
                        before=before, after=after, position=position,
                        col_hdr=col_hdr)

    def add_column(self, name, data, unit='', null='',
                   description='', format=None, dtype=None,
                   before=None, after=None, position=None, fill=None,
                   col_hdr=None ):
        """
        Add a column to the table

        INPUTS:
            name    [ string ]  The name of the column to add
            data    [ np.ndarray ]  The column data

        KEYWORDS:
            unit    [ string ]  The unit of the values in the column
            null    [ datatype ]    The values corresponding to 'null', if not NaN
            description [ string ]  A description of the content of the column
            format  [ string ]  ASCII printing format
            dtype   [ np.dtype ]    Numpy type of the column.

            before  [ string ]  Column before which the new column should be inserted
            after   [ string ]  Column after which the new column should be inserted
            position    [ integer ] Position at which the new column should be inserted (0 = first)

            col_hdr     [ ColumnHeader ] The metadata from an existing column to copy over.
            Metadata includes the unit, null value, description,
            format, and dtype.  For example, to create a column 'b'
            with identical metadata to column 'a' in table 't', use::

                t.add_column('b', column_header=t.columns[a])

        """

        _data = np.asarray(data)

        if col_hdr is not None:
            """ Priority to col_hdr """
            dtype = column_header.dtype
            unit = column_header.unit
            null = column_header.null
            description = column_header.description
            format = column_header.format
        else:
            dtype = _data.dtype

        # unknown type is converted to text
        if dtype.type == np.object_:
            if len(data) == 0:
                longest = 0
            else:
                longest = len(max(data, key=len))
                data = np.asarray(data, dtype='|%iS' % longest)

        dtype = data.dtype

        if data.ndim > 1:
            newdtype = (str(name), data.dtype, (data.shape[1],))
        else:
            newdtype = (str(name), data.dtype)

        # get position
        if before:
            try:
                position = list(self.colnames).index(before)
            except:
                raise Exception("Column %s does not exist" % before)
        elif after:
            try:
                position = list(self.colnames).index(after) + 1
            except:
                raise Exception("Column %s does not exist" % before)

        if len(self.columns) > 0:
            self.data = append_field(self.data, data, dtype=newdtype, position=position)
        else:
            self.data = np.array(data, dtype=[newdtype])

        if not format or format in ['e', 'g', 'f']:
            format = default_format[dtype.type]

        # Backward compatibility with tuple-style format
        if type(format) in [tuple, list]:
            format = string.join([str(x) for x in format], "")

        if format == 's':
            format = '%is' % data.itemsize

        column = ColumnHeader(dtype, unit=unit, description=description, null=null, format=format)

        if not np.equal(position, None):
            self.columns.insert(position, name, column)
        else:
            self.columns[name] = column

    def sort(self, keys):
        """
        Sort the table according to one or more keys. This operates
        on the existing table (and does not return a new table).

        INPUTS:
            keys: [ string | list of strings ]   The key(s) to order by
        """
        if not type(keys) == list:
            keys = [keys]
        self.data.sort(order=keys)

    @property
    def empty_row(self):
        """ Return an empty row array respecting the table format """
        return np.rec.recarray(shape=(1,), dtype=self.data.dtype)

    def append_row(self, iterable):
        """
        Append set of rows in this table.
        """
        assert( len(iterable) == self.ncols ), 'Expecting as many items as columns'
        r = self.empty_row
        for k, v in enumerate(iterable):
            r[0][k] = v
        self.stack(r)

    def addLine(self, iterable):
        """
        Append set of rows in this table.
        """
        self.append_row(iterable)

    def stack(self, r, defaults=None):
        """
        Superposes arrays fields by fields
        """
        self.data = recfunctions.stack_arrays( [self.data, r], defaults, usemask=False, asrecarray=True)

    def addCol(self, name, data, unit='', null='',
               description='', format=None, dtype=None,
               before=None, after=None, position=None, fill=None,
               col_hdr=None ):
        """ Add individual column to the table
            See add_column
        """
        self.add_column(name, data, unit=unit, null=null,
                        description=description, format=format, dtype=dtype,
                        before=before, after=after, position=position, fill=fill,
                        col_hdr=col_hdr)

    def delCol(self, name):
        """ Delete Table column
        INPUTS:
            name [ string ] Column to delete
        """
        self.remove_columns([name])

    def remove_column(self, name):
        """ Delete Table column
        INPUTS:
            name [ string ] Column to delete
        """
        self.remove_columns([name])

    def remove_columns(self, names):
        """
        Remove several columns from the table
        INPUTS:
            names   [ list ] A list containing the names of the columns to remove
        """
        self.pop_columns(names)

    def pop_columns(self, names):
        """
        Pop several columns from the table
        INPUTS:
            names   [ list ] A list containing the names of the columns to remove
        """

        if type(names) == str:
            names = [names]

        _names = []
        for k in names:
            _k = self.resolve_alias(k)
            if (_k in self._aliases):
                self._aliases.pop(_k)
            else:
                list(map(self._aliases.pop, self.reverse_alias(_k)))
                _names.append(_k)

        p = [self.columns.pop(name) for name in _names]
        self.data = drop_fields(self.data, _names)

        # Remove primary key if needed
        if self._primary_key in _names:
            self._primary_key = None

        return p

    def find_duplicate(self, index_only=False, values_only=False):
        """Find duplication in the table entries, return a list of duplicated elements
            Only works at this time is 2 lines are *the same entry*
            not if 2 lines have *the same values*
        """
        dup = []
        idd = []
        for i in range(len(self.data)):
            if (self.data[i] in self.data[i + 1:]):
                if (self.data[i] not in dup):
                    dup.append(self.data[i])
                    idd.append(i)
        if index_only:
            return idd
        elif values_only:
            return dup
        else:
            return list(zip(idd, dup))

    def __getitem__(self, v):
        return np.asarray(self.data.__getitem__(self.resolve_alias(v)))

    def __setitem__(self, v):
        return self.data.__setitem__(self.resolve_alias(v))

    def __getattr__(self, k):
        try:
            return object.__getattribute__(self, k)
        except:
            try:
                return self[k]
            except:
                raise AttributeError('Attribute {} not found'.format(k))

    def __deepcopy__(self, memo):
        return copyTable(self)

    def __copy__(self, memo):
        return self

    def __pretty_print__(self, idx=None, fields=None, ret=False):
        """ Pretty print the table content
            you can select the table parts to display using idx to
            select the rows and fields to only display some columns
            (ret is only for insternal use)"""

        if fields is None:
            fields = list(self.keys())
        if isinstance(fields, str):
            fields = fields.split(',')

        nfields = len(fields)

        fields = list(map( self.resolve_alias, fields ))

        if idx is None:
            if self.nrows < 10:
                rows = [ [ str(self[k][rk]) for k in fields ] for rk in range(self.nrows)]
            else:
                _idx = list(range(6))
                rows = [ [ str(self[k][rk]) for k in fields ] for rk in range(5) ]
                if nfields > 1:
                    rows += [ ['...' for k in range(len(fields)) ] ]
                else:
                    rows += [ ['...' for k in range(len(fields)) ] ]
                rows += [ [ str(self[k][rk]) for k in fields ] for rk in range(-5, 0)]
        elif isinstance(idx, slice):
            _idx = list(range(idx.start, idx.stop, idx.step or 1))
            rows = [ [ str(self[k][rk]) for k in fields ] for rk in _idx]
        else:
            rows = [ [ str(self[k][rk]) for k in fields ] for rk in idx]
        units = [ '(' + str( self.columns[k].unit or '') + ')' for k in fields ]
        #fmt   = [ '%' + self.columns[k].format for k in self.keys() ]

        #if nfields == 1:
        #   fields = [fields]
        if (''.join(units) == ''.join(['()'] * len(fields)) ):
            out = __indent__([fields] + rows, hasHeader=True, hasUnits=False)
        else:
            out = __indent__([fields] + [units] + rows, hasHeader=True, hasUnits=True)
        if ret is True:
            return out
        else:
            print(out)

    def __str__(self):
        return self.__pretty_print__(ret=True)

    def __repr__(self):
        s = object.__repr__(self)
        s += '\nTable: %s,\n  nrows=%i, ncols=%i (%s)' % (self.header['NAME'], self.nrows, self.ncols, pretty_size_print(self.nbytes))
        return s

    def __getslice__(self, i, j):
        return self.data.__getslice__(i, j)

    def __contains__(self, k):
        return (k in list(self.keys())) or (k in self._aliases)

    def __iter__(self):
        return self.data.__iter__()

    def iterkeys(self):
        return iter(self.columns.keys())

    def itervalues(self):
        return iter(self.data.values())

    def info(self):
        print(self.header)
        print("Table contains: %i row(s) in %i column(s)\n" % (self.nrows, self.ncols))
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

    def evalexpr(self, expr, exprvars=None, start=None, stop=None, step=None, dtype=float):
        """ evaluate expression based on the data and external variables
            all np function can be used (log, exp, pi...)
        """
        _globals = {}
        for k in ( list(self.keys()) + list(self._aliases.keys()) ):
            _globals[k] = self[k]

        if exprvars is not None:
            assert(hasattr(exprvars, 'keys') & hasattr(exprvars, '__getitem__' )), "Expecting a dictionary-like as condvars"
            for k, v in ( list(exprvars.items()) ):
                _globals[k] = v

        # evaluate expression, to obtain the final filter
        r    = np.empty( self.nrows, dtype=dtype)
        r[:] = eval(expr, _globals, np.__dict__)

        return r

    def where(self, condition, condvars=None, start=None, stop=None, step=None, *args):
        """ Read table data fulfilling the given `condition`.
            Only the rows fulfilling the `condition` are included in the result.
            INPUTS:
                condition   ndarray[dtype=bool]

            OUTPUTS:
                out     ndarray/ tuple of ndarrays

            Additional arguments are forwarded to np.where
        """
        ind = np.where(self.evalexpr(condition, condvars, start=start, stop=stop, step=step, dtype=bool ), *args)
        return ind

    def selectWhere(self, fields, condition, condvars=None, **kwargs):
        """ Read table data fulfilling the given `condition`.
            Only the rows fulfilling the `condition` are included in the result.
        """
        # make a copy without the data itself (memory gentle)
        tab = self.__class__()
        for k in list(self.__dict__.keys()):
            if k != 'data':
                setattr(tab, k, deepcopy(self.__dict__[k]))

        if fields.count(',') > 0:
            _fields = fields.split(',')
        elif fields.count(' ') > 0:
            _fields = fields.split()
        else:
            _fields = fields

        if condition in [True, 'True', None]:
            ind = None
        else:
            ind = self.where(condition, condvars, **kwargs)

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

        tab.header['COMMENT'] = 'SELECT %s FROM %s WHERE %s' % (','.join(_fields), self.header['NAME'], condition)
        return tab

    def setUnit(self, colName, unit):
        """ Set the unit of a column referenced by its name """
        self.columns[self.resolve_alias(colName)].unit = unit

    def setComment(self, colName, comment):
        """ Set the comment of a column referenced by its name """
        self.columns[self.resolve_alias(colName)].description = comment

    def setNull(self, colName, null):
        """ Set the comment of a column referenced by its name """
        self.columns[self.resolve_alias(colName)].null = null

    def setFormat(self, colName, fmt):
        """ Set the comment of a column referenced by its name """
        self.columns[self.resolve_alias(colName)].format = fmt

    def has_key(self, k):
        return k in self


def __indent__(rows, hasHeader=False, hasUnits=False, headerChar='-', delim=' | ', justify='left',
               separateRows=False, prefix='', postfix='', wrapfunc=lambda x: x):
    """Indents a table by column.

    INPUTS:
    rows    [sequence or sequences of items]    one sequence per row.

    KEYWORDS:
    hasHeader   [bool]      True if the first row consists of the columns' names
    headerChar  [char]      Character to be used for the row separator line
    delim       [char]      The column delimiter.
    justify     [str]       { center | right | left } data justified in their column.
    separateRows    [bool]      True if rows are to be separated by a line of 'headerChar's.
    prefix      [str]       A string prepended to each printed row.
    postfix     [str]       A string appended to each printed row.
    wrapfunc    [function]  A function f(text) for wrapping text; each element in the
                    table is first wrapped by this function.
    """
    # closure for breaking logical rows to physical, using wrapfunc
    def rowWrapper(row):
        newRows = [wrapfunc(item).split('\n') for item in row]
        return [[substr or '' for substr in item] for item in map(None, *newRows)]
    # break each logical row into one or more physical ones
    logicalRows = [rowWrapper(row) for row in rows]
    # columns of physical rows
    columns = map(None, *reduce(operator.add, logicalRows))
    # get the maximum of each column by the string length of its items
    maxWidths = [max([len(str(item)) for item in column]) for column in columns]
    rowSeparator = headerChar * (len(prefix) + len(postfix) + sum(maxWidths) + len(delim) * (len(maxWidths) - 1))

    # select the appropriate justify method
    justify = {'center': str.center, 'right': str.rjust, 'left': str.ljust}[justify.lower()]
    output = io.StringIO()
    if separateRows:
        print(rowSeparator, file=output)
    for physicalRows in logicalRows:
        for row in physicalRows:
            print(prefix \
                + delim.join([justify(str(item), width) for (item, width) in zip(row, maxWidths)]) \
                + postfix, file=output)
        if separateRows:
            print(rowSeparator, file=output)
        elif hasHeader & (not hasUnits):
            print(rowSeparator, file=output)
        elif (not hasHeader) & hasUnits:
            print(rowSeparator, file=output)
            hasUnits = False
        hasHeader = False

    return output.getvalue()


#==============================================================================
def from_dict(obj, **kwargs):
    """ Generate a table from a recArray or numpy array """
#==============================================================================
    assert( hasattr(obj, 'iteritems') or hasattr(obj, 'items') ), "expecting obj has iteritem attribute (dict-like)"

    tab = Table()
    for k, v in obj.items():
            _v = np.asarray(v)
            tab.add_column( k, _v, dtype=_v.dtype )
    return tab


#==============================================================================
def from_ndArray(obj, **kwargs):
    """ Generate a table from a recArray or numpy array """
#==============================================================================
    supported_types = [ np.core.records.recarray, np.ndarray ]
    assert( type(obj) in supported_types ), "expecting recArray, got %s " % type(obj)
    tab = Table()

    if obj.dtype.names is not None:
        for k in obj.dtype.names:
            tab.add_column( k, obj[k], dtype=obj.dtype[k] )
    else:
        for i in range(obj.shape[1]):
            tab.add_column( 'f%d' % i, obj[:, i], dtype=obj.dtype )

    return tab


#==============================================================================
def from_Table(obj, **kwargs):
    """ Generate a new table fr om a Table obj """
#==============================================================================
    assert( isinstance(obj, Table)), "Expecting Table object, got %s" % type(obj)
    return copyTable(obj)


#==============================================================================
def copyTable(obj, **kwargs):
    """ Copy a Table """
#==============================================================================
    t = obj.__class__()
    for k in list(obj.__dict__.keys()):
        setattr(t, k, deepcopy(obj.__dict__[k]))
    return (t)


@warning(message='Deprecated function. Use Table class constructor instead.')
def load(*args, **kwargs):
    return Table(*args, **kwargs)
