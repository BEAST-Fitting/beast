""" Backend for HDF5 format based on (py)tables """


import os
import inspect
localpath = '/'.join(os.path.abspath(inspect.getfile(inspect.currentframe())).split('/')[:-1])
import numpy as np
from .basebackend import BaseBackend
#from ..core.tableheader import TableHeader
from ..core import helpers
from ..table import Table

import numexpr
numexpr.set_num_threads(1)
import tables

typeDict = {}
typeDict[np.bool_]      = tables.BoolCol
typeDict[np.short]      = tables.Int8Col
typeDict[np.int8]       = tables.Int8Col
typeDict[np.uint8]      = tables.UInt8Col
typeDict[np.int16]      = tables.Int16Col
typeDict[np.uint16]     = tables.UInt16Col
typeDict[np.int32]      = tables.Int32Col
typeDict[np.uint32]     = tables.UInt32Col
typeDict[np.int64]      = tables.Int64Col
typeDict[np.uint64]     = tables.UInt64Col
typeDict[np.float16]    = tables.Float32Col
typeDict[np.float32]    = tables.Float32Col
typeDict[np.float64]    = tables.Float64Col
typeDict[np.complex64]  = tables.Complex64Col
typeDict[np.complex128] = tables.Complex128Col
typeDict[np.str]        = tables.StringCol
typeDict[np.string_]    = tables.StringCol
typeDict[str]           = tables.StringCol
typeDict[np.unicode_]   = tables.StringCol


#==============================================================================
class hdf5Backend(BaseBackend):
#==============================================================================
    def __init__(self):
        """ constructor """
        BaseBackend.__init__(self, tableType='hdf5')

    def descr_from_dtype(self, dtype_, shape_={}):
        """
        Get a description instance and byteorder from a (nested) NumPy dtype.

        INPUTS:

            dtype_  dtype   dtype of a ndarray
            shape_  dict    a dictionary of column's shapes

        OUTPUTS:
            Description of a new table
        """
        fields = {}
        fbyteorder = '|'
        for (name, (dtype, pos)) in list(dtype_.fields.items()):
            kind = dtype.base.kind
            shape = shape_.get(name, ())
            byteorder = dtype.base.byteorder
            if byteorder in '><=':
                if fbyteorder not in ['|', byteorder]:
                    raise NotImplementedError(
                        "record arrays with mixed byteorders "
                        "are not supported yet, sorry" )
                fbyteorder = byteorder
            # Non-nested column
            if kind in 'biufSc':
                #col = tables.Col.from_dtype(dtype, pos=pos)
                if len(shape) > 1:
                    col = typeDict[dtype.type](shape=shape)
                else:
                    col = typeDict[dtype.type]()
            # Nested column
            elif kind == 'V' and dtype.shape in [(), (1,)]:
                col, _ = self.descr_from_dtype(dtype)
                col._v_pos = pos
            else:
                raise NotImplementedError(
                    "record arrays with columns with type description ``%s`` "
                    "are not supported yet, sorry" % dtype )
            fields[name] = col

        return tables.Description(fields)

    def __createTable__(self, hd5, tab, group, tablename):
            shapes = {}
            for k in tab.colnames:
                shapes[k] = tab[k].shape
            desc = self.descr_from_dtype(tab.dtype, shapes)
            if group[-1] == '/':
                group = group[:-1]

            t = hd5.create_table(group, tablename, desc,
                    expectedrows=tab.nrows, createparents=True)
            hd5.flush()
            return t

    def write(self, tab, output='exportedData.hd5', group='/',
            mode='w', tablename=None,
            append=False, silent=False, keep_open=False):

        if hasattr(output, 'read'):
            raise Exception("HDF backend does not implement stream")

        tablename = tablename or tab.header['NAME']
        if tablename in ['', None, 'Noname', 'None']:
            tablename = 'data'

        if append:
            mode = 'a'

        hd5 = tables.open_file(output, mode=mode)

        #generate the empty table
        if not append:
            t = self.__createTable__(hd5, tab, group, tablename)
            try:
                t = self.__createTable__(hd5, tab, group, tablename)
            except:
                if not silent:
                    print("Warning: Table creation exception. Table may already exist.")
                t = hd5.get_node(group + tablename)
        else:
            try:
                t = hd5.get_node(group + tablename)
            except tables.NoSuchNodeError:
                if not silent:
                    print("Warning: Table does not exists.  New table will be created")
                t = self.__createTable__(hd5, tab, group, tablename)

        #fill the table
        dtype_wanted = t.description._v_dtype

        # Reorder the columns so that the data set can be appended correctly
        data_in_correct_order = np.zeros(len(tab.data), dtype_wanted)
        for k in dtype_wanted.names:
            data_in_correct_order[k] = tab.data[k]

        t.append(data_in_correct_order)
        hd5.flush()

        #update the header
        for k, v in tab.header.items():
            if k == 'FILTERS' and float(t.attrs['VERSION']) >= 2.0:
                t.attrs['filters'] = v
            else:
                t.attrs[k] = v
        if not 'TITLE' in tab.header:
            t.attrs['TITLE'] = tablename

        #add Column descriptions
        for colname, col in tab.columns.items():
            key = [ k for k in t.attrs._v_attrnames if t.attrs[k] == colname ]

            key = key[0].replace('NAME', '')
            t.attrs[key + 'UNIT'] = col.unit
            t.attrs[key + 'DESC'] = col.description
            t.attrs[key + 'FMT']  = col.format
            t.attrs[key + 'NULL'] = col.null

        #add aliases
        for i, (k, v) in enumerate(tab._aliases.items()):
            t.attrs['ALIAS%d' % i ] = '%s=%s' % (k, v)

        print(hd5)
        hd5.close()

    def readColDesc(self, tab, name):
        key = [ k for k in tab.attrs._v_attrnames if tab.attrs[k] == name ]

        key = key[0].replace('NAME', '')
        attrnames = 'desc unit fmt null'.split()
        attrmap   = 'description unit format null'.split()
        meta = {}
        for il, l in enumerate(attrnames):
            meta[ attrmap[il] ] = tab.attrs.__dict__.get( (key + l).upper(), '')

        return meta

    def read(self, filename, tableName=None, silent=False, *args, **kwargs):

        source = tables.open_file(filename, *args, **kwargs)

        if 'tablename' in kwargs:
            tableName = kwargs['tablename']
        if tableName is None:
            node = source.list_nodes('/')[0]
            tableName = node.name
        else:
            if tableName[0] != '/':
                tableName = '/' + tableName
            node = source.get_node(tableName)
        if silent is True:
            print("\tLoading table: %s" % tableName)

        #read data
        data = node[:]

        #generate a Table
        tab = Table(data)

        #update column meta
        for colname, colhdr in tab.columns.items():
            meta = self.readColDesc(node, colname)
            for mk, mv in meta.items():
                if mk != 'null':
                    if (mk == 'format') & (mv in ['', None, 'None']):
                        mv = helpers.default_format.get(tab.dtype[colname], '')
                    colhdr.__setattr__(mk, mv)

        #update header & aliases
        exclude = ['NROWS', 'VERSION', 'CLASS', 'EXTNAME']
        for k in node.attrs._v_attrnames:
            if (not k in exclude) & (k[:5] != 'FIELD') & (k[:5] != 'ALIAS'):
                tab.header[k] = node.attrs[k]
            if (k[:5] == 'ALIAS'):
                c0, c1 = node.attrs[k].split('=')
                tab.set_alias(c0, c1)

        empty_name = ['', 'None', 'Noname', None]
        if (tab.header['NAME'] in empty_name) & (tab.header.__dict__.get('TITLE', None) not in empty_name):
            tab.header['NAME'] = tab.header['TITLE']

        source.close()

        return tab
