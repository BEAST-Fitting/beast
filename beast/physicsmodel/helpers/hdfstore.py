""" HDFStore
High level interface to PyTables for reading and writing structures to disk
Could work either as a (key, value) storage or as a classic tables.file.File object.

Example Usage::
    import numpy as np
    #make a store
    with HDFStore('tmp.hd5', mode='w') as hd:
        #make some variables
        d = {}
        d['a'] = np.arange(10, dtype=float)
        d['b'] = np.arange(10, dtype='int')
        c = np.random.normal(0, 1, (10, 10))
        d['c'] = c
        #put values into the store
        hd['/subdir/table'] = d
        hd['/subdir1/array'] = c
        #check values
        print hd.keys()
        print hd['subdir/table']
        print hd['subdir1/array']
        hd.remove_node('/subdir1', recursive=True)
"""


import tables
import numpy as np

if int(tables.__version__[0]) < 3:
    _future = True
    from future.tables_file import create_earray
else:
    _future = False

__all__ = ['HDFStore', 'convert_dict_to_structured_ndarray']

# Type conversion mapping between numpy and pytables
_typeDict = {}
_typeDict[np.bool_]      = tables.BoolCol
_typeDict[np.short]      = tables.Int8Col
_typeDict[np.int8]       = tables.Int8Col
_typeDict[np.uint8]      = tables.UInt8Col
_typeDict[np.int16]      = tables.Int16Col
_typeDict[np.uint16]     = tables.UInt16Col
_typeDict[np.int32]      = tables.Int32Col
_typeDict[np.uint32]     = tables.UInt32Col
_typeDict[np.int64]      = tables.Int64Col
_typeDict[np.uint64]     = tables.UInt64Col
_typeDict[np.float16]    = tables.Float32Col
_typeDict[np.float32]    = tables.Float32Col
_typeDict[np.float64]    = tables.Float64Col
_typeDict[np.complex64]  = tables.Complex64Col
_typeDict[np.complex128] = tables.Complex128Col
_typeDict[np.str]        = tables.StringCol
_typeDict[np.string_]    = tables.StringCol
_typeDict[str]           = tables.StringCol
_typeDict[np.unicode_]   = tables.StringCol


def __descr_from_dtype__(dtype_, shape_={}):
    """
    Get a description instance and byteorder from a (nested) NumPy dtype.

    Parameters
    ----------
    dtype_:  dtype
        dtype of a ndarray

    shape_:  dict
        a dictionary of column's shapes

    returns
    -------
    r: tables.Description
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
                col = _typeDict[dtype.base.type](shape=shape)
            else:
                col = _typeDict[dtype.type]()
        # Nested column
        elif kind == 'V' and dtype.shape in [(), (1,)]:
            col, _ = __descr_from_dtype__(dtype)
            col._v_pos = pos
        else:
            raise NotImplementedError(
                "record arrays with columns with type description ``%s`` "
                "are not supported yet, sorry" % dtype )
        fields[name] = col

    return tables.Description(fields)


def __createTable__(hd5, tab, group, tablename):
    """__createTable__ -- make a new table in the hdf file

    Parameters
    ----------

    hd5: tables.file.File
        opened hdf file

    tab: structured ndarray
        structure that contains individual columns

    group: str
        path to where the table will be created

    tablename: str
        name of the node for the table

    returns
    -------
    t: tables.table.Table
        table node instance already set into the file

    """
    shapes = {}
    for k in tab.dtype.names:
        shapes[k] = tab[k].shape
    desc = __descr_from_dtype__(tab.dtype, shapes)
    if group[-1] == '/':
        group = group[:-1]

    t = hd5.create_table(group, tablename, desc, expectedrows=len(tab), createparents=True)
    hd5.flush()
    return t


def convert_dict_to_structured_ndarray(data):
    """convert_dict_to_structured_ndarray

    Parameters
    ----------

    data: dictionary like object
        data structure which provides iteritems and itervalues

    returns
    -------
    tab: structured ndarray
        structured numpy array
    """
    newdtype = []
    for key, dk in data.items():
        _dk = np.asarray(dk)
        dtype = _dk.dtype
        # unknown type is converted to text
        if dtype.type == np.object_:
            if len(data) == 0:
                longest = 0
            else:
                longest = len(max(_dk, key=len))
                _dk = _dk.astype('|%iS' % longest)
        if _dk.ndim > 1:
            newdtype.append((str(key), _dk.dtype, (_dk.shape[1],)))
        else:
            newdtype.append((str(key), _dk.dtype))
    tab = np.rec.fromarrays(iter(data.values()), dtype=newdtype)
    return tab


def __write__(data, hd5, group='/', tablename=None, append=False, silent=False, header={}, **kwargs):
    """__write__

    Parameters
    ----------
    data: dict like object
        structure that contains individual columns

    hd5: tables.file.File
        opened hdf file

    group: str
        path to where the table will be created

    tablename: str
        name of the node for the table

    append: bool
        if set, attempts to add content at the end of an existing node
        Only works for Tables (not Arrays)

    silent: bool
        if set, do not print information content

    header: dictionary like
        attributes to add to the node
    """
    if tablename in ['', None, 'Noname', 'None']:
        tablename = 'data'

    #fill the table
    if not hasattr(data, 'dtype'):
        #assumes dict like object
        tab = convert_dict_to_structured_ndarray(data)
    else:
        tab = data

    if tab.dtype.names is None:
        if not append:
            if _future:
                create_earray(hd5, group, tablename, obj=tab, createparents=True)
            else:
                hd5.create_earray(group, tablename, obj=tab, createparents=True)
            hd5.flush()
            t = hd5.get_node(group + '/' + tablename)
        else:
            t = hd5.get_node(group + '/' + tablename)
            t.append(tab)
    else:
        #generate the empty table
        if not append:
            t = __createTable__(hd5, tab, group, tablename)
        else:
            try:
                t = hd5.get_node(group + '/' + tablename)
            except tables.NoSuchNodeError:
                if not silent:
                    print("Warning: Table does not exists.  New table will be created")
                t = __createTable__(hd5, tab, group, tablename)
                try:
                    t = __createTable__(hd5, tab, group, tablename)
                except:
                    if not silent:
                        print("Warning: Table creation exception. Table may already exist.")
                    t = hd5.get_node(group + '/' + tablename)

        t.append( tab.astype(t.description._v_dtype) )
        hd5.flush()

    #update the header
    for k, v in header.items():
        if k == 'FILTERS' and float(t.attrs['VERSION']) >= 2.0:
            t.attrs['filters'] = v
        else:
            t.attrs[k] = v
    if not 'TITLE' in header:
        t.attrs['TITLE'] = tablename


class HDFStore(object):
    """ Handling quick in and out of the HDF5 file
        This class can be used as a context manager

        Any attribute of the HDF will be directly available transparently if
        the source if opened
    """
    def __init__(self, lnpfile, mode='r', **kwargs):
        """__init__

        Parameters
        ----------

        lnpfile: str or tables.file.File
            storage to use

        mode: str ('r', 'w', 'a', 'r+')
            The mode to open the file.
            (see set_mode)

        **kwargs can be use to set options to tables.open_file
        """
        self.lnpfile = lnpfile
        self.source = None
        self._keep_open = False
        self._open_kwargs = {}
        self._mode = None
        self.set_mode(mode, **kwargs)

    def set_mode(self, val=None, **kwargs):
        """set_mode - set a flag to open the file in a given mode operations
        if the mode changes, and the storage already opened, it will close the
        storage and reopen it in the new mode

        Parameters
        ----------

        val: str (optional)
            The mode to open the file.  It can be one of the following:
            'r'  -- Read-only; no data can be modified.
            'w'  -- Write; a new file is created (existing file would be deleted)
            'a'  -- Append; an existing file is opened for reading and writing,
                    and if the file does not exist it is created.
            'r+' -- It is similar to 'a', but the file must already exist.

        **kwargs is forwarded to tables.open_file

        returns
        -------
        status: str
            return the status only if None was provided
        """
        if val is None:
            return self._mode
        else:
            title = kwargs.get('title', None) or self._open_kwargs.get('title', '')
            root_uep = kwargs.get('root_uep', None) or self._open_kwargs.get('root_uep', '/')
            filters = kwargs.get('filters', None) or self._open_kwargs.get('filters', None)
            self._open_kwargs.update(dict(title=title, root_uep=root_uep, filters=filters))
            if (val != self._mode):
                self._mode = str(val)
                if self.source is not None:
                    self.close_source(force=True)
                    self.open_source()

    def keep_open(self, val=None):
        """keep_open - set a flag to keep the storage file open for multiple operations

        Parameters
        ----------

        val: bool (optional)
            if set, set the flags to apply

        returns
        -------
        status: bool
            return the status only if None was provided
        """
        if val is None:
            return self._keep_open
        else:
            self._keep_open = bool(val)

    def open_source(self):
        """open_source
        Open the file if needed, handles that the object was initialized with an opened file

        returns
        -------
        v: bool
            returns if any operation was actually done (True) else the file was already opened (False)
        """
        if self.source is None:
            if type(self.lnpfile) == str:
                self.source = tables.open_file(self.lnpfile, mode=self._mode, **self._open_kwargs)
            elif isinstance(self.lnpfile, tables.file.File):
                self.source = self.lnpfile
            return True
        return False

    def close_source(self, force=False):
        """close_source
        close the file only if the object opened any file

        Parameters
        ----------
        force: bool (default: False)
            if set, bypass keep_open status
        """
        if (force is True) or ((self.source is not None) and (not self.keep_open())):
            if isinstance(self.source, tables.file.File):
                self.source.close()
                self.source = None

    def keys(self):
        """
        Return a (potentially unordered) list of the keys corresponding to the
        objects stored in the HDFStore. These are ABSOLUTE path-names (e.g. have the leading '/'
        """
        with self as s:
            r = list(k._v_pathname for k in s.source.walk_nodes())
        return r

    def items(self):
        """
        key->group
        """
        return list(self.items())

    def iteritems(self):
        """
        iterate on key->group
        """
        with self as s:
            for g in s.groups():
                yield g._v_pathname, g

    def groups(self, where='/'):
        """ return a list of all the top-level nodes (that are not themselves a pandas storage object) """
        with self as s:
            g = list(s.source.walk_groups(where=where))
        return g

    def write(self, data, group='/', tablename=None, append=False, silent=False, header={}, **kwargs):
        """write -- write a data structure into the source file

        Parameters
        ----------
        data: dict like object
            structure that contains individual columns

        hd5: tables.file.File
            opened hdf file

        group: str
            path to where the table will be created

        tablename: str
            name of the node for the table

        append: bool
            if set, attempts to add content at the end of an existing node
            Only works for Tables (not Arrays)

        silent: bool
            if set, do not print information content

        header: dictionary like
            attributes to add to the node
        """
        with self as s:
            if s._mode[0] == 'r':
                s.set_mode('a')
            __write__(data, s.source, group, tablename, append, silent, header, **kwargs)

    def __getitem__(self, key):
        """__getitem__
        Returns the node corresponding to the key
        """
        with self as s:
            if key[0] != '/':
                p = '/' + key
            else:
                p = key
            if p not in list(s.keys()):
                raise KeyError('No object named {0} in this store'.format(p))
            return s.source.get_node(p)

    def __setitem__(self, key, value):
        """__setitem__
        create a node with the key path/name and set value as its content
        """
        with self as s:
            if key in list(s.keys()):
                raise KeyError('Object named {0} already exists in this store'.format(key))
            p = key.split('/')
            name = p[-1]
            p = '/'.join(p[:-1])
            if (p is None) or (p == ''):
                p = '/'
            if p[0] != '/':
                p = '/' + p
            self.write(value, group=p, tablename=name, silent=True)

    def __getattr__(self, name):
        """__getattr__
        Give direct access to any function from self.source (tables.file.File)
        """
        if name in self.__dict__:
            return self.__dict__[name]
        elif hasattr(self.source, name):
            return getattr(self.source, name)
        else:
            msg = "'{0}' object has no attribute '{1}'"
            raise AttributeError(msg.format(type(self).__name__, name))

    def __repr__(self):
        return object.__repr__(self)

    def __enter__(self):
        """ enter context """
        if self.open_source():
            self._context_close = [self.close_source]
        else:
            self._context_close = []
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """ close context """
        for k in self._context_close:
            k()

    def __del__(self):
        self.close_source()

    def __contains__(self, key):
        """ check for existance of this key
        can match the exact pathname or the pathnm w/o the leading '/'
        """
        return key in list(self.keys())

    def __len__(self):
        return len(self.groups())

    def get_Q_from_node(self, nodename, expr, condvars={}, coordinates=None):
        """ returns a quantity from a HDF5 Node given its math expression.
            Assuming that all quantities are either from the node or in condvars

            attempt to be clever and optimizing speed when coordinates are provided

            all np function can be used (log, exp, pi...)

            method
            ------
            Pure eval expression was too slow when using disk access and HD5
            nodes.  Instead, we use the tables.expression feature that parse
            "expr" and extracts the variable names for interest.  Based on this
            list we build the context of evaluation by adding missing values
            from the node.

            Parameters
            ----------
            nodename: str or ables.table.Table
                the Table node to get data from (need named columns)

            expr: str
                the mathematical expression which could include numpy functions or
                simply be a column name

            condvars: dict
                The condvars mapping may be used to define the variable names
                appearing in the condition.

            coordinates: iterable of int
                reduces the computations to only set of rows given their indexes
                Only works with nodes of tables.table.Table type

            returns
            -------
            q: ndarray like
                the column values requested with type and shape defined in the table.
        """
        with self as storage:
            node = storage.source.get_node(nodename)
            colnames = node.colnames
            if expr in colnames:
                if hasattr(node, 'read_coordinates') & (coordinates is not None):
                    q = node.read_coordinates(coordinates, field=expr)
                else:
                    q = node.read(field=expr)
            else:
                #parse qname, using tables.Expr does not help with log10 etc
                names = tables.expression.getExprNames(expr, {})[0]
                if hasattr(node, 'read_coordinates') & (coordinates is not None):
                    for k in names:
                        if k in colnames:
                            condvars[k] = node.read_coordinates(coordinates, field=k)
                elif (coordinates is not None):
                    _tmp = node[coordinates]
                    for k in names:
                        if k in colnames:
                            condvars[k] = _tmp[k]
                else:
                    for k in names:
                        if k in colnames:
                            condvars[k] = node.read(field=k)
                q = eval(expr, condvars, np.__dict__)
        return q


def unittest():
    """unittest
    Example usage
    """
    import numpy as np

    #make a store
    with HDFStore('tmp.hd5', mode='w') as hd:
        #make some variables
        d = {}
        d['a'] = np.arange(10, dtype=float)
        d['b'] = np.arange(10, dtype='int')
        c1 = np.random.normal(0, 1, (10, 10))
        c2 = np.random.normal(0, 1, (10, 10))
        d['c'] = c1

        #put values into the store
        hd['/subdir/table'] = d
        hd['/subdir1/array'] = c1
        hd.write(c2, '/subdir1', 'array', append=True)

        #check values
        print(list(hd.keys()))
        print(hd['subdir/table'])
        print(hd['subdir1/array'])
        hd.remove_node('/subdir1', recursive=True)
        print(list(hd.keys()))
