""" HDFStore
High level interface to PyTables for reading and writing structures to disk

TODO: easy writing is not their yet
TODO: get node from a __getitem__ call
TODO: set node from a __setitem__ call
"""

import tables
import numpy as np


class HDFStore(object):
    """ Handling quick in and out of the HDF5 file
        This class can be used as a context manager

        Any attribute of the HDF will be directly available transparently if
        the source if opened
    """
    def __init__(self, lnpfile, mode='r', **kwargs):
        """__init__

        keywords
        --------

        lnpfile: str or tables.file.File
            storage to use

        mode: str ('r', 'w', 'a', 'r+')
            The mode to open the file.
            (see set_mode)

        **kwargs can be use to set options to tables.openFile
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

        keywords
        --------

        val: str (optional)
            The mode to open the file.  It can be one of the following:
            'r'  -- Read-only; no data can be modified.
            'w'  -- Write; a new file is created (existing file would be deleted)
            'a'  -- Append; an existing file is opened for reading and writing,
                    and if the file does not exist it is created.
            'r+' -- It is similar to 'a', but the file must already exist.

        **kwargs is forwarded to tables.openFile

        returns
        -------
        status: str
            return the status only if None was provided
        """
        if val is None:
            return self._mode
        else:
            title = kwargs.get('title', None) or self._open_kwargs.get('title', '')
            rootUEP = kwargs.get('rootUEP', None) or self._open_kwargs.get('rootUEP', '/')
            filters = kwargs.get('filters', None) or self._open_kwargs.get('filters', None)
            self._open_kwargs.update(dict(title=title, rootUEP=rootUEP, filters=filters))
            if (val != self._mode) & (self._mode is not None):
                self._mode = str(val)
                if self.source is not None:
                    self.close_source(force=True)
                    self.open_source()

    def keep_open(self, val=None):
        """keep_open - set a flag to keep the storage file open for multiple operations

        keywords
        --------

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
                self.source = tables.openFile(self.lnpfile, **self._open_kwargs)
            elif isinstance(self.lnpfile, tables.file.File):
                self.source = self.lnpfile
            return True
        return False

    def close_source(self, force=False):
        """close_source
        close the file only if the object opened any file

        keywords
        --------
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
            r = list(k._v_pathname for k in s.source.walkNodes())
        return r

    def items(self):
        """
        key->group
        """
        return list(self.iteritems())

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
            g = list(s.source.walkGroups(where=where))
        return g

    def __getattr__(self, name):
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
        return key in self.keys()

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

            keywords
            --------
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
            node = storage.source.getNode(nodename)
            colnames = node.colnames
            if expr in colnames:
                if hasattr(node, 'readCoordinates') & (coordinates is not None):
                    q = node.readCoordinates(coordinates, field=expr)
                else:
                    q = node.read(field=expr)
            else:
                #parse qname, using tables.Expr does not help with log10 etc
                names = tables.expression.getExprNames(expr, {})[0]
                if hasattr(node, 'readCoordinates') & (coordinates is not None):
                    for k in names:
                        if k in colnames:
                            condvars[k] = node.readCoordinates(coordinates, field=k)
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
