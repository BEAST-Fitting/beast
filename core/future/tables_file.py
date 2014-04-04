from tables.flavor import flavor_of, array_as_internal
from tables.atom import Atom
from tables.earray import EArray
from tables.filters import Filters


def _checkfilters(filters):
    if not (filters is None or
            isinstance(filters, Filters)):
        raise TypeError("filter parameter has to be None or a Filter "
                        "instance and the passed type is: '%s'" %
                        type(filters))


def create_earray(fileobj, where, name, atom=None, shape=None, title="",
                  filters=None, expectedrows=1000, chunkshape=None,
                  byteorder=None, createparents=False, obj=None):
    """Create a new enlargeable array.

    Parameters
    ----------
    fileobj: File instance
        file in which creating the table
    where : str or Group
        The parent group from which the new array will hang. It can be a
        path string (for example '/level1/leaf5'), or a Group instance
        (see :ref:`GroupClassDescr`).
    name : str
        The name of the new array
    atom : Atom
        An Atom (see :ref:`AtomClassDescr`) instance representing the
        *type* and *shape* of the atomic objects to be saved.

        .. versionchanged:: 3.0
            The *atom* parameter can be None (default) if *obj* is
            provided.

    shape : tuple
        The shape of the new array.  One (and only one) of the shape
        dimensions *must* be 0.  The dimension being 0 means that the
        resulting EArray object can be extended along it.  Multiple
        enlargeable dimensions are not supported right now.

        .. versionchanged:: 3.0
            The *shape* parameter can be None (default) if *obj* is
            provided.

    title : str, optional
        A description for this node (it sets the TITLE HDF5 attribute on
        disk).
    expectedrows : int, optional
        A user estimate about the number of row elements that will be added
        to the growable dimension in the EArray node.  If not provided, the
        default value is EXPECTED_ROWS_EARRAY (see tables/parameters.py).
        If you plan to create either a much smaller or a much bigger array
        try providing a guess; this will optimize the HDF5 B-Tree creation
        and management process time and the amount of memory used.
    chunkshape : tuple, numeric, or None, optional
        The shape of the data chunk to be read or written in a single HDF5
        I/O operation.  Filters are applied to those chunks of data.  The
        dimensionality of chunkshape must be the same as that of shape
        (beware: no dimension should be 0 this time!).  If None, a sensible
        value is calculated based on the expectedrows parameter (which is
        recommended).
    byteorder : str, optional
        The byteorder of the data *on disk*, specified as 'little' or
        'big'. If this is not specified, the byteorder is that of the
        platform.
    createparents : bool, optional
        Whether to create the needed groups for the parent path to exist
        (not done by default).
    obj : python object
        The array or scalar to be saved.  Accepted types are NumPy
        arrays and scalars, as well as native Python sequences and
        scalars, provided that values are regular (i.e. they are
        not like ``[[1,2],2]``) and homogeneous (i.e. all the
        elements are of the same type).

        The *obj* parameter is optional and it can be provided in
        alternative to the *atom* and *shape* parameters.
        If both *obj* and *atom* and/or *shape* are provided they must
        be consistent with each other.

        .. versionadded:: 3.0

    See Also
    --------
    EArray : for more information on enlargeable arrays

    """
    if obj is not None:
        flavor = flavor_of(obj)
        obj = array_as_internal(obj, flavor)

        earray_shape = (0,) + obj.shape[1:]

        if shape is not None and shape != earray_shape:
            raise TypeError('the shape parameter is not compatible '
                            'with obj.shape.')
        else:
            shape = earray_shape

        if atom is not None and atom.dtype != obj.dtype:
            raise TypeError('the atom parameter is not consistent with '
                            'the data type of the obj parameter')
        elif atom is None:
            atom = Atom.from_dtype(obj.dtype)

    parentnode = fileobj._getOrCreatePath(where, createparents)
    _checkfilters(filters)
    ptobj = EArray(parentnode, name,
                   atom=atom, shape=shape, title=title,
                   filters=filters, expectedrows=expectedrows,
                   chunkshape=chunkshape, byteorder=byteorder)

    if obj is not None:
        ptobj.append(obj)

    return ptobj

createEArray = create_earray
