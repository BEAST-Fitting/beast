from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .backends.basebackend import BaseBackend

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

__all__ = ['register_extension', 'determine_type']

global __extensions__
__extensions__ = {}

#==============================================================================
# COMMON functions
#==============================================================================


def register_extension(extensions, backend=None, readerFunction=None,
                       writerFunction=None, override=False):
    """
    Register storage supports from either a Backend object or individual
    read/write functions

    INPUTS:
        extensions  [str or seq]    register particular extensions

    KEYWORDS:
        backend     [Backend]   core.basebackend.BaseBackend Object
        readerFunction  [function]  function required for the reading part
                        this function should return a Table object
        writerFunction  [function]  function required for the writing part
                        this function should return a Table object
                        Will receive a Table object as input
                        and should also handle extra **kwargs
        override    [bool]      redefine the file manager if already registered.
    """
    if not hasattr(extensions, '__iter__'):
        extensions = [extensions]

    #if backend is not None:
    #   assert(isNestedInstance(backend, basebackend.BaseBackend)), "backend is expected to be a core.basebackend.BaseBackend instance"
    if backend is None:
        assert((readerFunction is not None) & (writerFunction is not None)), "read/write functions required"
        backend = BaseBackend(extensions[0], readerFunction, writerFunction  )

    if isinstance(extensions, str): extensions = [extensions]
    for k in extensions:
        if (not k in __extensions__) or override:
            __extensions__[k] = backend
        else:
            raise Exception("Type %s is already defined" % k)


def determine_type(mystring, verbose=True):
    """
    Determine the type of a table from its extension and try to give the
    point to the appropriate registered extension
    """
    if type(mystring) != str and type(mystring) != unicode:
        raise Exception('Could not determine input type (non-string input)')

    if len(__extensions__) == 0:
        set_defaults()

    s = mystring.lower()
    if not '.' in s:
        extension = s
    else:
        extension = s.split('.')[-1]

        if extension in ['gz', 'bzip2', 'Z']:
            raise Exception('Compressed files are not managed yet.')

        elif extension in ['sav', 'idl', 'idlsav']:
            raise Exception('Warning: IDL save files are not managed yet.')

    if extension in __extensions__:
        tableType = __extensions__[extension]
        if verbose:
            print("Auto-detected type: %s" % extension)
    else:
            raise Exception('Could not determine input type for extension %s' % extension)
    return tableType


def list_backends():
    return __extensions__


def set_defaults():
    from .backends.asciibackend import csvBackend, asciiBackend
    from .backends.latexbackend import LatexBackend
    from .backends.fitsbackend import fitsBackend
    from .backends.jsonbackend import jsonBackend
    from .backends.hdf5backend import hdf5Backend
    register_extension('txt dat tsv'.split(), asciiBackend)
    register_extension('csv', csvBackend)
    register_extension('tex', LatexBackend)
    register_extension('fits', fitsBackend)
    register_extension('json', jsonBackend)
    register_extension('hd5 hdf5'.split(), hdf5Backend)
