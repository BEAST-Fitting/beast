""" some common tools """
import sys

__all__ = ['isNestedInstance', 'path_of_module', 'pretty_size_print', 'nbytes']


def isNestedInstance(obj, cl):
    """ Test for sub-classes types
    I could not find a universal test

    Parameters
    ----------
    obj: object instance
        object to test

    cl: Class
        top level class to test

    returns
    -------
    r: bool
        True if obj is indeed an instance or subclass instance of cl
    """
    tree = []
    for k in cl.__subclasses__():
        tree += k.__subclasses__()
    tree += cl.__subclasses__() + [ cl ]
    return issubclass(obj.__class__, tuple(tree))


def path_of_module(mod=None):
    """ returns the definition code path of a given module, object or function
    If nothing is provided, the current frame will be query, i.e., the current
    working directory of the calling function.

    Parameters
    ----------
    mod: module, class, function
        object to find the defintion
        if None, inspect.currentframe is used

    returns
    -------
    path: str
        path of the definition
    """
    import os
    import inspect

    if mod is None:
        mod = inspect.currentframe()
    return '/'.join(os.path.abspath(inspect.getfile(mod)).split('/')[:-1])


def pretty_size_print(num_bytes):
    """
    Output number of bytes in a human readable format

    Parameters
    ----------
    num_bytes: int
        number of bytes to convert

    returns
    -------
    output: str
        string representation of the size with appropriate unit scale
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


def nbytes(obj, pprint=False):
    """ return the number of bytes of the object, which includes size of nested
    structures

    Parameters
    ----------
    obj: object
        object to find the size of

    pprint: bool, optional (default=False)
        if set, returns the result after calling pretty_size_print

    returns
    -------
    num_bytes: int or str
        total number of bytes or human readable corresponding string
    """
    num_bytes = sum(k.nbytes if hasattr(k, 'nbytes') else sys.getsizeof(k) for k in obj.__dict__.values())
    if pprint:
        return pretty_size_print(num_bytes)
    else:
        return num_bytes
