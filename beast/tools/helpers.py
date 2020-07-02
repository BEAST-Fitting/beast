"""
This is a first collection of tools making the design easier
"""
import sys
from functools import wraps
import itertools

# replace the common range by the generator
try:
    range = xrange
except NameError:
    pass


__all__ = [
    "generator",
    "chunks",
    "isNestedInstance",
    "type_checker",
    "nbytes",
]


def generator(func):
    """ A dummy decorator that only make codes mode readable.
    It allow to explicitly mark a function as generator (yielding values)
    and does nothing more than calling the initial function
    """

    @wraps(func)
    def deco(*args, **kwargs):
        return func(*args, **kwargs)

    return deco


@generator
def chunks(ll, n):
    """ Yield successive n-sized chunks from l.

    Parameters
    ----------
    ll: iterable
        object to iter over

    n: int
        number of elements per slice

    Returns
    -------
    chunk: tuple
        n values from ll
    """
    it = iter(ll)
    while True:
        chunk = tuple(itertools.islice(it, n))
        if chunk:
            yield chunk
        else:
            return
            # raise StopIteration


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
    tree = [cl]
    if hasattr(cl, "__subclasses"):
        for k in cl.__subclasses():
            if hasattr(k, "__subclasses"):
                tree += k.__subclasses__()
    return issubclass(obj.__class__, tuple(tree))


def type_checker(name, obj, tp):
    """ Check a given type and raise a type error if not correct
    Parameters
    ----------
    name: str
        name of the variable to show in the exception text
    obj: object
        object to check
    tp: type
        expected type of obj
    Raises
    ------
    :exc:TypeError:
        raises a TypeError if object is not of the correct type of a subclass of it
    """
    if not isNestedInstance(obj, tp):
        txt = 'Expected "{0:s}" of type {1:s}, got {2:s} instead.'
        raise TypeError(txt.format(name, str(tp.__name__), str(type(obj).__name__)))


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
        output = "%.3g YB" % (num_bytes / YiB)
    elif num_bytes > ZiB:
        output = "%.3g ZB" % (num_bytes / ZiB)
    elif num_bytes > EiB:
        output = "%.3g EB" % (num_bytes / EiB)
    elif num_bytes > PiB:
        output = "%.3g PB" % (num_bytes / PiB)
    elif num_bytes > TiB:
        output = "%.3g TB" % (num_bytes / TiB)
    elif num_bytes > GiB:
        output = "%.3g GB" % (num_bytes / GiB)
    elif num_bytes > MiB:
        output = "%.3g MB" % (num_bytes / MiB)
    elif num_bytes > KiB:
        output = "%.3g KB" % (num_bytes / KiB)
    else:
        output = "%.3g Bytes" % (num_bytes)

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
    num_bytes = sum(
        k.nbytes if hasattr(k, "nbytes") else sys.getsizeof(k)
        for k in list(obj.__dict__.values())
    )
    if pprint:
        return pretty_size_print(num_bytes)
    else:
        return num_bytes
