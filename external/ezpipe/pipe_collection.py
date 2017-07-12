from .helpers import Pipe
import itertools
import sys

try:
    import builtins
except ImportError:
    import builtins as builtins


@Pipe
def groupby(iterable, keyfunc):
    return itertools.groupby(sorted(iterable, key=keyfunc), keyfunc)


@Pipe
def sort(iterable, **kwargs):
    return sorted(iterable, **kwargs)


@Pipe
def reverse(iterable):
    return reversed(iterable)


@Pipe
def null(x):
    list(x)
    return


@Pipe
def mean(iterable):
    """
    Build the average for the given iterable, starting with 0.0 as seed
    Will try a division by 0 if the iterable is empty...
    """
    total = 0.0
    qte = 0
    for x in iterable:
        total += x
        qte += 1
    if qte == 0:
        return None
    else:
        return total / qte


@Pipe
def count(iterable):
    "Count the size of the given iterable, walking thrue it."
    count = 0
    for x in iterable:
        count += 1
    return count


@Pipe
def max(iterable, **kwargs):
    return builtins.max(iterable, **kwargs)


@Pipe
def min(iterable, **kwargs):
    return builtins.min(iterable, **kwargs)


@Pipe
def ptp(iterable, **kwargs):
    return (builtins.min(iterable, **kwargs), builtins.max(iterable, **kwargs))


@Pipe
def median(iterable, **kwargs):
    return 0.5 * (builtins.min(iterable, **kwargs) + builtins.max(iterable, **kwargs))


@Pipe
def sum(x):
    return builtins.sum(x)


@Pipe
def reduce(iterable, function, initial=None):
    if initial is not None:
        return builtins.reduce(function, iterable, initial)
    else:
        return builtins.reduce(function, iterable)


@Pipe
def where(iterable, predicate):
    return (x for x in iterable if (predicate(x)))


@Pipe
def to_dict(iterable):
    return dict(iterable)


@Pipe
def to_list(iterable):
    return list(iterable)


@Pipe
def to_type(x, t):
    return t(x)


@Pipe
def to_tuple(iterable):
    return tuple(iterable)


@Pipe
def bwrite(val, fmt="%s", buffer=sys.stdout):
    buffer.write(fmt % val)


@Pipe
def tee(iterable):
    for item in iterable:
        sys.stdout.write(str(item) + "\n")
        yield item


def _fmap(val, fnseq):
    return (fk(val) for fk in fnseq)


@Pipe
def apply(vals, args=(), iter=True):
    r = _fmap(vals, args)
    if not iter:
        return list(r)
    else:
        return r
