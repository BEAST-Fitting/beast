"""
This is a first collection of tools making the design easier
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys
from functools import partial, wraps, update_wrapper
from inspect import getargspec, ismethod
import warnings
import numpy as np
import itertools

#replace the common range by the generator
try:
    range = xrange
except NameError:
    pass


__all__ = ['NameSpace', 'Pipe', 'Pipeable', 'Pipegroup', 'chunks',
           'deprecated', 'generator', 'isNestedInstance', 'keywords_first',
           'kfpartial', 'merge_records', 'missing_units_warning', 'nbytes',
           'path_of_module', 'pretty_size_print', 'type_checker']


class NameSpace(dict):
    """A dict subclass that exposes its items as attributes.
    """
    def __init__(self, name, *args, **kwargs):
        self.__name__ = name
        dict.__init__(self, *args, **kwargs)

    def __dir__(self):
        return tuple(self)

    def __repr__(self):
        names = ', '.join([k for k in dir(self) if k[0] != '_'])
        return "{s.__name__:s}: {r:s}".format(s=self, r=names)

    def __getattribute__(self, name):
        try:
            return self[name]
        except KeyError:
            msg = "'{s.__name__:s}' has no attribute '{name:s}'"
            raise AttributeError(msg.format(s=self, name=name))

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        del self[name]

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        return False


def generator(func):
    """ A dummy decorator that only make codes mode readable.
    It allow to explicitly mark a function as generator (yielding values)
    and does nothing more than calling the initial function
    """
    @wraps(func)
    def deco(*args, **kwargs):
        return func(*args, **kwargs)
    return deco


def deprecated(func):
    """ A dummy decorator that warns against using a deprecated function """
    @wraps(func)
    def deco(*args, **kwargs):
        txt = 'Function {0:s} is deprecated. You should avoid its usage'
        warnings.warn(txt.format(func.__name__))
        return func(*args, **kwargs)
    return deco


@generator
def chunks(l, n):
    """ Yield successive n-sized chunks from l.

    Parameters
    ----------
    l: iterable
        object to iter over

    n: int
        number of elements per slice

    Returns
    -------
    chunk: tuple
        n values from l
    """
    it = iter(l)
    while True:
        chunk = tuple(itertools.islice(it, n))
        if chunk:
            yield chunk
        else:
            raise StopIteration


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
    tree = [ cl ]
    if hasattr(cl, '__subclasses'):
        for k in cl.__subclasses():
            if hasattr(k, '__subclasses'):
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


class Pipeable(object):
    """ Decorator overloading | operator (__ror__) such that you can pipe
    functions where the first argument is the variable on the left side of the
    | operator.
    This decorator allows you to use the decorated function normally and uses
    the provided values when using in pipes.

    Example::

        import pylab as plt
        _p = Pipeable(plt.plot, color='red', linestyle='--')
        _p(range(10), 'o-')  #  works
        range(10) | _p      #  will plot a red dashed line

    """

    def __init__(self, func, *args, **kwargs):
        self.func = func
        self.args = args
        self.kwargs = kwargs

    def __ror__(self, lhs):
        return self.func(lhs, *self.args, **self.kwargs)

    def __call__(self, *args, **kwargs):
        return self.func(*args, **kwargs)

    def __repr__(self):
        return self.func.__repr__()


class Pipe(object):
    """ Decorator overloading | operator (__ror__) such that you can pipe
    functions where the first argument is the variable on the left side of the
    | operator.
    The difference with Pipeable is that you cannot use decorated function
    outside of pipes but you gain the possibility to update the calling
    parameters

    Used with keywords_first make this a powerful Task
    """
    def __init__(self, func, *args, **kwargs):
        self.func = func
        self.args = args
        self.kwargs = kwargs
        update_wrapper(self, func)

    def __or__(self, other):
        if isinstance(other, Pipe):
            return Pipegroup( (self, other), mode='multi' )

    def __and__(self, other):
        if isinstance(other, Pipe):
            return Pipegroup( (self, other) )

    def __ror__(self, other):
            return self.func(other, *self.args, **self.kwargs)

    def __call__(self, *args, **kwargs):
        return Pipeable(self.func, *args, **kwargs)

    def __repr__(self):
        return 'Pipe: {}'.format(self.func.__repr__())

    def __str__(self):
        return '{}'.format(self.func.__name__)


class Pipegroup(object):
    def __init__(self, pipes, mode='sequential'):
        self.pipes = list(pipes)
        self.mode = mode
        self.func = self

    def __len__(self):
        return len(self.pipes)

    def seq_call(self, val, *args, **kwargs):
        r = self.pipes[0].func(val)
        if len(self) > 1:
            for pk in self.pipes[1:]:
                r = pk.func(r)
        return r

    def append(self, other):
        self.pipes.append(other)

    def __or__(self, other):
        if isinstance(other, Pipe):
            if self.mode in ['multi', '|']:
                self.append(other)
                return self
            else:
                return Pipegroup( (self, other), mode='multi')

    def __and__(self, other):
        if isinstance(other, Pipe):
            if self.mode in ['sequential', '&']:
                self.append(other)
                return self
            else:
                return Pipegroup( (self, other), mode='sequential')

    def multi_call(self, vals, iter=True):
        return [ pk.func(vals) for pk in self.pipes ]

    def __call__(self, val, *args, **kwargs):
        mode = kwargs.get('mode', self.mode).lower()
        if mode in ['sequential', '&']:
            return self.seq_call(val, *args, **kwargs)
        elif mode in ['multi', '|']:
            return self.multi_call(val, *args, **kwargs)
        else:
            raise NotImplemented

    def __ror__(self, other):
        return self(other)

    def __repr__(self):
        txt = 'Pipegroup: mode={},\n\t | '.format(self.mode)

        if self.mode == 'sequential':
            txt += ' & '.join([str(pk) for pk in self.pipes])
        else:
            txt += '\n\t | '.join([str(pk) for pk in self.pipes])
        return txt

    def __str__(self):
        if self.mode == 'sequential':
            delim = ' & '
        else:
            delim = ' | '
        return '({})'.format(delim.join([str(pk) for pk in self.pipes])).replace('Pipe: ', '')


def keywords_first(f):
    """ Decorator that enables to access any argument or keyword as a keyword """
    ## http://code.activestate.com/recipes/577922/ (r2)
    @wraps(f)
    def wrapper(*a, **k):
        a = list(a)
        #for idx, arg in enumerate(f.func_code.co_varnames[:f.func_code.co_argcount], -ismethod(f)):
        for idx, arg in enumerate(getargspec(f).args, -ismethod(f)):  # or [0] in 2.5
            if arg in k:
                if idx < len(a):
                    a.insert(idx, k.pop(arg))
                else:
                    break
        return f(*a, **k)
    return wrapper


def kfpartial(fun, *args, **kwargs):
    """ Allows to create partial functions with arbitrary arguments/keywords """
    return partial(keywords_first(fun), *args, **kwargs)


def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return " {0:s}:{1:d} {2:s}:{3:s}".format(filename, lineno, category.__name__, message)


def merge_records(lst):
    """ generates a stack of records even with slightly different but compatible dtypes

    Parameters
    ----------

    lst: sequence of np.recarray
        sequence of individual records

    Returns
    -------
    val: np.recarray
        array of stacked records
        Note if if lst is empty, returns an empty list
    """
    r = []
    for rk in lst:
        r.append(rk.tolist()[0])
    names = rk.dtype.names
    if len(r) > 0:
        return np.rec.fromrecords(r, names=names)
    else:
        return []


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
    num_bytes = sum(k.nbytes if hasattr(k, 'nbytes') else sys.getsizeof(k) for k in list(obj.__dict__.values()))
    if pprint:
        return pretty_size_print(num_bytes)
    else:
        return num_bytes
