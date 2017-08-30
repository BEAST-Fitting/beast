""" Some useful decorators """

import time
import math
import functools
from functools import wraps, partial
from multiprocessing import Pool
import warnings


class timeit(object):
    """ Time a block of your code.
    This can be used as a

    CONTEXT MANAGER:

        with timeit(text):
            <code>

        text    str Text to display

    as a FUNCTION DECORATOR
        @timeit(verbose=True)
        def myfunc(...):
           ...

        you can then call back the time values associated to your
        function: myfunc.time

    KEYWORDS:
        f   function
    """
    def __init__(self, f=None, verbose=True, text=None):
        self.f = f
        if not self.f is None:
            if type(self.f) != str:
                functools.update_wrapper(self, f)
                self.text = self.__name__
            else:
                self.text = f

        else:
            self.text = text or ''
        self.verbose = verbose

    def __enter__(self):
        print("Timing %s" % (self.text))
        self.start = time.time()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.stop = time.time()
        print(self.time)

    def __pretty_print(self, t):
        units = ["s", "ms", 'us', "ns"]
        scaling = [1, 1e3, 1e6, 1e9]
        if t > 0.0 and t < 1000.0:
            order = min(-int(math.floor(math.log10(t)) // 3), 3)
        elif t >= 1000.0:
            order = 0
        else:
            order = 3

        return "%s Execution time: %.3g %s" % (self.text, t * scaling[order], units[order])

    @property
    def time(self):
        return self.__pretty_print(self.stop - self.start)

    def __call__(self, *args, **kwargs):
        self.start = time.time()
        r = self.f(*args, **kwargs)
        self.stop = time.time()
        if self.verbose:
            print(self.time)
        return r


class memoize(dict):
    '''Decorator. Caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned
    (not reevaluated).
    Fastest implementation according to http://code.activestate.com/recipes/578231/
    '''
    def __init__(self, func):
        self.func = func
        functools.update_wrapper(self, func)

    def __getitem__(self, *key):
        return dict.__getitem__(self, key)

    def __missing__(self, key):
        ret = self[key] = self.func(*key)
        return ret

    __call__ = __getitem__

    def __repr__(self):
        '''Return the function's docstring.'''
        return self.func.__doc__


def run_async(func, pool=None):
    """
    run_async(func)
        function decorator, intended to make "func" run in a separate
        thread (asynchronously).
        Returns the created Thread object

        E.g.:
        @run_async
        def task1():
            do_something

        @run_async
        def task2():
            do_something_too

        t1 = task1()
        t2 = task2()
        ...
        t1.join()
        t2.join()
    """
    _pool = pool or Pool()

    @wraps(func)
    def async_func(*args, **kwargs):
        func_hl = _pool.Process(target=func, args=args, kwargs=kwargs)
        func_hl.start()
        return func_hl

    return async_func


def deprecated(message=None, stacklevel=2):
    """
    This is a decorator which can be used to mark functions as deprecated. It
    will result in a warning being emitted when the function is used.
    """
    def deco(func):
        @wraps(func)
        def new_func(*args, **kwargs):
            warnings.warn("Call to deprecated function {}. {}".format(func.__name__, message or ''),
                        category=DeprecationWarning, stacklevel=stacklevel)
            return func(*args, **kwargs)
        return new_func
    return deco


def elementwise(func):
    """ Quick and dirty elementwise function decorator it provides a quick way
    to apply a function either on one element or a sequence of elements """
    @wraps(func)
    def wrapper(it, **kwargs):
        if hasattr(it, '__iter__'):  # is a Sequence
            _f = partial(func, **kwargs)
            return list(map(_f, it))
        else:
            return func(it, **kwargs)
    return wrapper


def warning(message=None, category=UserWarning, stacklevel=2):
    """
    This is a decorator which can be used to mark functions with warnings. It
    will result in a warning being emitted when the function is used.
    """
    def deco(func):
        @wraps(func)
        def new_func(*args, **kwargs):
            warnings.warn("Call to function {}. {}".format(func.__name__, message or ''),
                        category=category, stacklevel=stacklevel)
            return func(*args, **kwargs)
        return new_func
    return deco
