"""
This package groups a few simple decorators but especially useful during
debugging phases

    memoize         (fast version) Caches a function's return value each time it is called.
    persistent_locals   Keeps access to the local variables accessible after execution.
    timeit          Time a block of your code (also usable as a context manager, 'with').
    mapreduce       in dev...
    run_async       intended to make decorated functions to run in a separate thread (asynchronously).
    trace           Trace a function call and logs args, kwargs, result and optionaly execution time

Requires dependencies:
    memtrace        This decorator uses <memory_profiler> to trace line by line the memory

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import functools
from functools import wraps
import time
import math
import sys
from multiprocessing import Pool, cpu_count


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

try:
    import memory_profiler

    class memtrace(object):
        """ This decorator uses memory_profiler to trace line by line the memory
        of a function
        Show the results by using func.prof
        """
        def __init__(self, f):
            self.f = f
            self.m = None
            functools.update_wrapper(self, f)

        def __call__(self, *args):
            self.m = memory_profiler.LineProfiler()
            self.c = self.m(self.f)
            return self.c(*args)

        def __repr__(self):
            '''Return the function's docstring.'''
            return self.f.__repr__()

        @property
        def prof(self):
            if not self.m is None:
                return memory_profiler.show_results(self.m)
            else:
                return None

    @memtrace
    def test_memtrace(*args):
        """ Test memtrace decorator """
        a = 2
        b = [1] * 5
        c = args
        print(args)
        return args
except ImportError:
    pass


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


class persistent_locals(object):
    """Decorator. Keeps the local variables accessible after execution.
    Uses python profiler to retrieve local variables and keep them into the
    func._locals for later usage

    @persistent_locals
    def func(...):
        ...

    func(...)

    func._locals
    func.clear_locals()
    """
    def __init__(self, func):
        self.func = func
        functools.update_wrapper(self, func)
        self._locals = {}

    def __call__(self, *args, **kwargs):
        def tracer(frame, event, arg):
            if event == 'return':
                self._locals = frame.f_locals.copy()

        # tracer is activated on next call, return or exception
        sys.setprofile(tracer)
        try:
            # trace the function call
            res = self.func(*args, **kwargs)
        finally:
            # disable tracer and replace with old one
            sys.setprofile(None)
        return res

    def clear_locals(self):
        self._locals = {}


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


def mapreduce(func, pool=None, ncpu=None, chunksize=None):

    if ncpu is None:
       ncpu = cpu_count()
    _pool = pool or Pool(processes=ncpu)

    @wraps(func)
    def map(*args, **kwargs):
        return _pool.map(func, args, kwargs, chunksize=None)


def trace(f, output=sys.stdout, time=True):
    """ Trace a function call
    @trace(output=sys.stdout)
    def func(...):
        ...

    KEYWORDS:
        output  file    default is stdout

    ex usage:
        > r = trace(fibonacci)
        > r(10)
        fibonacci: (10,), {} > 55
    """

    @wraps(f)
    def func(*args, **kwargs):
        if time:
            @timeit
            @wraps(f)
            def f1(*args, **kargs):
                return f(*args, **kwargs)
            f1.verbose = False
            r  = f1(*args, **kwargs)
            output.write("%s: %s, %s > %s , %s\n" % (f.__name__, args, kwargs, r, f1.time.split(':')[-1]) )
        else:
            r = f(*args, **kwargs)
            output.write("%s: %s, %s > %s \n" % (f.__name__, args, kwargs, r) )
        return r
    return func


if __name__ == '__main__': 

    @trace
    @memoize
    def fibonacci(n):
        """Return the nth fibonacci number."""
        if n in (0, 1):
            return n
        return fibonacci(n - 1) + fibonacci(n - 2)


    @persistent_locals
    def test_locals(*args, **kwargs):
        a = 2
        b = list(range(10))
        c = 'hello'
        return a + a


    @timeit
    def example():
        with timeit('Fibonacci cold start'):
            fibonacci(50)

        with timeit('Fibonacci second call'):
            fibonacci(50)

        test_memtrace([2] * 4)
        test_memtrace.prof

        test_locals()
        print(test_locals._locals)


    @run_async
    def test_async():
        from time import sleep
        print('starting print_somedata')
        sleep(2)
        print('print_somedata: 2 sec passed')
        sleep(2)
        print('print_somedata: 2 sec passed')
        sleep(2)
        print('finished print_somedata')
