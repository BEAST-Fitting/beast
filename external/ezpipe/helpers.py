""" This is a first collection of tools making the design of a pipeline easy.
Pipeline being a succesion of  entities that produce and consume data.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import sys
from time import ctime
from functools import partial, wraps, update_wrapper
from inspect import getargspec, ismethod
from .pickleshare import PickleShareDB


class RequiredFile(object):
    """ This class implement a first part of pipeline flow. It defines files
    that are required for the chain of function to exists. If the producer of
    the file (and calling arguments) are provided, the process will be called
    to generate the file

    this class can be used as a context manager (with)
    """
    def __init__(self, fname, producer=None, *args, **kwargs):
        self.fname = fname
        self._producer = producer
        self.args = args
        self.kwargs = kwargs

    @property
    def producer(self):
        if self._producer is None:
            return None
        else:
            #return '{}({}, {})'.format(self._producer.__name__, self.args, self.kwargs)
            return '{}({}, {})'.format(self._producer.__name__, '...', '...')

    def __repr__(self):
        txt = 'Required File {}'.format(self.fname)
        if self._producer is not None:
            txt += '\n    Produced by {}'.format(self.producer)
        return txt

    def __call__(self, *args, **kwargs):
        if not (os.path.isfile(self.fname) | os.path.islink(self.fname)):
            if self.producer is not None:
                print(("File {} is missing.\n Running {})".format(self.fname, self.producer)))
                return self._producer(*(self.args), **(self.kwargs))
            else:
                raise Exception('File {} is required but no producer is referenced'.format(self.fname))
        else:
            return self.fname

    def __str__(self):
        return self.fname

    def __enter__(self):
        self()
        return self.fname

    def __exit__(self, exc_type, exc_value, traceback):
        pass


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


class Task(object):
    """ A Task is a Pipe-like function that also propagates job information
        note: a task is pickleable
        TODO:
           * add memoization / conditional work
              This should not run the function if not necessary
    """

    def __init__(self, func, args=(), name=None, memoized=False, logger=None, **kwargs):
        if type(memoized) not in [ bool, str ]:
            raise ValueError('memoized value "{}" not understood. Expecting Boolean or str')
        self.memoized = memoized
        self.func = func
        self.name = name or func.__name__
        self.args = args
        self.kwargs = kwargs
        update_wrapper(self, func)
        self.__set_memoize_backend__()
        self.set_logger(logger)

    def set_logger(self, logger, **kwargs):
        """ set the logger backend """
        if isinstance(logger, Logger):
            self.logger = logger
        else:
            self.logger = Logger(out=logger, **kwargs)

    def memoize(self, val=None):
        if val is None:
            return self.memoized
        else:
            self.memoized = val
            self.__set_memoize_backend__()

    def clear_memoized(self):
        if self._cache is not None:
            self._cache.clear()

    def get_job_status(self, job_id, strict=True):
        """ Check if previous results are available
            if strict, unkown jobs will raise KeyError exceptions, otherwise
            return False
        """
        key = '{}/{}'.format(job_id, self.name)
        if self.memoized is False:
            return False
        elif (key in self._cache):
            return True
        else:
            if strict:
                raise KeyError('Job "{}" not found'.format(job_id))
            else:
                return False

    def get(self, job_id, strict=True):
        """ get the result from a previous memoized job
            See also get_job_status
        """
        if self.get_job_status(job_id, strict=strict):
            key = '{}/{}'.format(job_id, self.name)
            return self._cache[key]

    def __set_memoize_backend__(self):
        """ Defines the backend for memoization of the results """
        if self.memoized is False:
            self._cache = None
        elif self.memoized is True:
            self._cache = {}
        elif isinstance(self.memoized, str):
            self._cache = PickleShareDB(self.memoized)

    def run(self, key, val, *args, **kwargs):
        if hasattr(val, '__iter__') and type(val) not in ['string']:
            _args = tuple(val) + args + self.args
        else:
            _args = (val, ) + args + self.args
        _kwargs = self.kwargs.copy()
        _kwargs.update(kwargs)
        if self.memoized is False:
            self.logger.write('Job {} -- started\n'.format(key))
            ret = self.func(_args, **_kwargs)
            self.logger.write('Job {} -- ended\n'.format(key))
            return ret
        else:
            if (key not in list(self._cache.keys())):
                self.logger.write('Job {} -- started\n'.format(key))
                self._cache[key] = self.func(*_args, **_kwargs)
                self.logger.write('Job {} -- ended\n'.format(key))
                return self._cache[key]
            else:
                self.logger.write('Job {} -- existing result, skipped\n'.format(key))
                return self._cache[key]

    def __ror__(self, other):
        """ Used for the | operator """
        job_id = other[0]
        key = '{}/{}'.format(job_id, self.name)
        return job_id, self.run(key, other[1])

    def __call__(self, *args, **kwargs):
        if len(args) > 0:
            self.args = args
        if len(kwargs) > 0:
            self.memoized = kwargs.pop('memoized', self.memoized)
            self.logger = kwargs.pop('memoized', self.logger)
            self.kwargs = kwargs
        return self

    def __repr__(self):
        txt = 'Task: {}'.format(self.func.__repr__())
        if self.memoized:
            txt += '\n memoized {}'.format(object.__repr__(self._cache))
        if self.logger is not None:
            txt += '\n Logger {}'.format(self.logger.__repr__())
        return txt

    def __getitem__(self, job_id):
        return self.get(job_id)


def task_decorator(name=None, memoized=False, logger=None, *args, **kwargs):
    """ Decorator to create tasks """
    def deco(func):
        return Task(func, args=args, name=name, memoized=memoized, logger=logger, **kwargs)
    return deco


class Checkpoint(Task):
    """ A checkpoint is a simple Task where the function does lambda x: x,
    i.e., it does not update the data stream but can be used to log information
    or save the stream value between tasks This is especially useful with
    complex pipelines, to save the initial inputs or results from parallel
    processing
    """
    def __init__(self, name, memoized=True, logger=None):
        if type(memoized) not in [ bool, str ]:
            raise ValueError('memoized value "{}" not understood. Expecting Boolean or str')
        self.memoized = memoized
        self.args = ()
        self.kwargs = {}
        self.__set_memoize_backend__()
        self.set_logger(logger)
        self.name = name

    def func(self, val, *args, **kwargs):
        """ identity only """
        return val

    def __repr__(self):
        txt = 'Checkpoint: {}'.format(self.name)
        if self.memoized:
            txt += '\n memoized {}'.format(object.__repr__(self._cache))
        if self.logger is not None:
            txt += '\n Logger {}'.format(self.logger.__repr__())
        return txt


class Logger(object):
    """ Keep track of calling sequences
    This class is a context Manager that can be used as a decorator as well.
    This class is pickleable as long as its output is a file  that is not kept open (keep_opened property)
    and not a buffer.

    Example usage:
        def fn(*args, **kwargs):
            return args, kwargs

        with Logger() as _log:
            res = _log(fn)(root='./data', pattern='*fits')
            _log.message('info', 'hello')
    """
    def __init__(self, out=sys.stdout, mode='a', keep_opened=False):
        self.out = out
        self._buffer = None
        self.mode = mode
        self.keep_opened = keep_opened or (type(out) != str)

    @property
    def buffer(self):
        """ returns the appropriate buffer object, open it if necessary """
        if self._buffer is not None:
            return self._buffer
        else:
            self.open_buffer()
            return self._buffer

    def open_buffer(self):
        """ open the buffer """
        if type(self.out) == str:
            self._buffer = open(self.out, self.mode)
        else:
            self._buffer = self.out

    def close_buffer(self, force=False):
        """ close the buffer unless keep_opened is True """
        if (self.keep_opened is False) & (not force):
            if type(self.out) == str:
                self._buffer.close()
                self._buffer = None

    def open(self):
        """ open the buffer """
        self.open_buffer()

    def close(self):
        """ close the buffer """
        self.close_buffer(force=True)

    def __enter__(self):
        """ context manager """
        if self.keep_opened:
            self.open_buffer()
        return self

    def __call__(self, f):
        def caller(*args, **kwargs):
            self.open_buffer()
            if self._buffer is not None:
                self._buffer.write('[{}] {}: {}, {}\n'.format(ctime(), f.__name__, args.__repr__(), kwargs.__repr__()))
            self.close_buffer()
            return f(*args, **kwargs)
        return caller

    def __exit__(self, exc_type, exc_value, traceback):
        """ context manager """
        self.close()

    def message(self, tag, *args):
        self.open_buffer()
        if self._buffer is not None:
            self._buffer.write('[{}] {}: {}\n'.format(ctime(), tag, args))
        self.close_buffer()

    def write(self, string):
        self.open_buffer()
        if self._buffer is not None:
            self._buffer.write('[{}] {}'.format(ctime(), string))
        self.close_buffer()

    def writelines(self, sequence):
        self.open_buffer()
        kp = self.keep_opened
        self.keep_opened = True
        for string in sequence:
            self.write(string)
        self.close_buffer()
        self.keep_opened = kp

    def __repr__(self):
        if self.out is None:
            txt = 'NULL'
        elif type(self.out) == str:
            txt = self.out
        else:
            txt = self.out.name
        return 'Logger: {}, {}'.format(txt, object.__repr__(self))


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
