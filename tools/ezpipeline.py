""" This is a first collection of tools making the design of a pipeline easy.
Pipeline being a succesion of  entities that produce and consume data.
"""
import os


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
                print "File {} is missing.\n Running {})".format(self.fname, self.producer)
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
    | operator
    """
    def __init__(self, func, *args, **kwargs):
        self.func = func
        self.args = args
        self.kwargs = kwargs

    def __ror__(self, lhs):
        return self.func(lhs, *self.args, **self.kwargs)

    def __call__(self, *args, **kwargs):
        return self.func(*args, **kwargs)
