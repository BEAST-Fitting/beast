""" A simpler implementation than OrderedDict which keeps ordering but allows to
move the keys too"""

__all__ = ['odict']
__author__ = 'MF'
__version__ = '1.0'


try:
    unicode = unicode
except NameError:
    # 'unicode' is undefined, must be Python 3
    str = str
    unicode = str
    bytes = bytes
    basestring = (str,bytes)
    __strtypes__ = [str, unicode]
else:
    # 'unicode' exists, must be Python 2
    str = str
    unicode = unicode
    bytes = str
    basestring = basestring
    __strtypes__ = [str, unicode]

# Add Numpy str type if possible
try:
    import numpy as __np__
    __strtypes__.append(__np__.string_)
except ImportError:
    pass


class odict(object):
    """ A simpler implementation than OrderedDict
    This implementation aims at keeping ordering and also  at allowing to
    move the keys too
    """
    def __init__(self, **kwargs):
        self.__keys__ = []
        self.__values__ = []

        if ( len( kwargs) > 0 ):
            for k, v in list(kwargs.items()):
                self[k] = v

    def keys(self):
        return self.__keys__

    def __setitem__(self, key, value):
        if type(key) == int:
            if key > len(self.__keys__) - 1:
                raise Exception("Element %i does not exist" % key)
            else:
                self.__values__[key] = value
        elif type(key) in __strtypes__:
            if key in self.__keys__:
                index = self.__keys__.index(key)
                self.__values__[index] = value
            else:
                self.__keys__.append(key)
                self.__values__.append(value)
        else:
            print(type(key))
            print(__strtypes__)
            raise Exception("Wrong type for key: %s" % type(key))

    def __getitem__(self, key):
        if type(key) == int:
            return self.__values__[key]
        elif type(key) in __strtypes__:
            if not key in self.__keys__:
                raise KeyError
            index = self.__keys__.index(key)
            return self.__values__[index]
        else:
            raise Exception("Wrong type for key: %s" % type(key))

    def __repr__(self):
        string = "odict({"
        for i, key in enumerate(self.__keys__):
            if i > 0:
                string += ", "
            string += "'%s': %s" % (key, self.__values__[i])
        string += "})"
        return string

    def __contains__(self, key):
        return key in self.__keys__

    def pop(self, key):
        """remove specified key and return the corresponding value."""
        if not key in self.__keys__:
            raise KeyError
        index = self.__keys__.index(key)
        self.__keys__.pop(index)
        return self.__values__.pop(index)

    def __len__(self):
        return len(self.__keys__)

    def move(self, key, position):
        """ Move a key-value position """
        if not key in self.__keys__:
            raise KeyError
        v = self.pop(key)
        self.insert(position, key, v)

    def rename(self, oldkey, newkey):
        """ Rename a key """
        if not oldkey in self.__keys__:
            raise KeyError
        index = self.__keys__.index(oldkey)
        self.__keys__[index] = newkey
        return

    def insert(self, position, key, value):
        """ Insert a key-value pair at a given position """
        self.__keys__.insert(position, key)
        self.__values__.insert(position, value)

    def __iter__(self):
        return iter(self.__keys__)

    def iteritems(self):
        """an iterator over the (key, values) of D"""
        return iter(zip(self.__keys__, self.__values__))

    def items(self):
        """ (key, values) of D """
        return list(zip(self.__keys__, self.__values__))

    def iterkeys(self):
        """an iterator over the keys of D"""
        return iter(self.__keys__)

    def itervalues(self):
        """an iterator over the values of D"""
        return iter(self.__values__)

    def get(self, key, default=None):
        if key in list(self.keys()):
            return self[key]
        else:
            return default

    def update(self, other):
        if hasattr(other, 'keys'):
            for e, v in list(other.items()):
                self[e] = v
        elif hasattr(other, '__iter__'):
            for e, v in other:
                self[e] = v
        else:
            raise AttributeError('argument must be a dict-like object or (key, value) pairs list')
