import numpy as np


#==============================================================================
class TableHeader(object):
    """ this class defines the context of a Table """
#==============================================================================
    def __init__(self, dic=None, *args, **kwargs):
        """ constructor """
        for k, v in kwargs:
            self.__setattr__(k, v)
        if (isinstance(dic, dict)) or (isinstance(dic, TableHeader)):
            for k, v in dic.iteritems():
                self.__setattr__(k, v)
        if not 'NAME' in kwargs:
            self.__setattr__('NAME', 'Noname')

    def __setattr__(self, attribute, value):
        """ set attribute """
        if (attribute.lower() in ['description', 'comment', 'history']) & (attribute in self):
                self.__dict__[attribute] = str(self[attribute]) + '\n' + str(value)
        else:
            self.__dict__[attribute] = value

    def __setitem__(self, item, val):
        """get item"""
        self.__setattr__(item, val)

    def __getitem__(self, item):
        """get item"""
        return self.__dict__[item]

    def __contains__(self, item):
        return self.__dict__.__contains__(item)

    def __repr__(self):
        """ representation """
        s = 'Table Header\n'
        if len(self.__dict__) == 0:
            s += "(Empty)"
        else:
            keys = np.sort(self.keys())
            for k in keys:
                if self[k] is not None:
                    vals = str(self[k]).split('\n')
                    if len(vals) == 1:
                        s += '%20s\t%s\n' % (k.upper(), self[k])
                    else:
                        s += '%20s\t%s\n' % (k.upper(), vals[0])
                        for kval in vals[1:]:
                            s += '%20s\t%s\n' % ('', kval)
                else:
                    s += '%20s\t%s\n' % (k.upper(), None)

        return s

    def get(self, k, default=None):
        if k in self.keys():
            return self[k]
        else:
            return default

    def keys(self):
        """return registered keywords list"""
        return self.__dict__.keys()

    def copy(self):
        """return a copy of the header"""
        return TableHeader(self.__dict__)

    def __iter__(self):
        return self.__dict__.__iter__()

    def iterkeys(self):
        return self.__dict__.iterkeys()

    def itervalues(self):
        return self.__dict__.itervalues()

    def iteritems(self):
        return self.__dict__.iteritems()

    def items(self):
        return self.__dict__.items()

    def pop(self, name):
        return self.__dict__.pop(name)
