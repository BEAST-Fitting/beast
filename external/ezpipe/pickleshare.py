#!/usr/bin/env python

""" PickleShare - a small 'shelve' like datastore with concurrency support

Like shelve, a PickleShareDB object acts like a normal dictionary. Unlike
shelve, many processes can access the database simultaneously. Changing a
value in database is immediately visible to other processes accessing the
same database.

Concurrency is possible because the values are stored in separate files. Hence
the "database" is a directory where *all* files are governed by PickleShare.

Example usage::

    from pickleshare import *
    db = PickleShareDB('~/testpickleshare')
    db.clear()
    print "Should be empty:",db.items()
    db['hello'] = 15
    db['aku ankka'] = [1,2,313]
    db['paths/are/ok/key'] = [1,(5,46)]
    print db.keys()
    del db['aku ankka']

This module is certainly not ZODB, but can be used for low-load
(non-mission-critical) situations where tiny code size trumps the
advanced features of a "real" object database.

Installation guide: easy_install pickleshare

Author: Ville Vainio <vivainio@gmail.com>
License: MIT open source license.

"""
__version__ = '0.3.1'
from .pickleshare_path import path as Path
import stat
import time
import pickle as pickle
try:
    from collections import UserDict as TDict
except ImportError:
    from UserDict import DictMixin as TDict

#import UserDict
import collections
import glob

class PickleShareDB(TDict):
    """ The main 'connection' object for PickleShare database """
    def __init__(self, root, protocol=1):
        """ Return a db object that will manage the specied directory"""
        self.root = Path(root).expanduser().abspath()
        self.protocol = protocol
        if not self.root.isdir():
            self.root.makedirs()
        # cache has { 'key' : (obj, orig_mod_time) }
        self.cache = {}

    def __getitem__(self, key):
        """ db['key'] reading """
        fil = self.root / key
        try:
            mtime = (fil.stat()[stat.ST_MTIME])
        except OSError:
            raise KeyError(key)

        if fil in self.cache and mtime == self.cache[fil][1]:
            return self.cache[fil][0]
        try:
            # The cached item has expired, need to read
            obj = pickle.load(fil.open())
        except:
            raise KeyError(key)

        self.cache[fil] = (obj, mtime)
        return obj

    def __missing__(self, k):
        print('missing key', k)

    def __setitem__(self, key, value):
        """ db['key'] = 5 """
        fil = self.root / key
        parent = fil.parent
        if parent and not parent.isdir():
            parent.makedirs()
        pickle.dump(value, fil.open('w'), protocol=self.protocol)
        try:
            self.cache[fil] = (value, fil.mtime)
        except OSError as e:
            if e.errno != 2:
                raise

    def __delitem__(self, key):
        """ del db["key"] """
        fil = self.root / key
        self.cache.pop(fil, None)
        try:
            fil.remove()
        except OSError:
            # notfound and permission denied are ok - we
            # lost, the other process wins the conflict
            pass

    def _normalized(self, p):
        """ Make a key suitable for user's eyes """
        return str(self.root.relpathto(p)).replace('\\', '/')

    def keys(self, globpat=None):
        """ All keys in DB, or all keys matching a glob"""

        if globpat is None:
            files = self.root.walkfiles()
        else:
            files = [Path(p) for p in glob.glob(self.root / globpat)]
        return [self._normalized(p) for p in files if p.isfile()]

    def uncache(self, *items):
        """ Removes all, or specified items from cache

        Use this after reading a large amount of large objects
        to free up memory, when you won't be needing the objects
        for a while.

        """
        if not items:
            self.cache = {}
        for it in items:
            self.cache.pop(it, None)

    def waitget(self, key, maxwaittime=60 ):
        """ Wait (poll) for a key to get a value

        Will wait for `maxwaittime` seconds before raising a KeyError.
        The call exits normally if the `key` field in db gets a value
        within the timeout period.

        Use this for synchronizing different processes or for ensuring
        that an unfortunately timed "db['key'] = newvalue" operation
        in another process (which causes all 'get' operation to cause a
        KeyError for the duration of pickling) won't screw up your program
        logic.
        """

        wtimes = [0.2] * 3 + [0.5] * 2 + [1]
        tries = 0
        waited = 0
        while 1:
            try:
                val = self[key]
                return val
            except KeyError:
                pass

            if waited > maxwaittime:
                raise KeyError(key)

            time.sleep(wtimes[tries])
            waited += wtimes[tries]
            if (tries < len(wtimes) - 1):
                tries += 1

    def getlink(self, folder):
        """ Get a convenient link for accessing items  """
        return PickleShareLink(self, folder)

    def __repr__(self):
        return "PickleShareDB('%s')" % self.root


class PickleShareLink:
    """ A shortdand for accessing nested PickleShare data conveniently.

    Created through PickleShareDB.getlink(), example::

        lnk = db.getlink('myobjects/test')
        lnk.foo = 2
        lnk.bar = lnk.foo + 5

    """
    def __init__(self, db, keydir ):
        self.__dict__.update(locals())

    def __getattr__(self, key):
        return self.__dict__['db'][self.__dict__['keydir'] + '/' + key]

    def __setattr__(self, key, val):
        self.db[self.keydir + '/' + key] = val

    def __repr__(self):
        db = self.__dict__['db']
        keys = db.keys( self.__dict__['keydir'] + "/*")
        return "<PickleShareLink '%s': %s>" % (
            self.__dict__['keydir'],
            ";".join([Path(k).basename() for k in keys]))


def test():
    db = PickleShareDB('~/testpickleshare')
    db.clear()
    print("Should be empty:", list(db.items()))
    db['hello'] = 15
    db['aku ankka'] = [1, 2, 313]
    db['paths/nest/ok/keyname'] = [1, (5, 46)]
    print(list(db.keys()))
    print(db.keys('paths/nest/ok/k*'))
    print(dict(db))  # snapsot of whole db
    db.uncache()  # frees memory, causes re-reads later

    # shorthand for accessing deeply nested files
    lnk = db.getlink('myobjects/test')
    lnk.foo = 2
    lnk.bar = lnk.foo + 5
    print(lnk.bar)  # 7
