""" JSON Backend
	read/write handles: units, column comments, aliases, header keywords
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os, inspect, sys
localpath = '/'.join(os.path.abspath(inspect.getfile(inspect.currentframe())).split('/')[:-1])
import numpy as np
from .basebackend import BaseBackend
from ..core.tableheader import TableHeader
from ..table import Table
import json
from json.encoder import *
try:
    from collections import OrderedDict
except ImportError:
    from ..core.odict import odict as OrderedDict


ctypes = {      "S" : "string",
		"b" : "integer",
		"h" : "short",
		"i" : "integer",
		"l" : "long",
		"c" : "char",
		"f" : "float",
		"d" : "double",
		"D" : "complex" }

class Encoder(json.JSONEncoder):
   """ Rewrite the default Encoder to make sure we can parse NaNs and Infinite
   values to follow JSON's standards.
   This class is mainly intended to be used from the dumps or dump call also
   provided in this module. (Their distinction with the originals is from
   optional arguments allowing to set the values that will be used for special
   floats (NaN, Inf...)
   """

   float_specials = dict( nan = 'NaN', inf = 'Infinity', neginf = '-Infinity')

   def set_float_specials(self, nan = None, inf = None):
       """ Set how to parse special values """
       if not (nan is None):
           self.float_specials['nan'] = str(nan)
       if not (inf is None):
           self.float_specials['inf'] = str(inf)
           self.float_specials['neginf'] = '-'+str(inf)

   def iterencode(self, o, _one_shot=False):
        """Encode the given object and yield each string
        representation as available.

        For example::

            for chunk in JSONEncoder().iterencode(bigobject):
                mysocket.write(chunk)

        """
        _one_shot = False
        if self.check_circular:
            markers = {}
        else:
            markers = None
        if self.ensure_ascii:
            _encoder = encode_basestring_ascii
        else:
            _encoder = encode_basestring
        if self.encoding != 'utf-8':
            def _encoder(o, _orig_encoder=_encoder, _encoding=self.encoding):
                if isinstance(o, str):
                    o = o.decode(_encoding)
                return _orig_encoder(o)

        def floatstr(o, allow_nan=self.allow_nan,
                _repr=FLOAT_REPR, _inf=INFINITY, _neginf=-INFINITY):
            # Check for specials.  Note that this type of test is processor
            # and/or platform-specific, so do tests which don't depend on the
            # internals.

            if o != o:
                text = self.float_specials['nan']
            elif o == _inf:
                text = self.float_specials['inf']
            elif o == _neginf:
                text = self.float_specials['neginf']
            else:
                return _repr(o)

            if not allow_nan:
                raise ValueError( "Out of range float values are not JSON compliant: " + repr(o))

            return text

        if (_one_shot and c_make_encoder is not None
                and self.indent is None and not self.sort_keys):
            _iterencode = c_make_encoder(
                markers, self.default, _encoder, self.indent,
                self.key_separator, self.item_separator, self.sort_keys,
                self.skipkeys, self.allow_nan)
        else:
            _iterencode = json.encoder._make_iterencode(
                markers, self.default, _encoder, self.indent, floatstr,
                self.key_separator, self.item_separator, self.sort_keys,
                self.skipkeys, _one_shot)
        return _iterencode(o, 0)

def dumps(obj, nan=None, inf=None):
	""" Dumps obj using the above encoder """
	enc = Encoder()
	enc.set_float_specials(nan=nan, inf=inf)
	return enc.encode(obj)

def loads(f):
	""" Load a json stream into an OrderedDict """
	return json.JSONDecoder(object_pairs_hook=OrderedDict).decode(f)


#==============================================================================
class jsonBackend(BaseBackend):
#==============================================================================
	def __init__(self, col_meta_names = None):
		""" constructor """
		BaseBackend.__init__(self, tableType='json')
		self.col_meta_names = col_meta_names or ['name', 'datatype', 'format', 'unit', 'description', 'null']

	def read(self, filename):
		if hasattr(filename, 'read'):
			unit = filename
		else:
			unit = open(filename, 'r')

		d = loads( unit.read() )

		if not hasattr(filename, 'read'):
			unit.close()

		#generate an empty table and fill it
		tab = Table()
		tab.header = TableHeader(d['tablemeta'])

		aliases = d['aliases']
		colmeta = d['columnmeta']
		colnames = [ k['name'] for k in colmeta ]
		if 'dtype' in k:
			coldt   = [ k['dtype'] for k in colmeta ]
		elif 'datatype' in k:
			coldt   = [ k['datatype'] for k in colmeta ]
		if 'unit' in k:
			colunit = [ k['unit'] for k in colmeta ]
		if 'format' in k:
			colfmt  = [ k['format'] for k in colmeta ]
		if 'description' in k:
			coldesc = [ k['description'] for k in colmeta ]
		if 'null' in k:
			colnull = [ k['null'] for k in colmeta ]

		data = np.rec.fromrecords(d['data'], names=colnames)

		for i,colName in enumerate(colnames):
			tab.add_column(colName, data[colName],
					unit = colunit[i] or '',
					null = colnull[i] or '',
					description = coldesc[i] or '',
					format = colfmt[i],
					dtype = data.dtype[colName] )
		#set aliases
		for k,v in d['aliases'].items():
			tab.set_alias(k, v)


		return tab

	def writeColMeta(self, tab, k):
		col = tab.columns[k].__dict__
		d = OrderedDict()
		for i,v in enumerate(self.col_meta_names):
			if v.lower() in ['dtype', 'datatype']:
				d[v] = ctypes[col.get('dtype').kind]
			elif v.lower() == 'name':
				d[v] = k
			else:
				d[v] = col.get(v, '')
		return d

	def write(self, tab, filename, nan=None, inf=None, **kwargs):
		if hasattr(filename, 'write'):
			unit = filename
		else:
			unit = open(filename, 'w')

		d = OrderedDict()
		d['tablemeta']  = tab.header.__dict__
		d['columnmeta'] = [ self.writeColMeta(tab, k) for k in list(tab.keys()) ]
		d['aliases']    = tab._aliases
		d['data']       = tab.data.tolist()

		unit.write(dumps(d, nan=nan, inf=inf, **kwargs))

		if hasattr(filename, 'write'):
			return unit
		else:
			unit.close()


