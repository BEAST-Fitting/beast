""" Some tools for manipulating types """
import numpy as np

#==============================================================================
# SOME HELPERS
#==============================================================================
def isNestedInstance(obj, cl):
	""" test is an object is an instance of a deep subclass """
	tree = []
	for k in cl.__subclasses__():
		tree+=k.__subclasses__()
	tree += cl.__subclasses__() + [ cl ]
	print(tree)
	return  issubclass(obj.__class__, tuple(tree))

def name2dtype(name):
	""" returns the numpy type from a name """
	if (name, np.dtype):
		return name
	else:
		return np.typeDict[name]

def append_field(rec, data, dtype=None, position=None):
	""" Append field to a recarray """
	newdtype = rec.dtype.descr
	if position is None:
		newdtype.append(dtype)
	else:
		newdtype.insert(position, dtype)

	newdtype = np.dtype(newdtype)
	newrec = np.recarray(rec.shape, dtype=newdtype)
	for field in rec.dtype.fields:
		newrec[field] = rec[field]
	newrec[dtype[0]] = data
	return newrec

def drop_fields(rec, names):
	""" Drop a field from a recarray """
	names = set(names)
	newdtype = np.dtype([(name, rec.dtype[name]) for name in rec.dtype.names if name not in names])
	newrec = np.recarray(rec.shape, dtype=newdtype)
	for field in newdtype.fields:
		newrec[field] = rec[field]
	return newrec

def smart_dtype(dtype):
	if dtype.subdtype:
		return dtype.subdtype[0].type
	else:
		return dtype.type

def format_length(format):
	if '.' in format:
		return int(format.split('.')[0])
	else:
		return int(format[:-1])

default_format = {}
default_format[None.__class__] = '16.9e'
default_format[np.bool_] = '1i'
default_format[np.int8] = '3i'
default_format[np.uint8] = '3i'
default_format[np.int16] = '5i'
default_format[np.uint16] = '5i'
default_format[np.int32] = '12i'
default_format[np.uint32] = '12i'
default_format[np.int64] = '22i'
default_format[np.uint64] = '23i'
default_format[np.float16] = '16.8e'
default_format[np.float32] = '16.8e'
default_format[np.float64] = '25.17e'
default_format[np.str] = 's'
default_format[np.string_] = 's'
default_format[np.uint8] = 's'
default_format[str] = 's'
default_format[np.unicode_] = 's'
