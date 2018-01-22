""" Some helpers for sql related operations """
from .decorators import memoize
try:
    from collections import OrderedDict
except ImportError:
    from .odict import odict as OrderedDict
import numpy as np
import subprocess

# SQLite
try:
    import sqlite3
    sqlite3_installed = True
except:
    sqlite3_installed = False


# Type conversion dictionary
type_dict = {}
type_dict[np.bool_] = "BOOL"
type_dict[np.uint8] = "TINYINT"
type_dict[np.uint16] = "SMALLINT"
type_dict[np.uint32] = "INT"
type_dict[np.uint64] = "BIGINT"
type_dict[np.int8] = "TINYINT"
type_dict[np.int16] = "SMALLINT"
type_dict[np.int32] = "INT"
type_dict[np.int64] = "BIGINT"
type_dict[np.float16] = "FLOAT"
type_dict[np.float32] = "FLOAT"
type_dict[np.float64] = "DOUBLE PRECISION"
type_dict[np.str] = "TEXT"
type_dict[np.string_] = "TEXT"
type_dict[str] = "TEXT"


class remotePsqlConfig(object):
	""" Impose the conditions for the remote sql configuration """
	def __init__(self, gate, gate_user, db_name, db_user, db_host, db_cmd, **kwargs):
		self.gate        = gate
		self.gate_user   = gate_user
		self.db_name     = db_name
		self.db_user     = db_user
		self.db_host     = db_host
		self.db_cmd      = db_cmd
		for k,v in kwargs.items():
			self.__setattr__(k,v)

	def __repr__(self):
		return self.__dict__.__repr__()

	def iteritems(self):
		return iter(self.__dict__.items())

	def items(self):
		return list(self.__dict__.items())

	def values(self):
		return list(self.__dict__.values())

	def __getitem__(self, k):
		return self.__dict__[k]


class remotePsql(object):
	""" Very simple parsing of remote queries through networks
		Motivated by GreenPlum remote usage through ssh
	"""
	def __init__(self, config):
		for k,v in config.items():
			self.__dict__[k] = v
		self.__res__ = memoize(self.__runQuery__)

	@property
	def __protocol__(self):
		if self.gate is None:
			return ''

		a = self.gate.split('://')
		if (len(a) == 2):
			if a[0] == 'ssh':
				p = '%s %s@%s' % (a[0], self.gate_user, a[1])
			else:
				assert(False), 'Protocol %s not supported yet.' % a[0]
		else:
			p = ''
		return p

	@property
	def __precmd__(self):
		return (self.__protocol__, '%s -h %s -U %s %s -A -F, -c ' % (self.db_cmd, self.db_host, self.db_user, self.db_name) )

	@property
	def log(self):
		return list(self.__res__.keys())

	def clearlog(self):
		""" Clear memoized results """
		self.__res__.clear()

	def handle_special(self, q, p):
		if q[0] == '\\':
			NULL = p.stdout.readline()
			if len(q.split()) > 1:
				names = p.stdout.readline().split(',')
				r = np.recfromtxt(p.stdout, skip_footer=2, delimiter=',', names=names)
			else:
				names = p.stdout.readline().split(',')
				r = np.recfromtxt(p.stdout, skip_footer=1, delimiter=',', names=names)

		else:
			names = p.stdout.readline().split(',')
			r = np.recfromtxt(p.stdout, skip_footer=1, delimiter=',', names=names)
		return r

	def __runQuery__(self, txt):
		pcmd = self.__precmd__
		cmd = ' %s "%s" ' % (pcmd[0], pcmd[1]+"'"+txt+"'")
		p = subprocess.Popen(cmd , stdout=subprocess.PIPE, shell=True)
		r = self.handle_special(txt, p)
		return r

	def __call__(self, q):
		"""
		Execute the query q on the server and returns the associated
		recarray

		Calls are memoized so that identical queries return the first
		call result
		"""
		return self.__res__(q)

	__getitem__ = __call__
