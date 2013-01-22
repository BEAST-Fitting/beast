""" This file implements a Table class 
	that is designed to be the basis of any format
	It has been inspired by the atpy package

EXAMPLE USAGE:
	> import mytables
	> t = mytables.load('toto.fits')
	Auto-detected input type: fits
	## Any registered extention will work
	> t() or > t.info()
	## print information about the current table (col, part of the header...)
	> t[0] or t.getRow(0)
	## returns the first row of the table
	> t[0]['x']
	## rows are records, which means values keeps their names
	> t[2:9] or t[ [0,2,7] ]
	## returns a the [2:9] slice or the selected rows
	> t['x'] or t.getCol('x')
	## returns the column x
	> t['col'] = range(10)
	## this will add a new column named 'col' with the corresponding values
	> t.addCol(range(10), name='x', unit='myUnit', description='unit is arbitrary')
	## this also add a new column with more details	
	> t['col', unit='myUnit', description='a comment'] = range(10)
	## corresponds to the addCol expression
	> t.nrows
	## contains the number of rows in the table
	> t.delCol('col')
	> t.setComment('x', 'comment')
	> t.setunit('x', 'myunit')

DESCRIPTION:
    >t = Table()	 
    	t()		  perfom simple column operations and return results
	t.addCol()        add a new Column to the current table
	t.delCol()	  delete a given column from the table
	t.disp()       	  pretty print (part of) the table 
	t.getCol()	  returns a given column (tuple of columns)
	t.getRow()	  returns a given row (tuple of rows as record arrays)
	t.data		  contains a dictionary of columns
	t.evalexpr()   	  let you do some simple operations on the table using
		    	  column names as variables (incl. math symbols) and
		    	  external variables as well
		    	  TODO: need to find a way to extend math functions to
		    	  arrays
	t.extract()    	  Returns a sub-table based on a given selection of
		    	  lines or fields
	t.header	  header of the table (each column has its own header)
	t.info()	  print information about the table
	t.nrows		  contains the number of rows of the table
	t.pop()		  returns a given column and delete it from the table
	t.read()	  fill the table from a given file
			  (from either a regitered format or a given tableManager object)
	t.selectWhere()   is equivalent to Table[where()] and return a table
		    	  with only matching conditions (can also restrict 
		    	  selected	fields)
	t.setComment()	  add a comment to a given column
	t.setUnit()	  set the unit of a given column
	t.source	  if set, it contains the source file object related to
			  the source format of the table (e.g. in hd5 formats)
	t.stats()	  returns statistics of the different columns (requires mystats)
	t.where()      	  traditional where function based on evalexpr
	t.write()	  write the current table to a given file (including
			  format (from either a regitered format or a given tableManager
			  object)

     *Common functions for a transparent dictionary proxy of the columns*
	t.items(), t.iteritems, t.iterkeys(), t.itervalues, t.keys()
	
REQUIREMENTS:
	numpy		mandatory
	pyfits 		for using fits extensions
	pytables	for using hdf5 extensions
	mystats		for computing statistics on the table columns (t.stats())
	mypickleshare	for distributed tables
	sqlite3		for sqltables

HISTORY:
	0.1 -- 04/03/2011, MF
	Definition of the classes related to a Table object: 
			Table, 
			TableHeader,
			TableColumn,
			TableColumnHeader, 
			TableManager,
		Default registered formats:
			asciiManager (dat, txt, ascii) 
			csvManager   (csv)
			fitsManager  (fits) -- requires pyfits
		Import/Export manages also column description/comments & units
		in all the default registered formats
	0.2 -- 07/03/2011, MF
	New formats added:
			hd5Manager (hd5, hdf5) -- requires pytables
			latexManager (tex) only writing
	0.3 -- 07/03/2011, MF
	New feature added:
		Tables are now callable to produce simple operations
		you can now generate a column based on simple scalar operations
		or also +,-,*,/ operation.
		example:  > r = table('x+2*y')  or > r = table('a/b')
		but this is not tested to follow operation priorities, besides
		it could not perform complex operations sur as >table('x/(a+b)')
	0.3.1 -- 09/03/2011, MF
	Debug: Fits export, 
			order in the units, comments fixed
			string export fixed
			multiple comment & history keywords export fixed
	0.3.2 -- 01/04/2011, MF
	New features added
		Added join & match functions to the module
	0.3.3 -- 12/04/2011, MF
	Bug correction: opening HD5 tables that does not contain UNITS...
	0.3.4 -- 09/05/2011, MF
	Bug correction: opening HD5 tables and tablename issues leading
		to errors in reading table header
	0.3.5 -- 31/08/2011, MF
	Bug correction: ascii files issues
			delimiter was not always taken into account
			now you can force a given line to be the column
			description forceHealLine keyword
	0.3.6 -- 09/09/2011, MF
	Bug correction: ascii/csv files issues when mixing with text
			replacing loadfromtxt call by a recfrom<> call which
			does the job better
			if a column name already exists it will add a suffix
			instead of overwriting the column
	0.4 -- 09/09/2011, mf
	added features:
		disp        -- pretty print (part of) the table 
	  	__str__	    -- returns disp with default arguments
		evalexpr    -- let you do some simple operations on the table using
			       column names as variables (incl. math symbols) and
			       external variables as well
			       TODO: need to find a way to extend math functions to
			       arrays
		where       -- traditional where function based on evalexpr
		selectWhere -- is equivalent to Table[where()] and return a table
			       with only matching conditions (can also restrict 
			       selected	fields)
		extract     -- Returns a sub-table based on a given selection of
			       lines or fields
		__call__ updated to use evalexpr
		__getitem__ returns uses getRow for a single numerical value
				    uses getCol for a string
				    uses extract for an iterable
	0.4.1 -- 15/09/2011, mf
	added features:
		stats        -- returns a table with statistics over the
				selected columns (mean, std, hpd, quantiles...)
	0.5.0-- 13/10/2011, mf
	added features:
		distributedTable-- new table class that let you use disk space
				instead of RAM and allow multiple concurrent
				access
				(Still in development)
				Currently should work for reading
				Saving changes are currently not fully
				operational, data headers are not updated in
				realtime.
		Generating tables from recarray.
		Added a idl save file reader
	0.6.0 -- 30/03/2012, mf
	added features;
		sqltables are now managed through sqlManager nd sqlite3
		you can use it to load or generate sqlite tables. 
"""
import warnings
import numpy
import os
from math import *
import cStringIO,operator
import sys

__version__ = '0.6.0'
__author__  = 'M. Fouesneau'

global _extensions


#==============================================================================
# COMMON functions
#==============================================================================

def register_extension(manager=None,extName=None,  
			readerFunction=None, writerFunction=None, 
			override=False):
	""" Register a table type manager 
	    inputs: 
		can be either a TableManager object or individual functions
		manager --  TableManager Object
		extName --  register a particular extension if provided
		            (required if not giving a TableManager as input)
		readerFunction -- function that is required to do the reading
		                  part
				  this function should return a data argument
				  and a header argument
				  both should be at least dictionnaries but it
				  is recommanded to use the TableHeader object
				  and also to provide data with TableColumn
				  objects.
		writerFunction -- function that does the writing part
				  At least will receive a data and header
				  arguments (that can be accessed as dict
				  objects) and should also add **kwargs to
				  avoid unattended keyword arguments
	    
	    outputs:
		None
		It registers the manager in the module libs.
		
	    keywords:
		override -- redefine the file manager if already registered. 
	"""
	extName = extName or manager.tableType
	if not extName in _extensions or override:
		if manager != None:	
			assert(isinstance(manager, TableManager))
			_extensions[extName] = manager
		else:
			assert((readerFunction != None) & (writerFunction != None))
			_extensions[extName] = TableManager(extName, 
								readerFunction, 
								writerFunction  )
	else:
		raise Exception("Type %s is already defined" % ttype)

def _isiterable(val):
	""" determine if a given object is either iterable
		returns a boolean
	"""
	return hasattr(val, '__iter__')
	"""
	try:
		void = val.__getattribute__('__iter__')
		return True
	except:
		return False
	"""

def _determine_type(_extensions, string, verbose=True):
	""" 
	Determine the type of a table from its extension and try to give the
	point to the appropriate registered extension
	"""
	if type(string) <> str:
		raise Exception('Could not determine input type (non-string input)')
	s = string.lower()
	if not '.' in s:
		extension = s
	else:
		extension = s.split('.')[-1]
		if extension in ['gz', 'bzip2', 'Z']:
			raise Exception('Compressed files are not managed yet.')
		elif extension in ['sav', 'idl', 'idlsav']:
			raise Exception('Warning: IDL save files must be explicitely requested with "idlload" independent function.')
			
	if extension in _extensions:
		tableType = _extensions[extension]
		if verbose:
			print "Auto-detected type: %s" % extension
	else:
        	raise Exception('Could not determine input type for extension %s' % extension)
	return tableType


#==============================================================================
class Table(object): 
	""" This class implements a Table object which aims at being able to
	manage dataset table independently of its storage format.
	
	Tables allows row and column access using the __getitem__ intrinsic
	method. Columns can be added to or deleted from the table. 
	Any columns can be described by its own header which includes units.

	Besides Table object can be considered as dict objects and so are
	iterable over the columns.

	reading and writing methods are not directly defined in this class but
	apart using TableManager objects that are then registered with the
	"register_extension" function. 
	"""
#==============================================================================
	def __init__(self, iterable=None, colNames=None, header=None,
			verbose=True, name=None, *args,**kwargs):
		""" Constructor
			This can be called without any arguments.
			Without any arguments, you will generate an empty table
			that would mbe meant to accept new data later on.
			It can also convert iterable datasets (dict, array...)
			to Table objects
		inputs:
			iterable -- iterable object that stores the data
			            iteration over data dimensions (not rows)
			colNames -- [iterable] 
		                    if specified, defines the name that are used
			            for the columns
			header   -- [dict]
			            Any key of the dictionary will be stored
				    into the Table header.
			name 	 -- Table Name

		"""
		self.source = None
		self.nrows = 0
		self.ncols = 0
		#Table header
		self.header = TableHeader(header, **kwargs)
		if name != None:
			self.setName(name)
		#Gen Columns
		self.data = {}
		if iterable != None:
			assert(_isiterable(iterable))
			if isinstance(iterable, dict):
				for kCol in iterable:
					self.addCol( TableColumn(iterable[kCol], name=kCol) )
			elif isinstance(iterable, numpy.core.records.recarray):
				names = iterable.dtype.names
				for kName in names:
					self.addCol(TableColumn(iterable[kName], name=kName))
			else:	
				i = 0
				for kCol in iterable:
					i += 1
					if isinstance(kCol, TableColumn):
						self.addCol(TableColumn(kCol))
					else:
						colName = 'Col%s' % i
						if (colNames != None) & (numpy.size(colNames) >= i):
							colName = colNames[i-1]
							self.addCol(TableColumn(kCol,name=colName))
		elif colNames != None:
			assert(_isiterable(colNames))
			for kCol in colNames:
				self.addCol( TableColumn(None, name=kCol) )

	def __del__(self):
		if self.source != None:
			#print "Closing File %s" % self.source
			self.source.close()
		del self

	def disp(self, idx=None, fields=None, ret=False):
		""" Pretty print the table content
			you can select the table parts to display using idx to
			select the rows and fields to only display some columns
			(ret is only for insternal use)"""
		if fields==None:
			fields= numpy.sort(self.keys())


		if idx == None:
			if self.nrows < 10:
				rows = [ [ str(self[k][rk]) for k in fields ] for rk in range(self.nrows)] 
			else:
				_idx = range(5)	
				rows = [ [ str(self[k][rk]) for k in fields ] for rk in _idx] 
				rows += [ ['...' for k in range(len(fields)) ] ]
				rows += [ [ str(self[k][-rk]) for k in fields ] for rk in _idx] 
		elif isinstance(idx, slice):
			_idx = range(idx.start, idx.stop, idx.step or 1)
			rows = [ [ str(self[k][rk]) for k in fields ] for rk in _idx] 
		else:
			rows = [ [ str(self[k][rk]) for k in fields ] for rk in idx] 
		units = [ '('+str(self[k].header['unit'] or '')+')' for k in fields ]
		out = indent([fields]+[units]+rows, hasHeader=True, hasUnits=True)
		if ret == True :
			return out
		else:
			print out
	
	def addCol(self, array, name=None, **kwargs):
		""" Add individual column to the table
		    __setitem__ method is a shortcut to this function
			inputs:	
				array -- [ array | list | TableColumn ]
					 data to add, this can be multi-dimensional.
					 multi-dimentional columns storage is only 
					 managed by HDF5 format
			keywords:
				name -- [ string ]
					name of the column
				
				**kwargs refers to TableColumn arguments
		"""
		name = name or 'None'
		if len(self.data) > 0:
			shape = numpy.shape(array)
			sshape = self.data[self.data.keys()[0]].shape
			if shape[0] != sshape[0]:
				raise Exception("Column size mismatch, expecting %i rows, found %i" \
						% (sshape[0], shape[0]) )
		else:
			if numpy.size(numpy.shape(array)) > 0: 
				self.nrows = numpy.shape(array)[0]
		if isinstance(array, TableColumn):
			if name == 'None':
				name = array.header.name
			
				array.header['name'] = name
			self.data[name] = array
		else:
			self.data[name] = TableColumn(array, name=name, **kwargs)

	def delCol(self, name):
		""" Delete Table column 
			inputs:
				name -- [ string ]
					Column to delete
		"""
		cCol = self.data.pop(name)
		del cCol
	
	def pop(self, name):
		""" Pop Table column 
			inputs:
				name -- [ string ]
					Column to delete
			outputs:
				poped column
		"""
		return self.data.pop(name)

	def evalexpr(self, expr, exprvars=None, start=None, stop=None, step=None):
		""" evaluate expression based on the data and external variables
		    all numpy.math expression can be used (log, exp, pi...)
		"""
		_expr = expr
		keys = numpy.asarray(self.keys())
		lenkeys = numpy.asarray([ len(k) for k in keys ])
		keys = keys[ numpy.argsort(lenkeys)[::-1] ]
		for k in keys:
			_expr = _expr.replace(k, 'self.data["%s"]'%k)
		if exprvars != None:
			assert(isinstance(exprvars,dict)),"Expecting dictionary as condvars"
			for k in exprvars.keys():
				_expr = _expr.replace(k, 'exprvars["%s"]'%k)
		return eval(_expr)

	def extract(self, ind=None, fields=None):
		""" Generates a subtable containing lines defined by ind and
		fields. If fields is not provided, then the full set is
		retrieved. Header is copied into the new table. """
		if fields == None:
			fields = self.keys()
		if ind == None:
			ind = range(self.nrows)

		hdrs =[ TableColumn(self.getCol(k)).header         \
			for k in numpy.asarray(fields) if
			getattr(TableColumn(self.getCol(k)), 'header', False) ]
 
		data = [ TableColumn(self.getCol(k)[ind],          \
				name=self.getCol(k).header['name'],\
				unit=self.getCol(k).header['unit'],\
				description = self.getCol(k).header['description'],\
				dtype = self.getCol(k).header['dtype'], \
				null = self.getCol(k).header['null'],   \
				format=self.getCol(k).header['format']	)\
			for k in numpy.asarray(fields)                ]
		return Table(data, name=self.header['NAME']+'_extract')

	def selectWhere(self, condition, condvars=None, fields=None, start=None, stop=None, step=None):
		""" Read table data fulfilling the given `condition`.
			Only the rows fulfilling the `condition` are included in the result.
		"""
		ind = self.where(condition, condvars, start=start, stop=stop, step=step)
		return self.extract(ind, fields=fields)

	def where(self, condition, condvars=None, start=None, stop=None, step=None):
		""" Read table data fulfilling the given `condition`.
			Only the rows fulfilling the `condition` are included in the result.
		"""
		ind = numpy.where(self.evalexpr(condition, condvars, start=start, stop=stop, step=step ))
		return ind
			


	def getCol(self, name):
		""" Returns one/multiple Column(s)
			inputs:
				name -- [ list ]
					name of the columns to extract
		"""
		if _isiterable(name):
			return [ self.data[k] for k in name ]
		else:
			return self.data[name]	

	def getRow(self, idx, fields=None):
		""" Returns one/multiple row(s)
			inputs:
				idx -- [ list ]
					indexes of the rows to extract
		"""
		if fields == None:
			fields = self.keys()
		try:
			formats = [self.data[k].header.dtype for k in self]
		except:
			formats = [self.data[k].dtype for k in self]
		if _isiterable(idx):
			row = [ tuple([ self.data[k][kidx] for k in fields ]) for kidx in idx ] 
		else:
			row = tuple([ self.data[k][idx] for k in fields ])
		try:
			return numpy.array(row, {'names':fields, 'formats': formats})
		except ValueError:
			return row

	def __getslice__(self, i,j):
		return self.getRow(numpy.arange(i,j))

	def __iter__(self):
		return self.data.__iter__()
	
	def iterkeys(self):
		return self.data.iterkeys()

	def itervalues(self):
		return self.data.itervalues()

	def iteritems(self):
		return self.data.iteritems()

	def items(self):
		return self.data.items()

	def __getitem__(self, item, **kwargs):
		if _isiterable(item):
			ctype = item[0] 
			if isinstance(ctype, str):
				return self.extract(None, fields=item, **kwargs)
			elif isinstance(ctype, int):
				return self.extract(item, **kwargs)
		else:
			ctype = item
		if isinstance(ctype, str):
			return self.getCol(item)
		elif isinstance(ctype, int):
			return self.getRow(item, **kwargs)
		else:
			return self.extract(item, **kwargs)

	def __setitem__(self, key, val, name=None, **kwargs ):
		""" shortcut to addCol function """
		self.addCol(val, name=key, **kwargs)
	
	def __repr__(self):
		s = 'Table: %s,  nrows=%i, ncols=%i, ' % (self.header['NAME'], self.nrows, len(self.data))
		s += object.__repr__(self)
		return s

	def __str__(self):
		return self.disp(ret=True)

	def __call__(self, op=None):
		if op == None:
			self.info()
		else:
			return self.evalexpr(op)
			"""
			if ('(' in op) or (')' in op):
				print 'Parenthesis detected.\n Only simple operations are managed so far.'
				print 'You may want to compute each parentheses pair individually'
				raise Exception('Too complex operation requested')
			if op in self.keys():
				return self[op]
			elif len(op.split('+'))>1:
				txt = op.split('+')
				if len(txt[1].split('+'))>1:
					return TableColumn(self(txt[0])+self('+'.join(txt[1])), name=op)
				else: 
					return TableColumn(self(txt[0])+self(txt[1]), name=op)
			elif len(op.split('-'))>1:
				txt = op.split('-')
				if len(txt[1].split('-'))>1:
					return TableColumn(self(txt[0])-self('-'.join(txt[1])), name=op)
				else: 
					return TableColumn(self(txt[0])-self(txt[1]), name=op)
			elif len(op.split('*'))>1:
				txt = op.split('*')
				if len(txt[1].split('*'))>1:
					return TableColumn(self(txt[0])*self('*'.join(txt[1])), name=op)
				else: 
					return TableColumn(self(txt[0])*self(txt[1]), name=op)
			elif len(op.split('/'))>1:
				txt = op.split('/')
				if len(txt[1].split('/'))>1:
					return TableColumn(self(txt[0])/self('/'.join(txt[1])), name=op)
				else: 
					return TableColumn(self(txt[0])/self(txt[1]), name=op)
			else:
				try:
					return float(op)
				except:
					raise Exception("Operation is not understood")
			"""

	def info(self):
		print self.header
		print "Table contains: %i row(s) in %i column(s)" % (self.nrows, len(self.data))
		if self.source != None:
			print "Table linked to source file"
		keys = numpy.sort(self.keys())
		for k in keys:
			print self[k].__repr__().split('Table')[0][:-1].replace(',', '')
	
	def setUnit(self, colName, unit):
		""" Set the unit of a column referenced by its name """
		self.data[colName].header['unit'] = unit

	def setComment(self, colName, comment):
		""" Set the comment of a column referenced by its name """
		self.data[colName].header['description'] = comment

	def setName(self, name):
		""" Set Table name """
		self.header['NAME'] = name

	def keys(self):
		""" returns the data key names """
		return self.data.keys()
	def has_key(self, k):
		return self.data.has_key(k)

	def read(self, filename, type=None, manager=None, silent=False, **kwargs):
		""" This function is a general function aiming at reading files
		it uses the registered extensions to use the appropriate reading
		function.
		inputs:
			filename -- [ string ]
			            adress of the source file to be used
		keywords:
			type     -- [ string ]
				    if specified, this will force the function
				    to use considered the file to be of this
				    type.
			manager  -- [ TableManager ]
				    If specified, it will use this format
				    manager (even if not registered)

			**kwargs are sent to the TableManager.read function
		"""
		if manager == None:
			if type == None:
				manager = _determine_type(_extensions, filename, verbose=not silent)
			else:
				manager = _determine_type(_extensions, type, verbose=False)
		data, description = manager.read( filename, **kwargs)
		if isinstance(data, dict): 
			for k in data: self.addCol(data[k], name=k )
		else:
			assert(_isiterable(data))
			for k in data: self.addCol(k)
		if '_sourceFile' in description:
			self.source = description.pop('_sourceFile')
		self.header = description
		self.header['SOURCE'] = os.path.realpath(filename)

	def write(self, filename, type=None, manager=None, silent=False, **kwargs):
		""" This function is a general function aiming at writing files
		it uses the registered extensions to use the appropriate writing
		function.
		inputs:
			filename -- [ string ]
			            adress of the source file to be used
		outputs:
			optional returns from the TableManager write function

		keywords:
			type     -- [ string ]
				    if specified, this will force the function
				    to use considered the file to be of this
				    type.
			manager  -- [ TableManager ]
				    If specified, it will use this format
				    manager (even if not registered)

			**kwargs are sent to the TableManager.write function
		"""
		if manager == None:
			if type == None:
				manager = _determine_type(_extensions, filename, verbose=not silent)
			else:
				manager = _determine_type(_extensions, type, verbose=not silent)
		comments = [self[k].header['description'] for k in self]
		units    = [self[k].header['unit'] for k in self]
		return manager.write(self.data, header=self.header,
					output=filename, comments=comments, 
					units=units, silent=silent, **kwargs) 
	try:
		def stats(self, fields=None, val=['mean', 'min', 'max', 'std', 'q', 'n', 'hpd']):	
			"""returns a table with statistics over the selected columns (mean, std, hpd, quantiles...)
				inputs: 
					fields  -- selected columns to use
					val     -- list of statistics to return
							(see list below)
			
				output:
					table object containing the corresponding statistics (one line per field)

				Available statistics:
				
					hpd   -- highest probable density at 95%
					mean  -- mean value
					min   -- minimal value
					max   -- maximal value
					n     -- sample length
					q     -- quantiles: 2.5, 25, 50, 75, 97.5
					std   -- standard deviation
			"""
			import mystats
			from numpy import nan
			if fields==None:
				fields= self.keys()
			n = len(fields)

			def mean(a):
				try:
					return numpy.mean(a)
				except TypeError:
					return numpy.nan
			def min(a):
				try:
					return numpy.min(a)
				except TypeError:
					return numpy.nan
			def max(a):
				try:
					return numpy.max(a)
				except TypeError:
					return numpy.nan
			def std(a):
				try:
					return numpy.std(a)
				except TypeError:
					return numpy.nan
			def quantiles(a):
				try:
					q = mystats.quantiles(a)
					return [ q[k] for k in q]
				except:
					return [numpy.nan]*5
			def hpd(a):
				try:
					return mystats.hpd(a, 0.05)
				except:
					return [numpy.nan]*2

			s = Table(name='%s, Statistics' % self.header['NAME'])

			s.addCol(fields, name='Field')
			if 'n' in val:
				n = numpy.array([ len(self[k])  for k in fields ])
				s.addCol(n, name='n',  comment='sample length')
				del n
			if 'mean' in val:
				mean =  numpy.array([ mean(self[k]) for k in fields ])
				s.addCol(mean, name='mean',  comment='mean value')
				del mean
			if 'min' in val:
				_min =  numpy.array([ min(self[k]) for k in fields ])
				s.addCol(_min, name='min',  comment='minimum value')
				del _min
			if 'max' in val:
				_max =  numpy.array([ max(self[k]) for k in fields ])
				s.addCol(_max, name='max',  comment='maximum value')
				del _max
			if 'std' in val:
				_std =  numpy.array([ std(self[k]) for k in fields ])
				s.addCol(_std, name='std',  comment='Standard deviation')
				del _std
			if 'q' in val:
				q =  numpy.array([ quantiles(self[k]) for k in fields ]).T
				s.addCol(q[0], name='q2.5',  comment='2.5 quantile')
				s.addCol(q[1], name='q25',  comment='25 quantile')
				s.addCol(q[2], name='q50',  comment='50 quantile')
				s.addCol(q[3], name='q75',  comment='75 quantile')
				s.addCol(q[4], name='q97.5',  comment='97.5 quantile')
				del q
			if 'hpd' in val:
				h =  numpy.array([ hpd(self[k]) for k in fields ]).T
				s.addCol(h[0], name='hpdmin',  comment='min range of 95% highest probable density')
				s.addCol(h[1], name='hpdmax',  comment='min range of 95% highest probable density')

			return s
	
	except ImportError:
		pass

try:
	import mypickleshare as mps		
	#==============================================================================
	class distributedTable(Table): 
		""" This class implements a Table object that tends to have the minimum
		memory impact. Useful for huge tables"""

	#==============================================================================
		def __init__(self, storage, iterable=None, colNames=None, header=None,
				verbose=True, name=None, *args,**kwargs):
			""" Constructor
				This can be called without any arguments.
				Without any arguments, you will generate an empty table
				that would mbe meant to accept new data later on.
				It can also convert iterable datasets (dict, array...)
				to Table objects
			inputs:
				iterable -- iterable object that stores the data
					    iteration over data dimensions (not rows)
				colNames -- [iterable] 
					    if specified, defines the name that are used
					    for the columns
				header   -- [dict]
					    Any key of the dictionary will be stored
					    into the Table header.
				name 	 -- Table Name

			"""
			self.data = mps.PickleShareDB(storage)
			if len(self.data) > 0:
				self.source = None
				self.nrows = numpy.shape(self.data[self.keys()[0]])[0]
				self.ncols = len(self.data)
				if not '__HEADER__' in self:
					self.header = TableHeader(header, **kwargs)
				else:
					self.header = self.data['__HEADER__']
				self.header['STORAGE'] = storage
			else:
				self.source = None
				self.nrows = 0
				self.ncols = 0
				#Table header
				self.header = TableHeader(header, **kwargs)
				if name != None:
					self.setName(name)
				#Gen Columns
				if iterable != None:
					assert(_isiterable(iterable))
					if isinstance(iterable, dict):
						for kCol in iterable:
							self.addCol( TableColumn(iterable[kCol], name=kCol) )
					elif isinstance(iterable, numpy.core.records.recarray):
						print "Converting from recarray not implemented yet"
					else:	
						i = 0
						for kCol in iterable:
							i += 1
							if isinstance(kCol, TableColumn):
								self.addCol(TableColumn(kCol))
							else:
								colName = 'Col%s' % i
								if (colNames != None) & (numpy.size(colNames) >= i):
									colName = colNames[i-1]
									self.addCol(TableColumn(kCol,name=colName))
				elif colNames != None:
					assert(_isiterable(colNames))
					for kCol in colNames:
						self.addCol( TableColumn(None, name=kCol) )
		def __del__(self):
			self.saveheader()
			Table.__del__(self)

		def saveheader(self):
			if self.header != None:
				self.data['__HEADER__'] = self.header

		def keys(self):
			""" returns the data key names """
			k = self.data.keys()
			return [ l for l in k if (l[:2] != '__') ]

		def has_key(self, k):
			return self.data.has_key(k)

		def __iter__(self):
			return (l for l in self.data.__iter__() if l[:2] != '__')
except ImportError:
	print 'Distributed Table are not avalable'
#==============================================================================
class TableColumnHeader(object): 
	""" Manage how columns are described """
#==============================================================================
	def __init__(self, name, dtype, unit=None, description=None, null=None, format=None):
		self.__dict__['dtype'] = dtype
		self.__dict__['name'] = name
		self.unit = unit
		self.description = description
		self.__dict__['null'] = null
		self.format = format

	def __getitem__(self, k):
		return self.__dict__[k]

	def __setitem__(self, k, v):
		self.__dict__[k] = v

	def __setattr__(self, attribute, value):
		#if attribute in ['unit', 'description', 'format']:
		self.__dict__[attribute] = value
		#elif attribute in ['null', 'dtype']:
		#	raise Exception("Cannot change %s through the columns object" % attribute)
		#else:
		#	raise AttributeError(attribute)

	def __repr__(self):
		s = "type=%s" % str(self.dtype)
		if self.unit:
			s += ", unit=%s" % str(self.unit)
		if self.null:
			s += ", null=%s" % str(self.null)
		if self.description:
			s +=", description=%s" % self.description
		return s

	def __eq__(self, other):
		if other == None:
			return False
		if self.dtype <> other.dtype:
			return False
		if self.unit <> other.unit:
			return False
		if self.description <> other.description:
			return False
		if self.null <> other.null:
			if numpy.isnan(self.null):
				if not numpy.isnan(other.null):
					return False
			else:
				return False
		if self.format <> other.format:
			return False
		return True

	def __ne__(self, other):
		return not self.__eq__(other)

	def copy(self):
		"""return a copy of the header"""
		return TableColumnHeader(self.name, self.dtype, self.unit,
						self.description,
						self.null,
						self.format)
	def __iter__(self):
		return self.__dict__.__iter__()
	
	def keys(self):
		return self.__dict__.keys()

	def has_keys(self, k):
		return self.__dict__.has_key(k)

	def iterkeys(self):
		return self.__dict__.iterkeys()

	def itervalues(self):
		return self.__dict__.itervalues()

	def iteritems(self):
		return self.__dict__.iteritems()

	def items(self):
		return self.__dict__.items()
	
	def __getstate__(self):
		return self.__dict__.copy()

	def __setstate__(self, dic):
		for (name, value) in dic.iteritems():
			setattr(self, name, value)

#==============================================================================
class TableColumn(numpy.ndarray): 
	""" Manage one column in a table """
#==============================================================================
	def __new__(self, input_array, *args, **kwargs):
		# Input array is an already formed ndarray instance
		# We first cast to be our class type
		obj = numpy.array(input_array).view(self)
		# add the new attribute to the created instance
		self.header = None
		if isinstance(self, TableColumn):
			if obj.header != None:
				obj.header = self.header.copy()
		if isinstance(input_array, TableColumn):
			if obj.header != None:
				obj.header = input_array.header.copy()
		# Finally, we must return the newly created object:
		return obj

	def __init__(self, data, name='None', dtype=None, unit=None, description=None, null=None,
			format=None, *args, **kwargs):

		if self.header == None:
			self.header = TableColumnHeader(name, self.dtype, unit, description,
								null, format)
		if isinstance(data, TableColumn):
			if data.header != None:
				self.header = data.header.copy()
	
	def __repr__(self):
		s = ""	
		try:
			s = "Column: %s, " % self.header.name
			if self.header.unit != None:
				s+= '\tunit=%s, ' % self.header['unit']
			if self.header.description != None:
				s+= '\t(%s), ' %self.header['description']
		except AttributeError:
			pass
		s+= numpy.ndarray.__repr__(self)
		return s
	
	def copy(self):
		col = TableColumn(numpy.ndarray.copy(self))
		col.header = self.header.copy()
		return col

	
	def __getslice__(self, i,j):
		r = numpy.ndarray.__getslice__(self,i,j)
		if numpy.size(r) > 1:
			return numpy.asarray(r)
		else:
			return r
		"""
		col =  TableColumn(numpy.ndarray.__getslice__(self,i,j))
		if self.header != None:
			col.header = self.header.copy()
		if col.header['description'] != None:
			col.header['description'] += ' [slice]'
		else:
			col.header['description'] = '[slice]'
	
		return col
		"""
	def __getitem__(self, i):
		r = numpy.ndarray.__getitem__(self,i)
		if numpy.size(r) > 1:
			return numpy.asarray(r)
		else:
			return r
		"""
		col =  TableColumn(numpy.ndarray.__getitem__(self,i))
		if self.header != None:
			col.header = self.header.copy()
		if col.header['description'] != None:
			col.header['description'] += ' [slice]'
		else:
			col.header['description'] = '[slice]'
	
		return col
		"""
	
	def __le__(self, y):
		return numpy.asarray(self).__le__(numpy.asarray(y))

	def __lt__(self, y):
		return numpy.asarray(self).__lt__(numpy.asarray(y))

	def __ge__(self, y):
		return numpy.asarray(self).__ge__(numpy.asarray(y))

	def __gt__(self, y):
		return numpy.asarray(self).__gt__(numpy.asarray(y))

	def __eq__(self, y):
		return numpy.asarray(self).__eq__(numpy.asarray(y))

	def __ne__(self, y):
		return numpy.asarray(self).__ne__(numpy.asarray(y))


	def __reduce__(self):
		object_state = list(numpy.ndarray.__reduce__(self))
		subclass_state = (self.header,)
		object_state[2] = (object_state[2],subclass_state)
		return tuple(object_state)

	def __setstate__(self,state):
		nd_state, own_state = state
		numpy.ndarray.__setstate__(self,nd_state)

		header, = own_state
		self.header = header
	"""
	def __getstate__(self):
		return self.__dict__.copy()

	def __setstate__(self, dic):
		print dic
		for (name, value) in dic.iteritems():
			setattr(self, name, value)
	"""
#==============================================================================
class TableHeader(object): 
	""" this class defines the context of a Table """
#==============================================================================
	def __init__(self, dic=None, *args, **kwargs):
		""" constructor """
		for k, v in kwargs:
			self.__setattr__(k,v)
		if (isinstance(dic, dict)) or (isinstance(dic, TableHeader)):
			for k, v in dic.iteritems():
				self.__setattr__(k, v)
		if not 'NAME' in kwargs:
			self.__setattr__('NAME', 'Noname')

	def __setattr__(self, attribute, value):
		""" set attribute """
		if (attribute.lower() in ['description', 'comment', 'history']) & \
			(attribute in self) :
				self.__dict__[attribute] = str(self[attribute])+'\n'+str(value)
		else:
			self.__dict__[attribute] = str(value)
	
	def __setitem__(self, item, val):
		"""get item"""
		self.__setattr__(item,val)

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
			keys = numpy.sort(self.keys())
			for k in keys:
				if self[k] != None:
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

#==============================================================================
class TableManager(object): 
	""" Template class for managing access to table files
		The aim of this class is to propose various format managers by
		simply providing reading and writing methods

		Define any format manager by deriving this class and defining
		the two followin method:
		
		read  -- reading from the considered format
		write -- export to the considered format

		then you can register this format manager by using the
		register_extension function
	"""
#==============================================================================
	def __init__(self, tableType, readerFunction = None, writerFunction = None):
		""" constructor """
		self.tableType = tableType
		if readerFunction:
			del self.read
			self.read = readerFunction
		if writerFunction:
			del self.write
			self.write = writerFunction
	
	def read(self,*args, **kwargs):
		""" to be overwritten """
		pass

	def write(self, *args, **kwargs):
		""" to be overwritten """
		pass

#==============================================================================
# ASCII/CSV TABLE MANAGERS
#==============================================================================
import numpy

class csvManager(TableManager):
	def __init__(self):
		""" constructor """
		TableManager.__init__(self, tableType='csv')

	def read(self, filename, delimiter=',', noheader=False, skiprows=0,
	comment='#', *args, **kwargs):
		"""
		Read Csv file with header or not. Especially useful in association with
		exportdata module.
		So far it uses also the numpy.loadtxt method
		"""
		data = {}
		stream = open(filename, 'r')
		#get header
		description = TableHeader()
		colInfo = {}
		nHeadLines = 0
		header = None
		while header == None:
			line = stream.readline()[:-1]
			nHeadLines += 1
			if line[0] == comment:
				if line[1] != comment:
					k = line[1:].split('\t')
					key = k[0].split()[0] #remove trailing spaces
					for cv in k[1:]:
						description[key] = cv
				else:
					k = line[2:].split('\t')
					colName = k[0].split()[0]
					colUnit = k[1].split()[0]
					colComm = k[2] 
					if colUnit == 'None': colUnit = None
					if colComm == 'None': colComm = None
					colInfo[colName] = (colUnit,colComm)
			else:
				header = line.split(delimiter)
				if noheader:
					header = ['Col%s' % k \
							for k in range(len(header))]
		if not 'NAME' in description.keys():
			description['NAME'] = filename.split('/')[-1]
		stream.close()
		#get data
		skip = skiprows+(nHeadLines-int(noheader))
		d = numpy.recfromcsv(filename, skiprows=skip-1, *args, **kwargs) 
		dkeys = d.dtype.names

		for k in range(len(header)):
			colName = header[k]
			if colInfo.has_key(colName):
				colUnit, colComm = colInfo[colName]
				if colUnit == 'None': colUnit = None
				if colComm == 'None': colComm = None
			else:
				colUnit, colComm = (None,None)
			if colName in data:
				i=1
				while '%s_%d' % (colName, i) in data:
					i += 1
				colName = '%s_%d' % (colName, i)
			data[colName] = TableColumn(d[dkeys[k]], name=colName, unit=colUnit, description=colComm)
		return data, description

	def writeHeader(self, unit, header, comment='#'):
		""" Write File Header definition into the opened unit """
		keys = numpy.sort(header.keys())
		for key in keys:
			val = header[key]
			for kval in str(val).split('\n'):
				unit.write('%s %s\t%s\n' % (comment, key.upper(), kval) )
	
	def writeColHeader(self, data, unit, comment='#'):
		keys = numpy.sort(data.keys())

		unit.write('%s%s %20s\t%10s\t%s\n' % (comment, comment, 'Column', 'Unit', 'Comment') )
		for key in keys:
				txt = "%20s\t%10s\t%s" % (key, data[key].header.unit, 
							data[key].header.description)
				unit.write('%s%s %s\n' % (comment, comment, txt) )

	def writeColDef(self, data, unit, delimiter = ',', comment = ''):
		""" Write column definition into the opened unit """
		keys = numpy.sort(data.keys())
		txt = comment+keys[0]
		for k in keys[1:]: txt+=delimiter+k
		unit.write(txt+"\n")
	
	def writeData(self, data, unit, delimiter=','):
		""" Write data part into the opened unit """
		keys = numpy.sort(data.keys())
		size = numpy.size(data[data.keys()[0]])
		for ik in range(size):
			txt = str(data[keys[0]][ik])
			for k in keys[1:]: txt+=delimiter+str(data[k][ik])
			unit.write(txt+"\n")

	def write(self, data, header=None, output='exportedData.csv', 
		delimiter=',', comment='#', keep=False, unit = None,
		verbose=False, **kwargs):
		"""
		export data to a comma separated value file

		inputs:
			data -- data dictionnary to export
		
		outputs:
			output -- output file (def: exportedData.dat)

		keywords:
			header    -- header object describing the file
			delimiter -- delimiter to use (def: ',')
			comment   -- comment character for header (def: '#')
			keep      -- keeps unit opened
			unit	  -- uses opened stream, if provided 
		"""

		if unit != None:
			outputFile = unit
		else:
			outputFile = open(output, 'w')
			if header != None: 
				self.writeHeader(outputFile, header, comment='#')
			self.writeColHeader(data, outputFile, comment=comment)
			self.writeColDef(data, outputFile, delimiter=delimiter, 
				  		comment='')

		self.writeData(data, outputFile, delimiter=delimiter)	

		if keep == False:
			outputFile.close()
			if verbose: print "Data exported into %s" % output
		else:
			return outputFile
			
class asciiManager(TableManager):
	def __init__(self):
		""" constructor """
		TableManager.__init__(self, tableType='ascii')

	def read(self, filename, delimiter=None, noheader=False, skiprows=0, comment='#', forceHeadLine=0, *args, **kwargs):
		"""
		Read ascii file with header or not. Especially useful in association with
		exportdata module.
		So far it uses also the numpy.loadtxt method
		"""
		data = {}
		stream = open(filename, 'r')
		#get header
		description = TableHeader()
		colInfo = {}
		nHeadLines = 0
		header = None
		oldline = ''
		while header == None:
			line = stream.readline()[:-1]
			nHeadLines += 1
			if (line[0] == comment) | (nHeadLines == forceHeadLine):
				if line[1] != comment:
					k = line[1:].split(delimiter)
					key = k[0].split()[0] #remove trailing spaces
					for cv in k[1:]:
						description[key] = cv
				else:
					k = line[2:].split(delimiter)
					colName = k[0].split()[0]
					colUnit = k[1].split()[0]
					colComm = k[2] 
					if colUnit == 'None': colUnit = None
					if colComm == 'None': colComm = None
					colInfo[colName] = (colUnit,colComm)
				oldline = line
			else:
				header = oldline[1:].split(delimiter)
				if noheader:
					header = ['Col%s' % k \
							for k in range(len(header))]
		stream.close()
		if not 'NAME' in description.keys():
			description['NAME'] = filename.split('/')[-1]
		if 'Column' in colInfo.keys():
			colInfo.pop('Column')
		#get data
		nHeadLines -= 1
		skip = skiprows+(nHeadLines-int(noheader))
		d = numpy.recfromtxt(filename, skiprows=skip-1, delimiter=delimiter, *args, **kwargs) 
		dkeys = d.dtype.names
		for k in range(len(header)):
			colName = header[k]
			if colInfo.has_key(colName):
				colUnit, colComm = colInfo[colName]
				if colUnit == 'None': colUnit = None
				if colComm == 'None': colComm = None
			else:
				colUnit, colComm = (None,None)
			if colName in data:
				i=1
				while '%s_%d' % (colName, i) in data:
					i += 1
				colName = '%s_%d' % (colName, i)
			data[colName] = TableColumn(d[dkeys[k]], name=colName, unit=colUnit, description=colComm)
		return data, description

	def writeHeader(self, unit, header, comment='#'):
		""" Write File Header definition into the opened unit """
		keys = numpy.sort(header.keys())
		for key in keys:
			val = header[key]
			for kval in str(val).split('\n'):
				unit.write('%s %s\t%s\n' % (comment, key.upper(), kval) )

	def writeColHeader(self, data, unit, comment='#'):
		keys = numpy.sort(data.keys())

		unit.write('%s%s %20s\t%10s\t%s\n' % (comment, comment, 'Column', 'Unit', 'Comment') )
		for key in keys:
				txt = "%20s\t%10s\t%s" % (key, data[key].header.unit, 
							data[key].header.description)
				unit.write('%s%s %s\n' % (comment, comment, txt) )


	def writeColDef(self, data, unit, delimiter = None, comment = '#'):
		""" Write column definition into the opened unit """
		keys = numpy.sort(data.keys())
		txt = comment+keys[0]
		for k in keys[1:]: txt+=delimiter+k
		unit.write(txt+"\n")
	
	def writeData(self, data, unit, delimiter=None):
		""" Write data part into the opened unit """
		keys = numpy.sort(data.keys())
		size = numpy.size(data[data.keys()[0]])
		for ik in range(size):
			txt = str(data[keys[0]][ik])
			for k in keys[1:]: txt+=delimiter+str(data[k][ik])
			unit.write(txt+"\n")

	def write(self, data, header=None, output='exportedData.dat', 
		delimiter=None, comment='#', keep=False, unit = None,
		verbose=False, **kwargs):
		"""
		export data to a comma separated value file

		inputs:
			data -- data dictionnary to export
		
		outputs:
			output -- output file (def: exportedData.dat)

		keywords:
			header    -- header object describing the file
			delimiter -- delimiter to use (def: ',')
			comment   -- comment character for header (def: '#')
			keep      -- keeps unit opened
			unit	  -- uses opened stream, if provided 
		"""

		if unit != None:
			outputFile = unit
		else:
			outputFile = open(output, 'w')
			if header != None: 
				self.writeHeader(outputFile, header, comment=comment)
			self.writeColHeader(data, outputFile, comment=comment)
			self.writeColDef(data, outputFile, delimiter=' ', 
				  		comment=comment)

		self.writeData(data, outputFile, delimiter=' ')	

		if keep == False: 
			outputFile.close()
			if verbose: print "Data exported into %s" % output
		else:
			return outputFile


_extensions          = {}
register_extension(csvManager(),   'csv')
register_extension(asciiManager(), 'ascii')
register_extension(asciiManager(), 'dat')
register_extension(asciiManager(), 'txt')

class latexManager(TableManager):
	def __init__(self):
		""" constructor """
		TableManager.__init__(self, tableType='tex')

	def read(self, filename, *args, **kwargs):
		pass

	def writeData(self, data, unit, delimiter='&'):
		""" Write data part into the opened unit """
		keys = numpy.sort(data.keys())
		size = numpy.size(data[data.keys()[0]])
		for ik in range(size):
			txt = str(data[keys[0]][ik])
			for k in keys[1:]: txt+=delimiter+str(data[k][ik])
			unit.write(txt+"\\\\\n")

	def write(self, data, header=None, output='exportedData.tex',
			units=None, comments=None, verbose=False,
			**kwargs):
		"""
		export data to a latex Table

		inputs:
			data -- data dictionnary to export
		
		outputs:
			output -- output file (def: exportedData.dat)

		keywords:
			header    -- header object describing the file
			delimiter -- delimiter to use (def: ',')
			comment   -- comment character for header (def: '#')
			keep      -- keeps unit opened
			unit	  -- uses opened stream, if provided 
		"""
		
		unit = open(output, 'w')
		unit.write('\\begin{table}\n \\begin{center}\n')
		if 'NAME' in header:
			unit.write('\\caption{%s}\n' % header['NAME'])
		if isinstance(data, Table):
			aligntxt = ''.join(['c']*len(data.data))
		else:
			aligntxt = ''.join(['c']*len(data))
		unit.write('\\begin{tabular}{%s}\n' % aligntxt)
		names = numpy.sort(data.keys())
		colDefTxt = ''	
		_Notes = []
		for k in range(len(names)):
			colDefTxt += ' & %s' % names[k]
			if (comments != None) & (str(comments[k]) != 'None'):
				_Notes.append(comments[k])
				colDefTxt += "$^{\\rm{(%s)}}$" % len(_Notes)
		colDefTxt = colDefTxt[2:]
		unit.write('\\hline\\hline\n')
		unit.write('%s\\\\ \n' %colDefTxt )
		if (units == None) & isinstance(data, Table):
			units = [ str(k.header.unit) for k in data ]
		else:
			units = [ str(k) for k in units ]
		units = numpy.array(units)[numpy.argsort(data.keys())]
		if units != None:
			units = ' & '.join(units).replace('None', '')
			unit.write('%s \\\\\n' % units)
		unit.write('\\hline\n')

		self.writeData(data, unit, delimiter=' & ')	

		unit.write('\\hline\n')
		if len(_Notes) > 0:
			unit.write('\\begin{scriptsize}\n')
			for k in range(len(_Notes)):
				unit.write("$^{\\rm(%d)}$ %s\\\\\n" % (k+1, _Notes[k]))
			unit.write('\\end{scriptsize}\n')
		unit.write('\\end{tabular}\n\\end{center}\n\\end{table}\n')
		unit.close()

		if verbose: print "Data exported into %s" % output

register_extension(latexManager(), 'tex')

try:
	import pyfits

	class fitsManager(TableManager):

		def __init__(self):
			""" constructor """
			TableManager.__init__(self, tableType='fits')

		def _getFitsFmt(self, val):
			"""
			-- Internal use -- 
			return the format string to use while defining
			the fits table column.
			input:
				val -- column of values to define
			outputs:
				type -- the character definition associated
			"""
			t = type(val[0])
			if (t == numpy.int) | (t == numpy.int8)| (t == numpy.int16)| (t == numpy.int32):
				return('D')
			elif (t == numpy.float) | (t == numpy.float32)|(t == numpy.float64):
				return('F')
			elif (t == numpy.str) | (isinstance(val[0], str)):
				_len = numpy.max([len(k) for k in val])
				return(str(_len)+'A')
			elif t == numpy.short:
				return('I')
			elif (t == numpy.int32) | (t == numpy.int64):
				return('F')
			else:
				print "WARNING: type conversion not found", t

		def readHeader(self, hdu):
			header = TableHeader()
			genTerms = [ 'XTENSION', 'BITPIX', 'NAXIS', 'NAXIS1',
				     'NAXIS2', 'PCOUNT', 'GCOUNT', 'TFIELDS'  ]
			fieldTerms = ['TTYPE', 'TFORM', 'TUNIT']
			for k in hdu.header:
				if (not k in genTerms) & (not k[:5] in fieldTerms):
					header[k] = hdu.header[k]
			if 'EXTNAME' in hdu.header:
				header['NAME'] =  hdu.header['EXTNAME']
			return header

		def _readColComments(self, card):
			ttype =  str(card).split('/')
			if len(ttype) > 1:
				return ' '.join(ttype[-1].split())
			else:
				return None
			
		def readData(self, hdu):
			colDef = hdu.columns
			names = [ k.name for k in colDef ]
			units = [ k.unit for k in colDef ]
			#comms = [ k.comment for k in colDef ]
			comms = [ self._readColComments(k) for k in hdu.header.ascard['TTYPE*'] ]
			data  = [ TableColumn(hdu.data.field(names[k]), \
						name=names[k], \
						unit=units[k], \
						description=comms[k] )\
						for k in range(len(names)) ]
			return data	

		def read(self, filename, extension = 1, **kwargs):
			hdu = pyfits.open(filename)
			header = self.readHeader(hdu[extension])	
			data   = self.readData(hdu[extension])
			hdu.close()
			return data, header

		def writeColComment(self, header, colName, comment):
			cards = header.ascard['TTYPE*']
			refs = {}
			for k in cards: refs[k.value] = k.key
			header.update(refs[colName], colName, comment=comment) 

		def write(self, data, header=None, output='exportedData.fits', fmt=None,
			name=None, comments=None, units=None, clobber=False, append=True, global_attrs=None,
			silent=False, hdr_only = False, keep_open = False, hdr0 = None):
			"""
			export data to a FITS file

			inputs:
				data -- data dictionnary to export
			
			outputs:
				output -- output file (def: exportedData.dat)

			keywords:
				fmt     -- force a given column format 
				             (alphabetic data name order)
				name    -- extention name of the fits file (def: DATA)
				header  -- dictonnary of keywords and corresponding 
					   values
				comments-- list of column comments
				units   -- list of column units
				global_attrs   -- dictonnary of keywords to add to the
						  main container. (or update; No need
						  to include other keywords to do so.)
				clobber -- overwrite if set (def: False)
				append  -- add data to existing file (def: True)
			"""
			if (not os.path.isfile(output)) & append: 
				if not silent: 
					print "Warning: %s does not seem to exist" % output
					print "         Creating a new file."
				append=False	

			if data != None:
				keys = data.keys()
				indkeys = numpy.argsort(data.keys())
				if fmt == None:
					fmt = [ self._getFitsFmt(data[k]) for k in keys ]
				if (comments == None) | len(comments) != len(data) :
					comments = [None]*len(data)	
				if (units == None) | len(units) != len(data) :
					units = [None]*len(data)	
				cols = [ pyfits.Column( name=keys[ik],        \
							array=data[keys[ik]], \
							format=fmt[ik],       \
							unit=units[ik]      ) \
							for ik in indkeys  ]
				hdr = pyfits.new_table(cols)
				hdr.header.update('EXTNAME', name or header['NAME'])
				hdr.header.update('FILENAME', output.split('/')[-1])
				if header != None:
					for k in header:
						if (k != 'COMMENT') & (k != 'HISTORY'):
							hdr.header.update(k, header[k])	
						else:
							txt = header[k].split('\n')
							for j in txt:
								if k == 'COMMENT':
									hdr.header.add_comment(j)	
								elif k == 'HISTORY':
									hdr.header.add_history(j)	
				for ik in range(len(keys)):
					if comments[ik] != None:
						self.writeColComment(hdr.header, keys[ik], comments[ik])
				if hdr_only:
					return hdr

				if not append:
					hdr.writeto(output,clobber=clobber)
					if not silent: 
						print "Data exported into %s" % output
				else:
					if hdr0 == None:
						retHdr = True
						hdr0 = pyfits.open(output, mode='append')
					else:
						retHdr = False

					hdr0.append(hdr)
					if not keep_open:
						hdr0.flush()
						hdr0.close()
					else:
						if retHdr: return hdr0
					if not silent: 
						print "Data added into %s" % output

			if global_attrs != None: 
				if hdr0 == None:
					hdr0 = pyfits.open(output, mode='update')
				for k in global_attrs:
					hdr0[0].header.update(k, global_attrs[k])	
				hdr0.flush()
				hdr0.close()
				if not silent: 
					print "Keywords added to main table into %s" % output

	register_extension(fitsManager(), 'fits')

except ImportError:
	print "Warning: FITS files could not be managed."

try:
	import tables
	class hd5Manager(TableManager):

			def __init__(self):
				""" constructor """
				TableManager.__init__(self, tableType='hd5')


			def _getValClass(self, val):
				try:
					t = type(val[0])
				except:
					t = type(val)
				if (t == numpy.str) | (t==numpy.string_):
					return('StringCol')
				else:
					t = type(numpy.max(val))
					if t == numpy.int:
						return('IntCol')
					elif (t == numpy.float):
						return('FloatCol')
					elif (t == numpy.bool):
						return('BoolCol')
					elif (t == numpy.complex):
						return('ComplexCol')
					elif (t == numpy.float32):
						return('Float32Col')
					elif (t == numpy.float64):
						return('Float64Col')
					elif (t == numpy.short) | (t == numpy.int8):
						return('Int8Col')
					elif (t == numpy.int16):
						return('Int16Col')
					elif (t == numpy.int32):
						return('Int32Col')
					elif (t == numpy.int64):
						return('Int64Col')
					else:
						print "Warning: type unknown!", t

			def _getExtendedFmt(self, val, pos=0):
				"""
				return the format string to use while defining
				the pyTable description class.
				In this context, this has been extended to account for
				multidimentional values.
				input:
					val -- column of values to define
					pos -- position of the column (def: 0)
				outputs:
					type -- the string declaration
				"""
				if (not _isiterable(val[0])) :
					return self._getHD5Fmt(val,pos=pos)
				else:
					_val = val[0]
					s = numpy.shape(_val)
					_class = self._getValClass(_val)
					_args = ""
					if _class == 'StringCol':
						_len = numpy.max([len(k) for k in _val])
						#_args = str(_len)
					_args += 'shape='+str(s)+',pos='+str(pos)
					return(_class+'('+_args+')')

			def _getHD5Fmt(self, val, pos=0):
				"""
				-- Internal use -- 
				return the format string to use while defining
				the pyTable description class.
				input:
					val -- column of values to define
					pos -- position of the column (def: 0)
				outputs:
					type -- the string declaration
				"""
				try:
					t = type(val[0])
				except:
					t = type(val)
				if (t == numpy.str) | (t==numpy.string_):
					_len = numpy.max([len(k) for k in val])
					return('StringCol('+str(_len)+",pos="+str(pos)+')')
				else:
					t = type(numpy.max(val[:]))
					if t == numpy.int:
						return('IntCol(pos='+str(pos)+')')
					elif (t == numpy.float):
						return('FloatCol(pos='+str(pos)+')')
					elif (t == numpy.bool):
						return('BoolCol(pos='+str(pos)+')')
					elif (t == numpy.complex):
						return('ComplexCol(pos='+str(pos)+')')
					elif (t == numpy.float32):
						return('Float32Col(pos='+str(pos)+')')
					elif (t == numpy.float64):
						return('Float64Col(pos='+str(pos)+')')
					elif (t == numpy.short) | (t == numpy.int8):
						return('Int8Col(pos='+str(pos)+')')
					elif (t == numpy.int16):
						return('Int16Col(pos='+str(pos)+')')
					elif (t == numpy.int32):
						return('Int32Col(pos='+str(pos)+')')
					elif (t == numpy.int64):
						return('Int64Col(pos='+str(pos)+')')
					else:
						print "Warning: type unknown!", t
			
			def readColDesc(self, tab, name):
				key = [ k for k in tab.attrs._v_attrnames \
						if tab.attrs[k] == name ]

				key = key[0].replace('NAME', '')
				if key+'UNIT' in tab.attrs: 
					unit = tab.attrs[key+'UNIT']
				else: 
					unit = None
				if key+'DESC' in tab.attrs: 
					desc = tab.attrs[key+'DESC']
				else:
					desc = None
				return (unit, desc)

			def readCol(self, tab, colName):
					cunit, cdesc = self.readColDesc(tab, colName)	
					return TableColumn(tab.col(colName), 
								name = colName, 
								unit=cunit,
								description=cdesc)
			def readTabHeader(self, tab):
				head = TableHeader()
				exclude = ['NROWS', 'VERSION', 'CLASS', 'EXTNAME']
				for k in tab.attrs._v_attrnames:
					if (not k in exclude) & (k[:5] != 'FIELD'):
						head[k] = tab.attrs[k]
				return head

			def read(self, filename, tableName=None, silent=False, *args, **kwargs):
				source = tables.openFile(filename, *args, **kwargs)
				if 'tablename' in kwargs:
					tableName = kwargs['tablename']
				if tableName == None:
					node = source.listNodes('/')[0]
					tableName = node.name
				else:
					if tableName[0] != '/': tableName = '/'+tableName
					node = source.getNode(tableName)
				if silent != False:
					print "\tLoading table: %s" % tableName
				data = [self.readCol(node, k) for k in node.colnames]
				head = self.readTabHeader(node)
				head['_sourceFile'] = source
				if 'NAME' not in head or head['NAME']=='Noname' or head['NAME'] == None: 
					head['NAME'] = tableName
				if 'TITLE' not in head or head['TITLE']=='': 
					head['TITLE'] = node.title
				return data, head
			
			def writeColDesc(self, tab, name, unit=None, desc=None):
				key = [ k for k in tab.attrs._v_attrnames \
						if tab.attrs[k] == name ]

				key = key[0].replace('NAME', '')
				tab.attrs[key+'UNIT'] = unit
				tab.attrs[key+'DESC'] = desc
				
				
			def write(self, data, header=None,
					output='exportedData.hd5', tablename='data', 
					mode='w', group='/', silent=False,
					units=None, comments=None,
					appendTable=False, **kwargs):
				"""
				export data to a HDF5 file

				inputs:
					data -- data dictionnary to export
				
				outputs:
					output -- output file (def: exportedData.dat)

				keywords:
					tablename -- table name to create/modify (def: 'data')
					mode      -- 'w' will create a brand new file (default)
						     'a' will append an existing file (useful to add
						     tables in an existing file)
					group     -- path to the table (def: '/') 
					header    -- Dictionnary of attributes to add to the table
					silent    -- Do not print any message when set
				"""
				if header != None:
					if 'NAME' in header:
						tablename = header['NAME'].replace('.','_')
				if (comments == None) | len(comments) != len(data) :
					comments = [None]*len(data)	
				if (units == None) | len(units) != len(data) :
					units = [None]*len(data)	

				keys = data.keys()
				keys.sort()
				if appendTable == True:
					mode = 'a'
				if mode == 'w':
					hd5 = tables.openFile(output,mode=mode)
				else:
					hd5 = tables.openFile(output, mode=mode)

				#Gen the Table Colum def class
				if not appendTable:
					code = 'class newTable(tables.IsDescription):'
					pos = 0
					for k in keys:
						ktype = self._getExtendedFmt(data[k], pos)
						code +="\n\t"+ k + "= tables." + ktype 
						pos += 1
					# generate a table data class associated to data
					exec(code)	
					try:
						if group[-1] == '/':
							group = group[:-1]
						table = hd5.createTable(group, tablename, newTable, expectedrows = numpy.size(data[keys[0]]), createparents=True)
						hd5.flush()
					except:
						if not silent:
							print "Warning: Table creation exception. Table may already exist."
						table = hd5.getNode(group+tablename)
				else:
					try:
						table = hd5.getNode(group+tablename)
					except tables.NoSuchNodeError:
						code = 'class newTable(tables.IsDescription):'
						pos = 0
						for k in keys:
							ktype = self._getExtendedFmt(data[k], pos)
							code +="\n\t"+ k + "= tables." + ktype 
							pos += 1
						# generate a table data class associated to data
						exec(code)	
						if group[-1] == '/':
							group = group[:-1]
						table = hd5.createTable(group, tablename, newTable, expectedrows = numpy.size(data[keys[0]]), createparents=True)
						hd5.flush()
					
				row = table.row
				for i in range(numpy.size(data[keys[0]])):
					for k in keys:
						row[k] = data[k][i]
					row.append()
				table.flush()
				if header != None:
					for k in header:
						table.attrs[k] = header[k]
					if not 'TITLE' in header:
						table.attrs['TITLE'] = tablename
				
				for ik in range(len(keys)):
					self.writeColDesc(table, keys[ik], units[ik], comments[ik])
		
				hd5.close()
				if not silent: print "Data exported into %s" % output

	register_extension(hd5Manager(), 'hd5')
	register_extension(hd5Manager(), 'hdf5')

	
except:
	print "Warning: HD5 files could not be managed."


def load(filename, type=None, manager=None, distributed=False, storage=None, silent=False, **kwargs):
	""" Generates a Table object from a given file
		Basically it creates an empty Table object and call the
		Table.read() method to fill the table
		
	inputs:
		filename  -- [ string ]
			     file to read from

	keywords:
		type     -- [ string ]
			    if specified, this will force the function
			    to use considered the file to be of this
			    type.
		manager  -- [ TableManager ]
			    If specified, it will use this format
			    manager (even if not registered)
	     distributed -- generate a distributedTable object (BETA)
		**kwargs are sent to the TableManager.read function
	"""
	if distributed:
		t = distributedTable(storage)
	else:
		t = Table()
	t.read(filename, type=type, manager=manager, silent=silent, **kwargs)
	return t

try:
	from scipy import io
	def idlload(filename, distributed=False, storage=None, **kwargs):
		""" A shortcut to open an IDL save file into a list of recarrays 
			
		inputs:
			filename  -- [ string ]
				     file to read from

		keywords:
		
			**kwargs are sent to the scipy.io.idl.readsav function
		"""
		return io.idl.readsav(filename, **kwargs)
except:
	print "Warning: IDL savefiles could not be managed. (req: scipy.io.idl)"



try:
	import sqlite3
	
	class sqlManager(TableManager):

		def __init__(self):
			""" constructor """
			TableManager.__init__(self, tableType='sql')

		def getTableNames(self, *args, **kwargs):
			node = self.c.execute("SELECT name FROM sqlite_master WHERE type='table';")
			return [ t[0] for t in node if t[0][0] != '_' ]
		
		def execute(self, txt, *args, **kwargs):
			return self.c.execute(txt)

		def readCol(self, k):
			r = self.execute('select %s from %s;' % (k, self.tableName)).fetchall()
			return TableColumn(numpy.asarray([ rk[0] for rk in r]), name=k)
		
		def readTabHeader(self):
			""" assumes dr is a table with keyname, value in
			'_tableName_' """
			hdr = TableHeader()
			try:
				r = self.c.execute('select * from %s' % '_'+self.tableName+'_')
				names = [rk[0] for rk in r.description]
				vals  = [rv for rv in r.fetchall()[0]]
				for k in range(len(names)):
					hdr[names[k]] = vals[k]
				del r, names, vals
			except:
				pass
			return hdr
			
		def read(self, filename, tableName=None, silent=False, *args, **kwargs):
			self.source = sqlite3.connect(filename, *args, **kwargs)
			self.c = self.source.cursor()
			#TODO self.cache = False
			if 'tablename' in kwargs:
				tableName = kwargs['tablename']
			if tableName == None:
				tableName = self.getTableNames()[0]
			if silent != False:
				print "\tLoading table: %s" % tableName
			self.tableName = tableName

			r = self.c.execute('select * from %s' % tableName)
			colnames = [ t[0] for t in r.description ]
			data = [self.readCol(k) for k in colnames]

			head = self.readTabHeader()
			head['_sourceFile'] = self.source

			if 'NAME' not in head or head['NAME']=='Noname' or head['NAME'] == None:
				head['NAME'] = tableName
			return data, head

		def _getSQLFmt(self, val):
			"""
			-- Internal use -- 
			return the format string to use while defining
			the Table description.
			input:
				val -- column of values to define
			outputs:
				type -- the string declaration
			"""
			try:
				t = type(val[0])
			except:
				t = type(val)
			if (t == numpy.str) | (t==numpy.string_) | (t==numpy.unicode_) | (t==unicode):
				_len = numpy.max([len(k) for k in val])
				return('TEXT')
			else:
				t = type(numpy.max(val[:]))
				if t == numpy.int:
					return('INTEGER')
				elif (t == numpy.float):
					return('REAL')
				elif (t == numpy.float32):
					return('REAL')
				elif (t == numpy.float64):
					return('REAL')
				elif (t == numpy.short) | (t == numpy.int8) | (t == numpy.uint8):
					return('INTEGER')
				elif (t == numpy.int16):
					return('INTEGER')
				elif (t == numpy.int32):
					return('INTEGER')
				elif (t == numpy.int64):
					return('INTEGER')
				elif (t == None):
					return('NULL')
				else:
					print "Warning: type unknown!", t

		def write(self, data, header=None, output='exportedData.sqlite', 
					tablename=None, silent=False, units=None, 
					comments=None, append=False, **kwargs):
			
			if header != None:
				if 'NAME' in header:
					if tablename == None:
						tablename = header['NAME'].replace('.','_')
			if (comments == None) | len(comments) != len(data) :
				comments = [None]*len(data)	
			if (units == None) | len(units) != len(data) :
				units = [None]*len(data)	

			keys = data.keys()
			keys.sort()


			conn = sqlite3.connect(output)
			c = conn.cursor()
			if not append:
				txt = 'create table %s \n (' % tablename
				for k in keys:
					txt += '%s %s, ' % (k, self._getSQLFmt(data[k]) )
				txt = txt[:-2]+')'
				print txt
				c.execute(txt)
				conn.commit()
			nrows = len(data[keys[0]])
			rki = 1
			for rk in range(nrows):
				sys.stdout.write('Insert operation: %d/%d (%d%%)\r' %(rki, nrows,
							int(100.*float(rki)/float(nrows))))
				row = [ str(data[tk][rk]) != 'nan' and str(data[tk][rk]) or 'NULL' for tk in keys ]
				c.execute("INSERT INTO %s (%s) VALUES (%s);" % (tablename, ','.join(keys), ','.join(row)))
				rki += 1
			conn.commit()
			
			if header != None:
				keys = header.keys()
				keys.sort()
				if not append:
					txt = 'create table _'+tablename+'_ ('
					for k in keys:
						txt += '%s %s, ' % (k, self._getSQLFmt([header[k]]) )
					txt = txt[:-2]+')'
					c.execute(txt)
					conn.commit()
					txt = 'insert into _'+ tablename + '_ values ('+','.join(['?']*len(keys))+')'
					t = tuple([header[k] for k in keys])
					c.execute(txt, 	t )
				conn.commit()



			
	register_extension(sqlManager(),   'sql')
	register_extension(sqlManager(),   'sqlite')
	register_extension(sqlManager(),   'db')
except:
	print "Warning sqlite tables are not managed. (req: sqlite3)"

#==============================================================================
# Table manipulations
#==============================================================================

def match(c1,c2):
	""" Returns the indices at which the tables match 
	matching uses 2 columns (arrays) that are compared in values
	INPUTS:
		c1 -- array 1
		c2 -- array 2

	OUTPUS:
		tuple of both indices list where the two columns match.
	"""
	return numpy.where(numpy.equal.outer(c1,c2))

def join(t1,t2, n1, n2):
	""" returns a table containing matching lines 
	This will create a table containing both tables columns with matching
	lines only. (units and descriptions are preserved throught the join but
	not the headers)
	
	INPUTS:
		t1 -- Table 1 (Table Object)
		t2 -- Table 2 (Table Object)
		n1 -- Name of Column 1 (string)
		n2 -- Name of Column 2 (string)
	
	OUTPUTS:
		Joined Table 
	"""
	ind1, ind2 = match(t1[n1], t2[n2])
	tab = Table(name='Joined Table')
	# get matching lines
	r1 = t1[ind1]
	r2 = t2[ind2]
	# get column names
	k1 = t1.keys()
	k2 = t2.keys()

	# add columns
	for k in range(len(k1)):
		if t2.has_key(k1[k]):
			name = k1[k]+'_1'
		else:
			name = k1[k]
		tab.addCol(r1[k1[k]], name=name, 
				unit=t1[k1[k]].header['unit'],
				description=t1[k1[k]].header['description'] )	
	for k in range(len(k2)):
		if t1.has_key(k2[k]):
			name = k2[k]+'_2'
		else:
			name = k2[k]
		tab.addCol(r2[k2[k]], name=name,
				unit=t2[k2[k]].header['unit'],
				description=t2[k2[k]].header['description'])	
	#include header trace
	if n1 in k2:
		txt1 = n1+'_1'
	else:
		txt1 = n1
	if n2 in k1:
		txt2 = n2+'_2'
	else:
		txt2 = n2
	tab.header['COMMENT'] = 'Joined table using %s and %s' % (txt1, txt2)

	return tab


def indent(rows, hasHeader=False, hasUnits=False, headerChar='-', delim=' | ', justify='left',
           separateRows=False, prefix='', postfix='', wrapfunc=lambda x:x):
    """Indents a table by column.
       - rows: A sequence of sequences of items, one sequence per row.
       - hasHeader: True if the first row consists of the columns' names.
       - headerChar: Character to be used for the row separator line
         (if hasHeader==True or separateRows==True).
       - delim: The column delimiter.
       - justify: Determines how are data justified in their column. 
         Valid values are 'left','right' and 'center'.
       - separateRows: True if rows are to be separated by a line
         of 'headerChar's.
       - prefix: A string prepended to each printed row.
       - postfix: A string appended to each printed row.
       - wrapfunc: A function f(text) for wrapping text; each element in
         the table is first wrapped by this function."""
    # closure for breaking logical rows to physical, using wrapfunc
    def rowWrapper(row):
        newRows = [wrapfunc(item).split('\n') for item in row]
        return [[substr or '' for substr in item] for item in map(None,*newRows)]
    # break each logical row into one or more physical ones
    logicalRows = [rowWrapper(row) for row in rows]
    # columns of physical rows
    columns = map(None,*reduce(operator.add,logicalRows))
    # get the maximum of each column by the string length of its items
    maxWidths = [max([len(str(item)) for item in column]) for column in columns]
    rowSeparator = headerChar * (len(prefix) + len(postfix) + sum(maxWidths) + \
                                 len(delim)*(len(maxWidths)-1))
    # select the appropriate justify method
    justify = {'center':str.center, 'right':str.rjust, 'left':str.ljust}[justify.lower()]
    output=cStringIO.StringIO()
    if separateRows: print >> output, rowSeparator
    for physicalRows in logicalRows:
        for row in physicalRows:
            print >> output, \
                prefix \
                + delim.join([justify(str(item),width) for (item,width) in zip(row,maxWidths)]) \
                + postfix
        if separateRows:
		print >> output, rowSeparator
	elif hasHeader & (not hasUnits):
		print >> output, rowSeparator
	elif (not hasHeader) & hasUnits:
		print >> output, rowSeparator
		hasUnits=False
	hasHeader=False 

    return output.getvalue()


