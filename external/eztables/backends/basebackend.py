""" Defines the template of a format handling backend
Mainly a generic template to read and write functions.

The registery of backends is located in backends/__init__.py
"""

#==============================================================================
class BaseBackend(object): 
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


