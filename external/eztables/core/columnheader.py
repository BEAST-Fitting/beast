import numpy as np
#================================================================================
class ColumnHeader(object):
#================================================================================
	def __init__(self, dtype, unit=None, description=None, null=None, format=None):
		self.__dict__['dtype'] = dtype
		self.unit = unit
		self.description = description
		self.__dict__['null'] = null
		self.format = format

	def __setattr__(self, attribute, value):
		if attribute in ['unit', 'description', 'format']:
			self.__dict__[attribute] = value
		elif attribute in ['null', 'dtype']:
			raise Exception("Cannot change %s through the columns object" % attribute)
		else:
			raise AttributeError(attribute)

	def __repr__(self):
		s = "type=%s" % str(self.dtype)
		if self.unit:
			s += ", unit=%s" % str(self.unit)
		if self.null:
			s += ", null=%s" % str(self.null)
		if self.description:
			s += ", description=%s" % self.description
		return s

	def __eq__(self, other):
		if self.dtype != other.dtype:
			return False
		if self.unit != other.unit:
			return False
		if self.description != other.description:
			return False
		if self.null != other.null:
			if np.isnan(self.null):
				if not np.isnan(other.null):
					return False
			else:
				return False
		if self.format != other.format:
			return False
		return True

	def __ne__(self, other):
		return not self.__eq__(other)

