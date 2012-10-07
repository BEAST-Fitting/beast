""" Manage Various SED/spectral grids is a generic way """

import stellib, photometry, isochrone
import numpy
import pyfits
import mytables
from decorators import timeit
from itertools import chain

def isNestedInstance(obj, cl):
	tree = []
	for k in cl.__subclasses__():
		tree+=k.__subclasses__()
	tree += cl.__subclasses__() + [ cl ]
	return  issubclass(obj.__class__, tuple(tree))

class ModelGrid(object):
	""" Generic class for a minimum update of future codes """
	def __init__(self, *args, **kwargs):
		self.lamb = None
		self.seds = None
		self.grid = None

	def keys(self):
		""" returns the grid dimension names """
		return []

	def getGridPoints(self, *args, **kwargs):
		""" Returns age, mass, logg, logT, logL... """
		pass

	def getPDF(self, Qname, lnp, *args, **kwargs):
		assert (Qname in self.keys() ), "Cannot find %s in the grid description" % Qname

	def write(self, fname, *args, **kwargs):
		if ((self.lamb != None) & (self.seds != None) & (self.grid != None)):
			assert (isinstance(self.grid, mytables.Table)), 'Only mytables.Table are supported so far'
			r = numpy.vstack( [ numpy.copy(self.seds), self.lamb ])
			pyfits.writeto(fname, r, **kwargs)
			self.grid.write(fname, append=True)	

class MemoryGrid(ModelGrid):
	""" Instanciate an grid object that has no physical storage
		Helps to create new grids on the fly. Because it deriveds from
		ModelGrid, this can be exported on disk too.
	"""
	def __init__(self, lamb, seds=None, grid=None):
		""" MemoryGrid constructor
		INPUTS:
			lamb	ModelGrid or subclass	New ref to the given grid (debug purpose)
				lamb			wavelengths
				sed			seds
				grid			grid associated to seds
		"""
		if isNestedInstance(lamb, ModelGrid):
			self.lamb = lamb.lamb
			self.seds = lamb.seds
			self.grid = lamb.grid
		else:
			assert ((seds != None) & (grid != None)), 'Wrong number of arguments'
			self.lamb = lamb
			self.seds = seds
			self.grid = grid

class SpectralGrid(ModelGrid):

	def getSEDs(self, filter_names):
		"""
		Extract integrated fluxes through filters
		INPUTS:
			filter_names	list	list of filter names according to the filter lib
		"""
		lamb = self.lamb
		specs = self.seds

		flist = photometry.load_filters(filter_names)

		r = numpy.empty( (len(specs), len(flist) ), dtype=float)
		lf = numpy.empty( len(flist), dtype=float )


		for kf in range(len(flist)):
			for ks in range(len(specs)):
				r[ks,kf] = flist[kf].getFlux(lamb, specs[ks])
				lf[kf] = flist[kf].cl

		return MemoryGrid(lf, r, self.grid)

class FileSEDGrid(SpectralGrid):
	""" Generate a grid from a spectral library """

	def __init__(self, fname, *args, **kwargs):
		with pyfits.open(fname) as f:
			self.seds = f[0].data[:-1]
			self.lamb = f[0].data[-1] 
		self.grid = mytables.load(fname)
		#lamb, seds = self.getSEDs(filters, self.osl.wavelength, self.osl.spectra)
		for k in self.grid.keys():
			self.__dict__[k] = self.grid[k]

	def keys(self):
		""" returns the grid dimension names """
		return self.grid.keys()

class FileSpectralGrid(SpectralGrid):
	""" Generate a grid from a spectral library """

	def __init__(self, fname, *args, **kwargs):
		with pyfits.open(fname) as f:
			self.seds = f[0].data[:-1]
			self.lamb = f[0].data[-1] 
		self.grid = mytables.load(fname)
		#lamb, seds = self.getSEDs(filters, self.osl.wavelength, self.osl.spectra)
		for k in self.grid.keys():
			self.__dict__[k] = self.grid[k]

	def keys(self):
		""" returns the grid dimension names """
		return self.grid.keys()

class StellibGrid(SpectralGrid):
	""" Generate a grid from a spectral library """

	def __init__(self, filters, *args, **kwargs):
		self.osl = stellib.BaSeL()
		lamb, seds = self.getSEDs(filters, self.osl.wavelength, self.osl.spectra)
		self.lamb = lamb
		self.seds = seds
		self.filters = filters
		self.grid = self.osl.grid
		for k in self.grid.keys():
			self.__dict__[k] = self.grid[k]

	def keys(self):
		""" returns the grid dimension names """
		return self.grid.keys()





def generate_spectral_grid_from_isochrones(outfile, osl, oiso, Z=0.02):
	""" Reinterpolate a given stellar spectral library on to an Isochrone grid
	INPUTS:
		outfile		str		   	fits file to export to
		osl		stellib.stellib		a stellar library
		oiso		isochrone.Isochrone	an isochrone library
		Z		float			metallicity to use
	
	OUTPUTS:
		None

		only write into outfile
	"""
	
	assert(isNestedInstance(osl, stellib.Stellib) )
	assert(isNestedInstance(oiso, isochrone.Isochrone) )
	specs = numpy.empty( (oiso.data.nrows+1, len(osl.wavelength)), dtype=float )
	specs[-1] = osl.wavelength[:]

	progress = 0
	with timeit('interpolation'):
		for k in range(oiso.data.nrows):
			if progress < int(100*(k+1)/oiso.data.nrows):
				progress = int(100*(k+1)/oiso.data.nrows)
				print "progress... %d / 100" % progress
			r = numpy.array( osl.interp(oiso.data['logT'][k], oiso.data['logg'][k], Z, oiso.data['logL'][k]) ).T
			specs[k,:] = osl.genSpectrum(r)
	pyfits.writeto(outfile, specs)

	#copy pars
	data = {}
	for k in oiso.data.keys():
		data[k] = oiso.data[k]
	pars  = mytables.Table(data, name='Reinterpolated stellib grid')
	pars.header['stellib'] = osl.source
	pars.header['isoch'] = oiso.source

	pars.write(outfile, append=True)


