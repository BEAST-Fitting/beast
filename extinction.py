"""
# EXTINCTION MANAGEMENT
# Here are defined extinction functions that relies
# on given laws. 
# 
# ==================================================
"""
import os,sys
import numpy

__version__ = '0.0.2'

class ExtinctionLaw(object):
	""" Template class """
	def __init__(self):
		self.name = 'None'
	def function(self, lamb, *arg, **kwargs):
		""" expected to contain a function of lambda that return the
		extinction values
		"""
		pass
	def inFilter(self, s, *args, **kwargs):
		""" Exected to return the extinction value for a given filter
		band or filter color"""
		pass

	def __call__(self, *args, **kwargs):
		pass

class Calzetti(ExtinctionLaw):
	"""
	Calzetti et al.  (2000, ApJ 533, 682) developed a recipe for dereddening
	the spectra of galaxies where massive stars dominate the radiation
	output, valid between 0.12 to 2.2 microns.  
	Extrapolation down to 0.0912 microns
	
	Note that the supplied color excess should be that derived for the
	stellar  continuum, EBV(stars), which is related to the reddening
	derived from the gas, EBV(gas), via the Balmer decrement by 
	EBV(stars) = 0.44*EBV(gas)

	R_V - Ratio of total to selective extinction, default = 4.05.  Calzetti
	et al. (2000) estimate R_V = 4.05 +/- 0.80 from optical-IR observations
	of 4 starbursts.

	"""
	def __init__(self):
		self.name = 'Calzetti'
	
	def function(self, lamb, Av=1, Rv=4.05, Alambda=True):
		"""
		Returns Alambda or tau for a Calzetti law
		
		!!Expecting lamb in microns!!

		"""
		if type(lamb) == float:
			_lamb = numpy.asarray([lamb])
		else:
			_lamb = lamb[:]
	
		x = 1./_lamb
		k = numpy.zeros(numpy.size(x))

		ind = numpy.where( (_lamb >= 0.630 ) & (_lamb <= 2.2) )
		k[ind] = 2.659*(-1.857 + 1.040*x[ind]) + Rv

		ind = numpy.where((_lamb >= 0.0912 ) & (_lamb < 0.630) )
		k[ind] = 2.659*(-2.156 + 1.509*x[ind] - 0.198*x[ind]**2 + 0.011*x[ind]**3 ) + Rv

		if Alambda == True:
			return k
		else:	
			return 10**(0.4*k)
		
		
	

class Cardelli(ExtinctionLaw):
	def __init__(self):
		self.name = 'Cardelli'
	def function(self, lamb, Av=1., Rv=3.1, Alambda = True, debug=False):
		""" 
		Cardelli extinction Law
		Lamb as to be in microns!!!

		input:
		    lamb    <float>    wavelength of the extinction point !! in microns !!
		output:
		    tau        <float> returns tau as in redflux = flux*exp(-tau)
		keywords:
		    Alambda        <bool>  returns +2.5*1./log(10.)*tau
		    Av        <float>    extinction value (def: 1.0)
		    Rv        <float> extinction param. (def: 3.1)
		"""
		if type(lamb) == float:
			_lamb = numpy.asarray([lamb])
		else:
			_lamb = lamb[:]

		#init variables
		x = 1./(_lamb) #wavenumber in um^-1
		a = numpy.zeros(numpy.size(x))
		b = numpy.zeros(numpy.size(x))
		#Infrared (Eq 2a,2b)
		ind = numpy.where ((x >= 0.3) & (x < 1.1))
		a[ind] =  0.574*x[ind]**1.61
		b[ind] = -0.527*x[ind]**1.61
		#Optical & Near IR
		#Eq 3a, 3b
		ind = numpy.where ((x >= 1.1) & (x <= 3.3))
		y = x[ind]-1.82
		a[ind] = 1. + 0.17699*y   - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
		b[ind] =      1.41338*y   + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7
		#UV 
		#Eq 4a, 4b
		ind = numpy.where ((x >= 3.3) & (x <= 8.0))
		a[ind] =  1.752 - 0.316*x[ind] - 0.104/((x[ind]-4.67)**2+0.341)
		b[ind] = -3.090 + 1.825*x[ind] + 1.206/((x[ind]-4.62)**2+0.263)

		ind = numpy.where ((x >= 5.9) & (x <= 8.0))
		Fa     = -0.04473*(x[ind]-5.9)**2 - 0.009779*(x[ind]-5.9)**3
		Fb     =  0.21300*(x[ind]-5.9)**2 + 0.120700*(x[ind]-5.9)**3
		a[ind] = a[ind] + Fa
		b[ind] = b[ind] + Fb
		#Far UV 
		#Eq 5a, 5b
		ind = numpy.where ((x >= 8.0) & (x <= 10.0))
		#Fa = Fb = 0
		a[ind] = -1.073 - 0.628*(x[ind]-8.) + 0.137*((x[ind]-8.)**2) - 0.070*(x[ind]-8.)**3 
		b[ind] = 13.670 + 4.257*(x[ind]-8.) + 0.420*((x[ind]-8.)**2) + 0.374*(x[ind]-8.)**3

		# Case of -values x out of range [0.3,10.0]
		ind = numpy.where ((x > 10.0) | (x < 0.3))
		a[ind] = 0.0
		b[ind] = 0.0

		#Return Extinction vector 
		#Eq 1
		if (Alambda == True):
		    return( 2.5*1./numpy.log(10.)*( a + b/Rv ) * Av)
		else:
		    return( ( a + b/Rv ) * Av)
	
	def inFilter(self, s, cat=None, debug=False, *args, **kwargs):
		""" 
		getCardelliVect: returns the Extinction vector from usual fluxes or
		colors (e.g. U, U-B ...)

		input:
		    s    <string> the flux or color to determine
		output:
		    f    <float>  Extinction value for Av=1
		"""
		if debug:
		    print "debug: getCardelliVect(%s)" % s
		if cat == None:
		    catalog = mytables.load(cfg['filterCatalog'][0], tableName=cfg['filterCatalog'][1])
		else:
		    catalog = cat
		if (type(s) != str) & (getattr(s, '__iter__', False) != False):
		    return numpy.asarray([self.inFilter(sk, catalog, debug=debug, *args, **kwargs) for sk in s])
		else:
		    ss = s.upper().split('-')
		    if len(ss) == 1:
			cwave = catalog.selectWhere("A == ss", condvars={'A':catalog['TABLENAME'], 'ss':ss},fields=["CWAVE"])  
			cwave = float(cwave['CWAVE']) * qt.angstrom
			val =  self.function(cwave, *args, **kwargs)    # has to be in microns
		    else:
			c1 = catalog.selectWhere("A == ss", condvars={'A':catalog['TABLENAME'], 'ss':ss[0]},fields=["CWAVE"])  
			c1 = float(cwave['CWAVE']) * qt.angstrom
			c2 = catalog.selectWhere("A == ss", condvars={'A':catalog['TABLENAME'], 'ss':ss[1]},fields=["CWAVE"])  
			c2 = float(cwave['CWAVE']) * qt.angstrom
			val =  self.function(c1, *args, **kwargs)-self.function(c2,*args, **kwargs)
		    if cat == None:
			del catalog
		    return(val[0])

