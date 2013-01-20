"""
# EXTINCTION MANAGEMENT
# Here are defined extinction functions that relies
# on given laws. 
# 
# ==================================================
"""
import os,sys
import numpy
from scipy import interpolate

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
		""" Expected to return the extinction value for a given filter
		band or filter color"""
		pass

	def __call__(self, *args, **kwargs):
		return self.function(*args, **kwargs)

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
	
	def function(self, lamb, Av=1, Rv=4.05, Alambda=True, **kwargs):
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
	""" Cardelli, Clayton, and Mathis (1989, ApJ, 345, 245)"""
	def __init__(self):
		self.name = 'Cardelli'
	def function(self, lamb, Av=1., Rv=3.1, Alambda = True, debug=False, **kwargs):
		""" 
		Cardelli extinction Law
		Lamb has to be in microns!!!

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
		    return( ( a + b/Rv ) * Av)
		else:
		    return( 2.5*1./numpy.log(10.)*( a + b/Rv ) * Av)
	
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

class Fitzpatrick99(ExtinctionLaw):
	"""
	Fitzpatrick (1999, PASP, 111, 63)
	R(V) dependent extinction curve that explicitly deals with optical/NIR
	extinction being measured from broad/medium band photometry.
	Based on fm_unred.pro from the IDL astronomy library
	"""
	def __init__(self):
		self.name = 'Fitzpatrick99'

	def function(self, lamb, Av=1, Rv=3.1, Alambda=True, **kwargs):
		""" 
		Fitzparick99 extinction Law
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

		c2 = -0.824 + 4.717/Rv
		c1 = 2.030 - 3.007*c2
		c3 = 3.23
		c4 = 0.41
		x0 = 4.596
		gamma = 0.99

		x = 1./_lamb
		k = numpy.zeros(numpy.size(x))

		# compute the UV portion of A(lambda)/E(B-V)
		xcutuv = 10000.0/2700.
		xspluv = 10000.0/numpy.array([2700.,2600.])
		ind = numpy.where(x>=xcutuv)

		if numpy.size(ind) > 0:
			k[ind] = c1 + (c2*x[ind]) + c3*((x[ind])**2)/( ((x[ind])**2 - (x0**2))**2 + (gamma**2)*((x[ind])**2 ))
			yspluv = c1 + (c2*xspluv) + c3*((xspluv)**2)/( ((xspluv)**2 - (x0**2))**2 + (gamma**2)*((xspluv)**2 ))
		
			# FUV portion
			fuvind = numpy.where(x>=5.9)
			k[fuvind] += c4*(0.5392*((x[fuvind]-5.9)**2) + 0.05644*((x[fuvind]-5.9)**3))

			k[ind] += Rv
			yspluv += Rv

		# Optical/NIR portion

		ind = numpy.where(x<xcutuv)
		if numpy.size(ind) > 0:
			xsplopir = numpy.zeros(7)
			xsplopir[0] = 0.0
			xsplopir[1:7] = 10000.0/numpy.array([26500.0,12200.0,6000.0,5470.0,4670.0,4110.0])

			ysplopir = numpy.zeros(7)
			ysplopir[0:3] =  numpy.array([0.0,0.26469,0.82925])*Rv/3.1

			
			ysplopir[3:7] = numpy.array([numpy.poly1d([2.13572e-04 , 1.00270, -4.22809e-01])(Rv),
						     numpy.poly1d([-7.35778e-05, 1.00216, -5.13540e-02])(Rv),
						     numpy.poly1d([-3.32598e-05, 1.00184,  7.00127e-01])(Rv),
						     numpy.poly1d([ 1.19456, 1.01707, -5.46959e-03, 7.97809e-04,
								    -4.45636e-05][::-1])(Rv)])

			tck = interpolate.splrep(numpy.hstack([xsplopir,xspluv]), numpy.hstack([ysplopir,yspluv]), k=3)
			k[ind] = interpolate.splev(x[ind], tck)

		# convert from A(lambda)/E(B-V) to A(lambda)/A(V)
		k /= Rv
			
		if (Alambda == True):
			return(k*Av)
	        else:
			return(k*Av*(numpy.log(10.)/2.5))

class Gordon03_SMCBar(ExtinctionLaw):
	""" Gordon et al. 2003 (ApJ, 594:279-293)"""
	def __init__(self):
		self.name = 'Gordon03_SMCBar'
	def function(self, lamb, Av=1, Alambda=True, **kwargs):
		"""
		Lamb as to be in microns!!!

		input:
		    lamb    <float>    wavelength of the extinction point !! in microns !!
		output:
		    tau        <float> returns tau as in redflux = flux*exp(-tau)
		keywords:
		    Alambda        <bool>  returns +2.5*1./log(10.)*tau
		    Av        <float>    extinction value (def: 1.0)
		"""
		if type(lamb) == float:
			_lamb = numpy.asarray([lamb])
		else:
			_lamb = lamb[:]

		Rv = 2.74
		
#		c1= -3.618/Rv
#		c2 = 1.994/Rv
#		c3 = 0.239/Rv
#		c4 = 0.618/Rv
		c1= -4.959/Rv
		c2 = 2.264/Rv
		c3 = 0.389/Rv
		c4 = 0.461/Rv
		x0 = 4.6
		gamma = 1.0
		
		x = 1./_lamb
		k = numpy.zeros(numpy.size(x))
		
		# UV part
		xcutuv = 10000.0/2700.
		xspluv = 10000.0/numpy.array([2700.,2600.])

		ind = numpy.where(x>=xcutuv)
		if numpy.size(ind) > 0:
			k[ind] = 1.0 + c1 + (c2*x[ind]) + c3*((x[ind])**2)/( ((x[ind])**2 - (x0**2))**2 + (gamma**2)*((x[ind])**2 ))
			yspluv = 1.0 + c1 + (c2*xspluv) + c3*((xspluv)**2)/( ((xspluv)**2 - (x0**2))**2 + (gamma**2)*((xspluv)**2 ))
								   
			ind = numpy.where(x>=5.9)
			k[ind] += c4*(0.5392*((x[ind]-5.9)**2) + 0.05644*((x[ind]-5.9)**3))

		# Opt/NIR part
		ind = numpy.where(x<xcutuv)
		if numpy.size(ind) > 0:
			xsplopir = numpy.zeros(9)
			xsplopir[0] = 0.0
			xsplopir[1:10] = 1.0/numpy.array([2.198,1.65,1.25,0.81,0.65,0.55,0.44,0.37])

			# Values directly from Gordon et al. (2003)
#			ysplopir =  numpy.array([0.0,0.016,0.169,0.131,0.567,0.801,1.00,1.374,1.672])
			# K & J values adjusted to provide a smooth, non-negative cubic spline interpolation
			ysplopir =  numpy.array([0.0,0.11,0.169,0.25,0.567,0.801,1.00,1.374,1.672])

			tck = interpolate.splrep(numpy.hstack([xsplopir,xspluv]), numpy.hstack([ysplopir,yspluv]), k=3)
			k[ind] = interpolate.splev(x[ind], tck)
		
		if (Alambda == True):
			return( k*Av)
		else:
			return( k*Av*(numpy.log(10.)/2.5))

class RvFbumpLaw(ExtinctionLaw):
	""" Mixture of Fitzpatrick99 and Gordon03_SMCBar allowing to vary the
	bump amplitude in the extinction law
	"""
	def __init__(self):
		self.RvLaw = Fitzpatrick99()
		self.NoBumpLaw = Gordon03_SMCBar()

	def function(self, lamb, Av=1, Rv=3.1, Alambda=True, f_bump=0.5, **kwargs):
		"""
		Lamb as to be in microns!!!

		input:
		    lamb      <float>    wavelength of the extinction point !! in microns !!
		output:
		    tau       <float>    returns tau as in redflux = flux*exp(-tau)
		keywords:
		    Alambda   <bool>     returns +2.5*1./log(10.)*tau
		    Av        <float>    extinction value (def: 1.0)
		    Rv        <float>    extinction param. (def: 3.1)
		    f_bump    <float>    mixture fraction defining the bump amplitude
		"""
		return f_bump*self.RvLaw.function(lamb,Av=Av,Rv=Rv,Alambda=Alambda) + (1.-f_bump)*self.NoBumpLaw.function(lamb,Av=Av,Alambda=Alambda)

if __name__ == "__main__":
	# check that things look correct
	# -> make some plots
	import pylab

	x = (numpy.arange(100)/100.)*10. + 0.1
	lamb = 1./x

	ccm  = Cardelli()
	f99  = Fitzpatrick99()
	gsmc = Gordon03_SMCBar()

        fig = pylab.figure()
        plt = fig.add_subplot(1,1,1)
	
	Rv_vals = numpy.arange(2, 6, dtype=float)
	for Rv in Rv_vals:
		yccm = ccm.function(lamb, Rv=Rv)
		yf99 = f99.function(lamb, Rv=Rv)

		#pylab.plot(x,yccm)
		plt.plot(x, yf99, label='F99, Rv=%0.1f' % (Rv) )

	ygsmc = gsmc.function(lamb)
	plt.plot(x, ygsmc, label='G. SMC')

	mixlaw = RvFbumpLaw()
	ymix = mixlaw(lamb, Rv=3.1, f_bump=0.75)
	plt.plot(x, ymix, label='Mixture f(bump)=0.75')

	ymix = mixlaw(lamb, Rv=3.1, f_bump=0.5)
	plt.plot(x, ymix, label='Mixture f(bump)=0.5')

	ymix = mixlaw(lamb, Rv=3.1, f_bump=0.25)
	plt.plot(x, ymix, label='Mixture f(bump=0.25')

	plt.set_ylabel('A($\lambda$)/A(V)')
	plt.set_xlabel('1/x [$\mu$m$^{-1}$]')

	plt.legend(loc=0, frameon=False)

	pylab.show()
