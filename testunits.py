
__version__ = '0.1dev'

import numpy
from numpy import exp, pi 
import inspect
import itertools
from copy import deepcopy

from anased import *
import extinction
import grid
import output


#TODO: check the Alambda definition and update the calls is necessary

def getFake(g, idx, err=0.1, oAv = None, **kwargs):
	""" Generate a fake sed from a model grid 
	INPUTS:
		idx	int			index number on the grid
	OUTPUTS:
		fakein	int			the index of the model on the grid
		lamb	ndarray[float, ndim=1]	wavelength taken from the grid
		fakesed	ndarray[float, ndim=1]	resulting SED
		
	KEYWORDS:
		err	float			proportional error to consider on the fluxes
		oAv	ExtinctionLaw		if provided, apply extinction function using **kwargs

		**kwargs if provided, extra keywords are used to apply an extinction
	"""
	oAv      = extinction.Cardelli()
	lamb     = g.lamb
	fakein   = idx
	fakesed  = numpy.copy(g.seds[fakein,:])
	if (oAv is not None) & (len(kwargs)>0):
		tau      = oAv.function(g.lamb*1-4, Alambda=True, **kwargs) 
		fakesed *= exp(-tau)
	#magerr  = 0.05
	#fakeerr = fakesed * (1. - 10**(-0.4*magerr) ) 
	fakeerr =  err*fakesed
	return fakein, lamb, fakesed, fakeerr


def meshgrid(arrs, ravel=False ):
	""" Generate a n-dim grid from a list of n arrays containing the points
	to use. The gridding occurs by varying the last value first.

	INPUTS:
		args	list of arrays	a,b,c containing the grid points of each dimension
	OUTPUTS:
		ans	grid/array	output comparable to numpy.meshgrid or ndarray
	KEYWORDS:
		ravel	bool	return a 2d array if set, where each line is a n-d vector point
	"""
	arrs = tuple(arrs)
	lens = map(len, arrs)
	dim = len(arrs)
	sz = 1
	for s in lens:
		sz *= s
	ans = []
	for i, arr in enumerate(arrs):
		slc = [1]*dim
		slc[i] = lens[i]
		arr2 = numpy.asarray(arr).reshape(slc)
		for j, sz in enumerate(lens):
			if j != i:
				arr2 = arr2.repeat(sz, axis=j)
		ans.append(arr2)
	if ravel:
		return numpy.vstack( map( numpy.ravel, tuple(ans) ) ).T
	else:
		return tuple(ans)

def iter_Av_grid(g0, oAv, iterate=True, **kwargs):
	""" generate a grid by applying extinction values 
	INPUTS:
		g0	grid		initial spectral grid (unreddened)
		oAv	ExtinctionLaw	extinction law to apply

	OUTPUTS:
		git	iterator	iterator on SpectralGrid chunks
	
	KEYWORDS:
		**kwargs	grid point values.
				each keyword will be used to call oAv.function
	"""

	assert( isinstance(oAv, extinction.ExtinctionLaw)), 'Extinction must be an ExtinctionLaw instance, got %s'% type(oAv)
	assert( isinstance(g0,  grid.ModelGrid) ), 'Extinction must be an ExtinctionLaw instance, got %s'% type(oAv)

	#checking arguments of extinction.function
	av_args = [k for k in inspect.getargspec(extinction.RvFbumpLaw.function).args if k not in ['self', 'lamb', 'Alambda'] ]

	for k, v in kwargs.iteritems():
		assert( k in av_args ), 'Argument %s not in the extinction parameters' % k

	#generate the grid points
	if len(kwargs) > 1:
		gpts = meshgrid( kwargs.values() , ravel=True )
	else:
		gpts = kwargs.values()[0]

	#prepare output
	n0 = g0.grid.nrows
	npts = len(gpts)

	def gensubgrid( theta_av ):
		""" 
		**Partial function**
		Generate the chunk of grid corresponding to a given theta_av
		based on the initial grid g0
		"""
		_args = {}
		for e,k in enumerate(kwargs): 
			if hasattr(theta_av, '__iter__'):
				_args[k] = theta_av[e]
			else:
				_args[k] = theta_av
		tau        = oAv.function(g0.lamb*1-4, Alambda=True, **_args) 
		outputSEDs = g0.seds * numpy.exp(-tau) [None, :]
		t = deepcopy(g0.grid)
		for k, v in _args.iteritems():
			t.addCol( [ v ] * t.nrows, name=k) 
		g = grid.SpectralGrid()
		g.lamb = g0.lamb[:]
		g.seds = outputSEDs
		g.grid = t
		return g

	# If the grid is small enough, you can compute all at once and
	# concatenate the result. However, most of the time the grid is big
	# enough that the maximum of element per array is reached (~10^9 on
	# 64bits arch). So that the default is to return an iterator on the grid
	# chunks
	# This means each chunk is independent from the others as well, ergo
	# parallel jobs are possible.
	return itertools.imap( gensubgrid, gpts )

	

	
	


def test_seds(err=0.1):
	filters  = 'hst_wfc3_f225w hst_wfc3_f336w hst_acs_hrc_f475w hst_acs_hrc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split()


	#Load the initial model grid
	g0 = grid.FileSpectralGrid('/home/morgan/Work/anased/libs/kurucz2004.grid.fits')
	lamb = g0.lamb


	# define the extinction priors
	oAv = extinction.Cardelli()

	##define Av with a step size
	#Av = numpy.arange(0.,3., 0.1)	

	##define Av with a number of points
	Av = numpy.linspace(0, 1, 3)
	Rv = numpy.array([3.1])
	fb = numpy.array([0.5])	#only one value still needs to be array/list/tuple


	# Parameters from which we generate fake data
	## Number of fake SEDs to play with
	N   = 1
	## Extinction parameters are randomly drawn from the extinction space
	Av0    = numpy.random.uniform(0, len(Av), N)
	Rv0    = numpy.random.uniform(0, len(Rv), N)
	fb0    = numpy.random.uniform(0, len(fb), N)
	## Initial SEDs are randomly drawn from the model space
	fakein = numpy.random.randint(0, g0.grid.nrows, N)

	for tn in range(N):
		#fake DATA
		idx, l, fakesed, fakeerr = getFake(g, fakein[tn], Av=Av0[tn], Rv=Rv0[tn], f_bump=fb0[tn], err=err)
		mask = numpy.zeros(fakesed.shape, dtype=bool)
		#mask[3] = True
		#mask[2] = True


		r = numpy.empty( (g.seds.shape[0], len(Av)), dtype=float )
		with timeit('Likelihood Object %d' % tn):
			for k in range(len(Av)):
				r[:, k] = job(lamb[:], numpy.copy(fakesed), numpy.copy(fakeerr), mask, numpy.copy(g.seds), oAv, Av=Av[k], Rv=3.1) 

		output.fullTable('Tests/t%d' % tn, g, r, Av, lamb, fakesed, fakeerr, filters, Av0 = Av0[tn], orig = g.grid[idx])

