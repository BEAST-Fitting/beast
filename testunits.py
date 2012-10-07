
__version__ = '0.1dev'

import numpy
from numpy import log, log10, exp, pi 
#import numpy ufunc (c-coded for speed-up)
from numpy import multiply, subtract, add, divide, power

from anased import *
import grid
import output

def getFake(g, idx, Av0=1., Rv0=3.1, err=0.1):
	oAv      = extinction.Cardelli()
	lamb     = g.lamb
	fakein   = idx
	fakesed  = numpy.copy(g.seds[fakein,:])
	tau      = getFluxAttenuation(oAv, lamb, Av = Av0, Rv = Rv0)
	fakesed *= exp(-tau)
	#magerr  = 0.05
	#fakeerr = fakesed * (1. - 10**(-0.4*magerr) ) 
	fakeerr =  err*fakesed
	return fakein, lamb, fakesed, fakeerr


def test_seds(err=0.1):
	filters  = 'hst_wfc3_f225w hst_wfc3_f336w hst_acs_hrc_f475w hst_acs_hrc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split()


	#osl = stellib.BaSeL()
	oAv = extinction.Cardelli()
	g = grid.FileSpectralGrid('libs/SEDs_basel_padovaiso.fits')
	lamb = g.lamb

	Av = numpy.arange(0.,3., 0.1)	

	N   = 2
	Av0    = numpy.random.uniform(0,3,N)
	fakein = numpy.random.randint(0,g.grid.nrows,N)
	for tn in range(N):
		#fake DATA
		idx, l, fakesed, fakeerr = getFake(g, fakein[tn], Av0[tn], 3.1, err=err)
		mask = numpy.zeros(fakesed.shape, dtype=bool)
		#mask[3] = True
		#mask[2] = True


		r = numpy.empty( (g.seds.shape[0], len(Av)), dtype=float )
		with timeit('Likelihood Object %d' % tn):
			for k in range(len(Av)):
				r[:, k] = job(lamb[:], numpy.copy(fakesed), numpy.copy(fakeerr), mask, numpy.copy(g.seds), oAv, Av=Av[k], Rv=3.1) 

		output.fullTable('Tests/t%d' % tn, g, r, Av, lamb, fakesed, fakeerr, filters, Av0 = Av0[tn], orig = g.grid[idx])

