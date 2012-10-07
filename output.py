import mypickleshare
import mytables
import numpy

def weighted_percentile(data, wt, percentiles): 
	"""Compute weighted percentiles. 
	If the weights are equal, this is the same as normal percentiles. 
	Elements of the C{data} and C{wt} arrays correspond to 
	each other and must have equal length (unless C{wt} is C{None}). 

	@param data: The data. 
	@type data: A L{numpy.ndarray} array or a C{list} of numbers. 
	@param wt: How important is a given piece of data. 
	@type wt: C{None} or a L{numpy.ndarray} array or a C{list} of numbers. 
		All the weights must be non-negative and the sum must be 
		greater than zero. 
	@param percentiles: what percentiles to use.  (Not really percentiles, 
		as the range is 0-1 rather than 0-100.) 
	@type percentiles: a C{list} of numbers between 0 and 1. 
	@rtype: [ C{float}, ... ] 
	@return: the weighted percentiles of the data. 
	""" 
	assert numpy.greater_equal(percentiles, 0.0).all(), "Percentiles less than zero" 
	assert numpy.less_equal(percentiles, 1.0).all(), "Percentiles greater than one" 
	data = numpy.asarray(data) 
	assert len(data.shape) == 1 
	if wt is None: 
		wt = numpy.ones(data.shape, numpy.float) 
	else: 
		wt = numpy.asarray(wt, numpy.float) 
		assert wt.shape == data.shape 
		assert numpy.greater_equal(wt, 0.0).all(), "Not all weights are non-negative." 
	assert len(wt.shape) == 1 
	n = data.shape[0] 
	assert n > 0 
	i = numpy.argsort(data) 
	sd = numpy.take(data, i, axis=0) 
	sw = numpy.take(wt, i, axis=0) 
	aw = numpy.add.accumulate(sw) 
	if not aw[-1] > 0: 
		raise ValueError, "Nonpositive weight sum" 
	w = (aw-0.5*sw)/aw[-1] 
	spots = numpy.searchsorted(w, percentiles) 
	o = [] 
	for (s, p) in zip(spots, percentiles): 
		if s == 0: 
			o.append(sd[0]) 
		elif s == n: 
			o.append(sd[n-1]) 
		else: 
			f1 = (w[s] - p)/(w[s] - w[s-1]) 
			f2 = (p - w[s-1])/(w[s] - w[s-1]) 
			assert f1>=0 and f2>=0 and f1<=1 and f2<=1 
			assert abs(f1+f2-1.0) < 1e-6 
			o.append(sd[s-1]*f1 + sd[s]*f2) 
	return o 


def summuaryTable(fname, g, r, Av, lamb, fakesed, fakeerr, filters, Q = 'logg logT logL logM logA Av Z', **kwargs):

	_q = Q.split()
	_r = numpy.exp(r)
	_r /= _r.sum()

	d = {}
	r1 = _r.sum(1)
	rAv = _r.sum(0)

	for k in _q:
		if k.lower() != 'av': 
			m0, m1, m2 = weighted_percentile( g.grid[k], r1, [0.25,0.5,0.75] )
		else:
			m0, m1, m2 = weighted_percentile( Av, rAv, [0.25,0.5,0.75] )

		d['%s_p25' %k ] = m0
		d['%s_p50' %k ] = m1
		d['%s_p75' %k ] = m2

	for k in range(len(filters)):
		m0, m1, m2 = weighted_percentile( g.seds[:,k], r1, [0.25,0.5,0.75] )
		d['%s_p25' % (filters[k]) ] = m0
		d['%s_p50' % (filters[k]) ] = m1
		d['%s_p75' % (filters[k]) ] = m2

	t = mytables.Table(d, name='SED analysis summury table')
	for k in kwargs:
		t.header[k] = kwargs[k]
	
	t.write(fname)


	
def fullTable(fname, g, r, Av, lamb, fakesed, fakeerr, filters, **kwargs):
	d = mypickleshare.PickleShareDB(fname)
	d['lnp'] = r
	d['Av'] = Av
	d['filters'] = filters
	d['GRID'] = g.grid.header['SOURCE']
	d['LAMB'] = lamb
	d['INPUTSED'] = fakesed
	d['INPUTERR'] = fakeerr

	for k in kwargs:
		d[k] = kwargs[k]
	
	del d

#tests
filters  = 'hst_wfc3_f225w hst_wfc3_f336w hst_acs_hrc_f475w hst_acs_hrc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split()

