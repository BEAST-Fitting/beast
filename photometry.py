"""
This python module aims at generating the integrated photometry of spectra
"""

# TODO: replace trapz by simps

import pyfits, numpy, tables
import inspect, os
localpath = '/'.join(os.path.abspath(inspect.getfile(inspect.currentframe())).split('/')[:-1])

__default__      = localpath+'/libs/filters.hd5'
__default_vega__ = localpath+'/libs/vega.hd5'

class filter(object):
    """Class filter
    Define a filter by its name, wavelength and transmission
    """
    #----------------------------------------------------------------------
    def info(self):
    	""" display information about the current filter"""
	print "Filter object information:"
	print "   name: %s" % self.name
	print "   central wavelength: %f" % self.cl
	print "   norm: %f" % self.norm
	print "   pivot wavelength: %f" % self.lpivot  
	print "   definition contains %d points" % self.transmit.size

    def __repr__(self):
	    return "Filter: %s, %s" % (self.name, object.__repr__(self))


    def getFlux(self, slamb, sflux):
        """getFlux
        Integrate the flux within the filter and return the integrated energy
        INPUTS:
           slamb: spectrum wavelength definition domain
           sflux: associated flux
        OUTPUTS:
           <float>: Energy of the spectrum within the filter
        """
	if True in numpy.isinf(sflux):
		indinf = numpy.where(numpy.isinf(sflux))
		indfin = numpy.where(numpy.isfinite(sflux))
		sflux[indinf] = numpy.interp(slamb[indinf],
			slamb[indfin], sflux[indfin])
        ifT = numpy.interp(slamb, self.wavelength,self.transmit, 
		left=0.,right=0.)
	if True in (ifT > 0.):
		ind = numpy.where(ifT > 0.)
		a = numpy.trapz(slamb[ind]*ifT[ind]*sflux[ind], slamb[ind])
		b = numpy.trapz(slamb[ind]*ifT[ind], slamb[ind])
		if numpy.isinf(a) | numpy.isinf(b): print self.name, "Warn for inf value"
		return a/b 
	else:
		return 0.
    #----------------------------------------------------------------------        
    def __call__(self, slamb, sflux):
	return self.applyTo(slamb, sflux)

    def applyTo(self, slamb, sflux):
        """applyTo
        Apply filter to a spectrum
        INPUTS:
           slamb: spectrum wavelength definition domain
           sflux: associated flux
        OUTPUTS:
           [<float>]: new spectrum values accounting for the filter
        """
        ifT = numpy.interp(slamb, self.wavelength,self.transmit)
        return ifT*sflux
            
    #----------------------------------------------------------------------
    def __init__(self, wavelength, transmit, name=''):
        """Constructor"""
        self.name       = name
        self.wavelength = wavelength
        self.transmit   = transmit
        self.norm       = numpy.trapz(transmit,wavelength)
        self.cl         = numpy.trapz(wavelength*transmit, wavelength)/numpy.trapz(transmit,wavelength)
        self.lpivot     = numpy.sqrt(numpy.trapz(wavelength*transmit, wavelength)/numpy.trapz(transmit/wavelength,wavelength))
#                ---------------- end class filter ---------------


def __load__(fname, ftab):
	fnode = ftab.getNode('/filters/'+fname)
	return filter( fnode[:]['WAVELENGTH'], fnode[:]['THROUGHPUT'], name = fnode.name )
				       
def load_all_filters(filterLib = __default__):
	with tables.openFile(filterLib, 'r') as ftab:
		filters = [ __load__(fname, ftab) for fname in ftab.root.content.cols.TABLENAME ]
	return(filters)

def load_filters(names, filterLib = __default__):
	with tables.openFile(filterLib, 'r') as ftab:
		filters = [ __load__(fname, ftab) for fname in names ]
	return(filters)

class __newFilterTable__(tables.IsDescription):
	""" define table to store filter dataset """
	WAVELENGTH = tables.FloatCol(pos=0)
	THROUGHPUT = tables.FloatCol(pos=1)

def append_filter(lamb, flux, 
		   tablename, observatory, instrument, name, comment=None,
		   filterLib=__default__, 
		   updateVegaLib = True):
	""" 
	Edit the filter catalog and append a new one given by its transfer function
	"""
	ftab = tables.openFile(filterLib, 'a')
	contentTab = ftab.getNode('/content')
	if contentTab.readWhere('TABLENAME == tablename').size > 0:
		print '% '+sys.argv[0]+": Filter Table %s already exists.  Returning." % tablename
		return
		
	# Gen Filter object including relevant details
	filtInst = filter(lamb,flux, name=name)
	# Add a new line in the content table
	newRow = contentTab.row
	newRow['TABLENAME'] = tablename
	newRow['OBSERVATORY'] = observatory
	newRow['INSTRUMENT'] = instrument
	newRow['NAME'] = filtInst.name
	newRow['NORM'] = filtInst.norm
	newRow['CWAVE'] = filtInst.cl
	newRow['PWAVE'] = filtInst.lpivot
	if comment != None: newRow['COMMENT'] = comment
	newRow.append()
	contentTab.flush()
	# Create Table
	newTab = ftab.createTable('/filters', tablename, __newFilterTable__, title=filtInst.name,
		expectedrows = filtInst.wavelength.size)
	newRow = newTab.row
	for i in xrange(filtInst.wavelength.size):
		newRow["WAVELENGTH"] = filtInst.wavelength[i]
		newRow["THROUGHPUT"] = filtInst.transmit[i]
		newRow.append()
	newTab.flush()
	ftab.flush()
	ftab.close()
	print '% '+sys.argv[0]+": Filter %s added to %s." % (name, filterLib)
	if updateVegaLib == True:
		appendVegaFilter(filtInst)

#-------------------------------------------------------------------------------
# VEGA SPECTRUM and VEGA ZEROPOINTS
#-------------------------------------------------------------------------------
def analyseVegaSpectrum(w, f, filters):
	nFilters = len(filters)
	phot     = numpy.zeros((nFilters))
	cwave    = numpy.zeros((nFilters))
	fname    = []
	mag      = numpy.zeros((nFilters))
	for j in xrange(0,nFilters):
	    fname.append(filters[j].name)
	    cwave[j] = filters[j].cl
	    phot[j] = filters[j].getFlux( w, f )
	    mag[j] = -2.5 * numpy.log10(phot[j])
	return ({ 'fname':fname, 'cwave':cwave, 'lum':phot, 'mag':mag})
	
def appendVegaFilter(filtInst, VegaLib=__default_vega__):
	import tables
	vtab = tables.openFile(VegaLib, 'a')
	vl = vtab.root.spectrum[:]['WAVELENGTH']
	vf = vtab.root.spectrum[:]['FLUX']
	sedTab = vtab.getNode('/sed')
	fname = filtInst.name
	if sedTab.readWhere('FNAME == fname').size > 0:
		print '% '+sys.argv[0]+": Filter %s already exists.  Returning." % filtInst.name
		return
	
	data = analyseVegaSpectrum(vl, vf, [filtInst])
	newRow = sedTab.row
	newRow['FNAME'] = filtInst.name
	newRow['CWAVE'] = filtInst.cl
	newRow['LUM']   = data['lum'][0]
	newRow['MAG']   = data['mag'][0]
	newRow.append()
	sedTab.flush()
	vtab.close()
	print '% '+sys.argv[0]+": Filter %s added to %s." % (filtInst.name, VegaLib)


def extractPhotometry(l, f, filters, silent=True):
    """ Extract integrated fluxes from a given filter list """

    nFilters = len(filters)

    phot     = numpy.zeros(nFilters)

    #output
    for j in range(0,nFilters):
        phot[j] = filters[j].getFlux( l, f )

    return phot
