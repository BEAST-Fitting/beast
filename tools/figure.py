"""
This module aims to regroup helpful functions to easy create pretty figures.
	changeTickLength    - Rescale ticks
	cmap_discretize     - Return a discrete colormap from the continuous colormap cmap
	cmap_function       - Apply function to the given colormap
	cmap_fromList       - Generate a colormap from a given list of hexacolor values
	cplot               - Plot on the current graphe
	densityMap          - plot data points as a density map
	distPlot            - plot data points and marginal distributions aside axes
	getPercentileLevels - Determine the percentile values of an nd-array
	GridData            - object class gridding x,y,z data by interpolation
	Lasso		    - object class proxy to select points by drawing a region
	legend_font         - Change the legend font size
	paper1c             - Define params to make pretty 1 col. fig.
	plotDensity         - Density plot including scattered points
	plotMAP		    - Plot the distribution and MAP of a given sample
	reverseAxis	    - Reverse axis range
	second_axis         - Add a secondary axis
	set_axis_scale      - Set the scaling of the axis
	set_major_locator   - Set major tick positions
	set_minor_locator   - Set minor tick positions
	setMargins          - Tune subplot margins
	setNmajors          - set N major ticks on axis
	setNminors          - set N minor ticks between each major pairs

This module also gives access to matplotlib colormaps and pylab functions
	> from pylab import *"
	> import matplotlib.cm as cm"

HISTORY:
	04/10/11, MF	added plotDensity
	05/10/11, MF	optimized plotDensity single point selection
	06/10/11, MF	plotDensity includes percentile contouring options
			added function: getPercentileLevels
			updated doc
	10/10/11, MF	added function: plotMAP
			updated doc
	04/24/12, MF	added function Lasso (selection points from regions)

"""
__version__ = '1.1.4'
__author__ = 'MF'

from pylab import *
import numpy
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, MaxNLocator
from matplotlib.mlab import griddata

inches_per_cm = 0.393700787


def help():
	print __doc__()

def second_axis(axis=None, lim=None, title=None):
	"""
	second_axis - Add a secondary axis

	inputs:
		axis	{'x' | 'y'} to define the orientation
	output:
		ax	matplotlib.axis object
	keywords:
		lim	set this keyword to the range you want to fix
		titel	set this keyword to the desire title of the axis
	"""
	if axis != None:
		if axis == 'x':
			ax = twiny()
			if lim != None:
				ax.set_xlim(lim)
				ax.set_xlabel(title)
		if axis == 'y':
			ax = twinx()
			if lim != None:
				ax.set_ylim(lim)
				ax.set_ylabel(title)
    		draw_if_interactive()
		return ax

def cplot(*args,**kwargs):
	"""
	cplot - Plot on the current graphe
	This is an "alias" to gcf().gca().plot()
	"""
	return(gcf().gca().plot(*args,**kwargs))


def legend_font(size):
	"""
	legend_font - change the legend font size
	"""
	params = {'legend.fontsize': size}
	rcParams.update(params)

def reverseAxis(x=None, y=None, ax=None):
	"""
	reverseAxis - Reverse range of an axis
	"""
	if ax == None:
		ax = gca()
	if x != None:
		ax.set_xlim(ax.get_xlim()[::-1])
	if y != None:
		ax.set_ylim(ax.get_ylim()[::-1])
    	draw_if_interactive()

def set_axis_scale(xscale=None, yscale=None):
	"""
	set_axis_scale - Set the scaling of the x-axis: 'linear' | 'log' | 'symlog'
	"""
	ax = gca()
	if xscale != None:
		ax.set_xscale(xscale)
	if yscale != None:
		ax.set_yscale(yscale)
    	draw_if_interactive()

def set_major_locator(xval=None, yval=None, ax=None):
	"""
	set_major_locator - set major tick positions
	"""
	if ax == None:
		ax = gca()
	if xval != None:
		ax.xaxis.set_major_locator(MultipleLocator(xval))
	if yval != None:
		ax.yaxis.set_major_locator(MultipleLocator(yval))
    	draw_if_interactive()

def set_minor_locator(ax = None, xval=None, yval=None):
	"""
	set_minor_locator - set minor tick positions
	"""
	if ax == None:
		ax = gca()
	if xval != None:
		ax.xaxis.set_minor_locator(MultipleLocator(xval))
	if yval != None:
		ax.yaxis.set_minor_locator(MultipleLocator(yval))

    	draw_if_interactive()

def changeTickLength(size):
	"""
	changeTickLength - rescale ticks
	"""
	ax = gca()
	for line in ax.xaxis.get_ticklines(minor=False):
		s = line.get_markersize()
		line.set_markersize(size*s)
	for line in ax.xaxis.get_ticklines(minor=True):
		s = line.get_markersize()
		line.set_markersize(size*s)
	for line in ax.yaxis.get_ticklines(minor=False):
		s = line.get_markersize()
		line.set_markersize(size*s)
	for line in ax.yaxis.get_ticklines(minor=True):
		s = line.get_markersize()
		line.set_markersize(size*s)

    	draw_if_interactive()

def setNmajors(xval=None, yval=None, ax=None, mode='auto', **kwargs):
	"""
	setNmajors - set major tick number
	see figure.MaxNLocator for kwargs
	"""
	if ax == None:
		ax = gca()
	if (mode == 'fixed'):
		if xval != None:
			ax.xaxis.set_major_locator(MaxNLocator(xval, **kwargs))
		if yval != None:
			ax.yaxis.set_major_locator(MaxNLocator(yval, **kwargs))
	elif (mode == 'auto'):
		if xval != None:
			ax.xaxis.set_major_locator(AutoLocator(xval, **kwargs))
		if yval != None:
			ax.yaxis.set_major_locator(AutoLocator(yval, **kwargs))

    	draw_if_interactive()

def setNminors(nx=5, ny=5, ax=None, mode='auto'):
	"""
	setNminors - set N minor ticks between each major pairs
	mode 'auto' will determine how many ticks (4 or 5) to display
	mode 'fixed' will impose N ticks between the current majors.
		note that the fixed mode will set minor values until a tick
		reset.
	"""
	if ax == None:
		ax = gcf().gca()
	xMajorLoc = ax.xaxis.get_majorticklocs()
	yMajorLoc = ax.yaxis.get_majorticklocs()
	if nx != None:
		if (mode == 'fixed'): set_minor_locator(xval=(xMajorLoc[1]-xMajorLoc[0])/nx, ax=ax)
		if (mode == 'auto'): ax.xaxis.set_minor_locator(AutoMinorLocator())
	if ny != None:
		if (mode == 'fixed'): set_minor_locator(yval=(yMajorLoc[1]-yMajorLoc[0])/ny, ax=ax)
		if (mode == 'auto'): ax.yaxis.set_minor_locator(AutoMinorLocator())

    	draw_if_interactive()
"""
Define figure Size ========================================================
"""
def paper_map(scale=1.2):
    """
    paper1c - Define params to make pretty 1col fig publication
    """
    print 'define params'
    rescale  = scale/1.2
    fontSize = rescale*10.
    rc('figure', figsize=(3.*scale,2.3*scale))
    rc('figure.subplot', left=0.15, right=0.97, bottom=0.18, top=0.95)
    rc('lines', linewidth=0.5*rescale)
    rc('axes', linewidth=0.5*rescale)
    rc('font', size=fontSize, family='serif', weight='small')
    rc('xtick', labelsize='small')
    rc('ytick', labelsize='small')
    rc('legend', fontsize='x-small', borderpad=0.1, markerscale=1.,
    	fancybox=False)
    rc('text', usetex=True)
    rc('path', simplify=True)
    rc('image', aspect='auto')
    rc('ps', useafm=True, fonttype=3)
    rcParams['xtick.major.size']=4*rescale
    rcParams['xtick.minor.size']=2*rescale
    rcParams['ytick.major.size']=4*rescale
    rcParams['ytick.minor.size']=2*rescale
    rcParams['font.sans-serif']='Helvetica'
    rcParams['font.serif']='Helvetica'
    rcParams['text.latex.preamble']='\usepackage{pslatex}'
    #rc('ps', usedistiller="xpdf")
    #rcParams['ps.distiller.res']=80
    rcParams['image.resample']=True


def paper1c(sc=1.2):
    """
    paper1c - Define params to make pretty 1col fig publication
    """
    rescale  = sc/1.2
    fontSize = rescale*10.
    rc('figure', figsize=(3.*sc,2.3*sc))
    rc('figure.subplot', left=0.15, right=0.97, bottom=0.18, top=0.95)
    rc('lines', linewidth=0.5*rescale)
    rc('axes', linewidth=0.5*rescale)
    rc('font', size=fontSize, family='serif', weight='small')
    rc('xtick', labelsize='small')
    rc('ytick', labelsize='small')
    rc('legend', fontsize='x-small', borderpad=0.1, markerscale=1.,
 	   fancybox=False)
    rc('text', usetex=True)
    #rc('path', simplify=True)
    rc('image', aspect='auto')
    rc('ps', useafm=True, fonttype=3)
    rcParams['xtick.major.size']=8*rescale
    rcParams['xtick.minor.size']=4*rescale
    rcParams['ytick.major.size']=8*rescale
    rcParams['ytick.minor.size']=4*rescale
    rcParams['font.sans-serif']='Helvetica'
    rcParams['font.serif']='Helvetica'
    rcParams['text.latex.preamble']='\usepackage{pslatex}'

def screen(sc=2, figsize=[3.,2.3]):
    """
    paper1c - Define params to make pretty 1col fig publication
    """
    rescale  = sc
    fontSize = rescale*10
    rc('figure', figsize=(figsize[0]*sc,figsize[1]*sc))
    rc('figure.subplot', left=0.12, right=0.95, bottom=0.15, top=0.95)
    rc('lines', linewidth=1)
    rc('axes', linewidth=1.)
    rc('font', size=fontSize*sc/3., family='serif', weight='small')
    rc('xtick', labelsize='small')
    rc('ytick', labelsize='small')
    rc('legend', fontsize='x-small', borderpad=0.1, markerscale=1.,
 	   fancybox=False)
    rc('text', usetex=True)
    #rc('path', simplify=True)
    rc('image', aspect='auto')
    rc('ps', useafm=True, fonttype=3)
    rcParams['xtick.major.size']=4*rescale
    rcParams['xtick.minor.size']=2*rescale
    rcParams['ytick.major.size']=4*rescale
    rcParams['ytick.minor.size']=2*rescale
    rcParams['font.sans-serif']='Helvetica'
    rcParams['font.serif']='Helvetica'
    rcParams['text.latex.preamble']='\usepackage{pslatex}'

def setMargins(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None):
	"""
	Tune the subplot layout via the meanings (and suggested defaults) are::

	    left  = 0.125  # the left side of the subplots of the figure
	    right = 0.9    # the right side of the subplots of the figure
	    bottom = 0.1   # the bottom of the subplots of the figure
	    top = 0.9      # the top of the subplots of the figure
	    wspace = 0.2   # the amount of width reserved for blank space between subplots
	    hspace = 0.2   # the amount of height reserved for white space between subplots

	The actual defaults are controlled by the rc file

	"""
	subplots_adjust(left, bottom, right, top, wspace, hspace)
    	draw_if_interactive()

"""
Define ColorMap functions =====================================================
"""

def cmap_function(function = lambda x:x,cmap = None):
    """
    cmap_function - apply function to the given colormap
    Applies function (which should operate on vectors of shape 3:
    [r, g, b], on colormap cmap. This routine will break any discontinuous
    points in a colormap.
    """
    if cmap == None:
    	return
    cdict = cmap._segmentdata
    step_dict = {}
    # Firt get the list of points where the segments start or end
    for key in ('red','green','blue'):         step_dict[key] = map(lambda x: x[0], cdict[key])
    step_list = reduce(lambda x, y: x+y, step_dict.values())
    step_list = numpy.array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step : numpy.array(cmap(step)[0:3])
    old_LUT = numpy.array(map( reduced_cmap, step_list))
    new_LUT = numpy.array(map( function, old_LUT))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i,key in enumerate(('red','green','blue')):
        this_cdict = {}
        for j,step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j,i]
            elif new_LUT[j,i]!=old_LUT[j,i]:
                this_cdict[step] = new_LUT[j,i]
        colorvector=  map(lambda x: x + (x[1], ), this_cdict.items())
        colorvector.sort()
        cdict[key] = colorvector

    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)

def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet.
        N: Number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    cdict = cmap._segmentdata.copy()
    # N colors
    colors_i = linspace(0,1.,N)
    # N+1 indices
    indices = linspace(0,1.,N+1)
    for key in ('red','green','blue'):
        # Find the N colors
        D = array(cdict[key])
        I = interpolate.interp1d(D[:,0], D[:,1])
        colors = I(colors_i)
        # Place these colors at the correct indices.
        A = zeros((N+1,3), float)
        A[:,0] = indices
        A[1:,1] = colors
        A[:-1,2] = colors
        # Create a tuple for the dictionary.
        L = []
        for l in A:
            L.append(tuple(l))
        cdict[key] = tuple(L)
    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)

def cmap_fromList(colorList=None, name='mycm', discrete=False):
	"""
	cmap_fromList - Generate a colormap from a given list of hexacolor values
	inputs:
		colorList	list of hexadecimal color values
	outputs:
		cmap		matplotlib.colormap object
	keywords:
		name		name of this colormap
		discrete	if set, this will create a discrete table,
				avoiding the interpolation (ListedColormap)
	"""
	if colorList == None:
		colorList =['#6C9F7B','#783F68','#F4F4F2','#22F322','#F3F322','#0000F3']
	if discrete==False:
		cdic = {'red':[], 'green':[], 'blue':[]}
		x    = numpy.linspace(0.,1., len(colorList))
		for ic in range(len(colorList)):
			r,g,b = matplotlib.colors.hex2color(colorList[ic])
			cdic['red'].append((x[ic],  r, r))
			cdic['green'].append((x[ic],  g, g))
			cdic['blue'].append((x[ic],  b, b))
		return(matplotlib.colors.LinearSegmentedColormap(name,cdic,1024))
	else:
		return(matplotlib.colors.ListedColormap(colorList, name=name ))

def logScale(x):
	h = x.copy()
	ind = where(x > 0.)
	h[ind] = log10(h[ind])
	return(h)

def binaryScale(x):
	h = x.copy()
	h[where(x > 0.)] = 1.
	h[where(x <= 0.)] = 0.
	return(h)

def sinhScale(x):
	h = x.copy()
	ind = where(x > 0.)
	h[ind] = sinh(h[ind])
	return(h)

def sqrtScale(x):
    h = x.copy()
    ind = numpy.where(x>0)
    h[ind] = numpy.sqrt(h[ind])
    return(h)

try:
	import zscale
	def zScale(x,nsamples=1000, contrast=0.25, bpmask=None, zmask=None):

		from zscale import zscale
		zmin, zmax = zscale(x, nsamples=1000, contrast=0.25, bpmask=None, zmask=None)
		y = x.copy()
		# set all points less than zmin to zmin and points greater than
		# zmax to zmax
		y = numpy.where(x > zmin, y, zmin)
		y = numpy.where(x < zmax, y, zmax)
		return y
except ImportError:
	pass

def densityMap(x,y, weights=None, bins=200, cmap=None, interpolate=None,
		scaling=lambda x:x, axes=None, *args, **kargs):

	ind = numpy.where(numpy.isfinite(x) & numpy.isfinite(y))
	try:
		h, xe, ye = histogram2d(x[ind],y[ind], weights=weights,bins=bins,*args,**kargs )
	except:
		h, xe, ye = histogram2d(x,y, weights=weights,bins=bins,*args,**kargs )
	extent  = [xe[0], xe[-1], ye[0], ye[-1]]
	if cmap == None:
		cmap = cm.Blues

	if axes == None:
		imshow(numpy.transpose(scaling(h)),extent=extent, origin='lower',
				aspect='auto', cmap=cmap, interpolation=interpolate,
				*args, **kargs)
	else:
		axes.imshow(numpy.transpose(scaling(h)),extent=extent, origin='lower',
				aspect='auto', cmap=cmap,*args, **kargs)

	return(numpy.transpose(h), extent)

class distPlot:

	def __init__(self,x, y, marker='o', xlabel=None, ylabel = None,
		ms=2, lw=0.2, bins=20, mec='0.0',mfc='0.0',
		histtype='step', hec=None, hls='solid', normed=False,
		plottype='dots',scaling=lambda x:x, *args, **kwargs):

		self.nx, self.bx  = [],[]
		self.ny, self.by  = [],[]
		self.fig   = figure()
		self.main  = axes([0.15,0.15,0.6,0.6])
		self.xplot = axes([0.15,0.76,0.6,0.20], sharex=self.main)
		self.yplot = axes([0.75+2.3/3.*0.01,0.15,0.20*2.3/3.,0.6],
					sharey=self.main)
		if plottype=='dots':
			self.main.plot(x,y, marker, ms=ms, lw=lw, mec=mec, mfc=mfc,
						*args, **kwargs)
		if plottype=='density':
			densityMap(x,y, bins=bins, axes=self.main,
			scaling=scaling, *args,**kwargs)

		if (xlabel != None): self.main.set_xlabel(xlabel)
		if (ylabel != None): self.main.set_ylabel(ylabel)

		if (hec == None): hec = mec

		n,b,p = self.xplot.hist(x, bins=bins, histtype=histtype,
					ec=hec, ls=hls,normed=normed)
		self.nx.append(n)
		self.bx = b
		n,b,p = self.yplot.hist(y, bins=bins, histtype=histtype,
					ec=hec, ls=hls, normed=normed,
					orientation='horizontal')
		self.ny.append(n)
		self.by = b
		setp(self.xplot.get_xticklabels()+
				self.yplot.get_yticklabels(), visible=False)

		setNmajors(ax = self.xplot, yval=2)
		setNmajors(ax = self.yplot, xval=2)

	def add(self, x, y, marker='o', xlabel=None, ylabel = None,
		ms=2, lw=0.2, bins=None, mec='0.0',mfc='0.0',
		histtype='step', hec='0.0', hls='solid', normed=False,
		*args, **kwargs):

		self.main.plot(x,y, marker, ms=ms, lw=lw, mec=mec, mfc=mfc,
					*args, **kwargs)
		if bins == None:
			bx = self.bx
		else:
			bx = bins
		n,b, p = self.xplot.hist(x, bins=bx, histtype=histtype,
					ec=hec, ls=hls,normed=normed)
		self.nx.append(n)
		if bins == None:
			by = self.by
		else:
			by = bins
		n,b, p = self.yplot.hist(y, bins=by, histtype=histtype,
					ec=hec, ls=hls, normed=normed,
					orientation='horizontal')
		self.ny.append(n)
	def setNminors(self, n, mode='auto'):
		setNminors(ax = self.main, nx=n, ny=n, mode=mode)
		setNminors(ax = self.yplot, nx=n, mode=mode)
		setNminors(ax = self.xplot, ny=n, mode=mode)

	def legend(self, frame=False, loc=[1.0, 1.05],*args, **kwargs):
		legend = self.main.legend(loc=loc, *args,**kwargs)
		legend.draw_frame(frame)
		return legend


def demo():
	plot(arange(10))
	xlabel(r'$x-axis$')
	ylabel(r'y-axis')
	setNminors()
	show()



# ============================ beta libs ============================
from matplotlib.ticker import Locator

class AutoLocator(MaxNLocator):
    def __init__(self, nbins=9, steps=[1, 2, 5, 10], **kwargs):
        MaxNLocator.__init__(self, nbins=nbins, steps=steps, **kwargs )

class AutoMinorLocator(Locator):
    """
    Dynamically find minor tick positions based on the positions of
    major ticks. Assumes the scale is linear and major ticks are
    evenly spaced.
    """
    def __init__(self, n=None):
	    	self.ndivs = n

    def __call__(self):
        'Return the locations of the ticks'
        majorlocs = self.axis.get_majorticklocs()
        try:
            majorstep = majorlocs[1] - majorlocs[0]
        except IndexError:
            raise ValueError('Need at least two major ticks to find minor '
                             'tick locations')
        # see whether major step should be divided by 5, 4 or 2. This
        # should cover most cases.
	if self.ndivs is None:
            x = int(round(10 ** (np.log10(majorstep) % 1)))
            if x in [1, 5, 10]:
                ndivs = 5
            else:
                ndivs = 4
        else:
            ndivs = self.ndivs

        minorstep = majorstep / ndivs

	vmin, vmax = self.axis.get_view_interval()
        if vmin > vmax:
            vmin,vmax = vmax,vmin

        t0 = majorlocs[0]
        tmin = np.ceil((vmin - t0) / minorstep) * minorstep
        tmax = np.floor((vmax - t0) / minorstep) * minorstep
        locs = np.arange(tmin, tmax, minorstep) + t0
        cond = np.abs((locs - t0) % majorstep) > minorstep/10.0
        locs = locs.compress(cond)

        return self.raise_if_exceeds(np.array(locs))

	"""
        tmin = majorlocs[0] - majorstep
        tmax = majorlocs[-1] + majorstep
        locs = np.arange(tmin, tmax, minorstep)
        vmin, vmax = self.axis.get_view_interval()
        if vmin > vmax:
            vmin,vmax = vmax,vmin

        return locs[(vmin < locs) & (locs < vmax)]
	"""

def minorticks_on():
    """
    Display minor ticks on the current plot.

    Displaying minor ticks reduces performance; turn them off using
    minorticks_off() if drawing speed is a problem.
    """
    gca().minorticks_on()
    draw_if_interactive()

def minorticks_off():
    """
    Remove minor ticks from the current plot.
    """
    gca().minorticks_off()
    draw_if_interactive()



# ========================================== BETA LIBS END =====================

def theme(type='beauty', ax=None):
	if ax == None:
		ax = gca()
	setNmajors(5,5, ax=ax)
	setNminors(5,5, ax=ax)
	draw_if_interactive()

def plotCorr(l, pars, plotfunc=None, lbls=None, *args, **kwargs):
	""" Plot correlation matrix between variables
		inputs
			l  -- dictionary of variables (could be a Table)
		     pars  -- parameters to use

		*args, **kwargs are forwarded to the plot function
	"""

	fontmap = {1:10, 2:8, 3:6, 4:5, 5:4}
	if not len(pars)-1 in fontmap:
		fontmap[len(pars)-1] = 3

	if lbls == None:
		lbls = pars

	k = 1
	axes = numpy.empty((len(pars),len(pars)), dtype=object)
	for j in range(len(pars)):
		newrow=True
		for i in range(len(pars)):
			if i>j:
				if j>0:
					sharex = axes[j-1, i]
				else:
					sharex = None
				if i>0:
					sharey = axes[j, i-1]
				else:
					sharey = None
				ax = subplot(len(pars)-1,len(pars)-1,k, sharey=sharey, sharex=sharex)
				axes[j,i] = ax
				if plotfunc == None:
					plot(l[pars[i]],l[pars[j]],',',**kwargs)
				else:
					plotfunc(l[pars[i]],l[pars[j]],*args, **kwargs)

				theme(ax=ax)
				tlabels = gca().get_xticklabels()
				setp(tlabels, 'fontsize', 2*fontmap[len(pars)-1])
				tlabels = gca().get_yticklabels()
				setp(tlabels, 'fontsize', 2*fontmap[len(pars)-1])
				if not newrow:
					setp(ax.get_xticklabels()+ax.get_yticklabels(), visible=False)
				else:
					xlabel(lbls[i], fontsize=2.*fontmap[len(pars)-1])
					ylabel(lbls[j], fontsize=2.*fontmap[len(pars)-1])
					newrow=False
				N = int(0.5*fontmap[len(pars)-1])
				if N <2: N = 2
				setNmajors(N,N,ax=ax, prune='both')
			if i!=j:
				k+=1
	setMargins(hspace=0.0, wspace=0.0)

def getPercentileLevels(h, frac=[0.5, 0.65, 0.95, 0.975]):
	"""
	Return image levels that corresponds to given percentiles values
	Uses the cumulative distribution of the sorted image density values
	Hence this works also for any nd-arrays
	inputs:
		h	array
	outputs:
		res	array containing level values
	keywords:
		frac	sample fractions (percentiles)
			could be scalar or iterable
			default: 50%, 65%, 95%, and 97.5%

	"""
	if getattr(frac, '__iter__', False):
		return numpy.asarray( [getPercentileLevels(h, fk) for fk in frac])

	assert( (frac >= 0.) & (frac <1.)), "Expecting a sample fraction in 'frac' and got %f" %frac
	# flatten the array to a 1d list
	val = h.ravel()
	# inplace sort
	val.sort()
	#reverse order
	rval = val[::-1]
	#cumulative values
	cval = rval.cumsum()
	#retrieve the largest indice up to the fraction of the sample we want
	ind = numpy.where(cval <= cval[-1]*float(frac))[0].max()
	res = rval[ind]
	del val, cval, ind, rval
	return res

def plotMAP(x, ax = None, error=0.01, frac =[0.65,0.95, 0.975], hpd=True,
	hist={'histtype':'step'}, vlines={}, fill={}, optbins={'method':'freedman'}, *args, **kwargs):
	""" Plot the MAP of a given sample and add statistical info
	If not specified, binning is assumed from the error value or using
	mystats.optbins if available.
	if mystats module is not available, hpd keyword has no effect

	inputs:
		x 	dataset
	keywords
		ax	axe object to use during plotting
		error	error to consider on the estimations
		frac	fractions of sample to highlight (def 65%, 95%, 97.5%)
		hpd	if set, uses mystats.hpd to estimate the confidence intervals

		hist	keywords forwarded to hist command
		optbins	keywords forwarded to mystats.optbins command
		vlines	keywords forwarded to vlines command
		fill	keywords forwarded to fill command
		"""
	try:
		import mystats
		statsMod = True
	except ImportError:
		hpd = False
		statsMod = False

	_x = numpy.ravel(x)
	if ax == None:
		ax = gca()
	if not ('bins' in hist):
		if statsMod == True:
			bins = mystats.optbins(x,method=optbins['method'], ret='N')
		else:
			bins = numpy.linspace(numpy.min(_x), numpy.max(_x), numpy.round((numpy.max(_x)-numpy.min(_x))/error))
		n, b, p = ax.hist(_x, bins=bins, *args, **hist)
	else:
		n, b, p = ax.hist(_x, *args, **hist)
	c = 0.5*(b[:-1]+b[1:])
	dc = 0.5*(b[:-1]-b[1:])
	ind = n.argmax()
	_ylim = ylim() #(n.min(), n.max())
	if hpd == True:
		ax.vlines(mystats.hpd(_x,1-0.01), _ylim[0], _ylim[1], **vlines)
		for k in frac:
			nx = mystats.hpd(_x, 1.-k)
			ax.fill_between(nx, _ylim[0], _ylim[1], alpha=0.4/float(len(frac)), zorder=-1, **fill)
	elif hpd == False:
		ax.vlines(c[ind], _ylim[0], _ylim[1], **vlines)
		cx = c[ n.argsort() ][::-1]
		cn = n[ n.argsort() ][::-1].cumsum()
		for k in frac:
			sx = cx[numpy.where(cn <= cn[-1]*float(k))]
			sx = [sx.min(), sx.max()]
			ax.fill_between(sx, _ylim[0], _ylim[1], alpha=0.4/float(len(frac)), zorder=-1, **fill)
	theme()
	xlabel(r'Values')
	ylabel(r'Counts')

def plotDensity(x,y, bins=100, ax = None,
	Nlevels = None, levels=None, frac=None,
	contour = {'colors':'0.0', 'linewidths':0.5},
	contourf= {'cmap':cm.Greys_r},
	scatter = {'c':'0.0', 's':0.5, 'edgecolor':'None'},
	*args, **kwargs	):
	"""
	Plot a the density of x,y given certain contour paramters and includes
	individual points (not represented by contours)

	inputs:
		x,y	data to plot

	keywords:
		bins	bin definition for the density histogram
		ax	use a specific axis
		Nlevels	the number of levels to use with contour
		levels	levels
		frac	percentiles to contour if specified

		Extra keywords:
		*args, **kwargs forwarded to histogram2d
		**contour       forwarded to contour function
		**contourf      forwarded to contourf function
		**plot          forwarded to contourf function

	"""
	if ax == None:
		ax = gca()

	h, xe, ye = numpy.histogram2d(x,y, bins=bins, *args, **kwargs)

	if (Nlevels == None) & (levels == None) & (frac == None):
		levels = numpy.sort(getPercentileLevels(h))
	elif (Nlevels != None) & (levels == None) & (frac == None):
		levels = numpy.linspace(2., h.max(), Nlevels)[1:].tolist()+[h.max()]
	elif (frac != None):
		levels = getPercentileLevels(h, frac=frac)


	assert( getattr(levels, '__iter__', False) ), "Expecting levels variable to be iterable"

	if levels[-1] != h.max():
		levels = list(levels)+[h.max()]

	if isinstance(contourf, dict):
		cont = ax.contourf(h.T, extent=[xe[0],xe[-1], ye[0],ye[-1]],
					levels=levels, **contourf)
	else:
		cont = None
	if isinstance(contour, dict):
		ax.contour(h.T, extent=[xe[0],xe[-1], ye[0],ye[-1]],
					levels=levels,
					**contour)

	ind = numpy.asarray([False]*len(x))

	if cont != None:
		nx = numpy.ceil(numpy.interp(x,0.5*(xe[:-1]+xe[1:]),range(len(xe)-1)))
		ny = numpy.ceil(numpy.interp(y,0.5*(ye[:-1]+ye[1:]),range(len(ye)-1)))
		nh = [ h[nx[k],ny[k]] for k in range(len(x)) ]
		ind = numpy.where(nh < numpy.min(levels))
		ax.scatter(x[ind], y[ind], **scatter)
	else:
		ax.plot(x, y, **scatter)


class GridData():
	""" Implement a grid data interpolation
	This class is able to generate a regular grid from sparse data points
	"""
	def __init__(self, x,y,z, bins=100,
			xlim=None, ylim=None, zlim=None,
			xlabel=None, ylabel=None,zlabel=None):

		assert(numpy.size(x) == numpy.size(y) == numpy.size(z))

		self.x = numpy.copy(x)
		self.y = numpy.copy(y)
		self.z = numpy.copy(z)

		self.xi = None
		self.yi = None
		self.zi = None

		self.set_bins(bins, rebuild=False)

		self.xlabel = xlabel
		self.ylabel = ylabel
		self.zlabel = zlabel

		if xlim==None: xlim = [numpy.min(x), numpy.max(x)]
		if ylim==None: ylim = [numpy.min(y), numpy.max(y)]
		if zlim==None: zlim = [numpy.min(z), numpy.max(z)]

		self.xlim = xlim
		self.ylim = ylim
		self.zlim = zlim


	def __del__(self):
		del self.x, self.y, self.z, self.xi, self.yi, self.zi
		del self.xbin, self.ybin
		del self.xlim, self.ylim, self.zlim

	def set_bins(self, bins=None, xbin=None, ybin=None, rebuild=True):

		if bins:
			if numpy.shape(bins) == ():
				self.xbin = self.ybin = bins
			else:
				assert(numpy.size(bins) == 2)
				self.xbin = bins[0]
				self.ybin = bins[1]
		if xbin: self.xbin=xbin
		if ybin: self.ybin=ybin

		if rebuild: self.make(force=True)

	def make(self, force=False):
		if ((force == True) | (self.xi == None) | (self.yi==None) ):
			self.makeGrid()
			return
		else:
			return

	def makeGrid(self):
		self.xi = numpy.linspace(self.xlim[0],self.xlim[1],self.xbin)
		self.yi = numpy.linspace(self.ylim[0],self.ylim[1],self.ybin)
		self.zi = griddata(self.x,self.y,self.z,self.xi,self.yi)

	def getExtent(self):
		return [self.xi.min(),self.xi.max(),self.yi.min(), self.yi.max()]

	def imshow(self, ax=None, origin='lower', aspect='auto',
		cmap=cm.jet, addColorbar=True, *args, **kwargs):

		self.make()
		if ax==None: ax = gca()
		im = imshow(self.zi,
				extent=self.getExtent(), origin=origin,
				aspect=aspect,	cmap=cmap, *args, **kwargs)
		if addColorbar:
			cb = colorbar()
			if self.zlabel: cb.set_label(self.zlabel)

		if self.xlabel: ax.set_xlabel(self.xlabel)
		if self.ylabel: ax.set_ylabel(self.ylabel)
		ax.set_xlim(self.xlim)
		ax.set_ylim(self.ylim)
		theme(ax=ax)
		return im, cb

	def contour(self, fill=None, *args, **kwargs):

		self.make()

		if fill:
			cs = contourf(self.xi,self.yi,self.zi,*args, **kwargs)
		else:
			cs = contour(self.xi,self.yi,self.zi,*args, **kwargs)
		theme()
		return cs

	def scatter(self, marker='o', c=None, s=30, cmap=cm.jet, *args, **kwargs):
		if c == None:
			ret = scatter(self.x,self.y,marker='o',c=self.z, s = s,
				cmap=cmap, *args, **kwargs)
		else:
			ret = scatter(self.x,self.y,marker='o',c=c, cmap=cmap, s=s,
				*args, **kwargs)
		theme()
		return ret

def gridTest():
	""" Demo of the GridData Class """
	# make up some randomly distributed data
	npts = 200
	x = numpy.random.uniform(-2,2,npts)
	y = numpy.random.uniform(-2,2,npts)
	z = x*numpy.exp(-x**2-y**2)

	grid = GridData(x,y,z, bins=100)
	grid.imshow()
	grid.contour(15,linewidths=0.5,colors='k')
	#grid.contour(15,cmap=cm.jet, fill=True)
	scatter(x,y,marker='o',c=z,s=30, cmap=cm.jet)
	title('griddata test (%d points)' % npts)
	show()
	return grid

#============================ LASSO PLOT ==============================
from matplotlib.widgets import Lasso
from matplotlib.nxutils import points_inside_poly
from numpy import nonzero,array

class LassoPlot(object):
    """ This class is designed to select the datapoints by drawing the region
    around it.

    Example usage:

    > plot(xs, ys, ps=3 ) # first plot the data
    > las = lasso_plot.lasso_plot(xs,ys)

    Now click on the plot and do not release the mouse button till
    you draw your region
    After that the variable las.ind will contain the indices of those
    points inside the region and las.verts will contain the vertices of the
    polygon you've just drawn """

    def __init__(self, xs, ys, ax=None):
        self.axes = ax or gca()
	if len(self.axes.lines+self.axes.collections) == 0:
		self.axes.plot(xs,ys, 'o')
        self.canvas = self.axes.figure.canvas
        self.xys = array([xs,ys]).T#[d for d in zip(xs,ys)]
        fig = self.axes.figure
        self.cid = self.canvas.mpl_connect('button_press_event', self.onpress)
        self.ind = None
        self.mask = None
        self.verts = None
	show()

    def callback(self, verts):
        mask = points_inside_poly(self.xys, verts)
        ind = nonzero(mask)[0]
        self.canvas.draw_idle()
        self.canvas.widgetlock.release(self.lasso)
        self.canvas.mpl_disconnect(self.cid)
        del self.lasso
        del self.xys
        self.verts = verts
        self.ind = ind
        self.mask = mask

    def inside(self, xs,ys):
        tmpxys = zip(xs,ys)
        return points_inside_poly(tmpxys, self.verts)

    def onpress(self, event):
        if self.canvas.widgetlock.locked(): return
        if event.inaxes is None: return
        self.lasso = Lasso(event.inaxes, (event.xdata, event.ydata), self.callback)
        # acquire a lock on the widget drawing
        self.canvas.widgetlock(self.lasso)

import matplotlib.transforms as mtransforms

def on_draw(event):
   bboxes = []
   for label in labels:
       bbox = label.get_window_extent()
       # the figure transform goes from relative coords->pixels and we
       # want the inverse of that
       bboxi = bbox.inverse_transformed(fig.transFigure)
       bboxes.append(bboxi)

   # this is the bbox that bounds all the bboxes, again in relative
   # figure coords
   bbox = mtransforms.Bbox.union(bboxes)
   if fig.subplotpars.left < bbox.width:
       # we need to move it over
       fig.subplots_adjust(left=1.1*bbox.width) # pad a little
       fig.canvas.draw()
   return False

def autobbox(fig):
	fig.canvas.mpl_connect('draw_event', on_draw)





try:
	__figure_loaded
except NameError:
	__figure_loaded = True
	#print "Figure Package Loaded."
	#print "Figures are optimized for screen display."
	#screen()

