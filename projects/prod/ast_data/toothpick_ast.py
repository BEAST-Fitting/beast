import numpy as np
from scipy.spatial import cKDTree
from beast.external.eztables import Table
from beast.core.grid import FileSEDGrid
from beast.core.vega import Vega
from knn import KDTreeInterpolator
import tables

def _prepare_x(x):
        """ make sure the data is correctly presented """
        xi = np.asarray(x)
        shape = xi.shape

        if len(shape) < 1:
            raise ValueError('single value passed')

        if len(shape) == 1:
            return xi[..., np.newaxis]
        else:
            return xi

def nearest_neighbors(x,k=10,eps=0.):
    tree = cKDTree(x)
    d,ind = tree.query(x,k=k,eps=eps)
    return ind



def get_input_table(f):
    """get_input_table

    keywords
    --------

    f: str
        file to open

    returns
    -------

    tab: eztable.Table
        table instance
    """
    return Table(f)


def get_info(tab):
    """get_info - retrieve useful columns from the AST files

    Output information contains:
        MAG1IN      Input mag in filter 1
        MAG1OUT     Output mag in filter 1
        MAG2IN      Input mag in filter 1
        MAG2OUT     Output mag in filter 2
        MAG1_STD    Output uncertainty in filter 1
        MAG2_STD    Output uncertainty in filter 2
        BIAS1       Output - Input in filter 1
        BIAS2       Output - Input in filter 2

    keywords
    --------

    tab: eztable.Table instance
        table to extract from
        (can be replaced by a direct access at some point)

    returns
    -------
    res: ndarray
        array that contains the output array

    fields: sequence of strings
        string corresponding to the ordered content of the array

    names: sequence of strings
        strings used to make nice plot labels
    """

    #extract useful info
    fields = 'RA DEC MAG1IN MAG1OUT MAG2IN MAG2OUT MAG1_STD MAG2_STD'.split()
    res = np.asarray([ tab[k] for k in fields ])

    #add sub products
    sub = np.asarray([ res[2] - res[3], res[4] - res[5]])
    res = np.vstack([res, sub])

    fields += 'BIAS1 BIAS2'.split()
    names = 'RA DEC MAG1$_{IN}$ MAG1$_{OUT}$ MAG2$_{IN}$ MAG2$_{OUT}$ $\sigma_{MAG1}$ $\sigma_{MAG2}$ $\mu_{MAG1}$ $\mu_{MAG_2}$'.split()
    return res, fields, names

def get_AST_dict(files):
    """
    Returns a separate dictionnary containing AST data for each camera
    """

    tab_vis = get_input_table(files[1])
    tab_ir = get_input_table(files[2])
    tab_uv = get_input_table(files[0])
    res1,fields1,names1 = get_info(tab_uv)
    res2,fields2,names2 = get_info(tab_vis)
    res3,fields3,names3 = get_info(tab_ir)

    d_uv = dict([(nk,res1[ek]) for ek,nk in enumerate(fields1)])
    d_vis = dict([(nk,res2[ek]) for ek,nk in enumerate(fields2)])
    d_ir = dict([(nk,res3[ek]) for ek,nk in enumerate(fields3)])
    
    return d_uv,d_vis,d_ir

def NDinterp(x, y):
    """ generate an interpolator but taking care of non-finite values """
    ind = np.isfinite(y)
    inter = KDTreeInterpolator(x[ind], y[ind])
    return inter


#lazy functions

def toFlux(mag):
	return 10 ** (-0.4 * mag)

def toFlux2(mag,band):
	distanceModulus = 24.47
	filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W', 'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
	if band == 'f275w':
		ii  = 0
	elif band == 'f336w':
		ii = 1
	elif band == 'f475w':
		ii = 2
	elif band == 'f814w':
		ii = 3
	elif band == 'f110w':
		ii = 4
	elif band == 'f160w':
		ii = 5
	else :
		print 'band should be between 0 and 5'
	with Vega() as v:
		aa,bb,cc = v.getFlux(filters)
	return bb[ii] * 10 ** ((distanceModulus - mag)/2.5)

def toMag(flux):
    return -2.5 * np.log10(flux)

def compute_std(d,camera,k=10,eps=0):
    """
    Computes standard deviation and store the result in a dictionary
    INPUTS:
           d: dictionary
              AST data obtained with get_AST_dict()
           camera: string (uv, opt or ir)
              camera name to handle more easily filter names
           k: Integer
              Number of nearest neighbors taken in the standard deviation computation
           eps: non-negative float
    OUTPUT:
           d: dictionary
    """
    ii1,=np.where(d['MAG1OUT'] < 80) # removing unrecovered points
    ii2,=np.where(d['MAG2OUT'] < 80) # removing unrecovered points
    mag1in = d['MAG1IN'][ii1]
    mag1out = d['MAG1OUT'][ii1]
    mag2in = d['MAG2IN'][ii2]
    mag2out = d['MAG2OUT'][ii2]    
    
    mag1in_prep = _prepare_x(mag1in)
    mag2in_prep = _prepare_x(mag2in)
    NN1_mag = nearest_neighbors(mag1in_prep,k=k,eps=eps)
    NN2_mag = nearest_neighbors(mag2in_prep,k=k,eps=eps)
    std1_mag = np.zeros(len(mag1in))
    std2_mag = np.zeros(len(mag2in))
    for i in range(len(mag2in)):
        std1_mag[i]= np.sqrt((1./k)*((mag1out[NN1_mag[i]]-mag1in[i])**2).sum())
        std2_mag[i]= np.sqrt((1./k)*((mag2out[NN2_mag[i]]-mag2in[i])**2).sum())

    if camera == 'uv':
	    band1 = 'f275w'
	    band2 = 'f336w'
    elif camera == 'opt':
	    band1 = 'f475w'
	    band2 = 'f814w'
    elif camera == 'ir':
	    band1 = 'f110w'
	    band2 = 'f160w'
    f1in = toFlux2(mag1in, band1)
    f2in = toFlux2(mag2in, band2)
    f1out = toFlux2(mag1out, band1)
    f2out = toFlux2(mag2out, band2)

    mag1inup = mag1in + std1_mag
    mag1indown = mag1in - std1_mag
    mag2inup = mag2in + std2_mag
    mag2indown = mag2in - std2_mag

    f1up = toFlux2(mag1indown, band1)
    f1down = toFlux2(mag1inup, band1)
    std1_flux = (f1up-f1down)/2.
    f2up = toFlux2(mag2indown, band2)
    f2down = toFlux2(mag2inup, band2)
    std2_flux = (f2up-f2down)/2.
    
    B1F = f1out - f1in
    B2F = f2out - f2in

    D = {'STD1_MAG': std1_mag, 'STD2_MAG': std2_mag, 'STD1_FLUX':std1_flux, 'STD2_FLUX':std2_flux, 'BIAS1_FLUX': B1F, 'BIAS2_FLUX': B2F, 'BIAS1_MAG':d['BIAS1'][ii1], 'BIAS2_MAG':d['BIAS2'][ii2], 'MAG1IN':mag1in, 'MAG2IN':mag2in, 'MAG1OUT':mag1out, 'MAG2OUT':mag2out ,'RA1':d['RA'][ii1], 'F1IN':f1in, 'F1OUT':f1out,'F2IN':f2in, 'F2OUT':f2out, 'DEC1':d['DEC'][ii1], 'RA2':d['RA'][ii2], 'DEC2':d['DEC'][ii2] }
    
    return D

def extract_bias_std(brick,subdivision):
    """
    computes bias and noise in a given region from the ASTs
    and write the result in an output file.
    INPUTS:
           brick: string
           subdivision: string
    OUTPUT:
           D: dictionary
              contains AST results over a region -- written in a file
    """

    if brick > 10: 
	    direct = '/astro/dust_kg/harab/beast_data/b%s_bis/' % brick
	    file_vis = direct + 'fake_stars_b%s_%s_opt.fits' % (brick,subdivision)
	    file_ir = direct + 'fake_stars_b%s_%s_ir.fits' % (brick,subdivision)
	    file_uv = direct + 'fake_stars_b%s_%s_uv.fits' % (brick,subdivision)
    else:
	    direct = '/astro/dust_kg/harab/beast_data/b0%s_bis/' % brick
	    file_vis = direct + 'fake_stars_b0%s_%s_opt.fits' % (brick,subdivision)
	    file_ir = direct + 'fake_stars_b0%s_%s_ir.fits' % (brick,subdivision)
	    file_uv = direct + 'fake_stars_b0%s_%s_uv.fits' % (brick,subdivision)

    files = [file_uv,file_vis,file_ir]
    d1,d2,d3 = get_AST_dict(files)

    d1 = compute_std(d1,'uv',k=10,eps=0)
    d2 = compute_std(d2,'opt',k=10,eps=0)
    d3 = compute_std(d3,'ir',k=10,eps=0)
    t1 = Table(d1)
    t2 = Table(d2)
    t3 = Table(d3)

    D = {'uv':Table(d1), 'opt':Table(d2), 'ir':Table(d3)}

    outdir = 'PHAT_camera_AST/'
    outfile1 =  outdir + 'bias_std_uv_b%s_reg%s.dat' % (brick,subdivision)
    outfile2 = outdir +'bias_std_opt_b%s_reg%s.dat' % (brick,subdivision)
    outfile3 = outdir +'bias_std_ir_b%s_reg%s.dat' % (brick,subdivision)
    t1.write(outfile1)
    t2.write(outfile2)
    t3.write(outfile3)
    return D


def ast_interp(D,sedgrid):
    """
    Interpolate the results of the ASTs on a model grid
    INPUTS:
           D: dictionary
             output of extract_bias_std
           sedgrid: beast.core.grid type
             model grid to interpolate AST results on
    OUTPUT:
            (bias,sigma): tuple
             Two 1d numpy arrays
    """

    N = len(sedgrid.grid)
    filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W', 'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
    M = len(filters)
    bias = np.empty((N,M))
    sigma = np.empty((N,M))
    flux = sedgrid.seds
    with Vega() as v:
        names,vega_flux,lamb = v.getFlux(filters)
        
    for i in range(M):
        if i in [0,1]:
            key = 'uv'
            if i==0:
                var1 = D[key]['F1IN']
                var2 = D[key]['BIAS1_FLUX']
                var3 = D[key]['STD1_FLUX']
            else:
                var1 = D[key]['F2IN']
                var2 = D[key]['BIAS2_FLUX']
                var3 = D[key]['STD2_FLUX']
        elif i in [2,3]:
            key = 'opt'
            if i==2:
                var1 = D[key]['F1IN']
                var2 = D[key]['BIAS1_FLUX']
                var3 = D[key]['STD1_FLUX']
            else:
                var1 = D[key]['F2IN']
                var2 = D[key]['BIAS2_FLUX']
                var3 = D[key]['STD2_FLUX']
        else:
            key = 'ir'
            key = 'opt'
            if i==4:
                var1 = D[key]['F1IN']
                var2 = D[key]['BIAS1_FLUX']
                var3 = D[key]['STD1_FLUX']
            else:
                var1 = D[key]['F2IN']
                var2 = D[key]['BIAS2_FLUX']
                var3 = D[key]['STD2_FLUX']

        bias[:,i] = np.interp(flux[:,i],np.sort(var1),var2[np.argsort(var1)])
        sigma[:,i] = np.interp(flux[:,i],np.sort(var1),var3[np.argsort(var1)])

    return (bias,sigma)


brick = '15'
subdivision ='15'

project = 'b%s_%s' % (brick,subdivision)
dir_project = '/astro/dust_kg/harab/beast/projects/prod/%s/' % project
sedgrid = FileSEDGrid( dir_project + '/' + project + '_seds.grid.hd5')  

abs_calib = [1.19,1.04,0.71,0.74,0.84,0.88]
abs_calib_flux = sedgrid.seds*abs_calib*0.01

D = extract_bias_std(brick,subdivision)
(bias,sig) = ast_interp(D,sedgrid)

noise = np.sqrt(abs_calib_flux**2+sig**2)

outname = dir_project + 'noise_model_b%s_%s.hd5' % (brick,subdivision)
with tables.openFile(outname, 'w') as outfile:
	outfile.createArray(outfile.root,'bias',bias)
	outfile.createArray(outfile.root,'error',noise)

