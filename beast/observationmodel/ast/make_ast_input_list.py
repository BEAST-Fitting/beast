import sys, argparse
import numpy as np
import astropy.io.fits as pyfits
import string, tables
from matplotlib.colors import LogNorm
from vega import Vega
import datamodel


def flux_limits(seds,limits,Nfilter=1):
    """
    Selects models which have at least one filter above the limits

    INPUTS:
    -------
    seds:   np.array
            Magnitude array from BEAST grid 

    limits: list
            List of limit magnitudes 

    Nfilter: integer
             In how many filters, you want a fake star to be brighter than the limit
 
    OUTPUT:
    -------
    idx:    np.array
            Array of integers contining the indices of allowed models

    """
    flag = seds.copy() 

    # flag is True if the models are brigter (=smaller number in mag) than the limits
    for i,limit in enumerate(limits):
        flag[:,i] = seds[:,i] < limit 
    
    # Keep index where model is brighter than the limit in N filters
    s = np.sum(flag,axis=1)
    idx, = np.where(s >= Nfilter)            
        
    return idx


def pick_models(sedgrid, mag_cuts, Nfilter=1, N_stars= 625):
    """
    Creates a fake star catalog from a BEAST model grid and a list of gst cuts
    Uses the flux_limit function (see above) for PHAT ASTs (can be replaced with a simple where statement 
    to select the model that pass the gst cuts)
    Pick 4 Av values for each stellar models
    Writes the final catalog in a fits binTable

    INPUTS:
    -------
    sedgrid: beast.grid
               BEAST model grid in which the models are picked

    mag_cuts: list
               List of mag limits considered in the pick

    Nfilter: Integer
             In how many filters, you want a fake star to be brighter than the limit

    N_stars: Integer
               Number of stellar models picked (default=625)

    outfile: String
               Name of the output fake star catalog
               
    return:
    -------
    list: indices for the selected models
    """

    filters = [b'HST_WFC3_F225W', b'HST_WFC3_F336W', b'HST_WFC3_F438W',\
               b'HST_WFC3_F814W', b'HST_WFC3_F110W', b'HST_WFC3_F160W']
    #filters = sedgrid.filters       # Filters in the grid
    with Vega('vega.hd5') as v:               # Get the vega fluxes
        vega_f,vega_flux,lamb = v.getFlux(filters)

    sedsMags = -2.5 * np.log10(sedgrid.seds/vega_flux)  # Convert to Vega mags  

    # Select the models above the flux limits in N filters
    idx = flux_limits(sedsMags, mag_cuts-2.5*np.log10(0.2), Nfilter=Nfilter) 
    grid_cut = sedgrid.grid[idx] 
    mod_tot = len(idx)

    # Sample the model grid uniformly
    prime_params = np.column_stack((grid_cut['logA'],grid_cut['M_ini'],grid_cut['Av']))
    uniq_age = np.unique(prime_params[:,0])
    N_uniq_age = len(uniq_age)
    if N_uniq_age % 2 == 0: 
        search_age = np.append(uniq_age[np.arange(0,N_uniq_age,2)],uniq_age[N_uniq_age-1])
    else:
        search_age = uniq_age[np.arange(0,N_uniq_age+1,2)]

    N_sample = int(np.ceil(float(N_stars)/len(search_age)))
    models = []
    for iage in search_age:
        tmp, = np.where(prime_params[:,0] == iage)
        models.append(np.random.choice(tmp,N_sample))

         
    index = idx[np.array(models).reshape((-1))]       
    sedMags = sedsMags[index,:]


    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--ast", help="Generate an input AST file",
                        action="store_true")
    parser.add_argument("-s", "--split", help="Split real and fake stars into each SD bin",
                        action="store_true")
    parser.add_argument("-l", "--limits", dest='mag_cuts', default=[1.], help="Set custom magnitude limits",
                        action="store_true")
    args = parser.parse_args()

    if args.ast:
        sedFile = 'choi_ngc4214_seds.grid.hd5' # Input BEAST sed.grid file    
       
  
        mag_cuts = mag_cuts.astype('float') 
        if len(mag_cuts) == 1:
            tmp_cuts = mag_cuts
            obsdata = datamodel.get_obscat(datamodel.obsfile,
                                           datamodel.distanceModulus,
                                           datamodel.filters)

            min_mags = np.zeros(len(datamodel.filters))
            for k, filtername in enumerate(obsdata.filters):
                sfiltname = obsdata.data.resolve_alias(filtername)
                keep, = np.where(obsdata.data[:][sfiltername] < 99.)
                min_mags[k] = np.amax(obsdata.data[keep][sfiltername])
            
            mag_cuts = min_mags + tmp_cuts # max. mags from the gst observation cat. 
   
        N_models = 630
        Nfilter = 4

        sed = tables.open_file(inFile)
        sedgrid = sed.root
        index, sedsMags = pick_models(sedgrid, mag_cuts, Nfilter=Nfilter, N_stars= N_models)
        grid = sedgrid.grid[index]
        N_models_updated = len(index)
        Nast_per_model = 20
        astMags = np.repeat(sedsMags, Nast_per_model, 0)

        # Assign random x,y positions for ASTs based on source density 
        sd = pyfits.getdata('sd_20stars.fits')
        ref = pyfits.getdata('11360_NGC-4214_F438W_drz.chip1.fits')
        mask = ref>-7*1e7 
        sd = sd*mask
        levels = [sd.min(),10**(-1),10**(-0.5),10**(0),10**(0.5),10**(0.75),10**(1),10**(1.25),sd.max()]
        N_sd = len(levels) - 1
   
        N_ast = Nast_per_model*N_models_updated*N_sd
        N_ast_perSD = int(N_ast/N_sd)
        #x, y = np.zeros(N_ast), np.zeros(N_ast)
        add_noiseX = np.random.choice(np.arange(0.,1.,1e-7),N_ast)
        add_noiseY = np.random.choice(np.arange(0.,1.,1e-7),N_ast)

        outf = open('ngc4214_ASTs_trunchen_test.txt','a')
        outf2 = open('ngc4214_ASTs_trunchen_withSDbinNumber_test.txt','a')
        for i in range(len(levels)-1):
            good = (levels[i] < sd) & (sd <= levels[i+1])
            tmpy, tmpx = np.where(good==1)[0], np.where(good==1)[1]  
            add = np.random.choice(np.arange(len(tmpx)),N_ast_perSD)
            x = tmpx[add].astype('float64') + np.random.choice(add_noiseX,N_ast_perSD)
            y = tmpy[add].astype('float64') + np.random.choice(add_noiseY,N_ast_perSD)
      
            for j in range(N_ast_perSD):
                outf.write('0 1 %s %s %s %s %s %s %s %s\n' % (x[j],y[j],
                         astMags[j,0],astMags[j,1],astMags[j,2],astMags[j,3],astMags[j,4],astMags[j,5]))
                outf2.write('0 1 %s %s %s %s %s %s %s %s %d\n' % (x[j],y[j],
                         astMags[j,0],astMags[j,1],astMags[j,2],astMags[j,3],astMags[j,4],astMags[j,5],
                         i))

        outf.close()
        outf2.close()
   

    elif args.split:
        
        inReal = 'N4214_4band_detects_newPhotometry.fits'
        inFake = 'N4214_gst_fake_onlyNEW.fits'

        real, rhd = pyfits.getdata(inReal,1,header=True)
        realX = real.X; realY = real.Y
        fake, fhd = pyfits.getdata(inFake,1,header=True)
        fakeX = fake.XIN; fakeY = fake.YIN         
        sd = pyfits.getdata('sd_20stars.fits')
        ref = pyfits.getdata('11360_NGC-4214_F438W_drz.chip1.fits')
        mask = ref>-7*1e7 
        sd = sd*mask

        realX_ind = np.floor(realX).astype('int')
        realY_ind = np.floor(realY).astype('int')
        realSD = sd[realY_ind,realX_ind] # 2D array

        fakeX_ind = np.floor(fakeX).astype('int')
        fakeY_ind = np.floor(fakeY).astype('int')
        fakeSD = sd[fakeY_ind,fakeX_ind] # 2D array
        
        levels = [sd.min(),10**(-1),10**(-0.5),10**(0),10**(0.5),10**(0.75),10**(1),10**(1.25),sd.max()]
        N_sd = len(levels)-1

        for i in range(N_sd):
            realadd, = np.where((levels[i] < realSD) & (realSD <= levels[i+1]))
            fakeadd, = np.where((levels[i] < fakeSD) & (fakeSD <= levels[i+1]))

            pyfits.writeto('N4214_4band_SD%s.fits' % (i),real[realadd],header=rhd,clobber='True')
            pyfits.writeto('N4214_fake_SD%s.fits' % (i),fake[fakeadd],header=fhd,clobber='True')




    """
    # Create the pyfits columns
    newcols = pyfits.ColDefs([
              pyfits.Column(name='logA', format='D', array=grid['logA']),
              pyfits.Column(name='M_ini', format='D', array=grid['M_ini']),
              pyfits.Column(name='M_act', format='D', array=grid['M_act']),
              pyfits.Column(name='logT', format='D', array=grid['logT']),
              pyfits.Column(name='logL', format='D', array=grid['logL']),
              pyfits.Column(name='logg', format='D', array=grid['logg']),
              pyfits.Column(name='Av', format='D', array=grid['Av']),
              pyfits.Column(name='f225w_vega',format='D',array=sedsMags[:,0]),
              pyfits.Column(name='f336w_vega',format='D',array=sedsMags[:,1]),
              pyfits.Column(name='f438w_vega',format='D',array=sedsMags[:,2]),
              pyfits.Column(name='f814w_vega',format='D',array=sedsMags[:,3]),
              pyfits.Column(name='f110w_vega',format='D',array=sedsMags[:,4]),
              pyfits.Column(name='f160w_vega',format='D',array=sedsMags[:,5])])

    # Create the pyfits bin table
    tbhdu = pyfits.BinTableHDU.from_columns(newcols)

    # Write the output catalog
    tbhdu.writeto('ast_test.fits', clobber=True)
    """
