"""
Extract BEAST's best fits from lnp file
"""

## Imports
from beast.external.eztables.astro import AstroTable
from beast.core.observations import Observations
from beast.core.grid import FileSEDGrid
from beast.core.vega import Vega
from beast.tools import progressbar
from tables import *
from user_input_val import *
import numpy as np
import numexpr
import matplotlib.pyplot as plt

## Reading files and setting the SED grid
READ = True  #Set to TRUE for the first run
if READ == True:
    #path shortcut
    dir = "/astro/dust_kg/harab/beast/projects/real_data/"
    
    #BEAST grid
    BEASTgrid_fname = dir + "GridRvFbumpLaw.sed.grid.fits"  #tlusty and kurucz
    BEASTgrid = FileSEDGrid(BEASTgrid_fname)
    
    #Putting the grid from 10 pc to M31's distance
    distance = 0.75857757502918322 #from distanceModulus input
    #BEASTgrid.seds = BEASTgrid.seds * (1e-5/distance)**2

    #BEAST file  wanted to be analyzed - must matched the grid
    BEASTres = openFile(dir+"lnp.iso10.GridRvFbumpLaw_8004.hd5") #tlusty and kurucz
    
    # Best fit parameters will to written in this file      
    outname = 'test_BBF.iso10.GridRvFbumpLaw_8004.fits' 


#Getting Vega fluxes
with Vega() as v:
    vega_f, vega_flux, lamb = v.getMag(filters)


#Initialization
N_stars = 8004
#BEAST_chi2 = np.zeros(N_stars)
bestind = np.arange(N_stars)
bestmod = np.ndarray(shape=(N_stars,6), dtype=float)
expectSED = np.ndarray(shape=(N_stars,6), dtype=float)

#BEAST_Pmax = np.zeros(N_stars)

with progressbar.PBar(N_stars, txt="Extracting Best fits") as pbar:
    for tn in range(N_stars):
        star_node = BEASTres.getNode('/star_%d' % tn)
        idx = star_node.idx[:]
        lnp = star_node.lnp[:]
        #BEAST_chi2[tn] = chi2.min() # BEAST_chi2 contains the best chi2 for each star
        ind = np.where(lnp == lnp.max())
        bestind[tn] = idx[ind][0]
        #prob = numexpr.evaluate('exp(lnp)', local_dict={'lnp': np.float64(lnp)})
        #if prob.max()==0:
        #    expectSED[tn]=[0.,0.,0.,0.,0.,0.]
        #else:
        #    expectSED[tn] = np.average(BEASTgrid.seds[idx], axis=0, weights=prob)
        #norm = numexpr.evaluate('sum(prob)', local_dict={'prob': np.float64(lnp)})
        #BEAST_Pmax[tn] =  prob.max()/norm

        pbar.update(tn, force=True)


param = ['Av', 'Rv','f_bump', 'logg', 'logT', 'logM', 'logA', 'logL', 'Z']
d ={}
for e, k in enumerate(param):
    d[k] = BEASTgrid.grid[k][bestind]
    
#d['chi2'] = BEAST_chi2 #Computed assuming a constant 10% photometric error added to the catalog errors.
#d['Pmax'] = BEAST_Pmax #No normalization just max(e^[-chi2/2])
d['Bestind'] = bestind #Index of the best fit in the whole SED grid frame

tab = AstroTable(d)
tab.write(outname, clobber=True, append=False)

