"""
Priors as weights
================

The priors on age and mass are computed as weights to be used in the
likelihood computation.  This code was created by Heddy Arab and
integrated (likely badly) into the beast core by Karl Gordon.
"""
import numpy as np
from scipy.integrate import quad
import numexpr

from .grid import FileSEDGrid
from .grid import SpectralGrid
from ..external.eztables import Table

def add_priors_sedgrid(sedgrid_file, outname, filter_names):

    def interv(tab):
        """
        Computes the bin size in a grid
        """
        temp = tab[1:]-np.diff(tab)/2.
        tab2 = np.empty(len(tab)+1)
        tab2[0] = tab[0]-np.diff(tab)[0]/2.
        tab2[-1] = tab[-1]+np.diff(tab)[-1]/2.
        tab2[1:-1] = temp
        return tab2

    def integ (x): #Kroupa IMF
        m0 = 0.01
        m1 = 0.08
        m2 = 0.5
        alpha0 = -0.3
        alpha1 = -1.3
        alpha2 = -2.3
        if (x < m1):
            return x**alpha0
        elif (x >= m1) and (x < m2):
            return x**alpha1
        elif (x>=m2):
            return x**alpha2

    def integ_Salpeter(x):
        return x**(-2.35) # Salpeter IMF

    def stellar_extract(cols):
        """
        Extracts grid models for a given set of dust parameters
        """
        ind,=np.where((cols['Av'] == 0.) & (cols['Rv'] == 2.0) & (cols['f_bump']==1.))
        return cols[ind]
    
    def Z_extract(cols, z_val):
        """
        Extracts grid models for a constant metallicity
        """
        ind,=np.where(cols['Z'] == z_val)
        return cols[ind]
    
    def iso_extract(cols, iso):
        """
        Extract models of constant age
        """
        ind,=np.where(cols['logA'] == iso)
        return cols[ind]
    
    def weights(cols,z_val,iso_val,d_dages):
        """
        Computes weights to assign to each model in the full
        grid according to a given age-mass prior (in this case Kroupa IMF)
        """
        cols2 = stellar_extract(cols)                # extract stellar grid
        cols_z = Z_extract(cols2, z_val)             # extract models of cste metallicity
        cols2 = iso_extract(cols_z,iso_val)          # extract isochrone
        
        ind,=np.where(d_dages['logages']==iso_val)   # extract models for the given log(age)
        linage_binw = d_dages['weight'][ind]         # store the weights assign along this isochrone (linear age)
        
        d = {}
        
        isoc = np.sort(cols2['M_ini'])               # Get and sort the initial mass along this isochrone
        index_isoc = np.argsort(cols2['M_ini'])      
        isoc2 = interv(isoc)                         # Compute the initial mass bin width 
        res1 = quad(integ, isoc.min(), isoc.max())   # integrate according to the desired prior (integ function) along the isochrone
        denom = res1[0]
        I1 = np.empty(len(isoc))
        res = np.empty(len(isoc))
        for ik, uk in enumerate(isoc2[:-1]):
            res = quad(integ, isoc2[ik], isoc2[ik+1]) # integrate according to the prior on the mass bin
            I1[ik] = res[0]/denom*linage_binw         # Compute the final weight
        d['weight'] = I1                              # Store the result in a dict
        #d['index'] = index_isoc
        d['M_ini'] = isoc
        return d

    def weights_age(logages):
        logages2 = interv(logages)                  # Computes the bin size in log
        d_dages = {}
        d_dages['weight'] = np.diff(10**(logages2)) # Store the weight in linear 
        d_dages['logages'] = logages                
        return d_dages                              # Returns the dictionary with log(age) and corresponding weight
    
    def filled(D,grid):
        """
        Function to fill the full sedgrid with the computed weights
        """
        prior = np.empty(len(grid))
        for ek,a in enumerate(np.unique(grid['Z'])):
            for ik,b in enumerate(np.unique(grid['logA'])):
                ind,=np.where((grid['Z'] == a) & (grid['logA']==b))
                M_ini = grid['M_ini'][ind]
                for uk,c in enumerate(np.unique(M_ini)):
                    ind2,=np.where(M_ini == c) 
                    prior[ind[ind2]] = D[str(a)+','+str(b)]['weight'][uk]
                    
        return prior

    def prior_dict(cols):
        D = {}
        Zs = np.unique(cols['Z'])  # Values of metallicity
        for i in Zs:
            print('working on Z = ',i)
            cols2 = stellar_extract(sedgrid)       # Extract stellar grid from full grid
            cols = Z_extract(cols2,i)              # Extarct models of cste Z for each Z
            logages = np.unique(cols['logA'])      # Store the log(age) values
            d_dages = weights_age(logages)         # Computes the weights along an isochrone
            for j in logages:
                D[str(i)+','+str(j)] = weights(sedgrid,i,j,d_dages)  # Save the weights in a dict for each metallicity and age value
        prior = filled(D, sedgrid.grid)            # Fill the full sedgrid with the weights
        return prior

    ### work done here

    sedgrid = FileSEDGrid(sedgrid_file)

    prior = prior_dict(sedgrid)

    if 'weight' not in sedgrid.grid.keys():
        cols ={}
        for key in sedgrid.grid.keys():
            cols[key] = sedgrid.grid[key]
        cols['weight'] = prior
        _lamb = sedgrid.lamb
        _seds = sedgrid.seds
        
        g = SpectralGrid(_lamb, seds=_seds, grid = Table(cols), backend='memory')
        g.grid.header['filters'] = ' '.join(filter_names)
            
        if hasattr(g, 'writeHDF'):     #write to disk
            g.writeHDF(outname)
        else:
            for gk in g:
                gk.writeHDF(outname, append=True)
            
