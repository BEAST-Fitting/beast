"""
Priors as weights
================

The priors on age and mass are computed as weights to be used in the
likelihood computation.  This code was created by Heddy Arab and
integrated into the beast core by Karl Gordon.
"""
from __future__ import print_function

import numpy as np
from scipy.integrate import quad
import numexpr

from .grid import FileSEDGrid
from .grid import SpectralGrid
from ..external.eztables import Table

# compute the bin edges
# approximate the edge bins by adding 1/2 the adjacent bin width
def compute_bin_boundaries(tab):
    """
    Computes the bin boundaries
    """
    temp = tab[1:]-np.diff(tab)/2.
    tab2 = np.empty(len(tab)+1)
    tab2[0] = tab[0]-np.diff(tab)[0]/2.
    tab2[-1] = tab[-1]+np.diff(tab)[-1]/2.
    tab2[1:-1] = temp
    return tab2

# Kroupa IMF
def imf_kroupa(x):
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
    
# Salpeter IMF
def imf_salpeter(x):
    return x**(-2.35) # Salpeter IMF

# compute the age weights for a constant SFR in linear age
def compute_age_weights(logages):
    aindxs = np.argsort(logages)   # ages need to be monotonically increasing
    logages2 = compute_bin_boundaries(logages[aindxs])    # Computes the bin boundaries in log
    age_weights = np.full(len(aindxs),0.0)
    age_weights[aindxs] = np.diff(10**(logages2))           # Returns the age weight as a numpy array
    return age_weights    # return in the order that logages was passed

# compute the mass weights at a constant age
# uses an assumed IMF to generate the weights
#  IMF norm is the integral of the assumed IMF over the full possible mass range
#  must be precomputed as this information is not usual available for a specific isochrone age
def compute_mass_weights(masses, full_imf_integral):
    d = np.zeros(len(masses))
        
    isoc = np.sort(masses)               # sort the initial mass along this isochrone
    index_isoc = np.argsort(masses)
    
    isoc2 = compute_bin_boundaries(isoc)      # Compute the initial mass bin boundaries

    I1 = np.empty(len(isoc))
    for ik, uk in enumerate(isoc2[:-1]):
        res = quad(imf_kroupa, isoc2[ik], isoc2[ik+1]) # integrate according to the prior on the mass bin
        I1[index_isoc[ik]] = res[0]/full_imf_integral      # compute the integrated weight
    return I1

# compute age-mass-metallicity prior weights
# age prior is a constant SFR in linear age
# mass prior is a Kroupa IMF (need to update code to allow user to pick the function)
# metallicity prior is flat
def compute_age_mass_metallicity_prior_weights(_tgrid):

    uniq_Zs = np.unique(_tgrid['Z'])  # get the unique metallicities
    total_z_weight = np.zeros(len(uniq_Zs))
    for az, z_val in enumerate(uniq_Zs):
        print('working computing the age-mass prior for Z = ', z_val)
        
        zindxs, = np.where(_tgrid['Z'] == z_val)   # get the grid for a single metallicity
        uniq_ages = np.unique(_tgrid[zindxs]['logA']) # get the unique ages for this metallicity
        age_weights = compute_age_weights(uniq_ages)  # compute the age weights for a constant SFR in linear age

        # get the integral over the IMF
        # assuming the min/max masses at any age for a single metallicity are the assumed min/max masses
        #min_mass = np.min(_tgrid[zindxs]['M_ini'])
        #max_mass = np.max(_tgrid[zindxs]['M_ini'])
        #imf_norm = quad(imf_kroupa, min_mass, max_mass)   # integrate according to the desired IMF along the isochrone
        imf_norm = [1.0]  # effectively ignore this term - assumes same mass range for all metallicities

        for ak, age_val in enumerate(uniq_ages):
            aindxs, = np.where((_tgrid['logA'] == age_val) & (_tgrid['Z'] == z_val))   # get the grid for a single age
            _tgrid_single_age = _tgrid[aindxs]
            if len(aindxs) > 1:
                mass_weights = compute_mass_weights(_tgrid_single_age['M_ini'],imf_norm[0])
            else:
                # must be a single mass for this age,z combination
                # set mass weight to zero to remove this point from the grid
                mass_weights = np.zeros(1)

            for i, k in enumerate(aindxs):
                _tgrid[k]['weight'] *= mass_weights[i]*age_weights[ak]

        total_z_weight[az] = np.sum(_tgrid[zindxs]['weight'])
        
    # ensure that the metallicity prior is uniform
    if len(uniq_Zs > 1):
        z_boundaries = compute_bin_boundaries(uniq_Zs)
        z_widths = np.diff(z_boundaries)

        total_z_weight *= z_widths   # very simple integration
        z_weights = total_z_weight/np.sum(total_z_weight)
        for az, z_val in enumerate(uniq_Zs):
            zindxs, = np.where(_tgrid['Z'] == z_val)   # get the grid for a single metallicity
            _tgrid[zindxs]['weight'] *= z_weights[az]

# previous version of the code that adds the weights to the sedgrid
# instead of the spectralgrid in the stelllib.py code
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
        #print(d_dages['logages'][ind])
        #print(linage_binw)
        
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
        
        #print(isoc)
        #print(I1)

        #exit()
                  
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

    if 'weight' in sedgrid.grid.keys():
        sedgrid.grid['weight'] = prior
        g = sedgrid
    elif 'weight' not in sedgrid.grid.keys():
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
            
