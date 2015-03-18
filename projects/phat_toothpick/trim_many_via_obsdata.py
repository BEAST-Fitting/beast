#!/usr/bin/env python2.7
"""
Code to create many trimmed model grids for batch runs
  Saves time by only reading the huge modelsed grid once
"""

# system imports
from __future__ import print_function
import sys
import os
import argparse
import time
import string

# BEAST imports
import datamodel_production as datamodel
import noisemodel 
from beast.core import trim_grid
from beast.core.grid import FileSEDGrid  

if __name__ == '__main__':
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("trimfile", help="file with modelgrid, astfiles, obsfiles to use")
    args = parser.parse_args()

    start_time = time.clock()

    f = open(args.trimfile, 'r')
    file_lines = list(f)

    modelfile = file_lines[0].rstrip()
    print('Reading the model grid files = ', modelfile)

    # get the modesedgrid on which to generate the noisemodel  
    modelsedgrid = FileSEDGrid(modelfile)  

    new_time = time.clock()
    print('time to read: ',(new_time - start_time)/60., ' min')
    
    old_noisefile = ''
    for k in range(1,len(file_lines)):
        line = file_lines[k]
        s1pos = string.find(line,' ')
        s2pos = string.rfind(line,' ')

        brick = line[0:s1pos].rstrip()
        source_density = line[s1pos+1:s2pos].rstrip()
        sub_source_density = line[s2pos+1:len(line)].rstrip()
    
        obsfile = 'BEAST_production/b'+brick+'/obscat/b' + brick + \
                  '-6filt-cut-4band-gst-bright-SD-' + string.replace(source_density,'_','-') + \
                  '-sub' + sub_source_density + '.fits'
        astfile = 'BEAST_production/merged_asts/PHAT_fake_stars_SD_' + \
                  string.replace(source_density,'-','_' ) + '.fits'
        noisefile = 'BEAST_production/BEAST_production_sd_' + \
                    string.replace(source_density,'_','-' ) + '_noisemodel.fits'
        
        stats_filebase = 'BEAST_production/b' + brick + '/b' + brick + \
                         '_sd' + string.replace(source_density,'_','-') + '_sub' + sub_source_density 
        sed_trimname = stats_filebase + '_sed_trim.grid.hd5'
        noisemodel_trimname = stats_filebase + '_noisemodel_trim.hd5'

        #print(k)
        #print(noisefile)
        #print(astfile)
        #print(obsfile)
        print('working on ' + sed_trimname)
        #print(noisemodel_trimname)
        
        start_time = time.clock()

        if noisefile == old_noisefile:
            print('not reading noisefile - same as last')
            #print(noisefile)
            #print(astfile)
        else:
            print('reading noisefile/astfile')
            # read in the noise model
            noisemodel_vals = noisemodel.get_noisemodelcat(noisefile)
            old_noisefile = noisefile

        # read in the observed data
        print('getting the observed data')
        obsdata = datamodel.get_obscat(obsfile, 24.0, modelsedgrid.filters)

        # trim the model sedgrid

        trim_grid.trim_models(modelsedgrid, noisemodel_vals, obsdata, sed_trimname, noisemodel_trimname, sigma_fac=3., n_detected=4)

        new_time = time.clock()
        print('time to trim: ',(new_time - start_time)/60., ' min')
