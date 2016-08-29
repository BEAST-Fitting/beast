""" Script to check the noise model by plotting the interpolated results on the model grid

.. history::
    Written 31Oct14 - KDG
"""
import numpy as np
import tables
import matplotlib.pyplot as pyplot

import noisemodel
from beast.core.grid import FileSEDGrid
from beast.external.eztables import Table

if __name__ == '__main__':

    g = FileSEDGrid('./mf_ngc4214_Aug2014_new/mf_ngc4214_Aug2014_new_seds.grid.hd5')

    nm = noisemodel.get_noisemodelcat('ngc4214_Oct2014/ngc4214_Oct2014_noisemodel.hd5')

    fig, ax = pyplot.subplots(6, sharex=True, figsize=(10,13))

    for i in range(6):
        ax[i].plot(g.seds[:,i],nm.root.bias[:,i]/g.seds[:,i],'bo')
        ax[i].set_ylabel('Bias/Flux(Input)')
        ax[i].set_xscale('log')
        ax[i].set_yscale('linear')
        ax[i].set_ylim(-1.1,5.1)
        ax[i].text(1e-16,4.,g.filters[i])

    ax[0].set_xlim(1e-20,1e-14)
    ax[5].set_xlabel('Flux')

    pyplot.show()
