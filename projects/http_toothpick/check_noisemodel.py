""" Script to check the noise model by plotting the interpolated results on the model grid

.. history::
    Written 31Oct14 - KDG
    Updated for HTTP 7Jul15 - KDG
"""
import numpy as np
import tables
import matplotlib.pyplot as pyplot
import datamodel_small as datamodel

import noisemodel
from beast.core.grid import FileSEDGrid
from beast.external.eztables import Table

if __name__ == '__main__':

    obs = datamodel.HTTPFluxCatalog(datamodel.obsfile,distanceModulus=datamodel.distanceModulus,
                                    filters=datamodel.filters)

    g = FileSEDGrid('{project:s}/{project:s}_seds.grid.hd5'.format(project=datamodel.project))  

    nm = noisemodel.get_noisemodelcat(datamodel.noisefile)

    fig, ax = pyplot.subplots(7, sharex=True, figsize=(10,13))

    for i in range(7):
        ax[i].plot(g.seds[:,i],nm.root.bias[:,i]/g.seds[:,i],'bo')
        ax[i].plot(g.seds[:,i],nm.root.error[:,i]/g.seds[:,i],'go')
        ax[i].set_ylabel('Bias/Flux(Input)')
        ax[i].set_xscale('log')
        ax[i].set_yscale('linear')
        ax[i].set_ylim(-1.1,5.1)
        ax[i].text(1e-16,4.,g.filters[i])

    #ax[0].set_xlim(1e-25,1e-15)
    ax[5].set_xlabel('Flux')

    for k in range(len(obs)):
        fluxes = obs.getFlux(k)
        print(fluxes)
        for i in range(7):
            if fluxes[i] > 0.0:
                ax[i].plot(fluxes[i],[1.0],'ro')

    pyplot.show()
