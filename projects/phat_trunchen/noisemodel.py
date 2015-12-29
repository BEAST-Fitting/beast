""" Script that generates the noise model for phat/m31 data,
based on the BEAST noise model package

This is the TRUNCHEN noise model version: the ASTs are assumed to be run
simultaneously allowing for the covariance to be compute.

The idea is to use artificial star tests (ASTs) to characterize the noise
introduced by crowding and selection function.

.. history::
   Started 18 Dec 2015 by Karl Gordon
"""

import numpy as np
import tables

from beast.core.noisemodel import trunchen

# define the object to hold the ASTs and processed versions
class PHAT_Trunchen_Noisemodel(trunchen.MultiFilterASTs):
    pass

def make_trunchen_noise_model(outname, astfile, basefilters, sedgrid):
    """ trunchen noise model with full covariance information

    Parameters
    ----------
    outname: str
        filename where the noise model will be saved

    astfile: str
        filename containing the FITS file with the AST results

    sedgrid: SEDGrid instance
        sed model grid needing a noise model

    returns
    -------
    noisefile: str
        noisemodel file name
    """
    # setup the noisemodel object
    #  including reading in the AST information
    model = PHAT_Trunchen_Noisemodel(astfile, sedgrid.filters)
    # compute the biases and covariance matrixies for all the
    # independent models in the AST file
    model.process_asts(basefilters)
    # evaluate the noise model for all the models in sedgrid
    results = model(sedgrid)

    # unpack the results
    bias = results[0]
    sigma = results[1]
    compl = results[2]
    q_norm = results[3]
    icov_diag = results[4]
    icov_offdiag = results[5]
    exit()
    
    # check if the noise model has been extrapolated at the faint flux levels
    # if so, then set the noise to a negative value (later may be used to
    # trim the model of "invalid" models)
    # we are assuming that extrapolation at high fluxes is ok as the noise
    # will be very small there
    for k in range(len(model.filters)):
        indxs, = np.where(sedgrid.seds[:,k] <= model._minmax_asts[0,k])
        if len(indxs) > 0:
            noise[indxs,k] *= -1.0

    print('Writting to disk into {0:s}'.format(outname))
    with tables.openFile(outname, 'w') as outfile:
        outfile.createArray(outfile.root,'bias', bias)
        outfile.createArray(outfile.root,'error', noise)
        outfile.createArray(outfile.root,'completeness', compl)

    return outname

if __name__ == '__main__':
    from beast.core.grid import FileSEDGrid

    # get the filters
    project = 'b15_nov15_test'
    
    modelseds = FileSEDGrid(project+'/'+project+'_seds.grid.hd5')

    uvastfile = 'data/fake_stars_b15_27_uv.fits'
    optastfile = 'data/fake_stars_b15_27_opt.fits'
    irastfile = 'data/fake_stars_b15_27_ir.fits'
    astfile = 'data/fake_stars_b15_27_all.hd5'
    Merge_PHAT_ASTs(uvastfile,optastfile,irastfile,astfile)

    noisefile = project + '/b15_27_noisemodel.hd5'
    make_toothpick_noise_model(noisemodel, astfile, modelseds)
