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

def make_trunchen_noise_model(outname, astfile, basefilters, sedgrid,
                              generic_absflux_a_matrix=None):
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
    results = model(sedgrid, generic_absflux_a_matrix=generic_absflux_a_matrix)

    # unpack the results
    bias = results[0]
    sigma = results[1]
    compl = results[2]
    q_norm = results[3]
    icov_diag = results[4]
    icov_offdiag = results[5]
    
    # check if the noise model has been extrapolated at the faint flux levels
    # if so, then set the noise to a negative value (later may be used to
    # trim the model of "invalid" models)
    # we are assuming that extrapolation at high fluxes is ok as the noise
    # will be very small there
    for k in range(len(model.filters)):
        indxs, = np.where(sedgrid.seds[:,k] <= model._minmax_asts[0,k])
        if len(indxs) > 0:
            sigma[indxs,k] *= -1.0

    print('Writting to disk into {0:s}'.format(outname))
    with tables.openFile(outname, 'w') as outfile:
        outfile.createArray(outfile.root,'bias', bias)
        outfile.createArray(outfile.root,'error', sigma)
        outfile.createArray(outfile.root,'completeness', compl)
        outfile.createArray(outfile.root,'q_norm', q_norm)
        outfile.createArray(outfile.root,'icov_diag', icov_diag)
        outfile.createArray(outfile.root,'icov_offdiag', icov_offdiag)

    return outname

def get_noisemodelcat(filename):
    """
    returns the noise model

    Parameters
    ----------
    filename: str
        file containing the outputs from OneD_ASTs_ModelGenerator

    Returns
    -------
    table: pytables.Table
        table containing the elements of the noise model
    """
    return tables.openFile(filename)

if __name__ == '__main__':

    pass
