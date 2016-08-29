""" Script that generates the noise model for phat/m31 data,
based on the BEAST noise model package

Currently TOOTHPICK noise model version: assumes independent filters
The idea is to use artificial star tests (ASTs) to characterize the noise
introduced by crowding and selection function.

.. history::
    Based on code from the ngc4214 products directory (18oct14, KDG)
       adapted for the PHAT AST files
       updated to sync up with HA's work on the PHAT ASTs
"""

import numpy as np
import tables

from beast.core.noisemodel import toothpick
from beast.external.ezpipe.helpers import RequiredFile, task_decorator

class PHAT_ToothPick_Noisemodel(toothpick.MultiFilterASTs):

    def set_data_mappings(self):
        """
        hard code mapping directly with the interface to ASTs format
        """
        for k in self.filters:
            try:
                self.data.set_alias(k + '_out',
                                    k.split('_')[-1].upper() + '_VEGA')
                self.data.set_alias(k + '_in',
                                    k.split('_')[-1].upper() + '_IN')
            except Exception as e:
                print(e)
                print('Warning: Mapping failed. This could lead to ' +
                      'wrong results')


def make_toothpick_noise_model(outname, astfile, sedgrid,
                               absflux_a_matrix=None, **kwargs):
    """ toothpick noise model assumes that every filter is independent with
    any other.

    Parameters
    ----------
    outname: str
        path and filename into which save the noise model

    astfile: str
        path to the file into which are ASTs results

    sedgrid: SEDGrid instance
        sed model grid for everyone of which we will evaluate the model

    absflux_a_matrix: ndarray
        absolute calibration a matrix giving the fractional uncertainties
        including correlated terms (off diagonals)

    returns
    -------
    noisefile: str
        noisemodel file name
    """

    # read mag_in, mag_out
    model = PHAT_ToothPick_Noisemodel(astfile, sedgrid.filters)
    # compute binned biases and uncertainties as a function of flux
    model.fit_bins(nbins=30, completeness_mag_cut=80)
    # evaluate the noise model for all the models in sedgrid
    bias, sigma, compl = model(sedgrid)

    # absolute flux calibration uncertainties
    #  currently we are ignoring the off-diagnonal terms
    if absflux_a_matrix is not None:
        if absflux_a_matrix.ndim == 1:
            abs_calib_2 = absflux_a_matrix[:] ** 2
        else:   # assumes a cov matrix
            abs_calib_2 = np.diag(absflux_a_matrix)

        noise = np.sqrt(abs_calib_2 * sedgrid.seds[:] ** 2 + sigma ** 2)
    else:
        noise = sigma

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

@task_decorator()
def t_gen_noise_model(project, sedgrid, astfile, outname=None,
                      absflux_a_matrix=None, *args, **kwargs):
    """
    Parameters
    ----------
    project: str
        token of the project this task belongs to

    sedgrid: SEDGrid instance
        sed model grid for everyone of which we will evaluate the model

    outname: str
        path and filename into which save the noise model

    astfile: str
        path to the file into which are ASTs results

    absflux_a_matrix: ndarray
        absolute calibration a matrix giving the fractional uncertainties
        including correlated terms (off diagonals)

    Returns
    -------

    project: str
        token of the project this task belongs to

    noisefile: str
        noisemodel file name

    sedgrid: grid.SEDGrid instance
        SED model grid instance
    """
    if outname is None:
        outname = '{0:s}_noisemodel.hd5'.format(project)

    noise_source = RequiredFile(outname, make_toothpick_noise_model, outname,
                                astfile, sedgrid,
                                absflux_a_matrix=absflux_a_matrix,
                                **kwargs)

    return project, noise_source(), sedgrid

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
