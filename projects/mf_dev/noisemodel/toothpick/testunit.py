import numpy as np
import tables

from beast.core.noisemodel import toothpick


def make_toothpick_noise_model(outname, astfile, sedgrid, covariance=None, **kwargs):
    """ toothpick noise model assumes that every filter is independent with
    any other.
    """
    #outname = dir_project + 'noise_model_b%s_%s.hd5' % (brick,subdivision)

    # read mag_in, mag_out
    model = toothpick.MultiFilterASTs(astfile, sedgrid.filters)
    # compute k-NN statistics: bias, stddev, completeness
    model.compute_knn_statistics(k=10, eps=0, completeness_mag_cut=80)
    # evaluate the noise model for all the models in sedgrid
    bias, sigma, compl = model(sedgrid)
    # save to disk/mem

    # absolute flux calibration uncertainties
    if covariance is not None:
        if covariance.ndim == 1:
            abs_calib_2 = covariance[:] ** 2
        else:   # assumes a cov matrix
            abs_calib_2 = np.diag(covariance)

        noise = np.sqrt(abs_calib_2 * sedgrid.seds ** 2 + sigma ** 2)
    else:
        noise = sigma

    with tables.openFile(outname, 'w') as outfile:
        outfile.createArray(outfile.root,'bias', bias)
        outfile.createArray(outfile.root,'error', noise)
        outfile.createArray(outfile.root,'completeness', compl)
