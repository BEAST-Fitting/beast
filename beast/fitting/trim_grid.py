import numpy as np
import tables

from beast.physicsmodel.grid import SEDGrid
from astropy.table import Table

__all__ = ["trim_models"]


def trim_models(
    sedgrid,
    sedgrid_noisemodel,
    obsdata,
    sed_outname,
    noisemodel_outname,
    sigma_fac=3.0,
    n_detected=4,
    inFlux=True,
    trunchen=False,
):
    """
    For a given set of observations, there will be models that are so
    bright or faint that they will always have ~0 probability of fitting
    the data.  This program trims those models out of the SED grid
    so that time is not spent calculating model points that are always
    zero probability.

    Parameters
    ----------
    sedgrid : grid.SEDgrid instance
        model grid
    sedgrid_noisemodel : beast noisemodel instance
        noise model data
    obsdata : Observation object instance
        observation catalog
    sed_outname : str
        name for output sed file
    noisemodel_outname : str
        name for output noisemodel file
    sigma_fac : float, optional
        factor for trimming the upper and lower range of grid so that
        the model range cuts off sigma_fac above and below the brightest
        and faintest models, respectively (default: 3.)
    n_detected : int, optional
        minimum number of bands where ASTs yielded a detection for
        a given model, if fewer detections than n_detected this model
        gets eliminated (default: 4)
    inFlux : bool, optional
        if true data are in fluxes (default: True)
    trunchen : bool, optional
        if true use the trunchen noise model (default: False)
    """
    # Store the brigtest and faintest fluxes in each band (for data and asts)
    n_filters = len(obsdata.filters)
    min_data = np.zeros(n_filters)
    max_data = np.zeros(n_filters)
    min_models = np.zeros(n_filters)
    max_models = np.zeros(n_filters)
    for k, filtername in enumerate(obsdata.filters):
        sfiltname = obsdata.filter_aliases[filtername]
        if inFlux:
            min_data[k] = np.amin(obsdata.data[sfiltname] * obsdata.vega_flux[k])
            max_data[k] = np.amax(obsdata.data[sfiltname] * obsdata.vega_flux[k])
        else:
            min_data[k] = np.amin(
                10 ** (-0.4 * obsdata.data[sfiltname]) * obsdata.vega_flux[k]
            )
            max_data[k] = np.amax(
                10 ** (-0.4 * obsdata.data[sfiltname]) * obsdata.vega_flux[k]
            )

        min_models[k] = np.amin(sedgrid.seds[:, k])
        max_models[k] = np.amax(sedgrid.seds[:, k])

    # link to the noisemodel values
    model_bias = sedgrid_noisemodel["bias"]
    model_unc = sedgrid_noisemodel["error"]
    model_compl = sedgrid_noisemodel["completeness"]
    if trunchen:
        model_q_norm = sedgrid_noisemodel["q_norm"]
        model_icov_diag = sedgrid_noisemodel["icov_diag"]
        model_icov_offdiag = sedgrid_noisemodel["icov_offdiag"]

    # has to be complete in all filters - otherwise observation model not defined
    # toothpick model means that if compl = 0, then bias = 0, and sigma = 0 from ASTs
    above_ast = model_compl > 0
    sum_above_ast = np.sum(above_ast, axis=1)
    (indxs,) = np.where(sum_above_ast >= n_filters)
    n_ast_indxs = len(indxs)

    print("number of original models = ", len(sedgrid.seds[:, 0]))
    print("number of ast trimmed models = ", n_ast_indxs)

    if n_ast_indxs <= 0:
        raise ValueError("no models are brighter than the minimum ASTs run")

    # Find models with fluxes (with margin) between faintest and brightest data
    for k in range(n_filters):
        print("working on filter # = ", k)

        # Get upper and lower values for the models given the noise model
        #  sigma_fac defaults to 3.
        model_val = sedgrid.seds[indxs, k] + model_bias[indxs, k]
        model_down = model_val - sigma_fac * model_unc[indxs, k]
        model_up = model_val + sigma_fac * model_unc[indxs, k]

        # print(k, min(model_val), max(model_val), min(model_bias[indxs, k]))

        (nindxs,) = np.where((model_up >= min_data[k]) & (model_down <= max_data[k]))
        if len(nindxs) > 0:
            indxs = indxs[nindxs]

    if len(indxs) == 0:
        raise ValueError("no models that are within the data range")

    print("number of original models = ", len(sedgrid.seds[:, 0]))
    print("number of ast trimmed models = ", n_ast_indxs)
    print("number of trimmed models = ", len(indxs))

    # Save the grid
    print("Writing trimmed sedgrid to disk into {0:s}".format(sed_outname))
    cols = {}
    for key in list(sedgrid.grid.keys()):
        cols[key] = sedgrid.grid[key][indxs]

    # New column to save the index of the model in the full grid
    cols["fullgrid_idx"] = indxs.astype(int)
    g = SEDGrid(
        sedgrid.lamb, seds=sedgrid.seds[indxs], grid=Table(cols), backend="memory"
    )
    filternames = obsdata.filters
    g.header["filters"] = " ".join(filternames)

    # trimmed grid name
    g.write(sed_outname)

    # save the trimmed noise model
    print("Writing trimmed noisemodel to disk into {0:s}".format(noisemodel_outname))
    with tables.open_file(noisemodel_outname, "w") as outfile:
        outfile.create_array(outfile.root, "bias", model_bias[indxs])
        outfile.create_array(outfile.root, "error", model_unc[indxs])
        outfile.create_array(outfile.root, "completeness", model_compl[indxs])
        if trunchen:
            outfile.create_array(outfile.root, "q_norm", model_q_norm[indxs])
            outfile.create_array(outfile.root, "icov_diag", model_icov_diag[indxs])
            outfile.create_array(
                outfile.root, "icov_offdiag", model_icov_offdiag[indxs]
            )
