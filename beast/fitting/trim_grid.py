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
        minimum factor for trimming the upper and lower range of grid so that
        the model range cuts off sigma_fac above and below the brightest
        and faintest models, respectively (default: 2.0)
    inFlux : bool, optional
        if true data are in fluxes (default: True)
    trunchen : bool, optional
        if true use the trunchen noise model (default: False)
    """

    n_models, n_filters = sedgrid.seds.shape
    
    # ------------------------------------------------------------
    # Calculate flux or magnitude limits - 0.1th and 99.9th percentiles
    # This ensures that models are not being kept solely to try to fit the largest outliers.
    # ------------------------------------------------------------
    
    min_data = np.empty(n_filters)
    max_data = np.empty(n_filters)
    
    for k, filtername in enumerate(obsdata.filters):
        sfiltname = obsdata.filter_aliases[filtername]
        data = obsdata.data[sfiltname]
    
        if inFlux:
            data_flux = data * obsdata.vega_flux[k]
        else:
            data_flux = 10.0 ** (-0.4 * data) * obsdata.vega_flux[k]

        min_data[k] = np.percentile(data_flux, 0.1)
        max_data[k] = np.percentile(data_flux, 99.9)

    # ------------------------------------------------------------
    # Import noisemodel arrays
    # ------------------------------------------------------------
    model_bias = sedgrid_noisemodel["bias"]
    model_unc = sedgrid_noisemodel["error"]
    model_compl = sedgrid_noisemodel["completeness"]
    if trunchen:
        model_q_norm = sedgrid_noisemodel["q_norm"]
        model_icov_diag = sedgrid_noisemodel["icov_diag"]
        model_icov_offdiag = sedgrid_noisemodel["icov_offdiag"]

    # ------------------------------------------------------------
    # Model-Trimming Criterion 1: require nonzero completeness in every filter.
    # ------------------------------------------------------------
    completeness_count = np.zeros(n_models, dtype=np.int8)

    for k in range(n_filters):
        completeness_count += model_compl[:, k] > 0

    completeness_mask = completeness_count >= n_filters
    del completeness_count

    n_compl_trim_models = np.count_nonzero(completeness_mask)

    if n_compl_trim_models == 0:
        raise ValueError("no models are complete")

    # ------------------------------------------------------------
    # Model-Trimming Criterion 2: trim models that have fluxes unlike anything in the source catalog
    #
    # Start by calculating sigma_fac for the upper and lower bounds
    # This ensures that no relevant models will be trimmed accidentally.
    
    # Faint-side test uses the lower sigma factor:
    #   model + lower_sigma * uncertainty >= observed minimum (0.1th percentile)
    #
    # Bright-side test uses the upper sigma factor:
    #   model - upper_sigma * uncertainty <= observed maximum (99.9th percentile)
    #
    # Evaluate one filter at a time to avoid full-grid temporary flux-bound arrays.
    # ------------------------------------------------------------
    
    n_complete = n_compl_trim_models
    sigma_pass = np.zeros(n_complete)

    for k in range(n_filters):
        flux_k = (sedgrid.seds[:, k][completeness_mask] + model_bias[:, k][completeness_mask])
        unc_k = model_unc[:, k][completeness_mask]

        sigma_lower = (min_data[k] - flux_k) / unc_k
        sigma_upper = (flux_k - max_data[k]) / unc_k

        # Minimum sigma factor for each model to overlap this filter.
        sigma_pass = np.maximum(sigma_pass,np.maximum(0.0, np.maximum(sigma_lower, sigma_upper)),)

    lower_sigma_per_filter = np.empty(n_filters)
    upper_sigma_per_filter = np.empty(n_filters)

    for k in range(n_filters):
        flux_k = (sedgrid.seds[:, k][completeness_mask] + model_bias[:, k][completeness_mask])
        unc_k = model_unc[:, k][completeness_mask]

        sigma_lower = np.maximum(0.0, (flux_k - min_data[k]) / unc_k)
        sigma_upper = np.maximum(0.0, (max_data[k] - flux_k) / unc_k)

        # Each candidate model must satisfy the overlap requirement
        # in every filter, represented by sigma_pass.
        lower_sigma_per_filter[k] = np.min(np.maximum(sigma_pass, sigma_lower))
        upper_sigma_per_filter[k] = np.min(np.maximum(sigma_pass, sigma_upper))

    del flux_k, unc_k, sigma_lower, sigma_upper, sigma_pass

    real_lower_sigma_fac = max(sigma_fac, np.max(lower_sigma_per_filter))
    real_upper_sigma_fac = max(sigma_fac, np.max(upper_sigma_per_filter))

    print(f"Lower-limit sigma_fac: {real_lower_sigma_fac:.3f}")
    print(f"Upper-limit sigma_fac: {real_upper_sigma_fac:.3f}")

    # ------------------------------------------------------------
    # Now do the flux trimming
    # ------------------------------------------------------------
    faint_mask = np.ones(n_models, dtype=bool)
    bright_mask = np.ones(n_models, dtype=bool)

    for k in range(n_filters):
        flux_k = sedgrid.seds[:, k] + model_bias[:, k]
        unc_k = model_unc[:, k]

        # Faint-side criterion
        bound = unc_k * real_lower_sigma_fac
        bound += flux_k
        faint_mask &= bound >= min_data[k]

        # Bright-side criterion
        bound = unc_k * real_upper_sigma_fac
        bound *= -1.0
        bound += flux_k
        bright_mask &= bound <= max_data[k]

    flux_mask = faint_mask & bright_mask
    
    del faint_mask, bright_mask, bound, flux_k, unc_k

    n_flux_trim_models = np.count_nonzero(flux_mask)

    # ------------------------------------------------------------
    # Combine the two criteria
    # ------------------------------------------------------------
    combined_mask = completeness_mask & flux_mask
    n_trim_models = np.count_nonzero(combined_mask)

    if n_trim_models == 0:
        raise ValueError("no models that are within the data range")

    # Convert the final Boolean mask to indices for saving.
    indxs = np.flatnonzero(combined_mask)

    print("number of original models =", n_models)
    print("number of models trimmed for completeness =", n_models - n_compl_trim_models)
    print("number of models trimmed for flux =", n_models - n_flux_trim_models)
    print("number of remaining models =", n_trim_models)

    # ------------------------------------------------------------
    # Save the grid
    # ------------------------------------------------------------

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

