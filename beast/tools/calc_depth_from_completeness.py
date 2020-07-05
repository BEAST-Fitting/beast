import numpy as np
from astropy.table import Table, vstack
from beast.physicsmodel.grid import SEDGrid
import beast.observationmodel.noisemodel.generic_noisemodel as noisemodel
from beast.observationmodel.vega import Vega


def calc_depth(
    physgrid_list,
    noise_model_list,
    completeness_value=0.5,
    vega_mag=True,
    vega_fname=None,
):
    """
    Calculate the observation depth of a field using the completeness.  Some
    fields have low completeness at both faint and bright fluxes; this finds
    the faintest flux at which the completeness exceeds the given value(s).

    Parameters
    ----------
    physgrid_list : string or list of strings
        Name of the physics model file.  If there are multiple physics model
        grids (i.e., if there are subgrids), list them all here.

    noise_model_list : string or list of strings
        Name of the noise model file.  If there are multiple files for
        physgrid_list (because of subgrids), list the noise model file
        associated with each physics model file.

    completeness_value : float or list of floats
        The completeness(es) at which to evaluate the depth. Completeness is
        defined in the range 0.0 to 1.0.

    vega_mag : boolean (default=True)
        If True, return results in Vega mags. Otherwise returns flux in
        erg/s/cm^2/A.

    vega_fname : string
        filename for the vega info (useful for testing)


    Returns
    -------
    depth_dict : dictionary
        keys are the filters present in the grid, each value is the flux or Vega
        mag for each of the requested completeness values

    """

    # ------ Reading in data

    # If there are subgrids, we can't read them all into memory.  Therefore,
    # we'll go through each one and just grab the relevant parts.
    compl_table_list = []

    # make a table for each physics model + noise model
    for physgrid, noise_model in zip(
        np.atleast_1d(physgrid_list), np.atleast_1d(noise_model_list)
    ):

        # get the physics model grid - includes priors
        modelsedgrid = SEDGrid(str(physgrid))
        if hasattr(modelsedgrid.seds, "read"):
            sed_grid = modelsedgrid.seds.read()
        else:
            sed_grid = modelsedgrid.seds
        # get list of filters
        filter_list = modelsedgrid.filters

        # read in the noise model
        noisegrid = noisemodel.get_noisemodelcat(str(noise_model))
        # get the completeness
        model_compl = noisegrid["completeness"]

        # put it all into a table
        table_dict = {filt: sed_grid[:, f] for f, filt in enumerate(filter_list)}
        table_dict.update(
            {filt + "_compl": model_compl[:, f] for f, filt in enumerate(filter_list)}
        )

        # append to the list
        compl_table_list.append(Table(table_dict))

    # stack all the tables into one
    compl_table = vstack(compl_table_list)

    # if chosen, get the vega fluxes for the filters
    if vega_mag:
        _, vega_flux, _ = Vega(source=vega_fname).getFlux(filter_list)

    # ------ Calculate depth

    # initialize dictionary to hold results
    depth_dict = {filt: [] for filt in filter_list}

    # grab numbers for each filter
    for f, filt in enumerate(filter_list):

        use_sed = compl_table[filt]
        use_comp = compl_table[filt + "_compl"]

        # get sorted versions of data
        sort_ind = np.argsort(use_sed)
        sort_sed = use_sed[sort_ind]
        sort_comp = use_comp[sort_ind]

        # grab depths
        for compl in np.atleast_1d(completeness_value):

            # first check whether the noise model even covers this completeness
            # (in case there weren't sufficient ASTs)
            if (compl < np.min(sort_comp)) or (compl > np.max(sort_comp)):
                depth_dict[filt].append(np.nan)
                continue

            # find first instance of completeness > N
            first_ind = np.where(sort_comp > compl)[0][0]
            # corresponding flux
            comp_flux = sort_sed[first_ind]
            # save it
            if vega_mag:
                depth_dict[filt].append(-2.5 * np.log10(comp_flux / vega_flux[f]))
            else:
                depth_dict[filt].append(comp_flux)

    # return the results
    return depth_dict
