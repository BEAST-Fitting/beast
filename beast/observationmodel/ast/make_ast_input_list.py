import os
import warnings
import numpy as np
from astropy.io import ascii
from astropy.table import Table
from astropy.table import Column

from beast.observationmodel.vega import Vega
from beast.physicsmodel.grid import SEDGrid

import h5py


def mag_limits(seds, faint_cut, Nfilter=1, bright_cut=None):
    """
    Selects models which have at least N filter above the limits

    Parameters
    ----------
    seds:   np.array
            Magnitude array from BEAST grid

    faint_cut: list
            List of limit magnitudes on the faint end

    Nfilter: integer
             In how many filters, you want a fake star to be brighter
             than the limit (or fainter than the upper limit)

    bright_cut: list
        List of limit magnitudes on the bright end. Useful for cutting
        out bright, nearby models, when no such bright nearby stars are
        present in the data

    Returns
    -------
    idx:    np.array
            Array of integers contining the indices of allowed models

    """
    flag = seds.copy()

    # flag is True if the models are brigter (=smaller number in mag)
    # than the limits
    for i, limit in enumerate(faint_cut):
        flag[:, i] = seds[:, i] < limit

    # flag is True if the models are frainter than the upper brightness
    # limits
    if bright_cut is not None:
        for i, limit in enumerate(bright_cut):
            flag[:, i] = np.logical_and(flag[:, i], seds[:, i] > limit)

    # Keep index where model is brighter than the limit in N filters
    s = np.sum(flag, axis=1)
    (idx,) = np.where(s >= Nfilter)

    return idx


def pick_models_toothpick_style(
    sedgrid_fname,
    filters,
    N_fluxes,
    min_N_per_flux,
    mag_cuts=None,
    Nfilters=3,
    outfile=None,
    outfile_params=None,
    bins_outfile=None,
):
    """
    Creates a fake star catalog from a BEAST model grid. The chosen seds
    are optimized for the toothpick model, by working with a given
    number of flux bins, and making sure that every flux bin is covered
    by at least a given number of models (for each filter individually,
    which is how the toothpick model works).

    Parameters
    ----------
    sedgrid_fname: string
        BEAST model grid from which the models are picked (hdf5 file)

    filters: list of string
        Names of the filters, to be used as columns of the output table

    N_fluxes: integer
        The number of flux bins into which the dynamic range of the
        model grid in each filter is divided

    min_N_per_flux: integer
        Minimum number of model seds that need to fall into each bin

    mag_cuts: dictionary (optional)
        Defines the magnitudes at which to cut the physics grid models
        for AST generation.

    Nfilters: integer (default=3)
        Defines the number of filters that must have fluxes inside the
        mag_cuts bounds to be accepted.

    outfile: string
        Output path for the models (optional). If this file already
        exists, the chosen seds are loaded from this file instead.

    outfile_params: string (default=None)
        If a file name is given, the physical parameters associated with
        each model will be written to disk

    bins_outfile: string
        Output path for a file containing the flux bin limits for each
        filter, and the number of samples for each (optional)

    bright_cut: list of float
        List of magnitude limits for each filter (won't sample model
        SEDs that are too bright)

    Returns
    -------
    sedsMags: astropy Table
        A table containing the selected model seds (columns are named
        after the filters)

    """
    if outfile is not None and os.path.isfile(outfile):
        print(
            "{} already exists. Will attempt to load SEDs for ASTs from there.".format(
                outfile
            )
        )
        t = Table.read(outfile, format="ascii")
        return t

    with Vega() as v:
        vega_f, vega_flux, lambd = v.getFlux(filters)

    modelsedgrid = SEDGrid(sedgrid_fname)

    sedsMags = -2.5 * np.log10(modelsedgrid.seds[:] / vega_flux)

    # Check if logL=-9.999 model points saliently sneak through
    if min(modelsedgrid.grid["logL"]) < -9:
        warnings.warn("There are logL=-9.999 model points in the SED grid!")
        print("Excluding those SED models from selecting input ASTs")
        idxs = np.where(modelsedgrid.grid["logL"] > -9)[0]
        sedsMags = sedsMags[idxs]

    # remove seds that have fluxes outside mag_cuts
    if mag_cuts is not None:
        if not isinstance(mag_cuts, dict):
            warnings.warn("ast_fluxbin_maglimits must be a dictionary like {'HST_ACS_WFC_F435W': [bright, faint]}. Skipping magnitude trimming.")
        else:
            faint_cuts = [mag_cuts[f][1] for f in filters]  # faint = upper (larger mag)
            bright_cuts = [mag_cuts[f][0] for f in filters] # bright = lower (smaller mag)
            idxs = mag_limits(sedsMags, faint_cuts, Nfilter=Nfilters, bright_cut=bright_cuts)
            print("Trimmed {} SEDs from physics model".format(len(sedsMags)-len(idxs)))
            sedsMags = sedsMags[idxs]

    Nseds = sedsMags.shape[0]
    Nf = sedsMags.shape[1]
    idxs = np.arange(Nseds)

    # Set up a number of flux bins for each filter
    maxes = np.amax(sedsMags, axis=0)
    mins = np.amin(sedsMags, axis=0)

    bin_edges = np.zeros((N_fluxes + 1, Nf))  # indexed on [fluxbin, nfilters]
    for f in range(Nf):
        bin_edges[:, f] = np.linspace(mins[f], maxes[f], N_fluxes + 1)
    bin_mins = bin_edges[:-1, :]
    bin_maxs = bin_edges[1:, :]
    if not len(bin_mins) == len(bin_maxs) == N_fluxes:
        raise AssertionError()

    bin_count = np.zeros((N_fluxes, Nf))
    include_mask = np.full(idxs.shape, True, dtype=bool)
    chosen_idxs = []
    counter = 0
    successes = 0
    chunksize = 100000
    
    while True:
        counter += 1
        # pick some random models
        rand_idx = np.random.choice(idxs[include_mask], size=chunksize)
        randomseds = sedsMags[rand_idx, :]

        # Find in which bin each model belongs, for each filter
        fluxbins = np.zeros(randomseds.shape, dtype=int)
        for fltr in range(Nf):
            fluxbins[:, fltr] = np.digitize(randomseds[:, fltr], bin_maxs[:, fltr])

        # Clip in place (models of which the flux is equal to the max
        # are assigned bin nr N_fluxes. Move these down to bin nr
        # N_fluxes - 1)
        np.clip(fluxbins, a_min=0, a_max=N_fluxes - 1, out=fluxbins)
        

        need = np.maximum(min_N_per_flux - bin_count, 0)   # how many still needed per (bin,filter)
        add_these = np.zeros(len(rand_idx), dtype=bool)   # which SEDs from this chunk to accept
        
        for f in range(Nf):
            bins_needed = np.nonzero(need[:, f] > 0)[0]   # integer bin indices that still need filling
            if bins_needed.size == 0:
                continue
        
            # For each bin that needs samples, choose up to `need` SEDs from this chunk that fall there
            for b in bins_needed:
                n_to_fill = int(need[b, f])
                if n_to_fill <= 0:
                    continue
        
                # robustly get integer indices of SEDs in this chunk that land in bin b for filter f
                sed_hits = np.flatnonzero(fluxbins[:, f] == b)  # always ndarray of ints (possibly empty)
                if sed_hits.size == 0:
                    continue
        
                # choose up to n_to_fill distinct SEDs (no replace if enough hits)
                n_select = min(n_to_fill, sed_hits.size)
                if sed_hits.size <= n_select:
                    chosen_local = sed_hits  # take them all
                else:
                    chosen_local = np.random.choice(sed_hits, size=n_select, replace=False)
        
                # mark them for addition, update counts, and mark as used in include_mask
                add_these[chosen_local] = True
                bin_count[b, f] += chosen_local.size


        # exclude the indices of the added models from being selected again
        include_mask[rand_idx[~add_these]] = False

        # update the number of successful model selections
        successes += add_these.sum()
        
        # Increment bin counts only for those models
        # We can vectorize updates using np.add.at (handles repeated indices correctly)
        np.add.at(bin_count, (fluxbins[add_these, :].ravel(), 
                              np.tile(np.arange(Nf), np.sum(add_these))), 1)    
        
        # Add the approved models
        chosen_idxs.extend(rand_idx[add_these])


        # If some of the randomly picked models were not added
        if not add_these.any():
            # ... check if we have enough samples everywhere, or if all
            # the models have been exhausted (and hence the bins are
            # impossible to fill).
            enough_samples = (bin_count.flatten() >= min_N_per_flux).all()
            still_models_left = include_mask.any()
            if enough_samples or not still_models_left:
                break

        if not counter % 10:
            print(
                "Sampled {} models. {} successfull seds. Ratio = {}".format(
                    counter * chunksize, successes, successes / counter / chunksize
                )
            )

    # Gather the selected model seds in a table
    sedsMags = Table(sedsMags[chosen_idxs, :], names=filters)

    if outfile is not None:
        ascii.write(
            sedsMags,
            outfile,
            overwrite=True,
            formats={k: "%.5f" for k in sedsMags.colnames},
        )

    # if chosen, save the corresponding model parameters
    if outfile_params is not None:
        grid_dict = {}
        for key in list(modelsedgrid.grid.keys()):
            grid_dict[key] = modelsedgrid.grid[key][chosen_idxs]
        grid_dict["sedgrid_indx"] = chosen_idxs
        ast_params = Table(grid_dict)
        ast_params.write(outfile_params, overwrite=True)

    if bins_outfile is not None:
        bin_info_table = Table()
        col_bigarrays = [bin_mins, bin_maxs, bin_count]
        col_basenames = ["bin_mins_", "bin_maxs_", "bin_count_"]
        for fltr, filter_name in enumerate(filters):
            for bigarray, basename in zip(col_bigarrays, col_basenames):
                bin_info_table.add_column(
                    Column(bigarray[:, fltr], name=basename + filter_name)
                )
        ascii.write(bin_info_table, bins_outfile, overwrite=True)

    return sedsMags


def pick_models(
    sedgrid_fname,
    filters,
    mag_cuts,
    Nfilter=3,
    N_stars=70,
    Nrealize=20,
    outfile=None,
    outfile_params=None,
    bright_cut=None,
    vega_fname=None,
    ranseed=None,
):
    """Creates a fake star catalog from a BEAST model grid

    Parameters
    ----------
    sedgrid_fname: string
        BEAST model grid from which the models are picked (hdf5 file)

    filters: list of string
        Names of the filters

    mag_cuts: list
        List of magnitude limits for each filter

    Nfilter: Integer
             In how many filters, you want a fake star to be brighter
             than the limit (mag_cut) (default = 3)

    N_stars: Integer
               Number of stellar models picked per a single log(age)
               (default=70)

    Nrealize: Integer
              Number of realization of each models (default = 20)

    outfile: str
        If a file name is given, the selected models will be written to
        disk

    outfile_params: str
        If a file name is given, the physical parameters associated with
        each model will be written to disk

    bright_cut: list of float
        Same as mag_cuts, but for the bright end

    vega_fname: str
        filename of vega file

    ranseed : int
        used to set the seed to make the results reproducable
        useful for testing

    Returns
    -------
    astropy Table of selected models
    - and optionally -
    ascii file: A list of selected models, written to 'outfile'
    fits file: the corresponding physical parameters, written to 'outfile_params'
    """

    with Vega(source=vega_fname) as v:  # Get the vega fluxes
        vega_f, vega_flux, lamb = v.getFlux(filters)

    modelsedgrid = SEDGrid(sedgrid_fname)

    # Convert to Vega mags
    sedsMags = -2.5 * np.log10(modelsedgrid.seds[:] / vega_flux)

    # make sure Nfilters isn't larger than the total number of filters
    if Nfilter > len(filters):
        Nfilter = len(filters)

    # Select the models above the magnitude limits in N filters
    idxs = mag_limits(sedsMags, mag_cuts, Nfilter=Nfilter, bright_cut=bright_cut)
    cols = {}
    for key in list(modelsedgrid.grid.keys()):
        cols[key] = modelsedgrid.grid[key][idxs]
    grid_cut = Table(cols)

    # Sample the model grid uniformly
    prime_params = np.column_stack(
        (grid_cut["logA"], grid_cut["M_ini"], grid_cut["Av"])
    )
    search_age = np.unique(prime_params[:, 0])

    N_sample = N_stars
    model_ind = []  # indices for the model grid
    ast_params = grid_cut[[]]  # the corresponding model parameters

    # set the random seed - mainly for testing
    if not None:
        np.random.seed(ranseed)

    for iage in search_age:
        (tmp,) = np.where(prime_params[:, 0] == iage)
        new_ind = np.random.choice(tmp, N_sample)
        model_ind.append(new_ind)
        [ast_params.add_row(grid_cut[new_ind[i]]) for i in range(len(new_ind))]

    index = np.repeat(idxs[np.array(model_ind).reshape((-1))], Nrealize)
    sedsMags = Table(sedsMags[index, :], names=filters)

    if outfile is not None:
        ascii.write(
            sedsMags,
            outfile,
            overwrite=True,
            formats={k: "%.5f" for k in sedsMags.colnames},
        )

    if outfile_params is not None:
        ast_params.write(outfile_params, overwrite=True)

    return sedsMags


def supplement_ast(
    sedgrid_fname,
    filters,
    nAST=1000,
    existingASTfile=None,
    outASTfile=None,
    outASTfile_params=None,
    mag_cuts=None,
    color_cuts=None,
):
    """
    Creates an additional fake star catalog from a BEAST model grid
    that fulfills the customized conditions to supplement input ASTs.
    If the existing input AST parameter file is given, already selected
    models will be excluded from this process. The input artificial
    stars are picked randomly from the remaining models.

    Parameters
    ----------
    sedgrid_fname: string
        BEAST model grid from which the models are picked (hdf5 file)

    filters: list of string
        Names of the filters

    nAST: int
        Number of unique additional ASTs per source density bin

    existingASTfile: string (optional, default=None)
        Name of the existing input AST parameter file. If not None,
        the models that were already listed in the existing list Will
        be removed by default

    outASTfile: string (optional, default=None)
        Output file name for the chosen models

    outASTfile_params: string (optional, default=None)
        If a file name is given, the physical parameters associated with
        each model will be written to disk

    mag_cut: dictionary (optional, default=None)
        Dictionary of bright and faint magnitude limits for given filters.
        The way to specify the cuts is by updating the "ast_suppl_maglimit" key
        in the beast_settings file. This is a dictionary that includes information
        for the magnitude cuts as a function of the filters included in observation.

        For example, for a field observed with HST_WFC3_F336W, HST_WFC3_F475W,
        and HST_WFC3_F814W, to set a magnitude range limit of 16<HST_WFC3_F475W<28 mag,
        and 15<HST_WFC3_F814W<27 mag you need to set the following within the beast_settings file:

        # specify that the ast_supplement mode should be on
        ast_supplement = True

        # initialize and populate the dictionary of desired magnitude limits
        ast_suppl_maglimits = {}
        # the magntidue limits are defined by the filter and a list of the limits in magnitudes
        ast_suppl_maglimits["HST_WFC3_F475W"] = [16,28]
        ast_suppl_maglimits["HST_WFC3_F814W"] = [15,27]

        # set the key word
        ast_suppl_maglimit = ast_suppl_maglimits

    color_cut: dictionary (optional, default=None)
        Dictionary of red color limits for given filters.
        The way to specify the cuts is by updating the "ast_suppl_colorlimit" key
        in the beast_settings file. This is a dictionary that includes information
        for the color cuts as a function of the filters included in observation.

        For example, for a field observed with HST_WFC3_F336W, HST_WFC3_F475W,
        and HST_WFC3_F814W, to set a color range limit of HST_WFC3_F475W-HST_WFC3_F814W<6,
        HST_WFC3_F336W-HST_WFC3_F475W<5 and HST_WFC3_F336W-HST_WFC3_F814W<4, you need
        to set the following within the beast_settings file:

        # specify that the ast_supplement mode should be on
        ast_supplement = True

        # initialize the dictionary of desired magnitude limits
        ast_suppl_colorlimits = {}

        # the color limits are defined by the first filter in the color (e.g, X for X-Y),
        # and the input is a list including the second filter (e.g., Y for X-Y) and the
        # color limit in magnitudes
        ast_suppl_colorlimits["HST_WFC3_F475W"] = [["HST_WFC3_F814W",6]]
        ast_suppl_colorlimits["HST_WFC3_F336W"] = [["HST_WFC3_F475W",5], ["HST_WFC3_F814W",4]]

        # set the key word
        ast_suppl_colorlimit =  ast_suppl_colorlimits

    Returns
    -------
    sedsMags: astropy Table
        A table containing the selected model seds (columns are named
        after the filters)

    """

    with Vega() as v:
        vega_f, vega_flux, lambd = v.getFlux(filters)

    modelsedgrid = SEDGrid(sedgrid_fname)

    # Convert to Vega mags
    sedsMags = -2.5 * np.log10(modelsedgrid.seds[:] / vega_flux)

    Nseds = sedsMags.shape[0]
    sedsIndx = np.arange(Nseds)

    if existingASTfile is not None and os.path.isfile(existingASTfile):
        print(
            "{} exists. Will attempt to load SEDs for ASTs from there \
            and remove those SEDs from the SED grid".format(
                existingASTfile
            )
        )
        print("existing AST file", existingASTfile)
        t = Table.read(existingASTfile, format="fits")
        sedsMags = np.delete(sedsMags, t["sedgrid_indx"], axis=0)
        sedsIndx = np.delete(sedsIndx, t["sedgrid_indx"])
        Nseds = sedsMags.shape[0]

    # Apply selection conditions if supplied
    # Just magnitude cuts
    print("mag_cuts", mag_cuts)
    print("color_cuts", color_cuts)
    if mag_cuts is not None:
        cond = np.ones(Nseds, dtype=bool)
        for key in list(mag_cuts.keys()):
            idx_filter = [i for i, iflt in enumerate(filters) if key in iflt]
            bright_cut = mag_cuts[key][0]
            faint_cut = mag_cuts[key][1]
            tmp_cond = np.logical_and(
                (sedsMags[:, idx_filter] >= bright_cut),
                (sedsMags[:, idx_filter] <= faint_cut),
            )

            if color_cuts is not None:
                if key in color_cuts:
                    for limit in color_cuts[key]:

                        idx_color_filter = [
                            i for i, iflt in enumerate(filters) if limit[0] in iflt
                        ]
                        tmp_cond = np.logical_and(
                            tmp_cond,
                            (
                                sedsMags[:, idx_filter] - sedsMags[:, idx_color_filter]
                                <= limit[1]
                            ),
                        )
            cond = np.logical_and(cond, tmp_cond.ravel())

        sedsMags = sedsMags[cond, :]
        sedsIndx = sedsIndx[cond]

    # Randomly select models
    # Supplementing ASTs does not need to follow
    # the toothpick-way selection
    chosen_idxs = np.random.choice(len(sedsIndx), nAST)
    sedsIndx = sedsIndx[chosen_idxs]

    # Gather the selected model seds in a table
    sedsMags = Table(sedsMags[chosen_idxs, :], names=filters)

    if outASTfile is not None:
        ascii.write(
            sedsMags,
            outASTfile,
            overwrite=True,
            formats={k: "%.5f" for k in sedsMags.colnames},
        )

    # if chosen, save the corresponding model parameters
    if outASTfile_params is not None:
        grid_dict = {}
        for key in list(modelsedgrid.grid.keys()):
            grid_dict[key] = modelsedgrid.grid[key][sedsIndx]
        grid_dict["sedgrid_indx"] = sedsIndx
        ast_params = Table(grid_dict)
        ast_params.write(outASTfile_params, overwrite=True)

    return sedsMags
