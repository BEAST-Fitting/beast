import os

import numpy as np
from astropy import units

from beast.physicsmodel.grid import SpectralGrid, SEDGrid
from beast.physicsmodel import creategrid
from beast.physicsmodel.stars import isochrone, stellib
from beast.physicsmodel.stars import evoltracks
from beast.physicsmodel.stars.isochrone import ezIsoch
from beast.physicsmodel.dust import extinction
from beast.physicsmodel.grid_and_prior_weights import (
    compute_distance_age_mass_metallicity_weights,
)
from beast.tools.beast_info import add_to_beast_info_file

__all__ = [
    "make_evoltrack_table",
    "make_spectral_grid",
    "add_stellar_priors",
    "make_extinguished_sed_grid",
]


def make_evoltrack_table(
    project,
    oet=None,
    age_info=[6.0, 10.13, 0.05],
    mass_info=[-1.0, 3.0, 0.05],
    z=[0.0152],
    et_fname=None,
    info_fname=None,
):
    """
    The evolutionary track table are loaded.
    Determined based on evolutionary tracks directly or isochrones derived
    from evolutionary tracks.

    Parameters
    ----------
    project : str
        project name

    oet : evoltracks.EvolTracks or isochrone.Isochrone object
        contains the full evolutionary track information

    age_info : list
        age information needed for isochrones [logtmin, logtmax, dlogt]

    mass_info : list
        mass information needed for evolutionary tracks [logmmin, logmmax, dlogm]

    z : float or sequence
        list of metalicity values, where default (Z=0.152) is adopted Z_sun
        for PARSEC/COLIBRI models

    et_fname : str
        Set to specify the filename to save the gridded evolutionary tracks
        to, otherwise saved to project/project_et.csv

    info_fname : str
        Set to specify the filename to save beast info to, otherwise
        saved to project/project_beast_info.asdf

    Returns
    -------
    fname: str
       name of saved file

    oet: evoltracks.EvolTracks or isochrone.Isochrone object
        contains the full evolutionary track information
    """
    if et_fname is None:
        et_fname = "%s/%s_et.csv" % (project, project)
    if not os.path.isfile(et_fname):
        if oet is None:
            oet = isochrone.PadovaWeb()

        if isinstance(oet, isochrone.Isochrone):
            logtmin, logtmax, dlogt = age_info
            t = oet._get_t_isochrones(max(5.0, logtmin), min(10.13, logtmax), dlogt, z)
            t.header["NAME"] = "{0} Isochrones".format("_".join(et_fname.split("_")[:-1]))
            print("{0} Isochrones".format("_".join(et_fname.split("_")[:-1])))
            info = {"project": project, "logt_input": age_info, "z_input": z}
            t.write(et_fname)
            # maybe needed as ezIsoch is a proxy for a Table
            # maybe we can just use a table????
            oet = ezIsoch(et_fname)
        elif isinstance(oet, evoltracks.EvolTracks):
            tab = oet.get_evoltracks(mass_info, z)
            tab.header["NAME"] = "{0} EvolTracks".format("_".join(et_fname.split("_")[:-1]))
            print("{0} EvolTracks".format("_".join(et_fname.split("_")[:-1])))
            info = {"project": project, "logm_input": mass_info, "z_input": z}
        else:
            print(f"Type {type(oet)} of evolutionary track not supported")

        # save info to the beast info file
        if info_fname is None:
            info_fname = f"{project}/{project}_beast_info.asdf"
        add_to_beast_info_file(info_fname, info)

    else:
        # read in the isochrone data from the file
        #   not sure why this is needed, but reproduces previous ezpipe method
        oet = ezIsoch(et_fname)

    return (et_fname, oet)


def make_spectral_grid(
    project,
    oiso,
    osl=None,
    bounds={},
    verbose=True,
    spec_fname=None,
    distance=10,
    distance_unit=units.pc,
    redshift=0.0,
    filterLib=None,
    add_spectral_properties_kwargs=None,
    extLaw=None,
    **kwargs,
):
    """
    The spectral grid is generated using the stellar parameters by
    interpolation of the isochrones and the generation of spectra into the
    physical units

    Parameters
    ----------
    project: str
        project name

    oiso: isochrone.Isochrone object
        set of isochrones to use

    osl: stellib.Stellib object
        Spectral library to use (default stellib.Kurucz)

    distance: float or list of float
        distances at which models should be shifted, specified as a
        single number or as [min, max, step]

        0 means absolute magnitude.

    distance_unit: astropy length unit or mag
        distances will be evenly spaced in this unit
        therefore, specifying a distance grid in mag units will lead to
        a log grid

    redshift: float
        Redshift to which wavelengths should be shifted
        Default is 0 (rest frame)

    spec_fname: str
        full filename to save the spectral grid into

    filterLib:  str
        full filename to the filter library hd5 file

    extLaw: extinction.ExtLaw (default=None)
        if set, only save the spectrum for the wavelengths over which the
        extinction law is valid

    add_spectral_properties_kwargs: dict
        keyword arguments to call :func:`add_spectral_properties`
        to add model properties from the spectra into the grid property table

    Returns
    -------
    fname: str
       name of saved file

    g: grid.SpectralGrid object
        spectral grid to transform
    """
    if spec_fname is None:
        spec_fname = "%s/%s_spec_grid.hd5" % (project, project)

    # remove the isochrone points with logL=-9.999
    oiso.data = oiso[oiso["logL"] > -9]

    if not os.path.isfile(spec_fname):
        osl = osl or stellib.Kurucz()

        # filter extrapolations of the grid with given sensitivities in
        # logg and logT
        if "dlogT" not in bounds:
            bounds["dlogT"] = 0.1
        if "dlogg" not in bounds:
            bounds["dlogg"] = 0.3

        # make the spectral grid
        if verbose:
            print("Make spectra")
        g = creategrid.gen_spectral_grid_from_stellib_given_points(
            osl, oiso.data, bounds=bounds
        )

        # Construct the distances array. Turn single value into
        # 1-element list if single distance is given.
        _distance = np.atleast_1d(distance)
        if len(_distance) == 3:
            mindist, maxdist, stepdist = _distance
            distances = np.arange(mindist, maxdist + stepdist, stepdist)
        elif len(_distance) == 1:
            distances = np.array(_distance)
        else:
            raise ValueError("distance needs to be (min, max, step) or single number")

        # calculate the distances in pc
        if distance_unit == units.mag:
            distances = np.power(10, distances / 5.0 + 1) * units.pc
        else:
            distances = (distances * distance_unit).to(units.pc)

        print("applying {} distances".format(len(distances)))

        if verbose:
            print(
                "Adding spectral properties:",
                add_spectral_properties_kwargs is not None,
            )
        if add_spectral_properties_kwargs is not None:
            nameformat = (
                add_spectral_properties_kwargs.pop("nameformat", "{0:s}") + "_nd"
            )

        # Apply the distances to the stars. Seds already at 10 pc, need
        # multiplication by the square of the ratio to this distance.
        # TODO: Applying the distances might have to happen in chunks
        # for larger grids.
        def apply_distance_and_spectral_props(g):
            # distance
            g = creategrid.apply_distance_grid(g, distances, redshift=redshift)

            # spectral props
            if add_spectral_properties_kwargs is not None:
                g = creategrid.add_spectral_properties(
                    g,
                    nameformat=nameformat,
                    filterLib=filterLib,
                    **add_spectral_properties_kwargs,
                )

            # extinction
            if extLaw is not None:

                ext_law_range_A = 1e4 / np.array(extLaw.x_range)
                valid_lambda = np.where(
                    (g.lamb > np.min(ext_law_range_A))
                    & (g.lamb < np.max(ext_law_range_A))
                )[0]

                g.lamb = g.lamb[valid_lambda]
                g.seds = g.seds[:, valid_lambda]

            return g

        # Perform the extensions defined above and Write to disk
        if hasattr(g, "write"):
            g = apply_distance_and_spectral_props(g)
            g.write(spec_fname)
        else:
            for gk in g:
                gk = apply_distance_and_spectral_props(gk)
                gk.write(spec_fname, append=True)

    g = SpectralGrid(spec_fname, backend="memory")

    return (spec_fname, g)


def add_stellar_priors(
    project,
    specgrid,
    distance_prior_model={"name": "flat"},
    age_prior_model={"name": "flat"},
    mass_prior_model={"name": "kroupa"},
    met_prior_model={"name": "flat"},
    verbose=True,
    priors_fname=None,
    info_fname=None,
    **kwargs,
):
    """
    make_priors -- compute the weights for the stellar priors

    Parameters
    ----------
    project: str
        project name

    specgrid: SpectralGrid object
        spectral grid to transform

    distance_prior_model: dict
        dict including prior model name and parameters

    age_prior_model: dict
        dict including prior model name and parameters

    mass_prior_model: dict
        dict including prior model name and parameters

    met_prior_model: dict
        dict including prior model name and parameters

    priors_fname: str
        full filename to which to save the spectral grid with priors

    info_fname : str
        Set to specify the filename to save beast info to, otherwise
        saved to project/project_beast_info.asdf

    Returns
    -------
    fname: str
       name of saved file

    g: SpectralGrid object
        spectral grid to transform
    """
    if priors_fname is None:
        priors_fname = "%s/%s_spec_w_priors.grid.hd5" % (project, project)
    if not os.path.isfile(priors_fname):

        if verbose:
            print("Make Prior Weights")

        compute_distance_age_mass_metallicity_weights(
            specgrid.grid,
            distance_prior_model=distance_prior_model,
            age_prior_model=age_prior_model,
            mass_prior_model=mass_prior_model,
            met_prior_model=met_prior_model,
            **kwargs,
        )

        # write to disk
        if hasattr(specgrid, "write"):
            specgrid.write(priors_fname)
        else:
            for gk in specgrid:
                gk.write(priors_fname, append=True)

    # save info to the beast info file
    info = {
        "distance_prior_model": distance_prior_model,
        "age_prior_model": age_prior_model,
        "mass_prior_model": mass_prior_model,
        "met_prior_model": met_prior_model,
    }
    if info_fname is None:
        info_fname = f"{project}/{project}_beast_info.asdf"
    add_to_beast_info_file(info_fname, info)

    # read in spectralgrid from file (possible not needed, need to check)
    g = SpectralGrid(priors_fname, backend="memory")

    return (priors_fname, g)


def make_extinguished_sed_grid(
    project,
    specgrid,
    filters,
    av=[0.0, 5, 0.1],
    rv=[0.0, 5, 0.2],
    fA=None,
    av_prior_model={"name": "flat"},
    rv_prior_model={"name": "flat"},
    fA_prior_model={"name": "flat"},
    extLaw=None,
    add_spectral_properties_kwargs=None,
    absflux_cov=False,
    verbose=True,
    seds_fname=None,
    filterLib=None,
    info_fname=None,
    **kwargs,
):

    """
    Create SED model grid integrated with filters and dust extinguished

    Parameters
    ----------
    project: str
        project name

    specgrid: SpectralGrid object
        spectral grid to transform

    filters: sequence
        ordered sequence of filters to use to extract the photometry
        filter names are the full names in core.filters

    av: sequence
        sequence of Av values to sample

    av_prior_model: dict
        dict including prior model name and parameters

    rv: sequence
        sequence of Rv values to sample

    rv_prior_model: dict
        dict including prior model name and parameters

    fA: sequence (optional)
        sequence of fA values to sample (depending on extLaw definition)

    fA_prior_model: dict
        dict including prior model name and parameters

    extLaw: extinction.ExtLaw
        extinction law to use during the process

    add_spectral_properties_kwargs: dict
        keyword arguments to call :func:`add_spectral_properties`
        to add model properties from the spectra into the grid property table

    asbflux_cov: boolean
        set to calculate the absflux covariance matrices for each model
        (can be very slow!!!  But it is the right thing to do)

    seds_fname: str
        full filename to save the sed grid into

    filterLib:  str
        full filename to the filter library hd5 file

    info_fname : str
        Set to specify the filename to save beast info to, otherwise
        saved to project/project_beast_info.asdf

    Returns
    -------
    fname: str
       name of saved file

    g: SpectralGrid object
        spectral grid to transform
    """

    # create the dust grid arrays
    if len(av) > 3:
        # check if a log grid is requested
        if av[3] == "log":
            print("generating a log av grid")
            avs = 10 ** np.arange(np.log10(av[0]), np.log10(av[1]), av[2])
        else:
            print("generating a linear av grid")
            avs = np.arange(av[0], av[1] + 0.5 * av[2], av[2])
    else:
        print("generating a linear av grid")
        avs = np.arange(av[0], av[1] + 0.5 * av[2], av[2])
    rvs = np.arange(rv[0], rv[1] + 0.5 * rv[2], rv[2])
    if fA is not None:
        fAs = np.arange(fA[0], fA[1] + 0.5 * fA[2], fA[2])
    else:
        fAs = [1.0]

    # create SED file name if needed
    if seds_fname is None:
        seds_fname = "%s/%s_seds.grid.hd5" % (project, project)

    # generate extinguished grids if SED file doesn't exist
    if not os.path.isfile(seds_fname):

        extLaw = extLaw or extinction.Cardelli()

        if verbose:
            print("Make SEDS")

        if fA is not None:
            g = creategrid.make_extinguished_grid(
                specgrid,
                filters,
                extLaw,
                avs,
                rvs,
                fAs,
                av_prior_model=av_prior_model,
                rv_prior_model=rv_prior_model,
                fA_prior_model=fA_prior_model,
                add_spectral_properties_kwargs=add_spectral_properties_kwargs,
                absflux_cov=absflux_cov,
                filterLib=filterLib,
            )
        else:
            g = creategrid.make_extinguished_grid(
                specgrid,
                filters,
                extLaw,
                avs,
                rvs,
                av_prior_model=av_prior_model,
                rv_prior_model=rv_prior_model,
                add_spectral_properties_kwargs=add_spectral_properties_kwargs,
                absflux_cov=absflux_cov,
            )

        # write to disk
        if hasattr(g, "write"):
            g.write(seds_fname)
        else:
            for gk in g:
                gk.write(seds_fname, append=True)

    # save info to the beast info file
    info = {
        "av_input": av,
        "rv_input": rv,
        "fA_input": fA,
        "avs": avs,
        "rvs": rvs,
        "fAs": fAs,
        "av_prior_model": av_prior_model,
        "rv_prior_model": rv_prior_model,
        "fA_prior_model": fA_prior_model,
    }
    if info_fname is None:
        info_fname = f"{project}/{project}_beast_info.asdf"
    add_to_beast_info_file(info_fname, info)

    g = SEDGrid(seds_fname, backend="memory")

    return (seds_fname, g)
