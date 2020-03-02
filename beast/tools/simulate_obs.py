import numpy as np
import argparse

from beast.physicsmodel.grid import FileSEDGrid
import beast.observationmodel.noisemodel.generic_noisemodel as noisemodel
from beast.observationmodel.observations import gen_SimObs_from_sedgrid

from astropy.table import vstack


def simulate_obs(
    physgrid_list,
    noise_model_list,
    output_catalog,
    nsim=100,
    compl_filter="F475W",
    weight_to_use='weight',
    ranseed=None,
):
    """
    Wrapper for creating a simulated photometry.

    Parameters
    ----------
    physgrid_list : list of strings
        Name of the physics model file.  If there are multiple physics model
        grids (i.e., if there are subgrids), list them all here, and they will
        each be sampled nsim/len(physgrid_list) times.

    noise_model_list : list of strings
        Name of the noise model file.  If there are multiple files for
        physgrid_list (because of subgrids), list the noise model file
        associated with each physics model file.

    output_catalog : string
        Name of the output simulated photometry catalog

    n_sim : int (default=100)
        Number of simulated objects to create.  If nsim/len(physgrid_list) isn't
        an integer, this will be increased so that each grid has the same
        number of samples.

    compl_filter : string (default='F475W')
        filter name to use for completeness

    weight_to_use : string (default='weight')
        Set to either 'weight' (prior+grid), 'prior_weight', or 'grid_weight' to
        choose the weighting for SED selection.

    ranseed : int
        seed for random number generator

    """

    # numbers of samples to do
    # (ensure there are enough for even sampling of multiple model grids)
    n_phys = len(physgrid_list)
    samples_per_grid = int(np.ceil(nsim / n_phys))

    # list to hold all simulation tables
    simtable_list = []

    # make a table for each physics model + noise model
    for physgrid, noise_model in zip(
        np.atleast_1d(physgrid_list), np.atleast_1d(noise_model_list)
    ):

        # get the physics model grid - includes priors
        modelsedgrid = FileSEDGrid(str(physgrid))

        # read in the noise model - includes bias, unc, and completeness
        noisegrid = noisemodel.get_noisemodelcat(str(noise_model))

        # generate the table
        simtable = gen_SimObs_from_sedgrid(
            modelsedgrid,
            noisegrid,
            nsim=samples_per_grid,
            compl_filter=compl_filter,
            weight_to_use=weight_to_use,
            ranseed=ranseed,
        )

        # append to the list
        simtable_list.append(simtable)

    # stack all the tables into one and write it out
    vstack(simtable_list).write(output_catalog, overwrite=True)


if __name__ == "__main__":  # pragma: no cover

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--physgrid_list",
        "-p",
        metavar="PHYS_MODEL",
        required=True,
        nargs="+",
        help="filename(s) of physics grid(s)",
    )
    parser.add_argument(
        "--noise_model_list",
        "-n",
        metavar="NOISE_MODEL",
        required=True,
        nargs="+",
        help="filename(s) of observation/noise grid(s)",
    )
    parser.add_argument(
        "--output_catalog",
        "-c",
        required=True,
        help="filename for simulated observations",
    )
    parser.add_argument(
        "--nsim", default=100, type=int, help="number of simulated objects"
    )
    parser.add_argument(
        "--compl_filter", default="F475W", help="filter name to use for completeness"
    )
    parser.add_argument(
        "--weight_to_use",
        default="weight",
        type=str,
        help="""Set to either 'weight' (prior+grid), 'prior_weight', or
        'grid_weight' to choose the weighting for SED selection."""
    )
    parser.add_argument(
        "--ranseed", default=None, type=int, help="seed for random number generator"
    )
    args = parser.parse_args()

    # run observation simulator
    simulate_obs(
        args.physgrid_list,
        args.noise_model_lists,
        args.output_catalog,
        nsim=args.nsim,
        compl_filter=args.compl_filter,
        weight_to_use=args.weight_to_use,
        ranseed=args.ranseed,
    )
