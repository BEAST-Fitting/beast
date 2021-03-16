import numpy as np
import argparse
import asdf

from beast.physicsmodel.grid import SEDGrid
import beast.observationmodel.noisemodel.generic_noisemodel as noisemodel
from beast.observationmodel.observations import gen_SimObs_from_sedgrid

from astropy.table import vstack


def simulate_obs(
    physgrid_list,
    noise_model_list,
    output_catalog,
    beastinfo_list=None,
    nsim=100,
    compl_filter="max",
    complcut=None,
    magcut=None,
    weight_to_use="weight",
    ranseed=None,
):
    """
    Simulate photometry based on a physicsmodel grid(s) and observation model(s).

    Parameters
    ----------
    physgrid_list : list of strings
        Name of the physics model file.  If there are multiple physics model
        grids (i.e., if there are subgrids), list them all here, and they will
        each be sampled nsim/len(physgrid_list) times.

    noise_model_list : list of strings
        Name of the noise model file.  If there are multiple files for
        physgrid_list (because of subgrids), list the noise model files
        associated with each physics model file.

    output_catalog : string
        Name of the output simulated photometry catalog

    beastinfo_list : list of strings (default=None)
        Name of the beast info file.  The mass and age prior models are read
        from this model to use to compute the number of stars to simulate. If
        there are multiple files for physgrid_list (because of subgrids), list
        the beast info files associated with each physics model file.

    n_sim : int (default=100)
        Number of simulated objects to create if beastinfo_list is not given. If
        nsim/len(physgrid_list) isn't an integer, this will be increased so that
        each grid has the same number of samples.

    compl_filter : str (default=max)
        filter to use for completeness (required for toothpick model)
        set to max to use the max value in all filters

    complcut : float (defualt=None)
        completeness cut for only including model seds above the given
        completeness, where the cut ranges from 0 to 1.

    magcut : float (defualt=None)
        faint-end magnitude cut for only including model seds brighter
        than the given magnitude in compl_filter.

    weight_to_use : string (default='weight')
        Set to either 'weight' (prior+grid), 'prior_weight', 'grid_weight',
        or 'uniform' (this option is valid only when nsim is supplied) to
        choose the weighting for SED selection.

    ranseed : int
        seed for random number generator
    """
    # numbers of samples to do
    # (ensure there are enough for even sampling of multiple model grids)
    n_phys = len(np.atleast_1d(physgrid_list))
    nsim = int(nsim)
    samples_per_grid = int(np.ceil(nsim / n_phys))

    if complcut is not None:
        complcut = float(complcut)

    if magcut is not None:
        magcut = float(magcut)

    if ranseed is not None:
        ranseed = int(ranseed)

    # list to hold all simulation tables
    simtable_list = []

    # make a table for each physics model + noise model
    for k, physgrid in enumerate(np.atleast_1d(physgrid_list)):
        noise_model = np.atleast_1d(noise_model_list)[k]

        # get the physics model grid - includes priors
        modelsedgrid = SEDGrid(str(physgrid))

        # read in the noise model - includes bias, unc, and completeness
        noisegrid = noisemodel.get_noisemodelcat(str(noise_model))

        if beastinfo_list is not None:
            with asdf.open(np.atleast_1d(beastinfo_list)[k]) as af:
                binfo = af.tree
                age_prior_model = binfo["age_prior_model"]
                mass_prior_model = binfo["mass_prior_model"]
        else:
            age_prior_model = None
            mass_prior_model = None

        # generate the table
        simtable = gen_SimObs_from_sedgrid(
            modelsedgrid,
            noisegrid,
            age_prior_model=age_prior_model,
            mass_prior_model=mass_prior_model,
            nsim=samples_per_grid,
            compl_filter=compl_filter,
            complcut=complcut,
            magcut=magcut,
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
        "--beastinfo_list",
        metavar="BEAST_INFO",
        nargs="+",
        help="filename(s) of beast info file(s)",
    )
    parser.add_argument(
        "--nsim", default=100, type=int, help="number of simulated objects"
    )
    parser.add_argument(
        "--compl_filter",
        default="F475W",
        help="filter for completeness, set to max for max of values from all filters",
    )
    parser.add_argument(
        "--complcut",
        default=None,
        type=float,
        help="completeness cut for selecting seds above the completeness cut"
    )
    parser.add_argument(
        "--magcut",
        default=None,
        type=float,
        help="magnitdue cut for selecting seds brighter than the magcut in compl_filter"
    )
    parser.add_argument(
        "--weight_to_use",
        default="weight",
        type=str,
        help="""Set to either 'weight' (prior+grid), 'prior_weight', or
        'grid_weight' to choose the weighting for SED selection.""",
    )
    parser.add_argument(
        "--ranseed",
        default=None,
        type=int,
        help="seed for random number generator"
    )
    args = parser.parse_args()

    # run observation simulator
    simulate_obs(
        args.physgrid_list,
        args.noise_model_lists,
        args.output_catalog,
        beastinfo_list=args.beastinfo_list,
        nsim=args.nsim,
        compl_filter=args.compl_filter,
        complcut=args.complcut,
        magcut=args.magcut,
        weight_to_use=args.weight_to_use,
        ranseed=args.ranseed,
    )
