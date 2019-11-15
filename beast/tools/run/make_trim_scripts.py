# system imports
import numpy as np
import os
import argparse

# BEAST imports
from beast.tools import verify_params, setup_batch_beast_trim
from beast.tools.run import create_filenames

from difflib import SequenceMatcher


import datamodel
import importlib


def make_trim_scripts(num_subtrim=1, nice=None, prefix=None):
    """
    `setup_batch_beast_trim.py` uses file names to create batch trim files.  This
    generates all of the file names for that function.

    NOTE: This assumes you're using source density or background dependent noise
    models.

    Parameters
    ----------
    num_subtrim : int (default = 1)
        number of trim batch jobs

    nice : int (default = None)
        set this to an integer (-20 to 20) to prepend a "nice" level
        to the trimming command

    prefix : string (default=None)
        Set this to a string (such as 'source activate astroconda') to prepend
        to each batch file (use '\n's to make multiple lines)

    Returns
    -------
    job_files : list of strings
        Names of the newly created job files
    """

    # before doing ANYTHING, force datamodel to re-import (otherwise, any
    # changes within this python session will not be loaded!)
    importlib.reload(datamodel)
    # check input parameters
    verify_params.verify_input_format(datamodel)

    # make lists of file names
    file_dict = create_filenames.create_filenames(
        use_sd=True, nsubs=datamodel.n_subgrid,
    )
    # extract some useful ones
    photometry_files = file_dict["photometry_files"]
    modelsedgrid_files = file_dict["modelsedgrid_files"]
    noise_files = file_dict["noise_files"]
    modelsedgrid_trim_files = file_dict["modelsedgrid_trim_files"]
    noise_trim_files = file_dict["noise_trim_files"]
    # the unique sets of things
    unique_sedgrid = [
        x for i, x in enumerate(modelsedgrid_files) if i == modelsedgrid_files.index(x)
    ]

    # save the list of job files
    job_file_list = []

    # iterate through each model grid
    for i in range(datamodel.n_subgrid):

        # indices for this model grid
        grid_ind = [
            ind
            for ind, mod in enumerate(modelsedgrid_files)
            if mod == unique_sedgrid[i]
        ]

        # create corresponding files for each of those
        input_noise = [noise_files[ind] for ind in grid_ind]
        input_phot = [photometry_files[ind] for ind in grid_ind]
        # to get the trim prefix, find the common string between trimmed noise
        # files and trimmed SED files
        input_trim_prefix = []
        for ind in grid_ind:
            str1 = modelsedgrid_trim_files[ind]
            str2 = noise_trim_files[ind]
            # find longest match
            match = SequenceMatcher(None, str1, str2).find_longest_match(
                0, len(str1), 0, len(str2)
            )
            # grab that substring (and remove trailing "_")
            input_trim_prefix.append(str1[match.a : match.a + match.size][:-1])

        # check if the trimmed grids exist before moving on
        check_trim = [os.path.isfile(noise_trim_files[ind]) for ind in grid_ind]

        # if any aren't trimmed for this model grid, set up trimming
        if np.sum(check_trim) < len(input_noise):

            job_path = "./{0}/trim_batch_jobs/".format(datamodel.project)
            if datamodel.n_subgrid > 1:
                file_prefix = "BEAST_gridsub" + str(i)
            if datamodel.n_subgrid == 1:
                file_prefix = "BEAST"

            # generate trimming at-queue commands
            setup_batch_beast_trim.generic_batch_trim(
                unique_sedgrid[i],
                input_noise,
                input_phot,
                input_trim_prefix,
                job_path=job_path,
                file_prefix=file_prefix,
                num_subtrim=num_subtrim,
                nice=nice,
                prefix=prefix,
            )

            job_file_list.append(job_path + file_prefix + "_batch_trim.joblist")

    return job_file_list


if __name__ == "__main__":  # pragma: no cover
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--num_subtrim", type=int, default=1, help="number of trim batch jobs",
    )
    parser.add_argument(
        "--nice",
        type=int,
        default=None,
        help="prepend a 'nice' level to the trimming command",
    )
    parser.add_argument(
        "--prefix", type=str, default=None, help="string to prepend to each batch file",
    )

    args = parser.parse_args()

    make_trim_scripts(
        num_subtrim=args.num_subtrim, nice=args.nice, prefix=args.prefix,
    )
