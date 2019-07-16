# system imports
from __future__ import (absolute_import, division, print_function)

# other imports
from multiprocessing import Pool



def subcatalog_fname(full_cat_fname, source_density, sub_source_density):
    """
    Return the name of a sub-catalog

    Parameters
    ----------
    full_cat_fname : string
        name of the photometry catalog

    source_density : string
        the current source density bin

    sub_source_density : string
        the current sub-file for the source density bin

    Returns
    -------
    string
        the file name of the sub-catalog

    """
    
    return full_cat_fname.replace('.fits', '_SD{}_sub{}.fits'.format(
            source_density.replace('_','-'), sub_source_density))



def parallel_wrapper(function, arg_tuples, nprocs=1):
    """
    A wrapper to automatically either run the function as-is or run it with parallel processes

    Parameters
    ----------
    function : function
        the function to be evaluated

    argument : list of tuples
        the input to the function (details of course depend on the function)

    nprocs : int (default=1)
        number of parallel processes (no parallelization if nprocs=1)


    Returns
    -------
    nothing

    """

    if (nprocs > 1):
        p = Pool(nprocs)
        for r in p.starmap(gen_subgrid, arg_tuples):
            print(r)
    else:
        for a in arg_tuples:
            r = function(*a)
            print(r)





def get_modelsubgridfiles(subgrid_names_file):
    """
    Read in the file that has the list of subgridded physicsmodel files

    Parameters
    ----------
    subgrid_names_file : string
        name of the file with the list of names


    Returns
    -------
    list of strings
        the names of the subgridded physicsmodel files

    """


    with open(subgrid_names_file, 'r') as f:
        modelsedgridfiles = f.read().split('\n')[:-1]

    return modelsedgridfiles
