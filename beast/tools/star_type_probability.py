import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt

from beast.tools import read_beast_data

from astropy.table import Table


def star_type_probability(
    sed_files,
    lnp_files,
    output_filebase,
):
    """
    Parameters
    ----------
    sed_files : string or list of strings
        Name of the file(s) with the SED grid(s).  If a list, it's each part of
        the subgrid.

    lnp_files : string or list of strings
        Name of the file(s) with the lnP values.  If a list, it's in the same
        order as the subgrids above.

    output_filebase : string
        prefix for saving the file of probabilities ("_startype.fits" will be
        appended)

    """

    # read in the data
    sed_data = defaultdict(list)
    lnp_data = defaultdict(list)
    for (sed_file, lnp_file) in zip(np.atleast_1d(sed_files), np.atleast_1d(lnp_files)):
        sed_temp = read_beast_data.read_sed_data(str(sed_file), param_list=['Av', 'M_ini', 'logA'])
        for key in sed_temp:
            sed_data[key].append(sed_temp[key])
        lnp_temp = read_beast_data.read_lnp_data(str(lnp_file), shift_lnp=False)
        for key in lnp_temp:
            lnp_data[key].append(lnp_temp[key])

    # combine arrays from each file
    for key in sed_data:
        sed_data[key] = np.concatenate(sed_data[key])
    for key in lnp_data:
        lnp_data[key] = np.concatenate(lnp_data[key])

    # get the physical parameters associated with each lnP
    lnp_grid_vals = read_beast_data.get_lnp_grid_vals(sed_data, lnp_data)

    # make an array of star probabilities (instead of log probabilities)
    # scale lnP for each star to have max of 1 to avoid numerical issues
    prob_vals = np.exp(lnp_data['vals'] - np.max(lnp_data['vals'], axis=0))


    # evaluate probabilities of things
    star_prob = {}

    # - extinguished O star
    star_prob['ext_O_star'] = ext_O_star(prob_vals, lnp_grid_vals)
    # test using subsample of lnPs
    test_with_subsample('ext_O_star', 0.85, prob_vals, lnp_grid_vals)

    # - other things



    # write out the table
    Table(star_prob).write(output_filebase+"_startype.fits", overwrite=True)



def ext_O_star(prob_vals, lnp_grid_vals, use_ind=True):
    """

    Calculate the probability that each star is an extinguished O star:
    * initial mass >= 10 Msun
    * A_V >= 0.5 mag

    Some useful references for O/B stars
    https://ui.adsabs.harvard.edu/abs/2019A%26A...625A.104R/abstract
    https://ui.adsabs.harvard.edu/abs/2018A%26A...615A..40R/abstract
    https://ui.adsabs.harvard.edu/abs/2018A%26A...609A...7R/abstract

    Parameters
    ----------
    prob_vals : np.array
        array with probability values (shape = n_lnp, n_stars).  This should be
        exp(lnP), with lnPs scaled to max=1 (if desired) to avoid numerical issues

    lnp_grid_vals : dict
        dictionary with corresponding physical parameters

    use_ind : array of booleans (default=True)
        use this to test what happens if you use a subset of the lnP values

    Returns
    -------
    star_prob : array
        probability for each star

    """

    return (
        np.ma.array(
            prob_vals,
            mask=np.invert(
                (lnp_grid_vals['Av'] >= 0.5) &
                (lnp_grid_vals['M_ini'] >= 10) &
                (use_ind)
            )
        ).sum(axis=0) / np.ma.array(prob_vals, mask=np.invert(use_ind)).sum(axis=0)
    ).data



def test_with_subsample(function_name, sample_frac, prob_vals, lnp_grid_vals):
    """
    Test what happens if we were to randomly sample fewer of the lnP values.
    Specifically, create a plot comparing the probabilities from function_name
    using the full lnP list to those derived from a partial sample.

    Parameters
    ----------
    function_name : string
        name of the function to evaluate

    sample_frac : float
        fraction of the original lnPs to sample (e.g., 0.25 for a quarter)

    prob_vals : np.array
        (see function input info)

    lnp_grid_vals : dict
        (see function input info)

    """

    # get indices of subset
    n_lnp, n_stars = prob_vals.shape
    use_ind = np.ones((n_lnp, n_stars)).astype(bool)
    for i in range(n_stars):
        tot_lnp = np.sum(prob_vals[:,i] > 0)
        n_to_remove = int(np.floor(tot_lnp * (1-sample_frac)))
        use_ind[np.random.choice(int(tot_lnp), n_to_remove, replace=False), i] = 0

    # evaluate probabilities of full/partial samples
    star_prob = {}
    star_prob['full'] = globals()[function_name](prob_vals, lnp_grid_vals)
    star_prob['partial'] = globals()[function_name](
        prob_vals,
        lnp_grid_vals,
        use_ind=use_ind,
    )

    fig = plt.figure(figsize=(5,4))

    plt.plot(star_prob['full'], star_prob['partial'],
        marker="o",
        mew=0,
        color="black",
        markersize=2,
        linestyle="None",
        alpha=0.2)

    plt.plot([0,1], [0,1], color='red', markersize=0, linestyle=':')
    ax = plt.gca()
    ax.set_xlabel('Prob using all lnP')
    ax.set_ylabel('Prob using {0:.3f} of lnP'.format(sample_frac))
    ax.set_title('Probability from '+function_name)
    plt.tight_layout()

    fig.savefig('test_lnp_sampling.pdf')
    plt.close(fig)
