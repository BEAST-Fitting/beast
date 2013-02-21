import numpy as np
import numexpr
import tables
import os
from sedfitter import grid
from sedfitter import eztables
from sedfitter import vega
from sedfitter import progressbar
from sedfitter import anased
from sedfitter import creategrid
from sedfitter import extinction
from sedfitter import stellib
from ezisoch import ezIsoch
from ezpipeline import RequiredFile
from figs import plot_keys


def create_sed_grid(sed_grid_fname, filter_names, iso_fname='iso.proposal.fits', spectral_grid_fname='spectral.iso.proposal.fits'):
    """
     a spectral grid will be generated using the stellar parameters by
     interpolation of the isochrones and the generation of spectra into the
     physical units

     a photometric grid precomputing attenuation values as well
        sed_grid_fname = spectral_grid_fname.replace('spectral', 'seds')
    """
    #select the stellar library and the isochrones to use
    osl = stellib.Kurucz()
    #oiso = isochrone.padova2010()
    oiso = ezIsoch(iso_fname)
    extLaw = extinction.RvFbumpLaw()

    #filter_names = 'hst_wfc3_f275w hst_wfc3_f336w hst_acs_wfc_f475w hst_acs_wfc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split()

    # grid spacing for stars
    __tiny_delta__ = 0.001
    # ages           = (1e7, 1e8, 1e9)
    # masses         = (1., 2., 3., 4., 50.)
    # Z              = (0.02, 0.004)
    ages           = 10 ** np.arange(6., 9. + __tiny_delta__, 0.1)
    masses         = 10 ** np.arange(0.5, 20 + __tiny_delta__, 0.1)
    Z              = [(0.02)]
    # grid spacing for dust
    # variable to ensure that range is fully covered in using np.arange
    avs            = np.arange(0.0, 5.0 + __tiny_delta__, 0.1)
    rvs            = np.arange(1.0, 6.0 + __tiny_delta__, 0.5)
    fbumps         = np.asarray([1.0])
    #fbumps         = np.arange(0.0, 1. + __tiny_delta__, 0.1)
    #avs            = np.arange(0.0, 5.0 + __tiny_delta__, 1.)
    #rvs            = np.arange(1.0, 6.0 + __tiny_delta__, 1.)
    #fbumps         = np.arange(0.0, 1. + __tiny_delta__, 0.5)

    #filter extrapolations of the grid with given sensitivities in logg and logT
    bounds = dict(dlogT=0.1, dlogg=0.3)

    #make the spectral grid
    creategrid.gen_spectral_grid_from_stellib(spectral_grid_fname, osl, oiso, ages=ages, masses=masses, Z=Z, bounds=bounds)

    # make the grid
    extgrid = creategrid.make_extinguished_grid(spectral_grid_fname, filter_names, extLaw, avs, rvs, fbumps)

    # save grid to file
    extgrid.write(sed_grid_fname, clobber=True)


def getFakeStar(g, idx, err=0.1, noise=True, err_in_mag=False):
    """ Generate a fake sed from a model grid with some Gaussian noise
    INPUTS:
        idx         int                     index number on the grid
    OUTPUTS:
        fakein      int                     the index of the model on the grid
        fakesed     ndarray[float, ndim=1]  resulting SED

    KEYWORDS:
        err         float                   proportional error to report on the fluxes
        noise       bool                    if set, add proportional gaussian noise to the fluxes according to err
        err_in_mag  bool                    if set, err will be considered in mag

    TODO: add an empirical law, err = function(flux,...)
    """
    fakein   = idx
    fakesed  = np.copy(g.seds[fakein, :])
    if err_in_mag:
        fakeerr = fakesed * ( 1. - np.power(10, -0.4 * err) )
        if noise:  # sampled noise (with a mean of 1 and standard deviation of err (input value)
            fakesed *= (1. - np.power(10, -0.4 * np.random.normal(1.0, err, len(fakesed)) ))
    else:
        fakeerr = err * fakesed
        if noise:  # sampled noise (with a mean of 1 and standard deviation of err (input value)
            fakesed *= np.random.normal(1.0, err, len(fakesed))
    return fakein, fakesed, fakeerr


def generate_fake_observations( Nstars, sedgrid,
                               distance_modulus=24.48, noise=0.05,
                               logM=1.0, logT=3.6, dlogM=0.1, dlogT=0.1,
                               Fcut=('F814W', 3e-18), idx_Fcut=3):
    """ Generate a fake set of observations at a given distance with some noise
    and selections """
    if (distance_modulus is None):
        print "Assuming absolute fluxes for the output seds"

    # set the models to M31's distance
    # MF: what is the point to use this for sensitivity tests?
    #     I keep it for the general template
    if (distance_modulus is not None):
        sedgrid.seds *= 10 ** ( -0.4 * distance_modulus )

    # select models from sedgrid where
    # |logT - grid.logT| < dlogT & |logM - grid.logM| < dlogM & grid.sed[:, idx_Fcut] > Fcut
    flux_sel = (sedgrid.seds[:, idx_Fcut] > Fcut[1])
    M_sel    = (np.abs(sedgrid.logM - logM) < dlogM)
    T_sel    = (np.abs(sedgrid.logT - logT) < dlogT)
    subsel   = np.where( flux_sel & M_sel & T_sel )[0]
    assert(len(subsel) > 0), 'nothing to select'
    # Keep only Nstars random draws from them
    fakein = subsel[np.random.randint(low=0, high=subsel.size, size=Nstars)]
    del subsel

    #TODO: return an observation class instance / write on disk the data
    #tab = eztables.Table(sedgrid.grid)
    #fake_obs = FakeObs(tab)

    return fakein


#==========================================================================================
#
#
#
#==========================================================================================
def fit_model_seds_pytables(fakein, n_test, ext_grid, filters=None, err=0.1, threshold=-40, outname='test1.hf5'):

    """
    Fit model seds with noise for sensitivity tests
    INPUTS:
        fakein            list      indices of fake stars in ext_grid
        n_test            int       number of noise realizations per model
        ext_grid          SEDgrid   stellar model SEDs (luminosities)
    KEYWORDS:
        filters           list      indices of columns to use in grid.sed
        err               float     fractional noise to assume
        threshold         float     toss out grid points where lnp - lnp_max < threshold
        outname           string    output file directory for results

    TODO: we can either compact the writing: instead of groups and arrays, do tables
    or use eztables to do so.
    HDF is not a bad choice since it compress and keep all in one file!
    """

    if filters is None:
        filters = ext_grid.filters

    #ext_grid.seds = ext_grid.seds[:, filters]
    #ext_grid.lamb = ext_grid.lamb[:, filters]

    # mask for non-detectors (currently none)
    #mask = np.zeros(len(ext_grid.lamb), dtype=bool)

    #Use the mask as filter selection
    mask = np.array([ (False if k in filters else True) for k in ext_grid.filters]).astype(bool)

    outfile = tables.openFile(outname, 'w')
    #Save wavelengths in root, remember #n_stars = root._v_nchildren -1
    outfile.createArray(outfile.root, 'waves', ext_grid.lamb)
    N = len(fakein)
    with progressbar.PBar(N, txt="Calculating lnp") as pbar:
        for tn in range(N):
            star_group = outfile.createGroup('/', 'fakeStar_%d'  % tn, title="Fake star %d" % tn)
            star_group._v_attrs.fakein = fakein[tn]
            for tt in range(n_test):
                #fake DATA
                idx, fakesed, fakeerr = getFakeStar(ext_grid, fakein[tn], err=err)
                lnp = anased.computeLogLikelihood(fakesed, fakeerr, ext_grid.seds, normed=False, mask=mask)
                fake_group = outfile.createGroup(star_group, 'fake_%d' % tt)
                fake_group._v_attrs.test_n = tt
                #Need ragged arrays rather than uniform table
                outfile.createArray(fake_group, 'fakesed', fakesed)
                outfile.createArray(fake_group, 'fakerr', fakeerr)
                indx = np.where((lnp - max(lnp)) > -40.)
                outfile.createArray(fake_group, 'idx', np.array(indx[0], dtype=np.int32))
                outfile.createArray(fake_group, 'lnp', np.array(lnp[indx[0]], dtype=np.float32))
            outfile.flush()
            pbar.update(tn)
    outfile.close()


def expectation_values_pytables(grid, keys=['Av', 'Rv', 'f_bump'], filepath='Tests/fake_many_0/fakestars_0.hf5', method='expectation'):
    """
    args:
         grid        Model parameters. (model_grid.grid)
    **kwargs:
         keys        Parameters to calculate expectation values for.
         filepath    File log likelihoods have been saved to.
    """
    #make sure keys are real keys
    for key in keys:
        if not (key in grid.keys()):
            raise KeyError('Key "%s" not recognized' % key)
    outfile = tables.openFile(filepath, 'r')
    n_stars = outfile.root._v_nchildren - 1  # nchildren - 1 since wavelength also saved there
    m_tests = outfile.root.fakeStar_0._v_nchildren

    # Mean of m_tests realizations for each star
    means = np.zeros((n_stars, len(keys)), dtype=np.float32)
    #Standard deviation of m_tests realizations for each star
    stds = np.zeros((n_stars, len(keys)), dtype=np.float32)
    fakeinds = np.zeros(n_stars, dtype=np.int32)
    temp_vals = np.zeros((m_tests, len(keys)), dtype=np.float32)

    with progressbar.PBar(n_stars, txt="E(key)") as pbar:
        for i in range(n_stars):
            for j in range(m_tests):
                if (method == 'expectation'):
                    #the usage of numexpr here is not necessary. (Not to mention tables already has a better version built in
                    prob = numexpr.evaluate('exp(lnp)', local_dict={'lnp': outfile.getNode('/fakeStar_%d/fake_%d/lnp' % (i, j)).read()})
                    indx = outfile.getNode('/fakeStar_%d/fake_%d/idx' % (i, j)).read()
                    for e, key in enumerate(keys):
                        temp_vals[j, e] = numexpr.evaluate('sum(prob*vals)', local_dict={'vals': grid.getCol(key)[indx], 'prob': prob})
                    temp_vals[j, :] /= prob.sum()
                elif (method == 'maxprob'):
                    lnp = outfile.getNode('/fakeStar_%d/fake_%d/lnp' % (i, j)).read()
                    indx = outfile.getNode('/fakeStar_%d/fake_%d/idx' % (i, j)).read()
                    sel = grid.getRow(indx[lnp.argmax()])
                    for e, key in enumerate(keys):
                        temp_vals[j, e] = sel[key]

            fakeinds[i] = outfile.getNode('/fakeStar_%d' % i)._v_attrs.fakein
            means[i] = temp_vals.mean(axis=0)
            stds[i] = temp_vals.std(axis=0)
            pbar.update(i)
    outfile.close()
    summary_tab = eztables.Table(name="Summary Table")
    for e, key in enumerate(keys):
        summary_tab.add_column(key + '_recovered', means[:, e])
        summary_tab.add_column(key + '_std', stds[:, e])
        summary_tab.add_column(key, grid.data[key][fakeinds])
    return summary_tab


#def gen_summary_tables(gext_grid, keys, outdir, outnames):
#    for i in range(len(outnames)):
#        summary_table = expectation_values_pytables(gext_grid, keys=keys, filepath=outdir + outnames[i] + '.hf5')
#        summary_table.write(outdir + 'summary_' + outnames[i] + '.fits', clobber=True, append=False)
#        del summary_table


def gen_summary_tables(gext_grid, keys, outdir, outname):
    _outdir = outdir
    if _outdir[-1] != '/':
        _outdir += '/'
    summary_table = expectation_values_pytables(gext_grid, keys=keys, filepath=_outdir + outname + '.hf5')
    summary_table.write(_outdir + 'summary_' + outname + '.fits', clobber=True, append=False)
    del summary_table


def main():
    """
    The current script creates a set of models from PADOVA isochrones and
    CK04 spectral library for a sampling in age, mass, Av, Rv, fbump
    then generates a set of fake observations by picking random star models
    for which Normal noise is added. A magnitude cut is applied at M_814W =
    -3 to the fake data to simulate observational limits.
    finally the script fits the fake data and do some plots.
    """

    #### Generate fake Nstars observations
    #TODO: No need for distance, or distance modulus, we work in abs mags!

    #Pick Nstars fake stars within delta_logParam of logParam
    Nstars = 300
    m31_distance_modulus = 24.48  # in magnitudes
    err = 0.05  # Photometric noise to apply

    filter_names = 'hst_wfc3_f275w hst_wfc3_f336w hst_acs_wfc_f475w hst_acs_wfc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split()

    #Star selection in Mass, Temperature and luminosity
    logM = 1.0
    logT = 3.6
    delt_logM = 0.1
    delt_logT = 0.1
    Fcut0 = ('HST_ACS_WFC_F814W', -3)  # Filter and absolute magnitude converted later

    iso_fname = 'iso.proposal.fits'
    spectral_grid_fname = 'spectral.iso.proposal.fits'
    stellar_filename = 'fbump_only.fits'

    #Convert cut into apparent flux and not in Vega mag
    with vega.Vega() as v:
        r = Fcut0[1] + float(v.getMag([Fcut0[0]])[1]) + m31_distance_modulus
    Fcut = ('HST_ACS_WFC_F814W', 10 ** (-0.4 * r ))  # filter and absolute

    ### Get model grid of SEDs (assuming all filters are present)
    ### If stellar_filename does not exists, run the grid generation function
    with RequiredFile(stellar_filename,
                      create_sed_grid,
                      stellar_filename, filter_names,
                      iso_fname=iso_fname,
                      spectral_grid_fname=spectral_grid_fname) as stellar_filename:

        gext = grid.FileSEDGrid(stellar_filename)

    # find the column in the seds that corresponds to F814W
    assert(hasattr(gext, 'filters')), "The grid {} does not have a filter definition attribute. \nPlease re-run the generation or add the attribute.".format(stellar_filename)
    idx_Fcut = int(np.where(np.asarray(gext.filters) == 'HST_ACS_WFC_F814W'.upper())[0])
    # replaced by the above assertion
    #if hasattr(gext, 'filters'):  # if new grid version it should contain this attr
    #    idx_Fcut = int(np.where(np.asarray(gext.filters) == 'HST_ACS_WFC_F814W'.upper())[0])
    #else:
    #    idx_Fcut = 3

    fakein = generate_fake_observations( Nstars, gext, distance_modulus=m31_distance_modulus,
                                        noise=err, logM=logM, logT=logT, dlogM=delt_logM,
                                        dlogT=delt_logT, Fcut=Fcut, idx_Fcut=idx_Fcut)

    ######### Generate the likelihoods of each fake star and set of filters
    ntests = 2  # number of MC per star (Noise sampling)
    keys = ['Av', 'Rv', 'logT', 'logM', 'logA', 'logL']
    outdir = 'Tests'
    # set filters for each tests
    #filters = [np.arange(6), [0, 1], [2, 3]]  # All six PHAT filters,  F275W and F336W,  F475W and F814W
    filters = [ 'hst_wfc3_f275w hst_wfc3_f336w hst_acs_wfc_f475w hst_acs_wfc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split(),
                'hst_wfc3_f275w hst_wfc3_f336w'.upper().split(),
                'hst_wfc3_f475w hst_wfc3_f814w'.upper().split() ]

    outnames = ['all_phat', 'F275_336', 'F475_814']

    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise Exception('Output directory "{}" already exists but is not a directory'.format(outdir))
    else:
        os.mkdir(outdir)

    for i in range(len(filters)):
        _outname = outdir + '/' + outnames[i] + '.hf5'
        fit_model_seds_pytables(fakein, ntests, gext, err=err, filters=filters[i], outname=_outname)

    with progressbar.PBar(len(outnames), txt="Generating sumary tables") as pbar:
        for e, ok in enumerate(outnames):
            pbar.update(ok)
            gen_summary_tables(gext.grid, keys, outdir, ok)  # generate and save summary tables

    del gext

    print "producing plots"
    plot_keys(keys, outdir, outnames)
