import numpy as np
import pylab as plt
import grid
import eztables
import tables
from observations import Observations
import vega
import progressbar
import anased
import numexpr


#===============================================================================
# HANDLE THE OBSERVATION DATA
#   *real* or *fake* we need to handle the inputs in a same manner
#   so that it becomes a minimal effort to change the dataset
#
#
# Class FakeObs is derived from observations.Observations and aims at being
# as general as possible
#
# Real data will only require to define something very similar (see phat.py)
#==============================================================================
class FakeObs(Observations):
    """ A quick insertion of the class that will be used eventually
        so that all the code can be used as a template for real data

        * This class reads either a file or a variable
          thanks to eztable.Table construtor flexibility.

        * Instances can be saved with calling .writeto()
          (it uses Table.write)
          TODO: save the mask and filter names
    """
    def __init__(self, fname):
        """ Generate a data interface object """
        self.filters   = None
        self.desc      = None
        self.badvalue  = None
        self.setDistanceModulus(0.)
        self.readData(fname)

    def getMags(self, num, filters):
        """ TODO update """
        return np.array([ self.data[tt][num] for tt in filters])

    def getErrors(self, num, filters):
        """ TODO update """
        return np.array([ self.data[tt + 'err'][num] for tt in filters])

    def getObs(self, num=0):
        #assert ( self.filters is not None), "No filter set."
        mags = self.getMags(num, self.filters)
        errs = np.ones(len(mags), dtype=float) * err
        if self.badvalue is not None:
            mask = (mags >= self.badvalue)
        else:
            mask = np.zeros(len(mags), dtype=bool)

        return mags, errs, mask

    def readData(self, fname):
        """ read the dataset from the original source file """
        from eztables import Table
        if type(fname) == str:
            self.inputFile = fname
        else:
            self.inputFile = None
        self.data = Table(fname)

    def writeto(self, *args, **kwargs):
        self.data.write(*args, **kwargs)


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


def generate_fake_fake_observations( Nstars, sedgrid,
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
    flux_sel = (sedgrid.seds[:, idx_Fcut] > Fcut)
    M_sel    = (np.abs(sedgrid.logM - logM) < dlogM)
    T_sel    = (np.abs(sedgrid.logT - logM) < dlogT)
    subsel   = np.where( flux_sel & M_sel & T_sel )[0]
    # Keep only Nstars random draws from them
    fakein = subsel[np.random.randint(low=0, high=subsel.size, size=Nstars)]
    del subsel

    #TODO: return an observation class instance / write on disk the data
    tab = eztables.Table(sedgrid.grid)
    fake_obs = FakeObs(tab)

    return fakein, fake_obs


#==========================================================================================
#
#
#
#==========================================================================================
def fit_model_seds_pytables(fakein, n_test, ext_grid, filters=None, err=0.1, threshold=-40, outname='Tests/fake_many_0/test1.hf5'):

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

    filters = filters or np.arange(6)

    ext_grid.seds = ext_grid.seds[:, filters]
    ext_grid.lamb = ext_grid.lamb[:, filters]

    # mask for non-detectors (currently none)
    mask = np.zeros(len(ext_grid.lamb), dtype=bool)
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


def gen_summary_tables(gext_grid, keys, outdir, outnames):
    for i in range(len(outnames)):
        summary_table = expectation_values_pytables(gext_grid, keys=keys, filepath=outdir + outnames[i] + '.hf5')
        summary_table.write(outdir + 'summary_' + outnames[i] + '.fits', clobber=True, append=False)
        del summary_table


def plot_keys(keys, outdir, outnames, show=True):
    """ Plot average offset per unique values in keys
    INPUTS:
        keys        list    columns to process
        outdir      str     directory where to find the summary fits tables
        outnames    list    list of names for tables

    TODO: replace inputs by a list of files and list of keys
    """
    for i in range(len(outnames)):
        summary_table = eztables.Table(outdir + 'summary_' + outnames[i] + '.fits')
        for j, key in enumerate(keys):
            #rec_vals = summary_table.data[key+'_recovered']  # recovered values
            true_vals = summary_table.data[key]              # true values
            #rec_vals = rec_vals - true_vals                  # offset
            rec_vals = summary_table.evalexpr('{}_recovered - {}'.format(key, key))  # offset
            uniq_vals = np.unique(true_vals)  # unique true values
            avg_offs = np.zeros(uniq_vals.size)  # Mean of recovered params for given input param
            for k in range(avg_offs.size):
                sel = np.where(true_vals == uniq_vals[k])
                avg_offs[k] = rec_vals[sel].mean()
            #figure
            plt.figure(j)
            ax = plt.subplot(111)
            ax.plot(uniq_vals, avg_offs, label=(outnames[i]).replace('_', ' '), marker='+', markersize=10)

    for j, key in enumerate(keys):
        plt.figure(j)
        ax = plt.subplot(111)
        ax.set_xlabel(key.replace('_', ''))
        ax.set_ylabel(key.replace('_', '') + ' (out - in)')
        ax.legend(loc='best')
        ax.plot(ax.get_xlim(), [0, 0], linestyle='dotted', color='black')
    if show:
        plt.show()


if __name__ == '__main__':
    #### Generate fake Nstars observations
    #TODO: No need for distance, or distance modulus, we work in abs mags!

    #Pick Nstars fake stars within delta_logParam of logParam
    Nstars = 300
    m31_distance_modulus = 24.48  # in magnitudes
    err = 0.05  # Photometric noise to apply

    #Star selection in Mass, Temperature and luminosity
    logM = 1.0
    logT = 3.6
    delt_logM = 0.1
    delt_logT = 0.1
    Fcut0 = ('HST_ACS_WFC_F814W', -3)  # Filter and absolute magnitude converted later

    stellar_filename = 'fbump_only.fits'

    #Convert cut into apparent flux and not in Vega mag
    with vega.Vega() as v:
        r = Fcut0[1] + float(v.getMag([Fcut0[0]])[1]) + m31_distance_modulus
    Fcut = ('HST_ACS_WFC_F814W', 10 ** (-0.4 * r ))  # filter and absolute
    idx_Fcut = 3  # TODO: get this automatically from the grid

    #Get model grid of SEDs (assuming all filters are present)
    gext = grid.FileSEDGrid(stellar_filename)

    fakein = generate_fake_fake_observations( Nstars, gext, distance_modulus=m31_distance_modulus,
                                             noise=err, logM=logM, logT=logT, delt_logM=delt_logM,
                                             delt_logT=delt_logT, Fcut=Fcut, idx_Fcut=idx_Fcut)

    ######### Generate the likelihoods of each fake star and set of filters
    ntests = 2  # number of MC per star (Noise sampling)
    keys = ['Av', 'Rv', 'logT', 'logM', 'logA', 'logL']
    filters = [np.arange(6), [0, 1], [2, 3]]  # All six PHAT filters,  F275W and F336W,  F475W and F814W
    outdir = 'Tests/'
    outnames = ['all_phat', 'F275_336', 'F475_814']

    for i in range(len(filters)):
        _outname = outdir + outnames[i] + '.hf5'
        fit_model_seds_pytables(fakein, ntests, gext, err=err, filters=filters[i], outname=_outname)

    gen_summary_tables(gext.grid, keys, outdir, outnames)  # generate and save summary tables
    del gext

    plot_keys(keys, outdir, outnames)
