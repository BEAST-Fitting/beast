"""
BEAST's running pipeline

CONTENT:

    Imports.................................[sec:imports]
        Needed packages incl. BEAST imports

    Globals.................................[sec:globals]
        User inputs and Data file

    Data interface.............................[sec:data]
        Deriving Observation class
            The class will take care of converting
            vegamags to fluxes, and bad value limits
        file column mappings

    Model definitions..........................[sec:model]

        Main ingredients.................[sec:models:main]
            define extinction, stellib, isochrones to use

        Grid sampling definition.........[sec:models:grid]
            define model grid sampling

        Grid generation................[sec:models:create]
            create the grid if it does not exist
            according to all the assumptions

     Fitting.....................................[sec:fit]
        fitting routines

     Pipeline................................[sec:pipeline]
        fast aliases checking the previous steps
        all in one call (assuming all params are ok)
"""

#---------------------------------------------------------
# Imports                                    [sec:imports]
#---------------------------------------------------------
import numpy as np
import tables
from beast.core.observations import Observations
from beast.core.isochrone import ezIsoch
from beast.core.vega import Vega, from_Vegamag_to_Flux
from beast.tools.ezpipeline import RequiredFile
from beast.core.anased import computeLogLikelihood
from beast.core import grid
from beast.core import creategrid
from beast.core import stellib
from beast.core import extinction
from beast.tools import progressbar

#---------------------------------------------------------
# Globals                                    [sec:globals]
#---------------------------------------------------------

# MF Comments on user input

# I don't think you need a ascii file. Especially when you impose an order of
# inputs. It was not obvious to me that \t was the separator!
#
# However, you should define a class that reads the whatever format you want to
# adopt. This way if anyone wants a different format, ascii, xml, python, yaml
# etc, you can change a few lines somewhere instead and it still work like a
# charm. This is Object oriented coding's magic!

# What I think is that python is the easiest input format, especially because
# the observation class needs to be rederived everytime you bring a new data
# input. You have to make distance, filters, mapping with the columns etc.

# also there are missing parameters such as the extinction law

#Read user's inputs
#cleaned up, reading the file only once
with open("input.txt","r") as infile:
    ## Put all user's inputs in InputList in the following order:
    ## project, subproject, observation file, filter list, distance modulus,
    ## grid parameters [Av, Rv, fb, age, mass, Z]
    # InputList = []
    # for e, tt in enumerate(infile):
    #   if tt[0] != '#':
    #        InputList.append(tt.split("\t")[0])

    # even more python
    InputList = [ tt.split("\t")[0] for tt in infile if tt[0] != '#' ]

project, subproject, obsfile, filters, distanceModulus = InputList[:4]

# Use a shortcut when using the 6 PHAT filters
if filters == 'All PHAT':
    filters = 'hst_wfc3_f275w hst_wfc3_f336w hst_acs_wfc_f475w hst_acs_wfc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split()

distanceModulus = float(InputList[4])
Av_min = float(InputList[5])
Av_max = float(InputList[6])
Av_step = float(InputList[7])
Rv_min = float(InputList[8])
Rv_max = float(InputList[9])
Rv_step = float(InputList[10])
fb_min = float(InputList[11])
fb_max = float(InputList[12])
fb_step = float(InputList[13])

logA_min = float(InputList[14])
logA_max = float(InputList[15])
logA_step = float(InputList[16])
logM_min = float(InputList[17])
logM_max = float(InputList[18])
logM_step = float(InputList[19])
Z = float(InputList[20])

#Isochrones, grid and result files
isofile = 'beast/libs/iso.proposal.fits'
lnp_outname = 'lnp.{0}.{1}.hd5'.format(project,subproject)
stat_outname = 'stat.{0}.{1}.hd5'.format(project,subproject)
res_outname = 'res.{0}.{1}.fits'.format(project,subproject)
spectral_grid_fname = 'coarse.spectral.grid.fits'
sed_grid_fname = 'medium.coarse.sed.grid.fits'

# MF: so I would only do a "from user_input import ..." to have the same effect
# while still providing easy and explicit inputs that the user cannot mess up

# However, at this point there is not check on the inputs or missing/wrong
# value which could be problematic

#---------------------------------------------------------
# Data interface                                [sec:data]
#---------------------------------------------------------

#Data are in Vega magnitudes
#  Need to use Vega
with Vega() as v:
    vega_f, vega_mag, lamb = v.getMag(filters)


# derive the global class and update what's needed
class Data(Observations):
    """ PHAT catalog for clusters in M31 """
    def __init__(self, inputFile, distanceModulus=distanceModulus):
        desc = 'PHAT star: %s' % inputFile
        Observations.__init__(self, inputFile, distanceModulus, desc=desc)
        self.setFilters( filters )
        self.setBadValue(50.0)  # some bad values smaller than expected
        self.minError = 0.001
        self.floorError = 0.05  # constant error term

    @from_Vegamag_to_Flux(lamb, vega_mag)
    def getObs(self, num):
        """ Using the decorator @from_Vegamag_to_Flux
            Hence, results are in flux (not in flux/flux_vega)
            Returns the fluxes, errors and mask of an observation.
        """
        return Observations.getObs(self, num)

    def getObsinMag(self, num):
        """ Returns the original catalog magnitudes """
        return Observations.getObs(self, num)

    def getErrors(self, num, filters):
        """ Redifined to impose a minimal error """
        err = np.array([ self.data[tt + 'err'][num] for tt in filters])
        if self.floorError > 0.:
            err = np.sqrt(err ** 2 + self.floorError ** 2)
        if self.minError > min(err):
            err[ err < self.minError ] = self.minError
        return err

obs = Data(obsfile, distanceModulus)
obs.setFilters(filters)
#Change the filters name to match
for k in filters:
    obs.data.set_alias(k, k.split('_')[-1] + '_mag')
    obs.data.set_alias(k + 'err', k.split('_')[-1] + '_magerr')

#---------------------------------------------------------
# Model definitions                            [sec:model]
#---------------------------------------------------------

# Main ingredients                       [sec:models:main]
#---------------------------------------------------------

extLaw = extinction.RvFbumpLaw()
osl = stellib.Kurucz()
oiso = ezIsoch(isofile)

# Grid sampling definition               [sec:models:grid]
#---------------------------------------------------------

# variable to ensure that range is fully covered in using np.arange
__tiny_delta__ = 0.001

ages = 10 ** np.arange(logA_min, logA_max + __tiny_delta__, logA_step)
masses = 10 ** np.arange(logM_min, logM_max + __tiny_delta__, logM_step)
Z = np.asarray([Z])
avs = np.arange(Av_min, Av_max + __tiny_delta__, Av_step)
rvs = np.arange(Rv_min, Rv_max + __tiny_delta__, Rv_step)
fbumps = np.arange(fb_min, fb_max + __tiny_delta__, fb_step)

griddef = (ages, masses, Z, avs, rvs, fbumps)


# Grid generation                      [sec:models:create]
#---------------------------------------------------------
def create_sed_grid(sed_grid_fname=sed_grid_fname,
                    filter_names=filters,
                    oiso=oiso,
                    osl=osl,
                    extLaw=extLaw,
                    spectral_grid_fname=spectral_grid_fname,
                    griddef=griddef):
    """
     a spectral grid will be generated using the stellar parameters by
     interpolation of the isochrones and the generation of spectra into the
     physical units

     a photometric grid precomputing attenuation values as well
        sed_grid_fname = spectral_grid_fname.replace('spectral', 'seds')
    """
    ages, masses, Z, avs, rvs, fbumps = griddef

    #filter extrapolations of the grid with given sensitivities in logg and logT
    bounds = dict(dlogT=0.1, dlogg=0.3)

    #make the spectral grid
    creategrid.gen_spectral_grid_from_stellib(spectral_grid_fname, osl, oiso, ages=ages, masses=masses, Z=Z, bounds=bounds)

    # make the grid
    extgrid = creategrid.make_extinguished_grid(spectral_grid_fname, filter_names, extLaw, avs, rvs, fbumps)

    # save grid to file
    extgrid.write(sed_grid_fname, clobber=True)


def get_sedgrid():

    with RequiredFile(sed_grid_fname, create_sed_grid, filter_names=filters,
                      oiso=oiso, osl=osl, extLaw=extLaw,
                      spectral_grid_fname=spectral_grid_fname, griddef=griddef) as __sed_grid_fname__:

            return grid.FileSEDGrid(__sed_grid_fname__)


#---------------------------------------------------------
# Fitting                                 [sec:fit]
#---------------------------------------------------------

def fit_model_seds_pytables(obs, sedgrid, threshold=-60, outname=lnp_outname):

    """
    Fit model seds with noise for sensitivity tests
    INPUTS:
        obs         Observtions     Observation object to analyze
        sedgrid     SEDgrid         stellar model SEDs (luminosities)

    KEYWORDS:
        threshold   float          toss out grid points where lnp - lnp_max < threshold
        outname     string         output file directory for results

    TODO: Clean up the disk output structures
    HDF is not a bad choice since it compress and keep all in one file!
    """
    filters = obs.getFilters()
    with tables.openFile(outname, 'w') as outfile:
        #Save wavelengths in root, remember #n_stars = root._v_nchildren -1
        outfile.createArray(outfile.root, 'grid_waves', sedgrid.lamb)
        outfile.createArray(outfile.root, 'obs_filters', filters)

        #loop over the obs and do the work
        #with progressbar.PBar(len(obs), txt="Calculating lnp") as pbar:
        #    for tn, (sed, err, mask) in obs.enumobs():
        # Restriction to 1000 obs
        with progressbar.PBar(1000, txt="Calculating lnp") as pbar:
            for tn, (sed, err, mask) in obs.enumobs():
                if tn > 1000:
                    break
                lnp = computeLogLikelihood(sed, err, sedgrid.seds, normed=False, mask=mask)

                #Need ragged arrays rather than uniform table
                star_group = outfile.createGroup('/', 'star_%d' % tn, title="star %d" % tn)
                indx = np.where((lnp - max(lnp[np.isfinite(lnp)])) > threshold)
                outfile.createArray(star_group, 'input', np.array([sed, err, mask]).T)
                outfile.createArray(star_group, 'idx', np.array(indx[0], dtype=np.int32))
                outfile.createArray(star_group, 'lnp', np.array(lnp[indx[0]], dtype=np.float32))
                #commit changes
                outfile.flush()

                pbar.update(tn, force=True)  # Forcing because it can be long to show the first ETA


#---------------------------------------------------------
# pipeline                                  [sec:pipeline]
#---------------------------------------------------------

def do_fit():
    sedgrid = get_sedgrid()
    fit_model_seds_pytables(obs, sedgrid, threshold=-20, outname=lnp_outname)


if __name__ == '__main__':
    do_fit()
