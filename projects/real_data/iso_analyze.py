import numpy as np
import numexpr
import tables
from beast.external.eztables.astro import AstroTable
from beast.core.vega import Vega, from_Vegamag_to_Flux
from beast.core.grid import FileSEDGrid
from beast.core.observations import Observations
from user_input_val import *
import matplotlib.pyplot as plt


#########  DATA ##########
#Data are in Vega magnitudes
#  Need to use Vega
with Vega() as v:
    vega_f, vega_mag, lamb = v.getMag(filters)
    vega_f, vega_flux, lamb = v.getFlux(filters)


# derive the global class and update what's needed
class Data2(Observations):
    """ PHAT catalog for clusters in M31 """
    def __init__(self, inputFile, distanceModulus=distanceModulus):
        desc = 'PHAT star: %s' % inputFile
        Observations.__init__(self, inputFile, distanceModulus, desc=desc)
        self.setFilters( filters )
        self.setBadValue(50.0)  # some bad values smaller than expected
        self.minError = 0.0001
        self.floorError = 0.0  # constant error term to match the IDL fitter
        
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


obs = Data2(obsfile, distanceModulus)
obs.setFilters(filters)

#Change the filters name to match
for k in filters:
    obs.data.set_alias(k, k.split('_')[-1] + '_MAG')
    obs.data.set_alias(k + 'err', k.split('_')[-1] + '_MAGERR')

###################################

### Grid 
sedgrid = FileSEDGrid("/astro/dust_kg/harab/beast/projects/real_data/GridRvFbumpLaw.sed.grid.fits")

### Expectation
#res = AstroTable(res_outname)

### Best Fits
bbf = AstroTable("test_BBF.iso10.GridRvFbumpLaw_8004.fits")

### Extracting Best models 
bestmod    = sedgrid.seds[bbf['Bestind']]
bestmodMag = -2.5 * np.log10(bestmod/vega_flux)
sedsMags=-2.5 * np.log10(sedgrid.seds/vega_flux)+distanceModulus
clrgrid=sedsMags[:,2]-sedsMags[:,3]


### Computing Expected models
#expectSED    = np.loadtxt("expectSED.iso10.test11_1000.txt")
#expectSEDmag = -2.5 * np.log10(expectSED/vega_flux) 


clr          = obs.data['F475W_MAG']-obs.data['F814W_MAG']
clerr        = obs.data['F475W_MAGERR']+obs.data['F814W_MAGERR']
bbfclr       = bestmodMag[:,2]-bestmodMag[:,3]
#expectSEDclr = expectSEDmag[:,2]-expectSEDmag[:,3]

ind = np.where(obs.data['AGE']==10.)
plt.ion()
plt.figure()
#plt.errorbar(clr[0:100],obs.data['F814W_MAG'][0:100],xerr=clerr[0:100],yerr=obs.data['F814W_MAGERR'][0:100],fmt='y.')
plt.plot(clrgrid,sedsMags[:,3],'go')
plt.plot(clr[0:8004],obs.data['F814W_MAG'][0:8004],'yo')
#plt.plot(bbfclr,bestmodMag[:,3]+distanceModulus,'ro')#,expectSEDclr,expectSEDmag[:,3]+distanceModulus,'bo')
#plt.plot(clr[indbad],obs.data['F814W_MAG'][indbad],'sy')
ax=plt.gca()
plt.xlabel("F475W - F 814W")
plt.ylabel("F814W")
ax.set_ylim(ax.get_ylim()[::-1])
#plt.xlim(0,5)
plt.figtext( 0.15,0.85, 'data', color='y')
plt.figtext( 0.15,0.825, 'grid', color='g')
#plt.figtext( 0.15,0.825, 'best fits', color='r')
#plt.figtext( 0.15,0.8, 'expected', color='b')
#plt.show()

plt.figure()
plt.errorbar(clr[0:8004],obs.data['F814W_MAG'][0:8004],xerr=clerr[0:8004],yerr=obs.data['F814W_MAGERR'][0:8004],fmt='y.')
plt.plot(bbfclr,bestmodMag[:,3]+distanceModulus,'r.')
ax=plt.gca()
plt.xlabel("F475W - F 814W")
plt.ylabel("F814W")
ax.set_ylim(ax.get_ylim()[::-1])
#plt.xlim(0,5)
plt.figtext( 0.15,0.85, 'data', color='y')
plt.figtext( 0.15,0.825, 'best fits', color='r')


print nn

plt.figure()
#plt.plot(np.log10(obs.data['MASS'][0:1000]),(res['logM'][0:1000]-np.log10(obs.data['MASS'][0:1000]))/np.log10(obs.data['MASS'][0:1000]),'bo')
plt.plot(np.log10(obs.data['MASS'][0:8004]),(bbf['logM']-np.log10(obs.data['MASS'][0:8004]))/np.log10(obs.data['MASS'][0:8004]),'ro')
plt.xlabel("log M")
plt.ylabel("Relative difference (Mass)")
#ax=plt.gca()
#ax.set_xlim(ax.get_xlim()[::-1])
#plt.show()

plt.figure()
#plt.plot(obs.data['LOGT'][0:1000],(res['logA'][0:1000] - np.log10(obs.data['AGE'][0:1000]*1e6)) / np.log10(1e6*obs.data['AGE'][0:1000]),'bo')
plt.plot(obs.data['LOGT'][0:8004],(bbf['logA'] - np.log10(obs.data['AGE'][0:8004]*1e6))/np.log10(obs.data['AGE'][0:8004]*1e6),'ro')
#plt.xlim(6.,7)
ax=plt.gca()
plt.xlabel("log T")
plt.ylabel("Relative difference (Age)")
ax.set_xlim(ax.get_xlim()[::-1])
#plt.ylim(-1.5.,2.5)
#plt.show()

plt.figure()
#plt.plot(obs.data['LOGT'][0:1000],(res['Av'][0:1000]-obs.data['AV'][0:1000])/obs.data['AV'][0:1000],'bo')
plt.plot(obs.data['LOGT'][0:8004],(bbf['Av']-obs.data['AV'][0:8004])/obs.data['AV'][0:8004],'ro')
ax=plt.gca()
plt.xlabel("log T")
plt.ylabel("Relative difference (Av)")
ax.set_xlim(ax.get_xlim()[::-1])
#plt.ylim([-1.5,2.5])
#plt.show()

plt.figure()
#plt.plot(np.log10(obs.data['LOGT'][0:1000]*1e6),(res['logT'][0:1000] - np.log10(obs.data['LOGT'][0:1000]*1e6)) / obs.data['LOGT'][0:1000],'bo')
plt.plot(obs.data['LOGT'][0:8004],(bbf['logT'] - obs.data['LOGT'][0:8004]) / obs.data['LOGT'][0:8004],'ro')
#plt.ylim(-1,1)
ax=plt.gca()
plt.xlabel("log T")
plt.ylabel("Relative difference (Temperature)")
ax.set_xlim(ax.get_xlim()[::-1])

plt.figure()
plt.plot(obs.data['LOGG'][0:8004],(bbf['logg'][0:8004] - obs.data['LOGG'][0:8004]) / obs.data['LOGG'][0:8004],'ro')
plt.xlabel("log g")
plt.ylabel("Relative difference (log g)")
#plt.show()
#print nn
import pylab as plb
for i in range(0,10):
    fig=plb.figure()
    plb.xscale('log')
    plb.xlabel('Wavelengths (A)')
    plb.ylabel('Flux (erg/s/cm^2/A)')
    plb.errorbar(lamb,obs.getObs(i)[0],yerr=obs.getObs(i)[1],fmt='go')
    plb.plot(lamb,bestmod[i],'r')
    #fig.savefig("figs/SED"+str(i)+"_GridRvFbumpLaw_fit.png")


    

