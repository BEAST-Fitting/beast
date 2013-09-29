"""
Make BEAST paper plots related to extinction
"""

from core import extinction
import pylab
import numpy

dir = './figs/'

x = (numpy.arange(100) / 100.) * 10. + 0.1
lamb = 1.e4/x    # wavelength in angstroms

mixlaw = extinction.RvFbumpLaw()

#BEAST paper extinction mixture law figure Rv_A=2.5 
fig = pylab.figure()
plt = fig.add_subplot(111)

ymix = mixlaw(lamb, Rv_A=2.5, f_bump=0.)
plt.plot(x, ymix,'-k', label='$f_A=0.0$')

ymix = mixlaw(lamb, Rv_A=2.5, f_bump=1./3.)
plt.plot(x, ymix,'--k', label='$f_A=0.333$')
    
ymix = mixlaw(lamb, Rv_A=2.5, f_bump=2./3.)
plt.plot(x, ymix, ':k', label='$f_A=0.667$')

ymix = mixlaw(lamb, Rv_A=2.5, f_bump=1.0)
plt.plot(x, ymix,'-.k', label='$f_A=1.0$')

pylab.figtext( 0.16,0.63, '$R_A(V)\  =\  2.5$', color='k',size=15)

pylab.xlim(0.,8.)
pylab.ylim(0.,7.)
pylab.rc('font', size=13, family='serif')
plt.set_ylabel('$A(\lambda)/A(V)$',size=15)
plt.set_xlabel('$1/\lambda$ [$\mu m^{-1}$]',size=15)

plt.legend(loc=0, frameon=False, fontsize = 15)
    
pylab.savefig(dir+'RvFbumpLaw_25.eps')
pylab.savefig(dir+'RvFbumpLaw_25.pdf')

    
#BEAST paper extinction mixture law figure Rv_A=5.5 
fig = pylab.figure()
plt = fig.add_subplot(111)

ymix = mixlaw(lamb, Rv_A=5.5, f_bump=0.)
plt.plot(x, ymix,'-k', label='$f_A=0.0$')

ymix = mixlaw(lamb, Rv_A=5.5, f_bump=1./3.)
plt.plot(x, ymix,'--k', label='$f_A=0.333$')
    
ymix = mixlaw(lamb, Rv_A=5.5, f_bump=2./3.)
plt.plot(x, ymix, ':k', label='$f_A=0.667$')

ymix = mixlaw(lamb, Rv_A=5.5, f_bump=1.0)
plt.plot(x, ymix,'-.k', label='$f_A=1.0$')

#pylab.figtext( 0.16,0.63, '$R_A(V)\  =\  5.5$', color='k',size=15)

pylab.xlim(0.,8.)
pylab.ylim(0.,7.)
pylab.rc('font', size=13, family='serif')
plt.set_ylabel('$A(\lambda)/A(V)$',size=15)
plt.set_xlabel('$1/\lambda$ [$\mu m^{-1}$]',size=15)

plt.legend(loc=0, frameon=False, fontsize = 15)
    
pylab.savefig(dir+'RvFbumpLaw_55.eps')
pylab.savefig(dir+'RvFbumpLaw_55.pdf')

#extinction mixture law figure 2 components 
fig = pylab.figure()
plt = fig.add_subplot(111)

f99  = extinction.Fitzpatrick99()
gsmc = extinction.Gordon03_SMCBar()

ymix = mixlaw(lamb, Rv_A=3.1, f_bump=0.5)
compA = 0.5*f99(lamb, Rv=3.1) 
compB = 0.5*gsmc(lamb)

plt.plot(x,ymix,'-k',label='Total')
plt.plot(x,compA,'--k',label='Bump Law')
plt.plot(x,compB,':k',label='No Bump')
#plt.plot(x,compA+compB,'-r')
plt.set_ylabel('$A(\lambda)/A(V)$',size=15)
plt.set_xlabel('$1/\lambda$ [$\mu m^{-1}$]',size=15)

pylab.figtext( 0.16,0.695, '$f_A\  =\  0.5$', color='k',size=15, fontweight='bold')

plt.legend(loc=0, frameon=False, fontsize = 15)

pylab.savefig(dir+'MixtLawComp.eps')
pylab.savefig(dir+'MixtLawComp.pdf')

#BEAST paper f_bump-Rv grid figure

Rv_vals=numpy.arange(2.0,6.01,0.08)
fb_vals=numpy.arange(0.0,1.02,0.02)
N_fb = len(fb_vals)
N_rv = len(Rv_vals)
Rv_res = numpy.ndarray(shape=(N_rv,N_fb), dtype=float)
ind = numpy.ndarray(shape=(N_rv,N_fb), dtype=float)
fig = pylab.figure()
plt = fig.add_subplot(1, 1, 1)
for i in range(N_fb):
    Rv_res[0:N_rv,i]=mixlaw.get_Rv_A(Rv_vals,fbump=fb_vals[i]) 
    ind = numpy.where((Rv_res[:,i] >= 2.) &  (Rv_res[:,i] <= 6.001))[0]
    ll=len(ind)
    plt.plot(numpy.ones(ll)*fb_vals[i],Rv_vals[ind],'ko',markersize=2)
        
pylab.xlim(-0.05,1.05)
pylab.ylim(1.9,6.1)
pylab.rc('font', size=13, family='serif')
plt.set_xlabel('$f_A$',size=15)
plt.set_ylabel('R($V$)',size=15)
    
pylab.savefig(dir+'RvFbump_grid.eps')
pylab.savefig(dir+'RvFbump_grid.pdf')

pylab.show()
