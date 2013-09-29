"""
Make BEAST paper plots related to extinction
"""
from core import extinction
import pylab
import numpy

pylab.rc('font', size=13, family='serif')


def savefig(name=None, fmt=['eps', 'pdf']):
    if name is not None:
        for fk in fmt:
            pylab.savefig('{}.{}'.format(savefig, fk), bbox_inches='tight')


def make_fig(lamb, law, Rv_A=2.5, ax=None, figname=None, save_fmt=['eps', 'pdf']):
    #BEAST paper extinction mixture law figure Rv_A=2.5
    if ax is None:
        fig = pylab.figure()
        ax = fig.add_subplot(111)

    ymix = mixlaw(lamb, Rv_A=2.5, f_bump=0.)
    ax.plot(x, ymix, '-k', label='$f_A=0$')

    ymix = mixlaw(lamb, Rv_A=2.5, f_bump=1. / 3.)
    ax.plot(x, ymix, '--k', label='$f_A=1/3$')

    ymix = mixlaw(lamb, Rv_A=2.5, f_bump=2. / 3.)
    ax.plot(x, ymix, ':k', label='$f_A=2/3$')

    ymix = mixlaw(lamb, Rv_A=2.5, f_bump=1.0)
    ax.plot(x, ymix, '-.k', label='$f_A=1$')

    ax.text(0.16, 0.63, '$R_A(V)$ = {:0.1f}'.format(Rv_A), color='k', size=15)
    ax.set_xlim(0., 8.)
    ax.set_ylim(0., 7.)
    ax.set_ylabel('$A(\lambda)/A(V)$', size=15)
    ax.set_xlabel('$1/\lambda$ [$\mu m^{-1}$]', size=15)
    ax.legend(loc=0, frameon=False, fontsize=15)

    savefig(figname, fmt=save_fmt)


def make_mixture_fig(f_A=0.5, ax=None, figname=None, save_fmt=['eps', 'pdf']):
    #extinction mixture law figure 2 components
    if ax is None:
        fig = pylab.figure()
        ax = fig.add_subplot(111)

    ymix = mixlaw(lamb, Rv_A=3.1, f_bump=f_A)
    compA = f_A * ymix.RvLaw(lamb, Rv=3.1)
    compB = f_A * ymix.NoBumpLaw(lamb)

    ax.plot(x,ymix, '-k', label='Mixture')
    ax.plot(x,compA, '--k', label='Bump Law')
    ax.plot(x,compB, ':k', label='No Bump')
    ax.set_ylabel(r'A($\lambda$)/A(V)', size=15)
    ax.set_xlabel(r'1/$\lambda$ [$\mu m^{-1}$]', size=15)

    ax.text(0.16, 0.695, '$f_A$ = {:0.1f}'.format(f_A), color='k', size=15, fontweight='bold')
    ax.legend(loc=0, frameon=False, fontsize=15)

    savefig(figname, fmt=save_fmt)


def Rv_fbump_grid(ax=None, figname=None, save_fmt=['eps', 'pdf']):
    if ax is None:
        fig = pylab.figure()
        ax = fig.add_subplot(111)

    mixlaw = extinction.RvFbumpLaw()
    Rv_vals = numpy.arange(2.0,6.01,0.08)
    fb_vals = numpy.arange(0.0,1.02,0.02)
    N_fb = len(fb_vals)
    N_rv = len(Rv_vals)

    Rv_res = numpy.ndarray(shape=(N_rv, N_fb), dtype=float)
    ind = numpy.ndarray(shape=(N_rv, N_fb), dtype=float)

    for i, fbk in enumerate(fb_vals):
        Rv_res[0: N_rv, i] = mixlaw.get_Rv_A(Rv_vals, fbump=fbk)
        ind = numpy.where((Rv_res[:, i] >= 2.) & (Rv_res[:, i] <= 6.001))[0]
        ll = len(ind)
        ax.plot(numpy.ones(ll) * fbk, Rv_vals[ind], 'ko', markersize=2)

    pylab.xlim(-0.05, 1.05)
    pylab.ylim(1.9, 6.1)
    ax.set_xlabel('$f_A$', size=15)
    ax.set_ylabel('R($V$)', size=15)

    savefig(figname, fmt=save_fmt)


dir = './figs/'
save_fmt = ['eps', 'pdf']

x = numpy.arange(0.1, 10, 0.1)
lamb = 1.e4 / x    # wavelength in angstroms

mixlaw = extinction.RvFbumpLaw()

# #BEAST paper extinction mixture law figure Rv_A=2.5
make_fig(lamb, mixlaw, Rv_A=2.5, figname=dir + 'RvFbumpLaw_25', save_fmt=save_fmt)

# #BEAST paper extinction mixture law figure Rv_A=5.5
make_fig(lamb, mixlaw, Rv_A=5.5, figname=dir + 'RvFbumpLaw_55', save_fmt=save_fmt)

#extinction mixture law figure 2 components
make_mixture_fig(f_A=0.5, figname=dir + 'MixtLawComp', save_fmt=save_fmt)

#BEAST paper f_bump-Rv grid figure
Rv_fbump_grid(figname=dir + 'RvFbump_grid', save_fmt=['eps', 'pdf'])

pylab.show()

# x = numpy.arange(0.1, 10, 0.1)
# lamb = 1.e4 / x    # wavelength in angstroms
#
# #BEAST paper extinction mixture law figure Rv_A=2.5
# fig = pylab.figure()
# plt = fig.add_subplot(111)
#
# ymix = mixlaw(lamb, Rv_A=2.5, f_bump=0.)
# plt.plot(x, ymix,'-k', label='$f_A=0.0$')
#
# ymix = mixlaw(lamb, Rv_A=2.5, f_bump=1. / 3.)
# plt.plot(x, ymix,'--k', label='$f_A=0.333$')
#
# ymix = mixlaw(lamb, Rv_A=2.5, f_bump=2. / 3.)
# plt.plot(x, ymix, ':k', label='$f_A=0.667$')
#
# ymix = mixlaw(lamb, Rv_A=2.5, f_bump=1.0)
# plt.plot(x, ymix,'-.k', label='$f_A=1.0$')
#
# pylab.figtext( 0.16,0.63, '$R_A(V)\  =\  2.5$', color='k',size=15)
#
# pylab.xlim(0.,8.)
# pylab.ylim(0.,7.)
# pylab.rc('font', size=13, family='serif')
# plt.set_ylabel('$A(\lambda)/A(V)$',size=15)
# plt.set_xlabel('$1/\lambda$ [$\mu m^{-1}$]',size=15)
#
# plt.legend(loc=0, frameon=False, fontsize = 15)
#
# pylab.savefig(dir+'RvFbumpLaw_25.eps')
# pylab.savefig(dir+'RvFbumpLaw_25.pdf')


# #BEAST paper extinction mixture law figure Rv_A=5.5
# fig = pylab.figure()
# plt = fig.add_subplot(111)
#
# ymix = mixlaw(lamb, Rv_A=5.5, f_bump=0.)
# plt.plot(x, ymix,'-k', label='$f_A=0.0$')
#
# ymix = mixlaw(lamb, Rv_A=5.5, f_bump=1./3.)
# plt.plot(x, ymix,'--k', label='$f_A=0.333$')
#
# ymix = mixlaw(lamb, Rv_A=5.5, f_bump=2./3.)
# plt.plot(x, ymix, ':k', label='$f_A=0.667$')
#
# ymix = mixlaw(lamb, Rv_A=5.5, f_bump=1.0)
# plt.plot(x, ymix,'-.k', label='$f_A=1.0$')
#
# #pylab.figtext( 0.16,0.63, '$R_A(V)\  =\  5.5$', color='k',size=15)
#
# pylab.xlim(0.,8.)
# pylab.ylim(0.,7.)
# pylab.rc('font', size=13, family='serif')
# plt.set_ylabel('$A(\lambda)/A(V)$',size=15)
# plt.set_xlabel('$1/\lambda$ [$\mu m^{-1}$]',size=15)
#
# plt.legend(loc=0, frameon=False, fontsize = 15)
#
# pylab.savefig(dir+'RvFbumpLaw_55.eps')
# pylab.savefig(dir+'RvFbumpLaw_55.pdf')

# #extinction mixture law figure 2 components
# fig = pylab.figure()
# plt = fig.add_subplot(111)
#
# f_A = 0.5
# ymix = mixlaw(lamb, Rv_A=3.1, f_bump=f_A)
# compA = f_A * ymix.RvLaw(lamb, Rv=3.1)
# compB = f_A * ymix.NoBumpLaw(lamb)
#
# # f99  = extinction.Fitzpatrick99()
# # gsmc = extinction.Gordon03_SMCBar()
# # compA = 0.5 * f99(lamb, Rv=3.1)
# # compB = 0.5 * gsmc(lamb)
#
# plt.plot(x,ymix,'-k',label='Total')
# plt.plot(x,compA,'--k',label='Bump Law')
# plt.plot(x,compB,':k',label='No Bump')
# #plt.plot(x,compA+compB,'-r')
# plt.set_ylabel('$A(\lambda)/A(V)$',size=15)
# plt.set_xlabel('$1/\lambda$ [$\mu m^{-1}$]',size=15)
#
# pylab.figtext( 0.16,0.695, '$f_A\  =\  0.5$', color='k',size=15, fontweight='bold')
#
# plt.legend(loc=0, frameon=False, fontsize = 15)
#
# pylab.savefig(dir+'MixtLawComp.eps')
# pylab.savefig(dir+'MixtLawComp.pdf')
#
# #BEAST paper f_bump-Rv grid figure
#
# Rv_vals=numpy.arange(2.0,6.01,0.08)
# fb_vals=numpy.arange(0.0,1.02,0.02)
# N_fb = len(fb_vals)
# N_rv = len(Rv_vals)
# Rv_res = numpy.ndarray(shape=(N_rv,N_fb), dtype=float)
# ind = numpy.ndarray(shape=(N_rv,N_fb), dtype=float)
# fig = pylab.figure()
# plt = fig.add_subplot(1, 1, 1)
# for i in range(N_fb):
#     Rv_res[0:N_rv,i]=mixlaw.get_Rv_A(Rv_vals,fbump=fb_vals[i])
#     ind = numpy.where((Rv_res[:,i] >= 2.) &  (Rv_res[:,i] <= 6.001))[0]
#     ll=len(ind)
#     plt.plot(numpy.ones(ll)*fb_vals[i],Rv_vals[ind],'ko',markersize=2)
#
# pylab.xlim(-0.05,1.05)
# pylab.ylim(1.9,6.1)
# pylab.rc('font', size=13, family='serif')
# plt.set_xlabel('$f_A$',size=15)
# plt.set_ylabel('R($V$)',size=15)
#
# pylab.savefig(dir+'RvFbump_grid.eps')
# pylab.savefig(dir+'RvFbump_grid.pdf')
#
# pylab.show()
