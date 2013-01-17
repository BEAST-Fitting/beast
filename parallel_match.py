import grid
from multiprocessing import Pool
import numexpr
import numpy as np

global _gext
global _keys
global _r_dof
global _uniques
global _uniques_len

grid_filename = 'libs/tlusty_ext_test.fits'

_gext = grid.FileSpectralGrid(grid_filename)
_uniques = {}
_r_dof = -0.5/len(_gext.lamb) #length of flux should be length of lamb

_keys = ['Rv','f_bump','T_eff','g'] #Which parameters do you want to match?
for i in range(len(_keys)):
    _uniques[_keys[i]] = np.unique(_gext.grid[_keys[i]]) 
_uniques_len = len(_uniques.keys())
_keys = _uniques.keys() #Re-sort the keys in case of _uniques.keys() calls further down the line

#For dealing with single argument requirement
def match_wrapper(args):
    return _match(*args)

#Takes observed flux, observed error, and matching method.
#Current matching methods are expectation and max_prob.
def _match(flux, err, match_by):
        flux[np.isinf(flux)] = 1.0e-5
        err[np.isinf(err)] = 1.
        chi2 = numexpr.evaluate('sum(((where((flux==0.),1.0e-5,flux) - fluxmod)/where((err==0.),1.,err))**2,axis=1)',local_dict={'flux':flux,'err':err,'fluxmod':_gext.seds})
        chi2 = numexpr.evaluate('where((chi2>50),50,chi2)',local_dict={'chi2':chi2})
        chi2 = numexpr.evaluate('where((chi2<-50),-50,chi2)',local_dict={'chi2':chi2})
        psum = np.log(numexpr.evaluate('sum(exp(chi2))',local_dict={'chi2':chi2}))
        lnp = numexpr.evaluate('chi2-psum',local_dict={'chi2':chi2,'psum':psum})

        result = []
        if match_by == 'max_prob':
            sel = lnp.argmax()
            for i in range(_uniques_len):
                result.append(_gext.grid[_keys[i]][sel])
            return result
        elif match_by == 'expectation':
            for i in range(_uniques_len):
                prob = numexpr.evaluate('exp(lnp)',local_dict={'lnp':lnp})
                cur_vars = np.array(_uniques[_keys[i]])
                edges = np.zeros(cur_vars.shape[0]+1)
                edges[0] = cur_vars[0]
                edges[1:-1] = 0.5*(cur_vars[1:]+cur_vars[:-1])
                edges[-1] = cur_vars[-1]
                hist,edge = np.histogram(_gext.grid[_keys[i]],weights=lnp,bins=edges)
                result.append(((cur_vars*hist).sum())/(hist.sum()))
            return result
        else:
            print "Matching method not recognized."
            return 0

def multiprocess_match(fluxes, errs, match_by='expectation', nproc=4):
    pool = Pool(processes=nproc)
    jobs = [(fluxes[i],errs[i],match_by) for i in range(fluxes.shape[0])]
    result = np.array(pool.map(match_wrapper,jobs))
    pool.close()
    pool.join()
    return result
