from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from beast.core import grid
import numpy as np
from matplotlib.ticker import MaxNLocator
#from beast.tools import figure
from beast.external import eztables
from beast.config import __NTHREADS__
from beast.config import __USE_NUMEXPR__
import numexpr
__USE_NUMEXPR__ = False
#if __USE_NUMEXPR__:
#    import numexpr
#    numexpr.set_num_threads(__NTHREADS__)


class binningAxis():
    """
    Binning axis object. Stores information necessary to properly make fixed and non-uniform 1-D pdfs.
    *args:
          g             sedfitter.core.grid.ModelGrid    Model grid binningAxis will be used to project from
          axis_key      string                           key for getting axis information out of g
          bin_edges     np.array, size = N_bins+1        Edges of desired bins. Can be uneven.
    **kwargs:
          bin_centers   np.array                         Supply if bin edges meant to be off center
          var_spacing   string                           Currently 'fixed' or 'rect'(angular)
          fixed_key     string                           Variable to use when calculating rectangular edges
          axis_name     string                           Alternative xlabel for axis, defaults to axis_key
                                                         without underscores
    
    """
    def __init__(self, g, axis_key, bin_edges, bin_centers=None, var_spacing='rect', fixed_key=None, axis_name=None):
        if (axis_name == None):
            self.name = axis_key.replace('_','')
        else:
            self.name = axis_name
        self.key = axis_key
        self.edges = bin_edges
        self.ngridrows = g.grid.nrows
        if (bin_centers==None):
            self.bins = (self.edges[1:]+self.edges[:-1])/2
        else:
            self.bins = bin_centers

        self.inds = np.empty(len(self.bins), dtype=object)
        self.weights = np.empty(len(self.bins), dtype=object)
            
        if (var_spacing == 'fixed'):
            self.fixed_edges(g)
        elif (var_spacing == 'rect'):
            assert fixed_key != None, "No fixed axis passed for rectangular case, axis not properly initialized."
            self.rect_edges(g, fixed_key)
        else:
            print("Not a recognized spacing, axis not properly initialized.")

        flat_proj = self.project(np.ones(self.ngridrows))
        self.flattening_weights = np.where(flat_proj != 0., 1./flat_proj, 0.)
        try:
            self.flattening_weights /= self.flattening_weights.sum()
        except:
            self.flattening_weights = 0
            print("Warning! No overlap between parameter grid and binningAxis!")


    def rect_edges(self, g, fixed_param):
        fixed_uniques = np.unique(g.grid.data[fixed_param])
        all_edges = np.zeros([g.grid.nrows, 2], dtype=np.float_)
        for i in range(len(fixed_uniques)):
            inds = np.where(g.grid.data[fixed_param] == fixed_uniques[i])[0]
            param_vals = g.grid.data[self.key][inds]
            params_unique = np.unique(param_vals)
            sub_edges = np.zeros(params_unique.size + 1,dtype=np.float_)
            sub_edges[1:-1] = numexpr.evaluate('(pars_left+pars_right)/2.', local_dict={'pars_left':params_unique[:-1], 'pars_right':params_unique[1:]})
            sub_edges[0] = params_unique[0] - (sub_edges[1]-params_unique[0])
            sub_edges[-1] = params_unique[-1] + (params_unique[-1] - sub_edges[-2])
            
            for e, val in enumerate(params_unique):
                all_edges[inds[np.where(param_vals==val)]] = sub_edges[e:e+2]

        
        for i in range(len(self.bins)):
            if __USE_NUMEXPR__:
                inds = np.where(numexpr.evaluate('where(ledges >= pledges, 1, 0)*where(ledges <= predges, 1, 0) + where(redges >= pledges, 1, 0)*where(redges <= predges, 1, 0) + where(ledges <= pledges, 1, 0)*where(redges >= predges, 1, 0) + where(ledges >= pledges, 1, 0)*where(redges <= predges, 1, 0)', local_dict={'ledges':self.edges[i],'redges':self.edges[i+1],'pledges':all_edges[:, 0],'predges':all_edges[:,1]}))[0]
            else:
                inds = np.where(((self.edges[i] >= all_edges[:,0])*(self.edges[i] <= all_edges[:,1])) + \
            ((self.edges[i+1] >= all_edges[:,0])*(self.edges[i+1] <= all_edges[:,1])) +\
            ((self.edges[i] <= all_edges[:,0])*(self.edges[i+1] >= all_edges[:,1])) +\
             ((self.edges[i] >= all_edges[:,0])*(self.edges[i+1] <= all_edges[:,1])))[0]
            if (len(inds)):
                self.weights[i] = (np.where(self.edges[i+1] <= all_edges[inds,1], self.edges[i+1], all_edges[inds,1]) - np.where(self.edges[i] >= all_edges[inds,0], self.edges[i], all_edges[inds,0]))/(all_edges[inds,1]-all_edges[inds,0])
                self.inds[i] = inds
            else:
                pass

    def fixed_edges(self, g):
        params = np.unique(g.grid.data[self.key])
        
        try:
            interval = params[1]-params[0]
        except:
            print("Parameter single-valued, binningAxis not properly initialized.")
            return False
        
        i_params = np.zeros(len(params),dtype=object)
        for i in range(len(params)):
            i_params[i] = np.where(g.grid.data[self.key] == params[i])[0]
        for i in range(len(self.bins)):
            if __USE_NUMEXPR__:
                param_inds = np.where(numexpr.evaluate('where(ledges >= params-interval, 1, 0)*where(ledges <= params+interval, 1, 0) + where(redges >= params-interval, 1, 0)*where(redges <= params+interval, 1, 0) + where(ledges <= params-interval, 1, 0)*where(redges >= params+interval, 1, 0) + where(ledges >= params-interval, 1, 0)*where(redges <= params+interval, 1, 0)', local_dict={'ledges':self.edges[i],'redges':self.edges[i+1],'params':params,'interval':interval/2}))[0]
            else:
                param_inds = np.where(((self.edges[i] >= (params-interval/2)) * (self.edges[i] <= (params+interval/2))) + ((self.edges[i+1] >= (params-interval/2)) * (self.edges[i+1] <= (params+interval/2))) + ((self.edges[i] <= (params-interval/2)) * (self.edges[i+1] >= (params+interval/2))) + ((self.edges[i] >= (params-interval/2)) * (self.edges[i+1] <= (params+interval/2))))[0]
                
            if (len(param_inds)):
                if __USE_NUMEXPR__:
                    weights = numexpr.evaluate('where(redges <= params+interval, redges, params+interval) - where(ledges >= params-interval, ledges, params-interval)', local_dict={'ledges':self.edges[i], 'redges':self.edges[i+1], 'params':params[param_inds], 'interval':interval/2})
                else:
                    weights = (np.where(self.edges[i+1] <= params[param_inds] + interval/2, self.edges[i+1], params[param_inds] + interval/2) - np.where(self.edges[i] >= params[param_inds] - interval/2, self.edges[i], params[param_inds] - interval/2))
                self.inds[i] = np.hstack(i_params[param_inds])
                self.weights[i] = np.hstack([(np.ones(len(i_params[param_inds[j]]))*weights[j]) for j in range(len(param_inds))])/interval
            else:
                pass

    def project(self, likelihood, inds=None, weights=1., to_return='pdf', flatten=False):
        pdf = np.zeros(len(self.bins))
        if (inds is not None):
            full_likelihood = np.zeros(self.ngridrows)
            full_likelihood[inds] = likelihood*weights
        else:
            full_likelihood = likelihood*weights
            
        for i in range(len(self.bins)):
            if (self.inds[i] == None):
                pdf[i] = 0
            else:
                if __USE_NUMEXPR__:
                    pdf[i] = numexpr.evaluate('sum(weights*full_likelihood)', local_dict={'weights':self.weights[i], 'full_likelihood':full_likelihood[self.inds[i]]})
                else:
                    pdf[i] = (self.weights[i]*full_likelihood[self.inds[i]]).sum()

        if (flatten):
            pdf *= self.flattening_weights
        
        if (to_return == 'pdf'):
            return pdf
        elif (to_return == 'max'):
            return self.bins[np.argmax(pdf)]
        elif (to_return == 'expectation'):
            pdf *= (self.edges[1:]-self.edges[:-1])
            return ((self.bins*pdf).sum())/(pdf.sum())
        elif (to_return == 'total'):
            return pdf.sum()
        else:
            print("to_return not recognized, returning nothing")
            return None


    def plot_projection(self, likelihood, inds=None, weights=1., flatten=False, align='b'):
        pdf = self.project(likelihood, inds=inds, weights=weights, to_return='pdf', flatten=flatten)
        if align=='b':
            figure.bar(self.bins, pdf, width=self.edges[1:]-self.edges[:-1], align='center', color='grey', bottom=0)           
            width = self.edges[-1]-self.edges[0]
            height = figure.ylim()[1]
            figure.xlim(self.edges[0]-0.03*width, self.edges[-1]+0.03*width)
            figure.ylim(-height*0.03, height)
            figure.xlabel(self.name)
        elif align=='r':
            figure.barh(self.bins, pdf, height=self.edges[1:]-self.edges[:-1], align='center', color='grey', left=0)
            height = self.edges[-1]-self.edges[0]
            width = figure.xlim()[1]
            figure.ylim(self.edges[0]-0.03*width, self.edges[-1]+0.03*width)
            figure.xlim(width, -width*0.03)
            figure.ylabel(self.name)

    def get_flattening_weights(self):
        return self.flattening_weights

def rectangular_param_weights(g, param, fixed_param):
    """
    For the log(M), log(A) case, param='logM'
    fixed_param='logA'. Calculates relative span of rectangular
    blocks in parameter space.

    *args:
      g             MemoryGrid      Grid to calculate weights for
      param         string          name of parameter to calculate weights for
      fixed_param   string          name of parameter to along which spacing of param varies 

    returns:
      param_weights numpy.ndarray   One dimensional array containing normalized widths along the
                                    param axis of rectangular grid cells
    """
    fixed_uniques = np.unique(g.grid.data[fixed_param])
    param_weights = np.zeros(g.grid.nrows, dtype=np.float_)
    fixed_interval = fixed_uniques[1]-fixed_uniques[0]
    for i in range(len(fixed_uniques)):
        inds = np.where(g.grid.data[fixed_param] == fixed_uniques[i])[0]
        param_vals = g.grid.data[param][inds]
        params_unique = np.unique(param_vals)
        edges = np.zeros(params_unique.size + 1,dtype=np.float_)
        edges[1:-1] = (params_unique[1:]+params_unique[:-1])/2.
        edges[0] = params_unique[0] - (edges[1]-params_unique[0])
        edges[-1] = params_unique[-1] + (params_unique[-1] - edges[-2])
        for val in params_unique:
            param_weights[inds[np.where(param_vals==val)]] = edges[1:] - edges[:-1]
    return param_weights/(param_weights.sum())

def fixed_param_weights(g, param, interval=0):
    """
    Assuming fixed intervals. Fixed interval can be given using
    the interval= keyword or calculated (crudely) from the parameter values.

    *args:
      g             MemoryGrid      Grid to calculate weights for
      param         string          name of parameter to calculate weights for
      
    **kwargs:
      interval      float           Spacing between grid points

    returns:
      param_weights numpy.ndarray   One dimensional array containing normalized widths of grid cells
                                    along the param axis
    """
    if (interval == 0):
        uniques = np.unique(g.grid.data[param])
        interval = uniques[1]-uniques[0]
    param_weights = interval * np.ones(g.grid.nrows, dtype=np.float_)
    return param_weights/(param_weights.sum())


