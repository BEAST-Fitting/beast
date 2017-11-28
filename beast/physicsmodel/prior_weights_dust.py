"""
Dust Prior Weights
==================
The priors on A(V), R(V), and f_A computed as weights to use
in the posterior calculations.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
#from scipy.integrate import quad

__all__ = ['PriorWeightsDust']

def _lognorm(x, max_pos, sigma=0.5, N=1.):
    """
    Lognormal distribution

    Parameters
    ----------
    xs: vector
       x values

    max_pos: float
       Position of the lognormal function's maximum

    sigma: float
       Sigma of the lognormal function

    N: floats
       Multiplicative factor

    Returns
    -------
    lognormal computed on the x grid
    """
    sqrt_2pi = 1. / np.sqrt(2 * np.pi)
    mu = np.log(max_pos) + sigma**2

    # avoid zero or negative due to log
    indxs, = np.where(x > 0)

    lnorm = np.zeros(len(x))

    log_x = np.log(x[indxs])
    normalization = sqrt_2pi / (x[indxs] * sigma)
    
    lnorm[indxs] = (N * normalization 
                    * np.exp(-0.5 * ((log_x - mu) / sigma)**2))
    return lnorm

def _two_lognorm(xs, 
                 max_pos1, max_pos2, 
                 sigma1=0.5, sigma2=0.5, 
                 N1=1., N2=1.):
    """
    Mixture of 2 lognormal functions

    Parameters
    ----------
    xs: vector
       x values

    max_pos1: float
       Position of the lognormal function's maximum for component 1
    max_pos2: float
       Position of the lognormal function's maximum for component 2

    sigma1: float
       Sigma of the lognormal function for component 1
    sigma2: float
       Sigma of the lognormal function for component 2

    N1: floats
       Multiplicative factor for component 1
    N2: floats
       Multiplicative factor for component 2

    Returns
    -------
    Mixture model: (LOGNORM1 + LOGNORM2) / INTEGRAL(LOGNORM1 + LOGNORM2)
    """
    pointwise = (_lognorm(xs, max_pos1, sigma=sigma1, N=N1)
                 + _lognorm(xs, max_pos2, sigma=sigma2, N=N2))
    normalization = np.trapz(pointwise, x=xs)
    return pointwise / normalization

def _exponential(x, a=2.0, N=1.):
    """
    Exponential distribution
    Parameters
    ----------
    x: vector
       x values
    a: float
       Decay Rate parameter in exp: N*e^-ax
       Distribution Mean = 1/a
    N: float
       Multiplicative factor
    Returns
    -------
    exponential computed on the x grid
    """
    return N * np.exp(-1. * a * x)

class PriorWeightsDust():
    """
    Compute the priors as weights given the input grid
    """

    def __init__(self, av_vals, av_model,
                 rv_vals, rv_model,
                 fA_vals, fA_model):
        """
        Initialize with basic information
        """

        # save the parameter grid values for later use
        self.av_vals = np.copy(av_vals)
        self.rv_vals = np.copy(rv_vals)
        self.fA_vals = np.copy(fA_vals)

        # initialize the prior_weights
        #   will use these variables to save the prior weights
        #   for use in adjusting priors in the megabeast
        self.set_av_weights(av_model)
        self.set_rv_weights(rv_model)
        self.set_fA_weights(fA_model)

    def get_av_weight(self, av):
        """
        Get the weight for one A(V)

        Parameters
        ----------
        av: float
           A(V) of point

        Returns
        -------
        weight: float
           weight fo the point
        """
        indx, = np.where(av == self.av_vals)
        
        return np.asscalar(self.av_priors[indx])

    def get_rv_weight(self, rv):
        """
        Get the weight for one R(V)

        Parameters
        ----------
        rv: float 
           R(V) of point

        Returns
        -------
        weight: float
           weight fo the point
        """
        indx, = np.where(rv == self.rv_vals)
        
        return np.asscalar(self.rv_priors[indx])

    def get_fA_weight(self, fA):
        """
        Get the weight for one f_A

        Parameters
        ----------
        fA: float
           f_A of point

        Returns
        -------
        weight: float
           weight fo the point
        """
        indx, = np.where(fA == self.fA_vals)
        
        return np.asscalar(self.fA_priors[indx])

    def get_weight(self, av, rv, fA):
        """
        Get the weight for one point in A(V), R(V), f_A space

        Parameters
        ----------
        av: float
           A(V) of point
        rv: float 
           R(V) of point
        fA: float
           f_A of point

        Returns
        -------
        weight: float
           weight fo the point
        """
        return(self.get_av_weight(av)*self.get_rv_weight(rv)
               *self.get_fA_weight(fA))

    def set_av_weights(self, model={'name': 'flat'}):
        """
        Weights on A(V) based on input model choice

        Parameters
        ----------
        model: string
          Choice of model type [default=flat]
          flat = flat prior on linear A(V)
          lognormal = lognormal prior on linear A(V)
          two_lognormal = two lognormal prior on linear A(V)
          exponential = exponential prior on linear A(V)
        """
        if model['name'] == 'flat':
            self.av_priors = np.full(self.av_vals.shape,1.0)
        elif model['name'] == 'lognormal':
            self.av_priors = _lognorm(self.av_vals,
                                      model['max_pos'],
                                      sigma=model['sigma'],
                                      N=model['N'])
        elif model['name'] == 'two_lognormal':
            self.av_priors = _two_lognorm(self.av_vals,
                                          model['max_pos1'],
                                          model['max_pos2'],
                                          sigma1=model['sigma1'],
                                          sigma2=model['sigma2'],
                                          N1=model['N1'],
                                          N2=model['N2'])
        elif model['name'] == 'exponential':
            self.av_priors = _exponential(self.av_vals,
                                          a=model['a'],
                                          N=model['N'])
        else:
            print('**error in setting the A(V) dust prior weights!**')
            print('**model ' + model['name'] + ' not supported**')
            exit()
    
    def set_rv_weights(self, model={'name': 'flat'}):
        """
        Weights on R(V) based on input model choice

        Parameters
        ----------
        model: string
          Choice of model type [default=flat]
          flat = flat prior on linear R(V)
          lognormal = lognormal prior on linear R(V)
          two_lognormal = two lognormal prior on linear R(V)
        """
        if model['name'] == 'flat':
            self.rv_priors = np.full(self.rv_vals.shape,1.0)
        elif model['name'] == 'lognormal':
            self.rv_priors = _lognorm(self.rv_vals,
                                      model['max_pos'],
                                      sigma=model['sigma'],
                                      N=model['N'])
        elif model['name'] == 'two_lognormal':
            self.rv_priors = _two_lognorm(self.rv_vals,
                                          model['max_pos1'],
                                          model['max_pos2'],
                                          sigma1=model['sigma1'],
                                          sigma2=model['sigma2'],
                                          N1=model['N1'],
                                          N2=model['N2'])
        else:
            print('**error in setting the R(V) dust prior weights!**')
            print('**model ' + model['name'] + ' not supported**')
            exit()
    
    def set_fA_weights(self, model={'name': 'flat'}):
        """
        Weights on f_A based on input model choice

        Parameters
        ----------
        model: string
          Choice of model type [default=flat]
          flat = flat prior on linear f_A
          lognormal = lognormal prior on linear f_A
          two_lognormal = two lognormal prior on linear f_A
        """
        if model['name'] == 'flat':
            self.fA_priors = np.full(self.fA_vals.shape,1.0)
        elif model['name'] == 'lognormal':
            self.fA_priors = _lognorm(self.fA_vals,
                                      model['max_pos'],
                                      sigma=model['sigma'],
                                      N=model['N'])
        elif model['name'] == 'two_lognormal':
            self.fA_priors = _two_lognorm(self.fA_vals,
                                          model['max_pos1'],
                                          model['max_pos2'],
                                          sigma1=model['sigma1'],
                                          sigma2=model['sigma2'],
                                          N1=model['N1'],
                                          N2=model['N2'])
        else:
            print('**error in setting the f_A dust prior weights!**')
            print('**model ' + model['name'] + ' not supported**')
            exit()
    
