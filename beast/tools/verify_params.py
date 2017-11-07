#!/usr/bin/env python
"""
	Script to verify wrong format of the BEAST input parameters 
           and terminate run_beast
	It must be called run_beast.py 
	Created by Maria Kapala on Feb 24th 2017
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from os.path import exists
from numpy import inf
from sys import exit


def verify_range(param, param_name, param_lim, noexit=False):
    # check if input param limits make sense
	
    if any(p < param_lim[0] for p in param):
        print((param_name + ' min value not physical.'))
        if not noexit: exit()

    if any(p > param_lim[1] for p in param):
        print((param_name + ' max value not physical.'))
        if not noexit: exit()

def check_grid(param, param_name, param_lim, noexit=False):
    # check if input param limits and grid initialisation make sense

    param_min, param_max, param_step = param
	
    if param_min < param_lim[0]:
        print((param_name + ' min value not physical.'))
        if not noexit: exit()

    if param_max > param_lim[1]:
        print((param_name + ' max value not physical.'))
        if not noexit: exit()

    if param_min > param_max:
        print((param_name + ' min value greater than max'))
        if not noexit: exit()

    if (param_max-param_min) < param_step:
        if param_max-param_min == 0.:
            print('Note: '+param_name+' grid is single-valued.')
        else:
            print((param_name + ' step value greater than (max-min)'))
            if not noexit: exit()

def verify_one_input_format(param, param_name, param_format, param_lim,
                            noexit=False):
    """
    Test for an input parameter correctness of format and limits.

    Parameters
    ----------
    param: str, float, or list(float)
        Input parameter for run_beast.py, defined in datamodel.py

    param_name: str
        Extact name of the param, used for printing purposes.

    param_format: str

    param_lim: list(float) or None
       pass information about any physical limitations of the param;
       [-inf, inf] when param is not constraint at all;
       None when concerns a str parameter

    noexit: boolean
       Override exit() commands if there is an error.
       Default = False (exit on error)
    """
	
    if 'list' in param_format:
        if type(param) is not list:
            if param is None:
                print('Warning: '+ param_name + ' is not defined.')
            else:
                print((param_name + ' is not in the right format - a list.'))
                if not noexit: exit()
        elif 'float' in param_format:
            is_list_of_floats = all(type(item) is float for item in param)
            if not is_list_of_floats:
                print((param_name +
                      ' is not in the right format - list of floats.'))
                if not noexit: exit()
            elif 'grid' in param_format:
                # when param is aranged from given [min, max, step],
                # instead of a specific list of values
                check_grid(param, param_name, param_lim, noexit=noexit)
            else:
                verify_range(param, param_name, param_lim, noexit=noexit)
	  			
        if 'str' in param_format:
            if type(param) is not str:
                print((param_name + ' is not in the right format - a string.')) 
                if not noexit: exit()
            elif 'file' in param_format:
                if not exists(param):
                    print((param_name +
                          ' does not exist. Please provide the file path.')) 
                    if not noexit: exit()

        if 'version' in param_format:
            if type(param) is not float:
                print((param_name + ' is not in the right format - a float'))
                if not noexit: exit()
            elif param not in param_lim:
                print((param_name +
                    ' is an invalid number, leading to version of the isochrone.'))
                if not noexit: exit()

def verify_input_format(datamodel, noexit=False):

    """
    Import BEAST input parameters from datamodel.
    Define relevant parameters, their correct names, format and limits.
    Call verify_one_input_format to test for correctness of format and 
    limits.

    Parameters
    ----------
    datamodel: module
        Input parameters are initialized in datamodel

    noexit: boolean
        Override exit() commands if there is an error.
        Default = False (exit on error)

    """

    parameters = [datamodel.z, datamodel.obsfile, \
                  datamodel.astfile, datamodel.logt, datamodel.avs,
                  datamodel.rvs, datamodel.fAs]
    parameters_names = ['z', 'obsfile', 'astfile', 'logt', 'avs', 'rvs',
                        'fAs']
    param_format = ['list_float', 'str_file', 'str_file',
                    'list_float_grid', 'list_float_grid',
                    'list_float_grid', 'list_float_grid']
    parameters_limits = [ [0., 0.1], None, None, [-inf, 10.15],
                          [0., inf], [1., 7.], [0., 1.]]

    for i, param_ in enumerate(parameters):
        verify_one_input_format(param_, parameters_names[i],
                                param_format[i], parameters_limits[i],
                                noexit=noexit)

if __name__ == '__main__':

    pass
