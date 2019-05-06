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
import numpy as np
from sys import exit
from astropy import units

from beast.physicsmodel.stars import isochrone
from beast.physicsmodel.stars import stellib
from beast.physicsmodel.dust import extinction
from beast.observationmodel.noisemodel import absflux_covmat

class datamodel():
    """
    All of the parameters from datamodel

    Attributes
    ----------
    datamodel_file : string
        name of the parameter file that was used for the input

    """


    def __init__(self, input_file):
        """
        Parameters
        ----------
        input_file : string
            Name of the input datamodel parameter file

        """

        self.datamodel_file = input_file
        #import pdb; pdb.set_trace()

        if self.datamodel_file is not None:
            # read in the file
            self.read_datamodel()
            # verify parameters
            self.verify_params()


    def read_datamodel(self):
        """
        Read in the datamodel file and set parameters
        """

        # read everything in as strings
        #import pdb; pdb.set_trace()
        #input_data = Table.read(self.datamodel_file, format='ascii.no_header',
        #                            delimiter='=', comment='#')
        #input_data = Table({'col1':['var','velocity','blah'],
        #                        'col2':['/path/to/somewhere','262.2 km/s','"some_file.txt"']})
        with open(self.datamodel_file,'r') as f:
            temp_data = f.readlines()
        input_data = [line.strip() for line in temp_data if line.strip() != '' and line.strip()[0] != '#']

        # list of parameters that have units
        # (and therefor need special treatment)
        params_with_units = ['distance_unit','velocity']

        # parse it into a dictionary
        datamodel_params = {}
        
        for i in range(len(input_data)):
            # extract parameter and value (as strings)
            param = input_data[i].split('=')[0].strip()
            val = input_data[i].split('=')[1].strip()
            # if it's a unit, put into Quantity
            if param in params_with_units:
                try:
                    datamodel_params[param] = units.Quantity(val)
                except:
                    datamodel_params[param] = eval(val)
            # if it's not a unit
            else:
                # exec the string to get parameter values accessible
                exec(input_data[i])
                datamodel_params[param] = eval(param)


        # turn dictionary into attributes
        for key in datamodel_params:
            setattr(self, key, datamodel_params[key])
                    
    

    def verify_params(self, noexit=False):

        """
        Import BEAST input parameters from datamodel.
        Define relevant parameters, their correct names, format and limits.
        Call verify_one_input_format to test for correctness of format and 
        limits.

        Parameters
        ----------
        noexit: boolean
            Override exit() commands if there is an error.
            Default = False (exit on error)

        """

        parameters_names = ['z', 'obsfile', 'astfile', 'logt', 'avs', 'rvs',
                            'fAs']
        parameters = [getattr(self, p) for p in parameters_names]
        param_format = ['list_float', 'str_file', 'str_file',
                        'list_float_grid', 'list_float_grid',
                        'list_float_grid', 'list_float_grid']
        parameters_limits = [ [0., 0.1], None, None, [-np.inf, 10.15],
                            [0., np.inf], [1., 7.], [0., 1.]]

        for i, param_ in enumerate(parameters):
           datamodel._verify_one_input_format(param_, parameters_names[i],
                                     param_format[i], parameters_limits[i],
                                     noexit=noexit)

    @staticmethod
    def _verify_one_input_format(param, param_name, param_format, param_lim,
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
           [-np.inf, np.inf] when param is not constraint at all;
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
                    datamodel._check_grid(param, param_name, param_lim, noexit=noexit)
                else:
                    datamodel._verify_range(param, param_name, param_lim, noexit=noexit)
	  			
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

    @staticmethod
    def _check_grid(param, param_name, param_lim, noexit=False):
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

    @staticmethod
    def _verify_range(param, param_name, param_lim, noexit=False):
        # check if input param limits make sense
	
        if any(p < param_lim[0] for p in param):
            print((param_name + ' min value not physical.'))
            if not noexit: exit()

        if any(p > param_lim[1] for p in param):
            print((param_name + ' max value not physical.'))
            if not noexit: exit()

