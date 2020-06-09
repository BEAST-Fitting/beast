#!/usr/bin/env python
"""
	Script to verify wrong format of the BEAST input parameters
           and terminate run_beast
	It must be called run_beast.py
	Created by Maria Kapala on Feb 24th 2017
"""

import numpy as np
from astropy import units

from beast.physicsmodel.stars import isochrone, stellib
from beast.physicsmodel.dust import extinction
from beast.observationmodel.noisemodel import absflux_covmat
from beast.tools import verify_params


class datamodel:
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

        if self.datamodel_file is not None:
            # read in the file
            self.read_datamodel()
            # verify parameters
            verify_params.verify_input_format(self)

    def read_datamodel(self):
        """
        Read in the datamodel file and set parameters
        """

        # read everything in as strings
        # import pdb; pdb.set_trace()
        # input_data = Table.read(self.datamodel_file, format='ascii.no_header',
        #                            delimiter='=', comment='#')
        # input_data = Table({'col1':['var','velocity','blah'],
        #                        'col2':['/path/to/somewhere','262.2 km/s','"some_file.txt"']})
        with open(self.datamodel_file, "r") as f:
            temp_data = f.readlines()
        input_data = [
            line.strip()
            for line in temp_data
            if line.strip() != "" and line.strip()[0] != "#"
        ]

        # list of parameters that have units
        # (and therefore need special treatment)
        params_with_units = ["distance_unit", "velocity"]

        # parse it into a dictionary
        datamodel_params = {}

        for i in range(len(input_data)):
            # extract parameter and value (as strings)
            param = input_data[i].split("=")[0].strip()
            val = input_data[i].split("=")[1].strip()
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
