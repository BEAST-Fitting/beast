#!/usr/bin/env python

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
        with open(self.datamodel_file, "r") as f:
            temp_data = f.readlines()
        # remove empty lines and comments
        input_data = [
            line.strip()
            for line in temp_data
            if line.strip() != "" and line.strip()[0] != "#"
        ]
        # if parameters are defined over multiple lines, combine lines
        for i in reversed(range(len(input_data))):
            if "=" not in input_data[i]:
                input_data[i - 1] += input_data[i]
                del input_data[i]

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
