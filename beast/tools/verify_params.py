#!/usr/bin/env python
"""
Script to verify and warn about wrong format of the BEAST input parameters.

Assumes that the datamodel.py file exists in the same directory as this script.
  And it must be called datamodel.py 
     used in make_model.py code in physicsmodel, more recoding needed to remove
     this dependency
"""

def verify_range(param, param_min, param_max, param_step):
  ''' Test for parameter 
  '''
  if param_min > param_max:
    print("WARNING: input min greater than max")
  if 