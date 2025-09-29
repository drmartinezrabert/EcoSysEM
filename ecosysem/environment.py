# -*- coding: utf-8 -*-
"""
Created on Fri Aug 29 12:00:06 2025

@author: emartinez

envdef - define and create environment(s)
=========================================

Summary
-------

Lorem ipsum...

Classes
-------
    
    :py:class:`ISA` - International Standard Atmosphere (ISA).

Functions
--------

    :py:fun:`lorem_ipsum` - Lorem ipsum.
"""

from thermodynamics import ThEq as eQ
from scipy.interpolate import RegularGridInterpolator
from molmass import Formula

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import earthaccess
import itertools
import xarray as xr
import datetime
import calendar
import os
import time
import cdsapi
import zipfile

import warnings
warnings.simplefilter(action = 'ignore')

class Environment:
    def _openNPZ():
        pass
    
    def _saveNPZ():
        pass
    
    def keys():
        pass
    
    def deleteKey():
        pass
    
    def _getDictionary():
        pass

# Atmosphere ------------------------------------------------------------------

class Atmosphere(Environment):
    def callAtmosphere():
        pass
    
    def plotAtmosphere():
        pass
    
    def _HfromP():
        pass
    
    def selectRegion():
        pass
    
    def _closestCoord():
        pass

class ISA(Atmosphere):
    pass

class MERRA2(Atmosphere):
    pass

class CAMS(Atmosphere):
    pass

class ISAMERRA2(ISA, MERRA2):
    pass

class CAMSMERRA2(CAMS, MERRA2):
    pass

# Hydrosphere -----------------------------------------------------------------

class Hydrosphere(Environment):
    pass

class Ocean(Hydrosphere):
    pass

class Sediments(Hydrosphere):
    pass

class River(Hydrosphere):
    pass

# Cryosphere ------------------------------------------------------------------

class Cryosphere(Environment):
    pass

class Glacier(Cryosphere):
    pass

# Lithosphere -----------------------------------------------------------------

class Lithosphere(Environment):
    pass

class Crust(Lithosphere):
    pass
