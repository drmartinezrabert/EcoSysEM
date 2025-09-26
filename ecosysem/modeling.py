# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 11:11:07 2025

@author: zLemaire
"""

# Import Python packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import importlib
import sys

# Import classes from modules (as abb.) or import with importlib
from envdef import ISA
from reactions import Reactions as Rxn
#from bioenergetics import ...
from thermodynamics import ThP, ThEq, ThSA
#TDmodule = importlib.import_module('thermodynamics')
#ThEq = TDmodule.ThEq
#ThP = TDmodule.ThP
#ThSA = TDmodule.ThSA


class MSMM:
    """
    Class for Multi-State Metabolic Model
    """
# ODE function
# plot function(s)
# -> identify useful functions from other scripts and import them
# -> respect presentation syntax
    def __init__(self, metabolism, Bini, K, mortality, DeltaGsynth):
        self.metabolism = metabolism #required reaction type
        self.Bini = Bini
        self.K = K #carrying capacity
        self.mortality = mortality
        self.DeltGsynth = DeltaGsynth
        #self.Alt = [method from ISA?]
        #self.time = [setTime method]
        
    def ODEsystem_MSMM(Bini):
        """
        Function for the differential equations system of the model.
        
        Parameters
        ----------
        
        Bini : LIST of INT
            Initial biomass in each metabolic state [cell per cubic meter]
        [...]
        
        Returns
        -------
        dB : LIST of FLOAT
            Biomass variation in each metabolic state
        """
        
        Bgrowth = Bini[0]
        Bmaintenance = Bini[1]
        Bsurvival = Bini[2]
        
        dB = [dBg, dBm, dBs, dBd]
        
        return dB
        
        
    def plotMSMM():
        """
        Function to plot solutions of the MSMM ODE system.
        
        Parameters
        ----------
        
        [...]
        
        Returns
        -------
        None
        """