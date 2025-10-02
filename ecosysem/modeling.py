# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 11:11:07 2025

@author: zLemaire
"""

# Import Python packages
import pandas as pd
import numpy as np
from scipy.integrate import odeint
#from IPython.display import display, Math
import matplotlib.pyplot as plt
#import importlib
import sys

# Import classes from modules (as abb.) or import with importlib
from envdef import ISA
from reactions import Reactions as Rxn, KinRates as KR
#from bioenergetics import ...
from thermodynamics import ThP, ThEq, ThSA
#ThModule = importlib.import_module('thermodynamics')
#ThEq = ThModule.ThEq
#ThP = ThModule.ThP
#ThSA = ThModule.ThSA
 

class MSMM:
    """
    Class for Multi-State Metabolic Model
    """
# ODE function
# plot function(s)
# -> identify useful functions from other scripts and import them
# -> respect presentation syntax
    def __init__(self, metabolism, Bini, K, mortality, DeltaGsynth):
        self.metabolism = metabolism    #reaction (STR), e.g. 'Mth'
        self.Bini = Bini                #initial biomass in each state (LIST)
        self.K = K                      #carrying capacity (FLOAT)
        self.mortality = mortality      #(FLOAT)
        self.DeltGsynth = DeltaGsynth   #cell synthesis required energy
        #self.Alt = [method from ISA?]
        #self.time = [setTime method]
        
    def _ODEsystem_MSMM(self, y, t):
        """ tooltip
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
        
        #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX definition
        """
        
        Bg = y[0]
        Bm = y[1]
        Bs = y[2]
        Btot = Bg + Bm + Bs
        
        mortality = self.mortality
        K = self.K
        
        typeMetabo = 'metabolisms'      #allow anaerobic types?
        rMetabo = self.metabolism       #only one community at a time?
        Cdict = {}                      #import dict from ISA.method?

        
        DGr, _ = ThSA.getDeltaGr(typeMetabo, rMetabo, 'L') #Ct ?
        cellRate = KR.getRs('MM-Arrhenius', ['MM_AtmMicr','ArrhCor_AtmMicr'], rMetabo, Cdict)
        
        """ bioenergetics
        R_MG = eta * theta_Req(1, h, aerobic_oxidation_type) * Bm                  
        R_GM = eta * (1 - theta_Req(1, h, aerobic_oxidation_type)) * Bg            
        R_SM = eta * theta_Req(2, h, aerobic_oxidation_type) * Bs                  
        R_MS = eta * (1 - theta_Req(2, h, aerobic_oxidation_type)) * Bm            
        R_SRIP = eta * (1 - theta_Req(3, h, aerobic_oxidation_type)) * Bs   
        
        eta value comes from metabolic shift rates dictionary where each key
        corresponds to pace among 'slow', 'moderate' & 'fast'.
        =>  create the dictionary in bioenergetics.py ?
            or set as self attribute in MSMM.__init__?
        => same goes for the steepness parameter
        
        """
        
        
        """
        dBg = DGr * cellRate * Bg * (1 - (Bg/K)) + R_MG - R_GM - mortality * Bg    
        dBm =  R_GM + R_SM - R_MG - R_MS - mortality * Bm                     
        dBs =  R_MS - R_SM - R_SRIP - mortality * Bs                          
        dBrip = mortality * Btot + R_SRIP  
        
        dB = [dBg, dBm, dBs, dBrip]
        
        return dB
        """
        
    def plotMSMM(self, time):
        
        """
        Function to plot solutions of the MSMM ODE system.
        
        Parameters
        ----------
        time : LIST or np.array
        [...]
        
        Returns
        -------
        None
        """
        
        
        Bini = self.Bini
        
        
        if not isinstance(time, np.ndarray): time = np.array(time)
        
        plot_ = odeint(MSMM.ODEsystem_MSMM, Bini, time) 
        