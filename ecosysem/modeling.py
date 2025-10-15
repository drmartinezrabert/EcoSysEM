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
from bioenergetics import CSP
 

class MSMM:
    """
    Class for Multi-State Metabolic Model
    """
# ODE function
# plot function(s)
# -> identify useful functions from other scripts and import them
# -> respect presentation syntax



    def __init__(self, typeMetabo, metabolism, eDonor, Wtype, K, mortality, envModel,
                 envArgs = None, DeltaGsynth = 9.54E-11, steepness = 0.2,
                 degradPace = 'Moderate', fluidType = 'ideal', actMethods = None):
        
        _metaboProperties = {}
        _metaboProperties['fast'] = {'protein turnover rate':1 ,'specific metabolic shift rates':1}     #resp. [h] & [1/h]
        _metaboProperties['moderate'] = {'protein turnover rate':5 ,'specific metabolic shift rates':0.2}
        _metaboProperties['slow'] = {'protein turnover rate':14 ,'specific metabolic shift rates':0.071}    
        dMtbRates = pd.DataFrame(data = _metaboProperties)
        self.envModel = envModel
        self.atmModels = ['ISA', 'MERRA2', 'CAMS', 'ISAMERRA2', 'CAMSMERRA2']
        #!!! other models?
        self.typeMtb = typeMetabo    #metabolism type (STR), e.g. 'AnMetabolisms'
        self.metabolism = metabolism    #reaction (STR), e.g. 'Mth' #??? only one community at a time?
        #self.phase = phase             #'L' as default
        self.Wtype = Wtype              # 'L_SW' or 'L_FW'
        self.K = K                      #carrying capacity (FLOAT)
        self.mortality = mortality      #(FLOAT) [1/h] #!!! and not 1/d
        self.DGsynth = DeltaGsynth      #cell synthesis required energy
        self.st = steepness
        self.mtbRates = degradPace      #'fast', 'moderate' or 'slow'
        self.eD = eDonor                #(specComp), e.g. 'CH4' for metabolism = 'Mth'
        self.fluidType = fluidType      #'ideal' or 'non-ideal'
        self.method = actMethods
        self._callEnvP(envModel = envModel, envArgs = envArgs)
        self._specMtbShiftRates = dMtbRates.loc['specific metabolic shift rates']
        #self.Bini = Bini                #initial biomass in each state (LIST) -> plot arg
        #self.time = time                #[setTime method] -> plot arg 
        #if plotMSMM == True:
        #self.plotMSMM(Bini, time) -> required args must be given when instance is created
    
    def _callEnvP(self, envModel, envArgs):
        """
        Function to import needed environment data (temperature, pH, etc.) 
        
        Parameters
        ----------
        
        envModel : STR
            Environment model from which data are extracted.
        envArgs : LIST or DICT
            Required arguments for the environment model.
            E.g. for envModel = 'ISA' :
                envArgs = {layers : 'All',      -> troposphere only ?
                           phase : 'All',       -> 'L-FW' / 'L_SW'
                           H2O : 0.0,           -> affects atmo composition
                           pH : 7.0,            -> affects bioenergetics
                           selCompounds : None, -> metabolism related only ?
                           selAlt = None,       -> subset of attributes related to portion of selected layers
                           resolution = 1000}
                or :
                envArgs = ['All', 'All', 0., 7., None, None, 1000]
            or: enArgs = None for other models
                
        Returns
        -------
        None (environment data set as MSMM attributes)
        """
        envModel = self.envModel
        atmModels = self.atmModels
        #!!! other models?
        if not isinstance(envModel, str):
            print(f'arg.error: environment model must be a string. current input: {envModel}')
            sys.exit()
        if envModel in atmModels :
            if envModel == 'ISA':
                if isinstance(envArgs, list):
                    ISAinst = ISA(*envArgs)
                elif isinstance(envArgs, dict):
                    ISAinst = ISA(**envArgs)
                else:
                    print('error in envArgs type')
                    sys.exit()
                self.pH = ISAinst.pH
                self.compositions = ISAinst.compositions
                self.compounds = ISAinst.compounds
                self.resolution = ISAinst.resolution
                self.temperature = ISAinst.temperature
                self.pressure = ISAinst.pressure
                self.altitude = ISAinst.altitude
                self.H2O = ISAinst.H2O
                self.Pi = ISAinst.Pi            #!!! could be useful for MSMMv2
                self.Ci_G = ISAinst.Ci_G        #!!! //
                if self.Wtype == 'L_FW':
                    self.Ct = ISAinst.Ci_LFW
                elif self.Wtype == 'L_SW':
                    self.Ct = ISAinst.Ci_LSW
                else: print('Liquid phase could not be recognized. Enter "L_FW" or "L_SW".')
                self.salinity = 0.0 #!!! for models other than atm, can be !=0
                self.typeKin = 'MM-Arrhenius'
                self.db = ['MM_AtmMicr', 'ArrhCor_AtmMicr']
                #!!! print('ISA attributes were set.')
            else: 
                if envArgs:
                    envArgs == None
                    print('arg.error: No arguments required for instances from other class than ISA.')
                #!!! import needed attributes from other models
        else:
            print('arg.error: environment model not found, see envModels')
            sys.exit()
        

    def _ODEsystem_MSMM(self, y, t):
        """
        Function for the differential equations system of the model.
        
        Parameters
        ----------
        
        y : LIST of INT
            Initial biomass in each metabolic state, e.g. [cell/m^3 air]
        t : LIST or np.array
            Time range over which biomass variation is computed
        
        Returns
        -------
        dB : LIST of FLOAT
            Biomass variation in each metabolic state [cell/h].
        
        """
        
        Bg = y[0]
        Bm = y[1]
        Bs = y[2]
        Blist = [Bg, Bm, Bs]
        Btot = sum(Blist)
        
        #import self.attributes
        mortality = self.mortality
        K = self.K
        # Set requested arguments for ThSA.getDeltaGr & KR.getRs
        DGr_args = {'typeRxn' : self.typeMtb,
                    'input_' : self.metabolism,
                    'phase' : 'L',
                    'specComp' : self.specComp,
                    'T' : self.temperature,
                    'pH' : self.pH,
                    'S' : self.salinity,
                    'Ct' : self.Ct,
                    'fluidType' : self.fluidType,
                    'molality' : True,
                    'methods' : None,
                    'solvent' : 'H2O',
                    'asm' : 'stoich', 
                    'warnings' : False,
                    'printDG0r' : False,
                    'printDH0r' : False
                    }       
        Rs_args = {'typeKin' : self.typeKin,
                   'paramDB' : self.db,
                   'reactions' : self.metabolism,
                   'T' : self.temperature,
                   'pH' : self.pH,
                   'Ct' : self.Ct, 
                   'sample' : 'All'
                   }
        # Compute the cell growth yield and cell-specific uptake rate    
        DGr = ThSA.getDeltaGr(**DGr_args)[0]
        Yx = DGr * 1000 * (0.5 / 1.04e-10)      # cell growth yield [cell/mol eD]
        cRate = KR.getRs(**Rs_args)   #cell-specific uptake rate [mol/cell.h]
        # Compute biomass transfer between metabolic states
        Rm_g, Rg_m, Rs_m, Rm_s, Rs_rip = MSMM._Bflux(self, Blist = Blist)
        # Compute biomass variation       
        dBg = Yx * cRate * Bg * (1 - (Bg/K)) + Rm_g - Rg_m - mortality * Bg    
        dBm =  Rg_m + Rs_m - Rm_g - Rm_s - mortality * Bm                     
        dBs =  Rm_s - Rs_m - Rs_rip - mortality * Bs                          
        dBrip = mortality * Btot + Rs_rip  
        dB = [dBg, dBm, dBs, dBrip]
        return dB

    def _Bflux(self, Blist):
        """
        Function to compute biomass transfer between metabolic states.
        
        Parameters
        ----------
        
        Blist : LIST
            List of 3 floats corresponding to biomass (e.g. [cell/m^3 air])
            in each state (growth, maintenance and survival) at time t.
            First element of the list must be for growth,
            second for maintenance and third for survival.
        
        Returns
        -------
        Rlist : LIST    #??? to be changed
             List of computed biomass transfer [cell/h] for each kind of metabolic shift:
                 - Rm_g : transfer from maintenance to growth
                 - Rg_m : transfer from growth to maintenance
                 - ...
                 - Rs_rip : transfer from survival to dead cells
        """
        #Biomass in each metabolic state
        Bg = Blist[0]
        Bm = Blist[1]
        Bs = Blist[2]
        
        eta = self._specMtbShiftRates[self.mtbRates] 
        Rm_g = Bm * eta * MSMM._stShifts(itheta = 'GxM')
        Rg_m = Bg * eta * (1 - MSMM._stShifts(itheta = 'GxM'))
        Rs_m = Bs * eta * MSMM._stShifts(itheta = 'MxS')
        Rm_s = Bm * eta * (1 - MSMM._stShifts(itheta = 'MxS'))
        Rs_rip = Bs * eta * (1 - MSMM._stShifts(itheta = 'S-RIP'))
        
        Rlist = [Rm_g, Rg_m, Rs_m, Rm_s, Rs_rip]
        
        return Rlist    #??? change return format 
        
    def _stShifts(self, itheta):
        """
        Function to compute shift control between two metabolic states.
        
        Parameters
        ----------
        
        itheta : STR
            'GxM' => shift from growth state to maintenance and conversely
            'MxS' => shift from maintenance state to survival and conversely
            'S-RIP' => shift from survival state to death
        
        Returns
        -------
        theta : TYPE
            Metabolic shift control [-]
        """
        # Set requested arguments for CSP.getAllCSP
        CSPargs = {'paramDB': self.db,
                   'typeKin': self.typeKin,
                   'typeMetabo': self.typeMtb,
                   'reaction': self.metabolism,
                   'specComp': self.eD,
                   'Ct': self.Ct,
                   'T': self.temperature,
                   'pH': self.pH,
                   'S': self.salinity,
                   'phase': 'L', 
                   'sample': 'All',
                   'fluidType': self.fluidType,
                   'molality': True,
                   'methods': None,
                   'solvent': 'H2O',
                   'asm': 'stoich',
                   'DGsynth': self.DGsynth}
        st = self.st
        Pcat = CSP.getAllCSP(**CSPargs)['Pcat']
        Pm = CSP.getAllCSP(**CSPargs)['Pm0']
        Ps = CSP.getAllCSP(**CSPargs)['Ps']
        Pcell = CSP.getAllCSP(**CSPargs)['Pcell']
        
        if itheta == 'GxM':
            theta = 1 / (np.exp((-Pcat + Pcell)/(st * Pcell)) +1)
        elif itheta == 'MxS':
            theta = 1 / (np.exp((-Pcat + Pm)/(st * Pm)) +1)
        elif itheta == 'S-RIP':
            theta = 1 / (np.exp((-Pcat + Ps)/(st * Ps)) +1)
        else: print('error in itheta value'), sys.exit() #!!!
        return theta
        
        
        
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
    
    """    
    def _runMSMM(self, envModel):
        print(f'MSMM was run for {envModel}.')
        #call plotMSMM ? useless ? non-systematical
    """



















        