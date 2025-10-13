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


    def __init__(self, typeMetabo, metabolism, Bini, K, mortality, envModel,
                 envArgs = None, sample = None, DeltaGsynth = 9.54E-11, steepness = 0.2,
                 degradation_pace = 'Moderate'):
        
              
        
        self.metaboType = typeMetabo    #metabolism type (STR), e.g. 'AnMetabolisms'
        self.metabolism = metabolism    #reaction (STR), e.g. 'Mth'
        self.Bini = Bini                #initial biomass in each state (LIST)
        self.K = K                      #carrying capacity (FLOAT)
        self.mortality = mortality      #(FLOAT) [1/h] #!!! and not 1/d
        self.DeltGsynth = DeltaGsynth   #cell synthesis required energy
        self.st = steepness
        self.metaboRates = degradation_pace
        #self.Alt = [method from ISA?]
        #self.time = [setTime method] -> plot arg 
        #other needed args : pH = 7., S = None, Ct = 1., T = 298.15, phase = 'L'
        
        
        self._callEnvP(envModel = envModel, envArgs = envArgs)
        #self._runMSMM -> useless? plot as callable method
        
          
    
    def _callEnvP(self, envModel, envArgs):
        """
        Function to import needed environment data (temperature, pH, etc.) 
        
        Parameters
        ----------
        
        envModel : STR
            {short description ; [unit] ; default or expected values ; error raise ; examples}
        envArgs : LIST or DICT
            E.g. for envModel = 'ISA' :
                envArgs = {layers : 'All',      -> troposphere only ?
                           phase : 'All',       -> 'L-FW' ?
                           H2O : 0.0,           -> affects atmo composition
                           pH : 7.0,            -> 
                           selCompounds : None, -> metabolism related only ?
                           selAlt = None,       -> subset of attributes related to portion of selected layers
                           resolution = 1000}
                or :
                envArgs = ['All', 'All', 0., 7., None, None, 1000]
        [...]
        
        Returns
        -------
        None (environment data set as MSMM attributes)
        """
        
        envModels = ['ISA', 'MERRA2', 'CAMS', 'ISAMERRA2', 'CAMSMERRA2'] #!!! other models?
        
        if not isinstance(envModel, str):
            print(f'arg.error: environment model must be a string. current input: {envModel}')
            sys.exit()
        if envModel in envModels :
            if envModel == 'ISA':
                if isinstance(envArgs, list):
                    ISAinst = ISA(*envArgs)
                elif isinstance(envArgs, dict):
                    ISAinst = ISA(**envArgs)
                else:
                    print('error in envArgs type')
                    sys.exit()
                self.temperature = ISAinst.temperature
                self.pressure = ISAinst.pressure
                self.altitude = ISAinst.altitude
                self.compounds = ISAinst.compounds
                [...] #!!! other requested attributes
                """
                #ISA args : 
                    self,
                    layers, 
                    phase = 'All', 
                    H2O = 0.0, 
                    pH = 7.0, 
                    selCompounds = None, -> all
                    selAlt = None, -> all
                    resolution = 1000
                    
                #ISA attributes : 
                    self.layers = layers
                    !self.pH = pH
                    self.resolution = resolution
                    self._computeTandP_ISA(layers, dISA)
                        !self.altitude = alt # [m]
                        !self.temperature = t + 273.15 # [K]
                        !self.pressure = p # [Pa]
                    self.compounds = dDC['Compounds']
                    self.compositions = pd.Series(dDC.Compositions.values, index = dDC.Compounds).to_dict()
                    self._computeWaterContent(H2O, dDC)
                        if H2O != 0 :
                            self.compositions = pd.Series(newComp, index = dDC.Compounds).to_dict()
                            self.H2O = H2O
                    if selAlt:
                        self._selectAltitude(selAlt)
                        self.altitude = prevAlt[imAlt:iMAlt]
                        self.temperature = prevT[imAlt:iMAlt]
                        self.pressure = prevP[imAlt:iMAlt]
                    !self._getConcISA(phase, selCompounds)
                        self.Pi = dict_Pi
                        self.Ci_G = dict_Ci_G
                        self.Ci_LFW = dict_Ci_LFW
                        self.Ci_LSW = dict_Ci_LSW
                """
                
                print('ISA attributes were set.')
            else: 
                if envArgs:
                    envArgs == None
                    print('arg.error: No arguments required for instances from other class than ISA.')
                #import needed attributes from other models through existing methods
        else:
            print('arg.error: environment model not found, see envModels')
            sys.exit()
        

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
        
        XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX definition
        """
        
        Bg = y[0]
        Bm = y[1]
        Bs = y[2]
        Blist = [Bg, Bm, Bs]
        Btot = sum(Blist)
        
        mortality = self.mortality
        K = self.K
        
        """
        typeMetabo = self.metaboType      #allow anaerobic types?
        metabo = self.metabolism       #only one community at a time?
        metaboRates = self.metaboRates
        st = self.st
        Cdict = {}                      #import dict from ISA.method?
        """
            
        DGr, _ = ThSA.getDeltaGr(typeMetabo, rMetabo, 'L') #Ct ?
        Yx = DGr * 1000 * (0.5 / 1.04e-10)      # cell growth yield [cell/mol eD]
        cRate = KR.getRs('MM-Arrhenius', ['MM_AtmMicr','ArrhCor_AtmMicr'], rMetabo, Cdict)   #cell-specific uptake rate [mol/cell.h]
        
        Rm_g, Rg_m, Rs_m, Rm_s, Rs_rip = MSMM._Bflux(self, Blist = Blist, metaboRates = self.metaboRates, st = st)
        
        """
        eta value comes from metabolic shift rates dictionary where each key
        corresponds to pace among 'slow', 'moderate' & 'fast'.
        =>  create the dictionary in bioenergetics.py ?
            or set as self attribute in MSMM.__init__?
        => same goes for the steepness parameter
        """

        dBg = Yx * cRate * Bg * (1 - (Bg/K)) + Rm_g - Rg_m - mortality * Bg    
        dBm =  Rg_m + Rs_m - Rm_g - Rm_s - mortality * Bm                     
        dBs =  Rm_s - Rs_m - Rs_rip - mortality * Bs                          
        dBrip = mortality * Btot + Rs_rip  
        
        dB = [dBg, dBm, dBs, dBrip]     # should be in [cell/m^3 air.h]
        
        return dB



    def _Bflux(self, Blist, st = 0.2):
        """
        Function to do blablabla [...].
        
        Parameters
        ----------
        
        input1 : TYPE(s)
            {short description ; [unit] ; default or expected values ; error raise ; examples}
        input2 : //
        [...]
        
        Returns
        -------
        output1 : TYPE
            {short description ; [unit] ; error raise ; interpretation ; examples}
        output2 : //
        [...]
        """
        
        metabo = self.metaboType
        rxn = self.metabolism
        metaboRates = self.metaboRates
        
        Bg = Blist[0]
        Bm = Blist[1]
        Bs = Blist[2]
        
        
        _metabo_properties = {}
        _metabo_properties['fast'] = {'protein turnover rate':1 ,'specific metabolic shift rates':1}     #!!! resp. [h] & [1/h]
        _metabo_properties['moderate'] = {'protein turnover rate':5 ,'specific metabolic shift rates':0.2}
        _metabo_properties['slow'] = {'protein turnover rate':14 ,'specific metabolic shift rates':0.071}
        #dMR = pd.DataFrame(data = _metabo_properties)
        
        eta = _metabo_properties[metaboRates]['specific metabolic shift rates']
                
        
        Rm_g = Bm * eta * MSMM._stShifts(itheta = 1, typeMetabo = metabo, reaction = rxn)
        Rg_m = Bg * eta * (1 - MSMM._stShifts(itheta = 1, typeMetabo = metabo, reaction = rxn))
        Rs_m = Bs * eta * MSMM._stShifts(itheta = 2, typeMetabo = metabo, reaction = rxn)
        Rm_s = Bm * eta * (1 - MSMM._stShifts(itheta = 2, typeMetabo = metabo, reaction = rxn))
        Rs_rip = Bs * eta * (1 - MSMM._stShifts(itheta = 3, typeMetabo = metabo, reaction = rxn))
        
        Rlist = [Rm_g, Rg_m, Rs_m, Rm_s, Rs_rip]    # should be in [cell/m^3 air.h]
        
        return Rlist    #!!! change return format 
        
        
        
    def _stShifts(self, itheta, typeMetabo, reaction,
                  pH = 7., S = None, Ct = 1., T = 298.15, #!!! Ct, pH, S, T, phase
                  phase = 'L', DGsynth = 9.54E-11, st = 0.2):
        """
        Function to do blablabla [...].
        
        Parameters
        ----------
        
        input1 : TYPE(s)
            {short description ; [unit] ; default or expected values ; error raise ; examples}
        input2 : //
        [...]
        
        Returns
        -------
        output1 : TYPE
            {short description ; [unit] ; error raise ; interpretation ; examples}
        output2 : //
        [...]
        """
        
        # !!! errors for args type
        
        Pcat = CSP.getPcat(typeMetabo, reaction, pH, S, Ct, T, phase)
        Pm = CSP.getPm0(T)
        Ps = CSP.getPs(T)
        Pcell = CSP.getAll(typeMetabo, reaction, Ct, T, phase, DGsynth)[-1]
        
        
        if itheta == 1:
            theta = 1 / (np.exp((-Pcat + Pcell)/(st * Pcell)) +1)
        elif itheta == 2:
            theta = 1 / (np.exp((-Pcat + Pm)/(st * Pm)) +1)
        elif itheta == 3:
            theta = 1 / (np.exp((-Pcat + Ps)/(st * Ps)) +1)
        else: print('error in itheta value'), sys.exit()
            
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



















        