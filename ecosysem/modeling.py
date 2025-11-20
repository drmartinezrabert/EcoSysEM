# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 11:11:07 2025

@author: zLemaire
"""

# Import Python packages
import pandas as pd
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import sys
from time import process_time

# Import classes from modules (as abb.) or import with importlib
from envdef import ISA
from reactions import KinRates as KR
from thermodynamics import ThSA
from bioenergetics import CSP

import copy


class MSMM:
    """
    Class for Multi-State Metabolic Model
    
    
    reminder for ISA args : #!!!
        self,
        'layers',   => troposphere index = 0
        'phase' : 'All',
        'H2O' : 0.0,
        'pH' : 7.0, 
        'selCompounds' : None, -> means all
        'selAlt' : None, -> means all
        'resolution' : 1000
    """

    def __init__(self, envModel, coord, typeMetabo, metabolism, eDonor, K, mortality, 
                 Wtype = 'L-FW', pH = 7.0, dataType = 'cyly', dataRange = [2020, 2024],
                 DeltaGsynth = 9.54E-11, steepness = 0.2, degradPace = 'moderate',
                 salinity = None, fluidType = 'ideal', actMethods = None):
        
        _metaboProperties = {}
        _metaboProperties['fast'] = {'protein turnover rate':1 ,'specific metabolic shift rates':1}     #resp. [h] & [1/h]
        _metaboProperties['moderate'] = {'protein turnover rate':5 ,'specific metabolic shift rates':0.2}
        _metaboProperties['slow'] = {'protein turnover rate':14 ,'specific metabolic shift rates':0.071}    
        dMtbRates = pd.DataFrame(data = _metaboProperties)
        
        self.envModel = envModel
        self.atmModels = ['ISA', 'ISAMERRA2', 'CAMSMERRA2']
        #!!! other models?
        if not isinstance(coord, (list, np.ndarray)): coord = [coord]
        self.coord = coord
        self.typeMtb = typeMetabo    #metabolism type (STR), e.g. 'AnMetabolisms'
        if not isinstance(metabolism, list): metabolism = [metabolism]
        if not len(metabolism) == 1:
            raise AttributeError(f'A single metabolism name must be given, current input: {metabolism}.')
        self.metabolism = metabolism    #reaction (STR), e.g. 'Mth' #??? only one community at a time?
        if not isinstance(Wtype, str):
            raise TypeError(f'Wtype must be a string, current type:{type(Wtype)}.')
        if not Wtype == 'L-FW' and not Wtype == 'L-SW':
            raise NameError(f'Given Wtype invalid ({Wtype}). Did you mean "L-FW" or "L-SW"?')
        self.Wtype = Wtype              # 'L_SW' or 'L_FW'
        if Wtype == 'L-SW':
            if not isinstance(salinity, list): salinity = [salinity]
            if not len(salinity) == 1:
                raise AttributeError(f'A single salinity value must be given, current salinity: {salinity}.')
            if not isinstance(salinity[0], (float, int)):
                raise TypeError(f'Given salinity value must be float or int, current type: {salinity[0]}.')
            self.salinity = salinity
        else: self.salinity = [0.0]
        if isinstance(pH, list):
            if not len(pH) == 1:
                raise AttributeError(f'A single pH value must be given, current pH: {pH}.')
            else: pH = pH[0]
        if not isinstance(pH, (float, int)):
            raise TypeError(f'Given pH value must be float or int, current type: {pH}.')
        self.pH = pH
        self.dataType = dataType
        self.dataRange = dataRange
        self.K = K                      #carrying capacity (FLOAT)
        self.mortality = mortality      #(FLOAT) [1/h] #!!! and not 1/d
        self.DGsynth = DeltaGsynth      #cell synthesis required energy
        self.st = steepness
        self.mtbRates = degradPace      #'fast', 'moderate' or 'slow'
        self.eD = eDonor                #(specComp), e.g. 'CH4' for metabolism = 'Mth'
        self.fluidType = fluidType      #'ideal' or 'non-ideal'
        self.method = actMethods
        self.typeKin = 'MM-Arrhenius'
        self.db = ['MM_AtmMicr', 'ArrhCor_AtmMicr']
        self._callEnvP(envModel)
        self._specMtbShiftRates = dMtbRates.loc['specific metabolic shift rates']
        self.Bsol = {}
    
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
                self.Pi = ISAinst.Pi        #!!! could be useful for MSMMv2
                self.Ci_G = ISAinst.Ci_G       #!!! //
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

    def _ODEsystem_MSMM(self, t, y, idP):
        """
        Function for the differential equations system of the model.
        
        Parameters
        ----------
        
        y : LIST of INT
            Initial biomass in each metabolic state, e.g. [cell/m^3 air]
        t : LIST or np.array
            Time range over which biomass variation is computed
        idP : INT
            Parameter index from the plot function (involved parameters depend on envModel)
            
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
        mortality = self.mortality  #??? no copy required for immutable data types like floats or strings
        K = self.K
        typeKin = self.typeKin
        paramDB = self.db.copy()
        typeRxn = self.typeMtb
        rxn = self.metabolism
        specComp = self.eD
        T = self.temperature.copy()
        pH = self.pH
        S = self.salinity
        C = self.Ct.copy()
        fluidType = self.fluidType
        # Set requested arguments for ThSA.getDeltaGr & KR.getRs
        if self.envModel in self.atmModels:
            DGr_args = {'typeRxn' : typeRxn,
                        'input_' : rxn,
                        'phase' : 'L',
                        'specComp' : specComp,
                        'T' : T,
                        'pH' : pH,
                        'S' : S,
                        'Ct' : C,
                        'fluidType' : fluidType,
                        'molality' : True,
                        'methods' : None,
                        'solvent' : 'H2O',
                        'asm' : 'stoich', 
                        'warnings' : False,
                        'printDG0r' : False,
                        'printDH0r' : False
                        }       
            Rs_args = {'typeKin' : typeKin,
                       'paramDB' : paramDB,
                       'reactions' : rxn,
                       'T' : T,
                       'pH' : pH,
                       'Ct' : C, 
                       'sample' : 'All'
                       }
            #!!! set getDeltaGr and getRs args for other models
        # Compute the cell growth yield and cell-specific uptake rate    
        DGr = (ThSA.getDeltaGr(**DGr_args)[0])[idP] * 1000  #[J/moleD]
        Yx = -(DGr * (0.5 / 1.04e-10))      # cell growth yield [cell/mol eD]
        cRate = pd.DataFrame((KR.getRs(**Rs_args)[0])[self.metabolism]).mean(axis=1)[idP]   #cell-specific uptake rate [mol/cell.h]
        # Compute biomass transfer between metabolic states
        Rm_g, Rg_m, Rs_m, Rm_s, Rs_rip = MSMM._Bflux(self, idP, Blist = Blist)
        # Compute biomass variation       
        dBg = Yx * cRate * Bg * (1 - (Bg/K)) + Rm_g - Rg_m - mortality * Bg    
        dBm =  Rg_m + Rs_m - Rm_g - Rm_s - mortality * Bm                     
        dBs =  Rm_s - Rs_m - Rs_rip - mortality * Bs                          
        dBrip = mortality * Btot + Rs_rip  
        dB = [dBg, dBm, dBs, dBrip]
        return dB

    def _Bflux(self, idP,  Blist):
        """
        Function to compute biomass transfer between metabolic states.
        
        Parameters
        ----------
        
        idP : INT
            Parameter index from the plot function (involved parameters depends on envModel)
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
        if len(Blist) != 3:
            print('error in initial biomass, Blist must contain 3 elements') #!!!
            sys.exit()
        #Biomass in each metabolic state
        Bg = Blist[0]
        Bm = Blist[1]
        Bs = Blist[2]
        
        thetaGM = MSMM._stShifts(self, shift = 'GxM')
        thetaMS = MSMM._stShifts(self, shift = 'MxS')
        thetaSRIP = MSMM._stShifts(self, shift = 'S-RIP')
        
        eta = self._specMtbShiftRates[self.mtbRates] 
        Rm_g = Bm * eta * thetaGM[idP]
        Rg_m = Bg * eta * (1 - thetaGM[idP])
        Rs_m = Bs * eta * thetaMS[idP]
        Rm_s = Bm * eta * (1 - thetaMS[idP])
        Rs_rip = Bs * eta * (1 - thetaSRIP[idP])
        Rlist = [Rm_g, Rg_m, Rs_m, Rm_s, Rs_rip]
        return Rlist    #??? change return format 
        
    def _stShifts(self, shift):
        """
        Function to compute shift control between two metabolic states.
        
        Parameters
        ----------
        
        shift : STR
            'GxM' => shift from growth state to maintenance and conversely
            'MxS' => shift from maintenance state to survival and conversely
            'S-RIP' => shift from survival state to death
        
        Returns
        -------
        theta : TYPE
            Metabolic shift control [-]
        """
        # Set requested arguments for CSP.getAllCSP
        if self.envModel in self.atmModels:
            paramDB = self.db.copy()
            typeKin = self.typeKin      # no copy required for immutable data types like floats or strings
            typeMetabo = self.typeMtb
            rxn = self.metabolism
            specComp = self.eD
            C = self.Ct.copy()
            T = self.temperature.copy()
            pH = self.pH
            S = self.salinity
            fluidType = self.fluidType
            DGsynth = self.DGsynth
        else: print('error in CSPargs'), sys.exit() #!!! set getAllCSP args for other models
        
        CSPargs = {'paramDB': paramDB, 
                   'typeKin': typeKin,
                   'typeMetabo': typeMetabo,
                   'reaction': rxn,
                   'specComp': specComp,
                   'Ct': C,
                   'T': T,
                   'pH': pH,
                   'S': S,
                   'phase': 'L', 
                   'sample': 'All',
                   'fluidType': fluidType,
                   'molality': True,
                   'methods': None,
                   'solvent': 'H2O',
                   'asm': 'stoich',
                   'DGsynth': DGsynth}
        
        st = self.st
        #Compute cell specific powers
        Pcat = CSP.getAllCSP(**CSPargs)['Pcat']
        Pm = CSP.getAllCSP(**CSPargs)['Pm0']
        Ps = CSP.getAllCSP(**CSPargs)['Ps']
        Pcell = CSP.getAllCSP(**CSPargs)['Pcell']
        #Compute shift controls (theta)
        if shift == 'GxM':
            theta = 1 / (np.exp((-Pcat + Pcell)/(st * Pcell)) +1)
        elif shift == 'MxS':
            theta = 1 / (np.exp((-Pcat + Pm)/(st * Pm)) +1)
        elif shift == 'S-RIP':
            theta = 1 / (np.exp((-Pcat + Ps)/(st * Ps)) +1)
        else: print('error in itheta value'), sys.exit() #!!!
        return theta

    def solveODE(self, Bini, tSpan, dt = 1, exportBint = False):
        print("debut")
        print(type(dt))
        print(tSpan)
        """ #!!! tooltip
        Function to plot solutions of the MSMM ODE system.
        
        Parameters
        ----------
        Bini : LIST of INT
            Initial biomass in each state (LIST)
        time : LIST or np.array
            Time range over which the microbial dynamic is computed
        exportBint : BOOL
            Command to export integrated biomass values as Excel document.
            Default is False. #!!! add export code
        
        Returns
        -------
        db : LIST
        [...]
        
        """
        if not isinstance(Bini, np.ndarray): Bini = np.array(Bini)
        if len(Bini) != 4:
            print('error in initial biomass, Bini must contain 4 elements') #!!!
            sys.exit()
        # set variables from self.attributes
        if self.envModel in self.atmModels:
            alt = self.altitude.copy()
            isolve = len(alt)
        else: print('envModel not found') #!!! set needed variables for other models
        if not isinstance(tSpan, int): tSpan = int(tSpan)
        #Initialize Bint matrix
        Bint = np.empty(Bini.shape)
        Bint = Bint[..., np.newaxis]
        Bint = np.repeat(Bint, tSpan+1, axis = -1)
        Bint = Bint[..., np.newaxis]
        Bint = np.repeat(Bint, isolve, axis = -1)
        #print('sol matrix shape:', Bint.shape)
        print("test")
        t0 = process_time()
        ODEsol = copy.copy(self._ODE_template)
        print(process_time() - t0)
        
        
        
        print("debut boucle for")
        t0 = process_time()
        for idP in range(isolve):
             ODEsol.set_f_params(idP)
             ODEsol.set_initial_value(Bini, 0)
             #print(vars(ODEsol))
             print(f'alt = {alt[idP]} m') # for envModel from atmModels only
             print("debut boucle while")
             while ODEsol.successful() and ODEsol.t < tSpan:
                 print("while...")
                 t0 = process_time()
                 print(ODEsol.t+dt)
                 sol = ODEsol.integrate(ODEsol.t+dt)
                 print(sol)
                 print(process_time() - t0)
                 #int_status = ODEsol.get_return_code()
                 #if int_status > 0: print('integration success')
                 #elif int_status < 0: print('integration stopped early or failed')
                 time = int(ODEsol.t)
                 #print(time, sol)
                 Bint[:,time,idP] = sol
        print("fin boucle for")
        print(process_time() - t0)
        return Bint
        
    def plotMSMM(self, Bini, idP, time, dt = 1):
        
        """ #!!! tooltip
        Function to plot solutions of the MSMM ODE system.
        
        Parameters
        ----------
        [...]
        
        Returns
        -------
        None (microbial dynamic is plotted, one microbial community at a time)
        """
        # call ODE solving function
        Bplot = self.solveODE(Bini, time, dt)

        if self.envModel in self.atmModels:
            datmMicr = {'CH4': 'Methanotrophs',
                        'H2': 'Hydrogen-oxidizing bacteria',
                        'CO': 'CO-oxidizing bacteria'}
            communityName = datmMicr[self.eD]
            plt.plot(time, Bplot[0,:,idP],'g-', linewidth=2.0)    #growth state curve
            plt.plot(time, Bplot[1,:,idP],'k-', linewidth=2.0)    #maintenance state curve
            plt.plot(time, Bplot[2,:,idP],'b-', linewidth=2.0)    #survival state curve
            plt.plot(time, Bplot[3,:,idP],'r--', linewidth=2.0)   #death state curve
            if self.envModel in self.atmModels:
                if self.envModel == 'ISA':
                    plt.xlabel('time (hours)')
                    plt.ylabel('Cell concentration (cell/m^3 air)')
                    plt.title(f'Dynamic of the {communityName} community at {idP} km altitude')
            else: sys.exit() #!!! labels for other models 
            plt.legend(['Growth', 'Maintenance', 'Survival', 'Dead cells'], bbox_to_anchor = (1.4, 1.0), borderaxespad = 1, title = 'Metabolic states:', title_fontproperties = {'size': 'large', 'weight': 'bold'})
            plt.grid() 
            plt.show()
        return

