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
    
    def _callEnvP(self, envModel):
        """
        Function to import needed environment data (temperature, pH, etc.) 
        
        Parameters
        ----------
        envModel : STR
            Environment model from which data are extracted.
                        
        Returns
        -------
        None (environment data set as MSMM attributes)
        """
        envModel = self.envModel
        atmModels = self.atmModels
        rxn = self.metabolism
        phase_ = self.Wtype
        pH_ = self.pH
        coord_ = self.coord.copy()
        dataType_ = self.dataType #??? if other than 'yly', what about dataRange ?
        dataRange_ = self.dataRange #??? dataRange other than years only ?
        #!!! other models?
        if not isinstance(envModel, str):
            raise ValueError(f'Environment model must be a string. current input: {envModel}')
        # check model and import its attributes:
        if envModel in atmModels :
            #import required ISA attributes:
            if envModel == 'ISA':
                if len(coord_) == 1: alt = coord_[0]
                else: raise AttributeError(f'A list of one element (selected altitude) should be given as coordinate if envModel is ISA, current coord: {coord_}.')    
                envArgs = {'layers' : 'All',
                           'phase' : phase_,
                           'H2O' : 0.0,
                           'pH' : pH_, 
                           'selCompounds' : None,
                           'selAlt' : [alt, alt],
                           'resolution' : 1000
                            }
                ISAinst = ISA(**envArgs)
                self.compositions = ISAinst.compositions
                self.compounds = ISAinst.compounds
                self.temperature = ISAinst.temperature
                self.pressure = ISAinst.pressure
                self.H2O = ISAinst.H2O
                if phase_ == 'L-FW':
                    self.Ct = ISAinst.Ci_LFW
                elif phase_ == 'L-SW':
                    self.Ct = ISAinst.Ci_LSW
                ISAinst.salinity = self.salinity     # set to None by default in ISA
                ISAinst.methods = self.method        # set to None by default in ISA
                ISAinst.getCSP(paramDB = self.db, typeKin = self.typeKin, typeMetabo = self.typeMtb,
                              reactions = self.metabolism, specComp = self.eD,
                              sample = 'All', DGsynth = 9.54E-11, EnvAttributes = True)
                self.DGr = ISAinst.DGr[f'{rxn[0]}_pH:{pH_}'] * 1000 #[J/moleD]
                self.Rs = ISAinst.Rs[f'{rxn[0]}'] / 3600 #[moleD/(cell.s)]
                self.CSP = ISAinst.CSP[f'{rxn[0]}_pH:{pH_}']
                self.ISAinst = ISAinst
            #import required ISAMERRA2 attributes:
            elif envModel == 'ISAMERRA2':
                if len(coord_) == 3:
                    alt = coord_[0]
                    lon = coord_[1]
                    lat = coord_[2]
                else: raise AttributeError(f'A list of 3 elements (selected altitude, longitude, latitude) should be given as coordinates for ISAMERRA2, current coord: {coord_}.')    
                envArgs = {'dataType': dataType_,
                           'y': dataRange_,
                           'm' : None,
                           'd' : None,
                           'pH' : pH_,
                           'bbox' : (lon, lat, lon, lat),
                           'compound' : None,
                           'phase' : phase_, 
                           'altArray' : [alt],
                           'numAlt' : 50,
                           'surftrop' : None,
                           'keysAsAttributes' : False, 
                           'showMessage' : True
                           }
                ISAMERRA2inst = ISAMERRA2(**envArgs)
                self.compositions = ISAMERRA2inst.compositions
                self.compounds = ISAMERRA2inst.compounds
                self.temperature = ISAMERRA2inst.temperature
                self.pressure = ISAMERRA2inst.pressure
                if phase_ == 'L-FW':
                    self.Ct = ISAMERRA2inst.Ci_LFW
                elif phase_ == 'L-SW':
                    self.Ct = ISAMERRA2inst.Ci_LSW
                ISAMERRA2inst.salinity = self.salinity      # set to None by default in ISAMERRA2
                ISAMERRA2inst.methods = self.method         # set to None by default in ISAMERRA2
                ISAMERRA2inst.getCSP(paramDB = self.db, typeKin = self.typeKin, typeMetabo = self.typeMtb,
                               reactions = self.metabolism, specComp = self.eD,
                               sample = 'All', DGsynth = 9.54E-11, EnvAttributes = True)
                self.DGr = ISAMERRA2inst.DGr[f'{rxn[0]}_pH:{pH_}'] * 1000 #[J/moleD]
                self.Rs = ISAMERRA2inst.Rs[f'{rxn[0]}'] / 3600 #[moleD/(cell.s)]
                self.CSP = ISAMERRA2inst.CSP[f'{rxn[0]}_pH:{pH_}']
                self.ISAMERRA2inst = ISAMERRA2inst
            elif envModel == 'CAMSMERRA2':
                if len(coord_) == 3:
                    alt = coord_[0]
                    lon = coord_[1]
                    lat = coord_[2]
                else: raise AttributeError(f'A list of 3 elements (selected altitude, longitude, latitude) should be given as coordinates for CAMSMERRA2, current coord: {coord_}.')
                envArgs = {'dataType': dataType_,
                           'y': dataRange_,
                           'm' : None,
                           'd' : None,
                           'pH' : pH_,
                           'bbox' : (lon, lat, lon, lat),
                           'keys' : 'All',
                           'phase' : phase_, 
                           'altArray' : [alt],
                           'numAlt' : 50,
                           'surftrop' : None,
                           'keysAsAttributes' : False, 
                           'showMessage' : True
                            }
                CAMSMERRA2inst = CAMSMERRA2(**envArgs)
                self.compounds = CAMSMERRA2inst.compounds
                self.temperature = CAMSMERRA2inst.temperature
                self.pressure = CAMSMERRA2inst.pressure
                if phase_ == 'L-FW':
                    self.Ct = CAMSMERRA2inst.Ci_LFW
                elif phase_ == 'L-SW':
                    self.Ct = CAMSMERRA2inst.Ci_LSW
                CAMSMERRA2inst.salinity = self.salinity     # set to None by default in CAMSMERRA2
                CAMSMERRA2inst.methods = self.method        # set to None by default in CAMSMERRA2
                CAMSMERRA2inst.getCSP(paramDB = self.db, typeKin = self.typeKin, typeMetabo = self.typeMtb,
                               reactions = self.metabolism, specComp = self.eD,
                               sample = 'All', DGsynth = 9.54E-11, EnvAttributes = True)
                self.DGr = CAMSMERRA2inst.DGr[f'{rxn[0]}_pH:{pH_}'] * 1000 #[J/moleD]
                self.Rs = CAMSMERRA2inst.Rs[f'{rxn[0]}'] / 3600  #[moleD/(cell.s)]
                self.CSP = CAMSMERRA2inst.CSP[f'{rxn[0]}_pH:{pH_}']
                self.CAMSMERRA2inst = CAMSMERRA2inst
        #elif envModel in ___ :    #!!! import needed attributes for other models (e.g. GWB)
        else:
            raise NameError(f'Environment model ({envModel}) not found, see valid models in the README.')

    def _ODEsystem_MSMM(self, t, y):
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
        mortality = self.mortality  #??? no copy required for immutable data types like floats or strings
        K = self.K
        DGr = self.DGr.copy()
        Rs = self.Rs.copy()
        Yx = -(DGr * (0.5 / 1.04e-10))      # cell growth yield [cell/mol eD]
        
        # Compute biomass transfer between metabolic states
        Rm_g, Rg_m, Rs_m, Rm_s, Rs_rip = MSMM._Bflux(self, Blist)
        # Compute biomass variation       
        dBg = Yx * Rs * Bg * (1 - (Bg / K)) + Rm_g - Rg_m - mortality * Bg    
        dBm = Rg_m + Rs_m - Rm_g - Rm_s - mortality * Bm                     
        dBs = Rm_s - Rs_m - Rs_rip - mortality * Bs                          
        dBrip = mortality * Btot + Rs_rip  
        dB = [dBg, dBm, dBs, dBrip]
        return dB

    def _stShifts(self):
        """
        Function to compute shift control between two metabolic states.
        
        Returns
        -------
        Save a dictionary of metabolic shift controls as attribute of MSMM.
            
        """
        CSPdict = self.CSP.copy()
        st = self.st
        #Compute cell specific powers
        Pcat = CSPdict['Pcat']
        Pm = CSPdict['Pm0']
        Ps = CSPdict['Ps']
        Pcell = CSPdict['Pcell']
        #Initialize thetaDict before shift controls (theta) calculations
        thetaDict = {}
        thetaDict['GxM'] = 1 / (np.exp((-Pcat + Pcell)/(st * Pcell)) +1)
        thetaDict['MxS'] = 1 / (np.exp((-Pcat + Pm)/(st * Pm)) +1)
        thetaDict['S-RIP'] = 1 / (np.exp((-Pcat + Ps)/(st * Ps)) +1)
        self.theta = thetaDict

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
        Rlist : LIST
             List of computed biomass transfer [cell/h] for each kind of metabolic shift:
                 - Rm_g : transfer from maintenance to growth
                 - Rg_m : transfer from growth to maintenance
                 - ...
                 - Rs_rip : transfer from survival to dead cells
        """
        if len(Blist) != 3:
            raise ValueError(f'Blist must contain 3 elements (Bg, Bm, Bs), current Blist: {Blist}.')
        #Create shift control dict in MSMM attributes
        self._stShifts()
        #Biomass in each metabolic state
        Bg = Blist[0] 
        Bm = Blist[1]
        Bs = Blist[2]
        #import metabolic shift controls and rates
        theta = self.theta.copy()
        eta = self._specMtbShiftRates[self.mtbRates]
        #compute biomass transfers
        Rm_g = Bm * eta * theta['GxM']
        Rg_m = Bg * eta * (1 - theta['GxM'])
        Rs_m = Bs * eta * theta['MxS']
        Rm_s = Bm * eta * (1 - theta['MxS'])
        Rs_rip = Bs * eta * (1 - theta['S-RIP'])
        Rlist = [Rm_g, Rg_m, Rs_m, Rm_s, Rs_rip]
        
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

