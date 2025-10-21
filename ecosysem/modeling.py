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

# Import classes from modules (as abb.) or import with importlib
from envdef import ISA
from reactions import KinRates as KR
from thermodynamics import ThSA
from bioenergetics import CSP


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



    def __init__(self, typeMetabo, metabolism, eDonor, Wtype, K, mortality, envModel,
                 envArgs = None, DeltaGsynth = 9.54E-11, steepness = 0.2,
                 degradPace = 'moderate', fluidType = 'ideal', actMethods = None):
        
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
                self.pH = ISAinst.pH.copy()
                self.compositions = ISAinst.compositions.copy()
                self.compounds = ISAinst.compounds.copy()
                self.resolution = ISAinst.resolution.copy()
                self.temperature = ISAinst.temperature.copy()
                self.pressure = ISAinst.pressure.copy()
                self.altitude = ISAinst.altitude.copy()
                self.H2O = ISAinst.H2O.copy()
                self.Pi = ISAinst.Pi.copy()           #!!! could be useful for MSMMv2
                self.Ci_G = ISAinst.Ci_G.copy()        #!!! //
                if self.Wtype == 'L_FW':
                    self.Ct = ISAinst.Ci_LFW.copy()
                elif self.Wtype == 'L_SW':
                    self.Ct = ISAinst.Ci_LSW.copy()
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

    def _ODEsystem_MSMM(self, y, t, idP):
        """
        Function for the differential equations system of the model.
        
        Parameters
        ----------
        
        y : LIST of INT
            Initial biomass in each metabolic state, e.g. [cell/m^3 air]
        t : LIST or np.array
            Time range over which biomass variation is computed
        idP : INT
            Parameter index from the plot function (involved parameters depends on envModel)
            
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
        mortality = self.mortality.copy()
        K = self.K.copy()
        typeKin = self.typeKin.copy()
        paramDB = self.db.copy()
        typeRxn = self.typeMtb.copy()
        rxn = self.metabolism.copy()
        specComp = self.eD.copy()
        T = self.temperature.copy()
        pH = self.pH.copy()
        S = self.salinity.copy()
        C = self.Ct.copy()
        fluidType = self.fluidType.copy()
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
        
        eta = self._specMtbShiftRates[self.mtbRates] 
        Rm_g = Bm * eta * MSMM._stShifts(self, shift = 'GxM')[idP]
        Rg_m = Bg * eta * (1 - MSMM._stShifts(self, shift = 'GxM')[idP])
        Rs_m = Bs * eta * MSMM._stShifts(self, shift = 'MxS')[idP]
        Rm_s = Bm * eta * (1 - MSMM._stShifts(self, shift = 'MxS')[idP])
        Rs_rip = Bs * eta * (1 - MSMM._stShifts(self, shift = 'S-RIP')[idP])
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
            typeKin = self.typeKin.copy()
            typeMetabo = self.typeMtb.copy()
            rxn = self.metabolism.copy()
            specComp = self.eD.copy()
            C = self.Ct.copy()
            T = self.temperature.copy()
            pH = self.pH.copy()
            S = self.salinity.copy()
            fluidType = self.fluidType.copy()
            DGsynth = self.DGsynth.copy()
        else: print('error in CSPargs'), sys.exit() #!!!
        
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
        #!!! set getAllCSP args for other models
        st = self.st.copy()
        Pcat = CSP.getAllCSP(**CSPargs)['Pcat']
        Pm = CSP.getAllCSP(**CSPargs)['Pm0']
        Ps = CSP.getAllCSP(**CSPargs)['Ps']
        Pcell = CSP.getAllCSP(**CSPargs)['Pcell']
        
        if shift == 'GxM':
            theta = 1 / (np.exp((-Pcat + Pcell)/(st * Pcell)) +1)
        elif shift == 'MxS':
            theta = 1 / (np.exp((-Pcat + Pm)/(st * Pm)) +1)
        elif shift == 'S-RIP':
            theta = 1 / (np.exp((-Pcat + Ps)/(st * Ps)) +1)
        else: print('error in itheta value'), sys.exit() #!!!
        return theta

    def solveODE(self, time, Bini, exportBint = False):
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
            Default is False.
        
        Returns
        -------
        db : LIST
        [...]
        
        """
        #ODEsol = ode(self._ODEsystem_MSMM)

        if len(Bini) != 4:
            print('error in initial biomass, Bini must contain 4 elements') #!!!
            sys.exit()
        # set variables from self.attributes
        if self.envModel in self.atmModels:
            alt = self.altitude.copy()
            isolve = len(alt)
        else: print('envModel not found') #!!! set needed variables for other models
        if not isinstance(time, np.ndarray): time = np.array(time)
        # Initialize Bint matrix
        Bint = np.empty(Bini.shape)
        Bint = Bint[..., np.newaxis]
        Bint = np.repeat(Bint, len(time), axis = -1)
        Bint = Bint[..., np.newaxis]
        Bint = np.repeat(Bint, isolve, axis = -1)

        for i in range(isolve):
            ode.set_f_params(i)
            Bint[..., i] = ode.integrate(time) #step = True ?
        return Bint
        
    def plotMSMM(self, time, Bini):
        
        """ #!!! tooltip
        Function to plot solutions of the MSMM ODE system.
        
        Parameters
        ----------
        [...]
        
        Returns
        -------
        None (microbial dynamic is plotted, one microbial community at a time)
        """
        # call ODEsolve fct
        # if self.envModel in self.atmModels:
        #     datmMicr = {'CH4': 'Methanotrophs',
        #                 'H2': 'Hydrogen-oxidizing bacteria',
        #                 'CO': 'CO-oxidizing bacteria'}
        #     communityName = datmMicr[self.eD]
            # plt.plot(time, Bplot[i][:,0],'g-', linewidth=2.0)    #growth state curve
            # plt.plot(time, Bplot[i][:,1],'k-', linewidth=2.0)    #maintenance state curve
            # plt.plot(time, Bplot[i][:,2],'b-', linewidth=2.0)    #survival state curve
            # plt.plot(time, Bplot[i][:,3],'r--', linewidth=2.0)   #death state curve
            # if self.envModel in self.atmModels:
            #     if self.envModel == 'ISA':
            #         plt.xlabel('time (hours)')
            #         plt.ylabel('Cell concentration (cell/m^3 air)')
            #         plt.title(f'Dynamic of the {communityName} community at {(alt[i]/1000)} km altitude')
            # else: sys.exit() #!!! labels for other models 
            # plt.legend(['Growth', 'Maintenance', 'Survival', 'Dead cells'], bbox_to_anchor = (1.4, 1.0), borderaxespad = 1, title = 'Metabolic states:', title_fontproperties = {'size': 'large', 'weight': 'bold'})
            # plt.grid() 
            # plt.show()


