# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 11:11:07 2025

@author: zLemaire
"""

# Import Python packages
import copy
import pandas as pd
import numpy as np
import os.path
import matplotlib.pyplot as plt
from scipy.integrate import ode

# Import environment classes
from environments import ISA, ISAMERRA2, CAMSMERRA2

class MSMM:
    """
    Class for Multi-State Metabolic Model
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
        return Rlist
    
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

    def solveODE(self, Bini, tSpan, dt = 1, solExport = False):
        """
        Function to solve the MSMM ODE system and export the results as .xlsx document.
        
        Parameters
        ----------
        Bini : LIST of INT
            Initial biomass in each state (LIST)
        tSpan : LIST or np.array
            Time range over which the microbial dynamic is computed, in hours
        dt : INT or FLOAT
            time step for the integration (default : 1 hour)
        solExport : BOOL
            Command to export integrated biomass values as Excel document.
            Default is False.
        
        Returns
        -------
        db : LIST
        [...]
        
        """
        if not isinstance(Bini, np.ndarray): Bini = np.array(Bini)
        if len(Bini) != 4:
            raise ValueError(f'Bini must contain 4 elements, current length: {len(Bini)}.')
        self.t_plot = np.linspace(0, tSpan, int(tSpan/dt)+1)
        #Initialize Bint matrix
        Bint = np.empty(4)
        Bint = Bint[..., np.newaxis]
        Bint = np.repeat(Bint, tSpan+1, axis = -1)
        Bint[:,0] = Bini
        #create ode instance & set initial values & integration method
        ODEsol = ode(self._ODEsystem_MSMM)
        ODEsol.set_initial_value(Bini, 0)
        ODEsol.set_integrator('vode', method='adams')
        while ODEsol.successful() and ODEsol.t < tSpan:
             sol = ODEsol.integrate(ODEsol.t+dt)
             time = int(ODEsol.t)
             Bint[:,time] = sol
        # save ODE solutions as MSMM attribute (rounded values)
        self.Bsol = np.round(Bint, 2)
        # export ODE solutions as .xlsx document
        if solExport == True:
            path = 'results/'
            nameDocument = input(' > Name of result document: ')
            self.fullPathSave = path + nameDocument + '.xlsx'
            if os.path.isfile(self.fullPathSave):
                val = input(' > '+ nameDocument + '.xlsx already exists in this directory. /!\ Make sure no instance of the file is currently open. Do you want to overwrite `' + nameDocument + '.xlsx`? [Y/N]: ')
                if val == 'Y' or val == 'y':       
                    os.remove(self.fullPathSave)
            MSMM._writeExcel(self)
    
    def _writeExcel(self):
        """
        Write calculated metabolic state biomass in Excel document.

        """
        nameSheet_B = 'MSMM biomass'
        # import solutions of the ODE and time array from MSMM attributes
        time = pd.DataFrame(self.t_plot, columns = ['time (h)| states :'])
        Bdf = pd.DataFrame({self.Bstates[i]: self.Bsol[i] for i in range(4)})
        # adapt header to environment model
        if self.envModel in self.atmModels:
            alt = self.coord[0]
            if self.envModel == 'ISA':
                introRowB = pd.DataFrame(np.array([f'States biomass [cell/m³ air] | Metabolism: {self.metabolism[0]} | Altitude: {alt}m']))
            else:
                lon = self.coord[1]
                lat = self.coord[2]
                introRowB = pd.DataFrame(np.array([f'States biomass [cell/m³ air] | Metabolism: {self.metabolism[0]} | Coordinates: {lon}LON;{lat}LAT | Altitude: {alt}m']))
        # write excel document 
        if not os.path.isfile(self.fullPathSave):
            with pd.ExcelWriter(self.fullPathSave) as writer:
                introRowB.to_excel(writer, sheet_name = nameSheet_B, index = False, header = False)
        #!!! elif: introRowB for other models 
        with pd.ExcelWriter(self.fullPathSave, engine='openpyxl', mode = 'a', if_sheet_exists='overlay') as writer:
            time.to_excel(writer, sheet_name = nameSheet_B, startrow = 2, startcol = 1, index = False, header = True)
            Bdf.to_excel(writer, sheet_name = nameSheet_B, startrow = 2, startcol = 2, index = False, header = True)    
        
    def plotMSMM(self):
        
        """ #!!! tooltip
        Function to plot MSMM microbial dynamic of a single point in the environment space.
        
        """
        Bplot = getattr(self, 'Bsol', None)
        if Bplot is None:
            raise AttributeError('MSMM attribute "Bsol" could not be found. Please first use MSMM.solveODE().')
        t_array = self.t_plot.copy()
        # plotting of metabolic state curves
        plt.plot(t_array, Bplot[0,:],'g-', linewidth=2.0)    #growth state curve
        plt.plot(t_array, Bplot[1,:],'k-', linewidth=2.0)    #maintenance state curve
        plt.plot(t_array, Bplot[2,:],'b-', linewidth=2.0)    #survival state curve
        plt.plot(t_array, Bplot[3,:],'r--', linewidth=2.0)   #death state curve
        # adapt labels to environment model
        if self.envModel in self.atmModels:
            datmMicr = {'CH4': 'Methanotrophs',
                        'H2': 'Hydrogen-oxidizing bacteria',
                        'CO': 'CO-oxidizing bacteria'}
            communityName = datmMicr[self.eD]
            plt.xlabel('time (hours)')
            plt.ylabel('Cell concentration (cell/m³ air)')
            plt.title(f"{communityName}'s dynamic at {self.coord[0]}m altitude ({self.envModel})")
            plt.legend(self.Bstates, bbox_to_anchor = (1.42, 1.0), borderaxespad = 1, title = 'Metabolic states:', title_fontproperties = {'size': 'large', 'weight': 'bold'})
        #!!! elif: labels for other models 
        plt.grid() 
        plt.show()
        return