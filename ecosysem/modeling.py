# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 11:11:07 2025

@author: zLemaire
"""

# Import Python packages
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

    def __init__(self, envModel, coord, typeMetabo, metabolism, K, mortality, 
                 Wtype = 'L-FW', pH = 7.0, dataType = 'cyly', dataRange = [2020, 2024],
                 DeltaGsynth = 9.54E-11, steepness = 0.2, degradPace = 'moderate',
                 salinity = None, fluidType = 'ideal', actMethods = None):
        
        validModels = {'ISA', 'ISAMERRA2', 'CAMSMERRA2', 'GWB'}
        validMetabo = ['Mth', 'HOB', 'COOB']
        if not isinstance(envModel, str):
            raise ValueError(f'Environment model must be a string. current input: {envModel}')
        if not envModel in validModels:
            raise NameError(f'Invalid model ({self.model}). Valid models: {validModels}.')
        if not isinstance(envModel, str):
            raise ValueError(f'Environment model must be a string. current input: {envModel}')
        self.envModel = envModel
        atmModels = ['ISA', 'ISAMERRA2', 'CAMSMERRA2']
        _metaboProperties = {'fast' : {'protein turnover rate':1 ,'specific metabolic shift rates':1},
                             'moderate': {'protein turnover rate':5 ,'specific metabolic shift rates':0.2},
                             'slow': {'protein turnover rate':14 ,'specific metabolic shift rates':0.071}
                             }  #resp. [h] & [1/h]
        dMtbRates = pd.DataFrame(data = _metaboProperties)
        _specMtbShiftRates = dMtbRates.loc['specific metabolic shift rates']
        _eDonor = {'Mth':'CH4',
                   'HOB': 'H2',
                   'COOB':'CO'}
        micrCommunity = {'CH4': 'Methanotrophs',
                         'H2': 'Hydrogen-oxidizing bacteria',
                         'CO': 'CO-oxidizing bacteria'}
        self.typeMtb = typeMetabo       #metabolism type (STR), e.g. 'AnMetabolisms'
        if not isinstance(metabolism, list): metabolism = [metabolism]
        if not len(metabolism) == 1:
            raise AttributeError(f'A single metabolism name must be given, current input: {metabolism}.')
        if not metabolism[0] in validMetabo:
            raise NameError(f'Invalid metabolism. Valid inputs: {validMetabo}')
        self.metabolism = metabolism    #reaction (STR), e.g. 'Mth'
        self.eD = _eDonor[metabolism[0]]    #(specComp) based on given metabolism
        self.communityName = micrCommunity[self.eD]
        if not isinstance(coord, (list, np.ndarray)): coord = [coord]
        self.coord = coord
        if envModel in atmModels:
            self.ylabel = 'Cell concentration (cell/m³ air)'
            if len(coord) == 1:
                self.title = f"{self.communityName}'s dynamic at {coord[0]}m altitude ({envModel})"
            elif len(coord) == 3:
                self.title = f"{self.communityName}'s dynamic at {coord[0]}m altitude, {coord[1]}LON ; {coord[2]}LAT ({envModel})"
            else : raise AttributeError('Invalid coordinates. Atmospheric models admit vertical (ISA) or 3D position. See README for more details.')   
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
        self.pH = pH                    #(FLOAT)
        self.dataType = dataType        #(STR)
        self.dataRange = dataRange      #(LIST)
        self.K = K                      #carrying capacity (FLOAT)
        if not isinstance(mortality, list): mortality = [mortality]
        if len(mortality) == 1: mortality *= 3
        elif len(mortality) != 3:
            raise AttributeError(f'Mortality rates must be either the same for all 3 states (Growth, Maintenance, Survival) or a list of 3 ordered Floats. Current input: {mortality}.')
        self.mortality = mortality      #(LIST) mortality rates of each metabolic state [1/h]
        self.DGsynth = DeltaGsynth      #cell synthesis required energy [J/cell]
        self.st = steepness             # [-]
        self.mtbRates = degradPace      #'fast', 'moderate' or 'slow'
        self.eta = _specMtbShiftRates[degradPace] #specific metabolic shift rate [1/h]
        self.fluidType = fluidType      #'ideal' or 'non-ideal'
        self.method = actMethods        #None or 'DH-ext' or 'SS'
        self.typeKin = 'MM-Arrhenius'
        self.db = ['MM_AtmMicr', 'ArrhCor_AtmMicr']
        self._callEnvP()
    
    def _callEnvP(self):
        """
        Function to import needed environment data (temperature, pH, etc.) as MSMM attributes.

        """
        #import required ISA attributes:
        if self.envModel == 'ISA':
            if len(self.coord) == 1: alt = self.coord[0]
            else: raise AttributeError(f'A list of one element (selected altitude) should be given as coordinate if envModel is ISA, current coord: {self.coord}.')    
            envArgs = {'layers' : 'All',
                       'phase' : self.Wtype,
                       'H2O' : 0.0,
                       'pH' : self.pH, 
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
            if self.Wtype == 'L-FW':
                self.Ct = ISAinst.Ci_LFW
            elif self.Wtype == 'L-SW':
                self.Ct = ISAinst.Ci_LSW
            ISAinst.salinity = self.salinity     # set to None by default in ISA
            ISAinst.methods = self.method        # set to None by default in ISA
            ISAinst.getCSP(paramDB = self.db, typeKin = self.typeKin, typeMetabo = self.typeMtb,
                          reactions = self.metabolism, specComp = self.eD,
                          sample = 'All', DGsynth = self.DGsynth, EnvAttributes = True)
            self.DGr = ISAinst.DGr[f'{self.metabolism[0]}_pH:{self.pH}'] * 1000 #[J/moleD]
            self.Rs = ISAinst.Rs[f'{self.metabolism[0]}'] / 3600 #[moleD/(cell.s)]
            self.CSP = ISAinst.CSP[f'{self.metabolism[0]}_pH:{self.pH}']
            self.ISAinst = ISAinst
        #import required ISAMERRA2 attributes:
        elif self.envModel == 'ISAMERRA2':
            if len(self.coord) == 3:
                alt = self.coord[0]
                lon = self.coord[1]
                lat = self.coord[2]
            else: raise AttributeError(f'A list of 3 elements (selected altitude, longitude, latitude) should be given as coordinates for ISAMERRA2, current coord: {self.coord}.')    
            envArgs = {'dataType': self.dataType,
                       'y': self.dataRange,
                       'm' : None,
                       'd' : None,
                       'pH' : self.pH,
                       'bbox' : (lon, lat, lon, lat),
                       'compound' : None,
                       'phase' : self.Wtype, 
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
            if self.Wtype == 'L-FW':
                self.Ct = ISAMERRA2inst.Ci_LFW
            elif self.Wtype == 'L-SW':
                self.Ct = ISAMERRA2inst.Ci_LSW
            ISAMERRA2inst.salinity = self.salinity      # set to None by default in ISAMERRA2
            ISAMERRA2inst.methods = self.method         # set to None by default in ISAMERRA2
            ISAMERRA2inst.getCSP(paramDB = self.db, typeKin = self.typeKin, typeMetabo = self.typeMtb,
                           reactions = self.metabolism, specComp = self.eD,
                           sample = 'All', DGsynth = self.DGsynth, EnvAttributes = True)
            self.DGr = ISAMERRA2inst.DGr[f'{self.metabolism[0]}_pH:{self.pH}'] * 1000 #[J/moleD]
            self.Rs = ISAMERRA2inst.Rs[f'{self.metabolism[0]}'] / 3600 #[moleD/(cell.s)]
            self.CSP = ISAMERRA2inst.CSP[f'{self.metabolism[0]}_pH:{self.pH}']
            self.ISAMERRA2inst = ISAMERRA2inst
        #import required CAMSMERRA2 attributes:   
        elif self.envModel == 'CAMSMERRA2':
            if len(self.coord) == 3:
                alt = self.coord[0]
                lon = self.coord[1]
                lat = self.coord[2]
            else: raise AttributeError(f'A list of 3 elements (selected altitude, longitude, latitude) should be given as coordinates for CAMSMERRA2, current coord: {self.coord}.')
            envArgs = {'dataType': self.dataType,
                       'y': self.dataRange,
                       'm' : None,
                       'd' : None,
                       'pH' : self.pH,
                       'bbox' : (lon, lat, lon, lat),
                       'keys' : 'All',
                       'phase' : self.Wtype, 
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
            if self.Wtype == 'L-FW':
                self.Ct = CAMSMERRA2inst.Ci_LFW
            elif self.Wtype == 'L-SW':
                self.Ct = CAMSMERRA2inst.Ci_LSW
            CAMSMERRA2inst.salinity = self.salinity     # set to None by default in CAMSMERRA2
            CAMSMERRA2inst.methods = self.method        # set to None by default in CAMSMERRA2
            CAMSMERRA2inst.getCSP(paramDB = self.db, typeKin = self.typeKin, typeMetabo = self.typeMtb,
                           reactions = self.metabolism, specComp = self.eD,
                           sample = 'All', DGsynth = self.DGsynth, EnvAttributes = True)
            self.DGr = CAMSMERRA2inst.DGr[f'{self.metabolism[0]}_pH:{self.pH}'] * 1000 #[J/moleD]
            self.Rs = CAMSMERRA2inst.Rs[f'{self.metabolism[0]}'] / 3600  #[moleD/(cell.s)]
            self.CSP = CAMSMERRA2inst.CSP[f'{self.metabolism[0]}_pH:{self.pH}']
            self.CAMSMERRA2inst = CAMSMERRA2inst
        if self.envModel == "GWB":
            raise AttributeError('Code part dedicated to general water body was not written yet.')

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
            Biomass variation [cell/h] in each metabolic state (Growth, Maintenance, Survival, Death).
        
        """
        Bg = y[0]
        Bm = y[1]
        Bs = y[2]
        Blist = [Bg, Bm, Bs]
        Btot = sum(Blist)
        
        #import self.attributes
        mortality = self.mortality.copy()
        mG = mortality[0]   #mortality in growth state
        mM = mortality[1]   #mortality in maintenance state
        mS = mortality[2]   #mortality in survival state
        K = self.K
        DGr = self.DGr.copy()
        Rs = self.Rs.copy()
        Yx = -(DGr * (0.5 / 1.04e-10))      # cell growth yield [cell/mol eD]
        
        # Compute biomass transfer between metabolic states
        Rm_g, Rg_m, Rs_m, Rm_s, Rs_rip = MSMM._Bflux(self, Blist)
        # Compute biomass variation       
        dBg = Yx * Rs * Bg * (1 - (Btot / K)) + Rm_g - Rg_m - mG * Bg    
        dBm = Rg_m + Rs_m - Rm_g - Rm_s - mM * Bm                     
        dBs = Rm_s - Rs_m - Rs_rip - mS * Bs                          
        dBrip = mG * Bg + mM * Bm + mS * Bs + Rs_rip  
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
        eta = self.eta
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
        Function to compute shift controls between metabolic states.
                   
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
        dt : INT or FLOAT, optional (default : 1h)
            time step for the integration
        solExport : BOOL, optional (default : False)
            Command to export integrated biomass values as Excel document if set to True.
        
        Returns
        -------
        
        None 
        ODE solutions (numpy.ndarray of shape [4, tSpan+1]) are saved as MSMM attribute ('Bsol')
        If solExport is set to True, creates an Excel document of the results.
        
        """
        # check Bini
        if not isinstance(Bini, np.ndarray): Bini = np.array(Bini)
        if len(Bini) != 4:
            raise ValueError(f'Bini must contain 4 elements, current length: {len(Bini)}.')
        # create time array for later plotting
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
        # compute solutions over given time range
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
        with pd.ExcelWriter(self.fullPathSave, engine='openpyxl', mode = 'a', if_sheet_exists='overlay') as writer:
            time.to_excel(writer, sheet_name = nameSheet_B, startrow = 2, startcol = 1, index = False, header = True)
            Bdf.to_excel(writer, sheet_name = nameSheet_B, startrow = 2, startcol = 2, index = False, header = True)    
        
    def plotMSMM(self):
        
        """
        Function to plot MSMM microbial dynamic of a single point in the environment space.
        
        """
        # import solutions of the ODE and time array from MSMM attributes
        Bplot = getattr(self, 'Bsol', None)
        if Bplot is None:
            raise AttributeError('MSMM attribute "Bsol" could not be found. Please first use MSMM.solveODE().')
        # plotting of metabolic state curves
        plt.plot(self.t_plot, Bplot[0,:],'g-', linewidth=2.0)    #growth state curve
        plt.plot(self.t_plot, Bplot[1,:],'k-', linewidth=2.0)    #maintenance state curve
        plt.plot(self.t_plot, Bplot[2,:],'b-', linewidth=2.0)    #survival state curve
        plt.plot(self.t_plot, Bplot[3,:],'r--', linewidth=2.0)   #death state curve
        plt.xlabel('time (hours)')
        plt.ylabel(self.ylabel)
        plt.title(self.title)
        plt.legend(['Growth', 'Maintenance', 'Survival', 'Dead cells'], bbox_to_anchor = (1.42, 1.0), borderaxespad = 1, title = 'Metabolic states:', title_fontproperties = {'size': 'large', 'weight': 'bold'})
        plt.grid() 
        plt.show()
        return