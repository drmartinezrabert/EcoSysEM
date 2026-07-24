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
                 DeltaGsynth = 9.54E-11, steepness = 0.2, salinity = None,
                 Wtype = 'L-FW', pH = 7.0, Wcontent = 0.0,  fluidType = 'ideal',
                 actMethods = None, molality = True, asm = 'stoich',
                 dataType = None, years = None, month = None, day = None,
                 turnoverRate = {'fast' : 1,'moderate': 5 ,'slow': 14}, degradPace = 'moderate',
                 kinDB = {'MM-Arrhenius': ['qs_FFAM', 'ArrhCor'], 'MM': ['qs_FFAM']},
                 typeKin = 'MM-Arrhenius', eD = {'Mth':'CH4', 'HOB': 'H2', 'COOB':'CO'},
                 microCommunity = {'Mth': 'Methanotrophs','HOB': 'Hydrogen-oxidizing bacteria','COOB': 'CO-oxidizing bacteria'}):
        validModels = {'ISA', 'ISAMERRA2', 'CAMSMERRA2', 'GWB'}
        validMetabo = ['Mth', 'HOB', 'COOB']
        if not isinstance(envModel, str):
            raise ValueError(f'Environment model must be a string. current input: {envModel}')
        if not envModel in validModels:
            raise NameError(f'Invalid model ({self.envModel}). Valid models: {validModels}.')
        self.envModel = envModel
        atmModels = ['ISA', 'ISAMERRA2', 'CAMSMERRA2']
        if not isinstance(degradPace, str):
            if isinstance(degradPace, (int, float)):
                self.specMSrate = 1 / degradPace #specific metabolic shift rate [1/h]
            else : raise TypeError(f'Degradation pace must be a float/int or str. Current type : {type(degradPace)}.')
        elif isinstance(degradPace, str):
            if not degradPace in ['fast', 'moderate', 'slow'] :
                raise NameError('Invalid str input for degradPace. Valid inputs : fast, moderate, slow.')
            else: 
                self.specMSrate = 1 / (turnoverRate[degradPace]) #specific metabolic shift rate [1/h]
        self.mtbRates = degradPace      #'fast', 'moderate', 'slow' or float/int.
        if not isinstance(typeMetabo, str):
            raise TypeError(f'typeMetabo must be a str. Current type: {type(typeMetabo)}.')
        self.typeMtb = typeMetabo       #metabolism type (STR), e.g. 'AnMetabolisms'
        if not isinstance(metabolisms, list): metabolisms = [metabolisms]
        for metabolism in metabolisms:
            if not metabolism in validMetabo:
                raise NameError(f'Invalid metabolism ({metabolism}). Valid inputs: {validMetabo}')
        self.metabolisms = metabolisms   #reaction (LIST), e.g. ['Mth']
        for metabolism in self.metabolisms:
            if eD.get(metabolism, None) == None:
                raise AttributeError(f'No {metabolism} key could be found in the eD dictionary. Please modify the corresponding argument.')
        self.eD = {metabolism: eD[metabolism] for metabolism in self.metabolisms}  #(specComp) based on given metabolisms
        self.compounds = [eD.get(metabolism) for metabolism in self.metabolisms] + ['O2','CO2']
        for metabolism in self.metabolisms:
            if microCommunity.get(metabolism, None) == None:
                raise AttributeError(f'No {metabolism} key could be found in the microCommunity dictionary. Please modify the corresponding argument.')
        self.communityNames = {metabolism: microCommunity[metabolism] for metabolism in self.metabolisms}
        if not isinstance(coord, (list, np.ndarray)): coord = [coord]
        self.coord = coord
        if envModel in atmModels:
            self.plotYlabel = 'Cell concentration (cell/m³ air)'
            if len(coord) == 1:
                self.plotTitle = f"Community dynamic at {coord[0]}m altitude ({envModel})"
            elif len(coord) == 3:
                self.plotTitle = f"Community dynamic at {coord[0]}m altitude, {coord[1]}LON ; {coord[2]}LAT ({envModel})"
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
                raise TypeError(f'Given salinity value must be float or int, current type: {type(salinity[0])}.')
        else: salinity = [0.0]
        if isinstance(pH, list):
            if not len(pH) == 1:
                raise AttributeError(f'A single pH value must be given, current pH: {pH}.')
            else: pH = pH[0]
        if not isinstance(pH, (float, int)):
            raise TypeError(f'Given pH value must be float or int, current type: {type(pH)}.')
        self.dataType = dataType        #(STR)
        self.dataYear = years #(INT or LIST)
        self.dataMonth = month #(INT)
        self.dataDay = day #(INT)
        if not isinstance(K, (float,int)):
            raise TypeError(f'Carrying capacity (K) must be a float or an int. Current type: {type(K)}.')
        self.K = K     #(INT or FLOAT), carrying capacity [cell/unit volume]
        if not isinstance(mortality, list): mortality = [mortality]
        if len(mortality) == 1: mortality *= 2
        elif len(mortality) != 2:
            raise AttributeError(f'Mortality rates must be either the same for both active states (Growth, Maintenance) or a list of 2 ordered Floats. Current input: {mortality}.')
        self.mortality = mortality      #(LIST) mortality rates of each metabolic state [1/h]
        if not isinstance(DeltaGsynth, (float,int)):
            raise TypeError(f'DeltaGsynth must be a float or an int. Current type: {type(DeltaGsynth)}.')
        self.DGsynth = DeltaGsynth      #cell synthesis required energy [J/cell]
        if not isinstance(steepness, list): steepness = [steepness]
        if not all(isinstance(k, (float,int)) for k in steepness):
            raise TypeError(f'Steepness must be a float or an int. Current type: {type(steepness)}.')
        if len(steepness) == 1: steepness *= 2
        elif len(steepness) != 2:
            raise AttributeError(f'Steepness parameter in the shift control functions must be either the same for both types of shift (GxM, M-RIP) or a list of 2 ordered Floats. Current input: {steepness}.')
        self.st = steepness             # [-]
        if not isinstance(fluidType, str):
            raise TypeError(f'fluidType must be a str. Current type: {type(fluidType)}.')
        if not fluidType == 'ideal' and not fluidType == 'non-ideal':
            raise AttributeError(f'Invalid input for fluidType. valid inputs: ideal, non-ideal. current input: {fluidType}.')
        self.fluidType = fluidType      #'ideal' or 'non-ideal'
        if not isinstance(typeKin, str):
            raise TypeError(f'typeKin must be a str. Current type: {type(typeKin)}.')
        if not isinstance(typeKin, str):
            raise TypeError(f'Argument typeKin must be a str. Current input: {type(typeKin)}.')
        self.typeKin = typeKin
        if kinDB.get(self.typeKin, None) == None:
            raise AttributeError(f'No {self.typeKin} key could be found in the kinDB dictionary. Please modify the corresponding argument.')
        self.kinDB = kinDB[self.typeKin]
        if not isinstance(Wcontent, (float, int)):
            raise TypeError(f'Wcontent must be an int or a float. current type: {type(Wcontent)}.')
        if isinstance(actMethods, str):
            if not actMethods == 'DH-ext' and not actMethods == 'SS':
                raise NameError(f'Str input for actMethod must be "DH-ext" or "SS". Current input: {actMethods}.')
        elif actMethods != None:
            raise TypeError(f'Argument actMethod must be a str ("DH-ext" or "SS") or None. Current type: {type(actMethods)}.')
        if not isinstance(molality, bool):
            raise TypeError(f'Argument molality must be a bool (True). Current type: {type(molality)}.')
        if not molality == True:
            raise AttributeError(f'Currently only admitted input for molality is True. invalid input: {molality}')
        if not asm == 'stoich':
            raise AttributeError(f'Currently only accepted input for asm is "stoich". Invalid input was given: {asm}.')
        self._callEnvP(salinity, pH, Wcontent, actMethods, molality, asm)
    
    def _callEnvP(self, salinity_, pH_, H2O_, method_, molality_, asm_):
        """
        Function to import needed environment data (temperature, pH, etc.) as MSMM attributes.

        """
        #import required ISA attributes:
        if self.envModel == 'ISA':
            if len(self.coord) == 1: alt = self.coord[0]
            else: raise AttributeError(f'A list of one element (selected altitude) should be given as coordinate if envModel is ISA, current coord: {self.coord}.')    
            ISAinst = ISA(layers = 'All',
                           phase = self.Wtype,
                           H2O = H2O_,
                           pH = pH_, 
                           selCompounds = None,
                           selAlt = [alt, alt],
                           resolution = 1000)
            ISAinst.salinity = salinity_     # set to None by default in ISA
            ISAinst.methods = method_        # set to None by default in ISA
            ISAinst.getCSP(paramDB = self.kinDB, typeKin = self.typeKin, 
                          typeMetabo = self.typeMtb, reactions = self.metabolisms, 
                          specComp = [self.eD.get(metabolism) for metabolism in self.metabolisms],
                          sample = 'All', DGsynth = self.DGsynth, solvent = 'H2O', 
                          molality = molality_, asm = asm_, phase = 'L-FW')
            self.DGr = {self.metabolisms[i] : ISAinst.DGr[f'{self.metabolisms[i]}_pH:{pH_}'] * 1000 
                        for i in range(len(self.metabolisms))} # Gibbs free energy [J/moleD]
            self.Rs = {self.metabolisms[i] : ISAinst.Rs[f'{self.metabolisms[i]}'] / 3600 
                        for i in range(len(self.metabolisms))} # reaction rates [moleD/(cell.s)]
            self.CSP = {self.metabolisms[i] : ISAinst.CSP[f'{self.metabolisms[i]}_pH:{pH_}'] 
                        for i in range(len(self.metabolisms))} # CSP [fW/cell]
            self.envConditions = ISAinst
        #import required ISAMERRA2 attributes:
        elif self.envModel == 'ISAMERRA2':
            if len(self.coord) == 3:
                alt = self.coord[0]
                lon = self.coord[1]
                lat = self.coord[2]
            else: raise AttributeError(f'A list of 3 elements (selected altitude, longitude, latitude) should be given as coordinates for ISAMERRA2, current coord: {self.coord}.')    
            ISAMERRA2inst = ISAMERRA2(dataType = self.dataType,
                                       y = self.dataYear,
                                       m = self.dataMonth,
                                       d = self.dataDay,
                                       pH = pH_,
                                       bbox = (lon, lat, lon, lat),
                                       compound = None,
                                       phase = self.Wtype, 
                                       altArray = [alt],
                                       numAlt = 50,
                                       surftrop = None,
                                       keysAsAttributes = False, 
                                       showMessage = True)
            ISAMERRA2inst.salinity = salinity_      # set to None by default in ISAMERRA2
            ISAMERRA2inst.methods = method_         # set to None by default in ISAMERRA2
            ISAMERRA2inst.getCSP(paramDB = self.kinDB, typeKin = self.typeKin, 
                            typeMetabo = self.typeMtb, reactions = self.metabolisms, 
                            specComp = [self.eD.get(metabolism) for metabolism in self.metabolisms],
                            sample = 'All', DGsynth = self.DGsynth, solvent = 'H2O', 
                            molality = molality_, asm = asm_, phase= 'L-FW')
            self.DGr = {self.metabolisms[i] : ISAMERRA2inst.DGr[f'{self.metabolisms[i]}_pH:{pH_}'] * 1000 
                        for i in range(len(self.metabolisms))} # Gibbs free energy [J/moleD]
            self.Rs = {self.metabolisms[i] : ISAMERRA2inst.Rs[f'{self.metabolisms[i]}'] / 3600 
                        for i in range(len(self.metabolisms))} # reaction rates [moleD/(cell.s)]
            self.CSP = {self.metabolisms[i] : ISAMERRA2inst.CSP[f'{self.metabolisms[i]}_pH:{pH_}'] 
                        for i in range(len(self.metabolisms))} # CSP [fW/cell]
            self.envConditions = ISAMERRA2inst
        #import required CAMSMERRA2 attributes:   
        elif self.envModel == 'CAMSMERRA2':
            if len(self.coord) == 3:
                alt = self.coord[0]
                lon = self.coord[1]
                lat = self.coord[2]
            else: raise AttributeError(f'A list of 3 elements (selected altitude, longitude, latitude) should be given as coordinates for CAMSMERRA2, current coord: {self.coord}.')
            CAMSMERRA2inst = CAMSMERRA2(dataType = self.dataType,
                                           y = self.dataYear,
                                           m = self.dataMonth,
                                           d = self.dataDay,
                                           pH = pH_,
                                           bbox = (lon, lat, lon, lat),
                                           keys = 'All',
                                           phase = self.Wtype, 
                                           altArray = [alt],
                                           numAlt = 50,
                                           surftrop = None,
                                           keysAsAttributes = False, 
                                           showMessage = True)
            CAMSMERRA2inst.salinity = salinity_     # set to None by default in CAMSMERRA2
            CAMSMERRA2inst.methods = method_        # set to None by default in CAMSMERRA2
            CAMSMERRA2inst.getCSP(paramDB = self.kinDB, typeKin = self.typeKin, 
                            typeMetabo = self.typeMtb, reactions = self.metabolisms, 
                            specComp = [self.eD.get(metabolism) for metabolism in self.metabolisms],
                            sample = 'All', DGsynth = self.DGsynth, solvent = 'H2O', 
                            molality = molality_, asm = asm_, phase = 'L-FW')
            self.DGr = {self.metabolisms[i] : CAMSMERRA2inst.DGr[f'{self.metabolisms[i]}_pH:{pH_}'] * 1000 
                        for i in range(len(self.metabolisms))} # Gibbs free energy [J/moleD]
            self.Rs = {self.metabolisms[i] :CAMSMERRA2inst.Rs[f'{self.metabolisms[i]}'] / 3600  
                        for i in range(len(self.metabolisms))} # reaction rates [moleD/(cell.s)]
            self.CSP = {self.metabolisms[i] :CAMSMERRA2inst.CSP[f'{self.metabolisms[i]}_pH:{pH_}'] 
                        for i in range(len(self.metabolisms))} # CSP [fW/cell]
            self.envConditions = CAMSMERRA2inst
        if self.envModel == "GWB":
            raise AttributeError('Code part dedicated to general water body was not written yet.')

    def _ODEsystem_MSMM(self, t, y):
        """
        Function for the model's system of differential equations.
        
        Parameters
        ----------
        
        y : LIST of INT or FLOAT
            Initial biomass in each metabolic state, e.g. [cell/m^3 air]
        t : LIST or np.array
            Time range over which biomass variation is computed
            
        Returns
        -------
        dB : LIST of FLOAT
            Biomass variation [cell/h] in each metabolic state (Growth, Maintenance, Death).
        
        """
        #import self.attributes
        mortality = self.mortality.copy()
        mG = mortality[0]   #mortality in growth state
        mM = mortality[1]   #mortality in maintenance state
        K = self.K
        DGr = self.DGr.copy()
        Rs = self.Rs.copy()
        dB = []
        
        for i in range(len(self.metabolisms)):
            Bg = y[0+3*i]
            Bm = y[1+3*i]
            Blist = [Bg, Bm]

            Yx = -(DGr[self.metabolisms[i]] * (0.5 / 1.04e-10))      # cell growth yield [cell/mol eD]
            Yx = CSP.estimate_yield(DGr[self.metabolisms[i]])      # cell growth yield [cell/mol eD]
            Btot = Bg + Bm
            # Compute biomass transfer between metabolic states
            Rm_g, Rg_m, Rm_rip = MSMM._Bflux(self, Blist, self.metabolisms[i])
            # Compute biomass variation       
            dBg = Yx * Rs[self.metabolisms[i]] * Bg * (1 - (Btot / K)) - mG * Bg - Rg_m + Rm_g 
            dBm =  - mM * Bm - Rm_g - Rm_rip + Rg_m               
            dBrip = mG * Bg + mM * Bm + Rm_rip  
            #print(f'cell growth yield for {self.metabolisms[i]}={Yx}')
            #print(f'uptake rate for {self.metabolisms[i]}={Rs[self.metabolisms[i]]}')
            dB = np.append(dB,[np.squeeze(dBg), 
                       np.squeeze(dBm), 
                       np.squeeze(dBrip)])
        return dB

    def _Bflux(self, Blist, metabolism):
        """
        Function to compute biomass transfer between metabolic states.
        
        Parameters
        ----------

        Blist : LIST
            List of 2 floats corresponding to biomass (e.g. [cell/m^3 air])
            in each state (growth and maintenance) at time t.
            First element of the list must be for growth,
            second for maintenance.
        
        Returns
        -------
        Rlist : LIST
             List of computed biomass transfer [cell/h] for each kind of metabolic shift:
                 - Rm_g : transfer from maintenance to growth
                 - Rg_m : transfer from growth to maintenance
                 - Rm_rip : transfer from maintenance to dead cells
        """
        if len(Blist) != 2:
            raise ValueError(f'Blist must contain 2 elements (Bg, Bm), current Blist: {Blist}.')
        #Create shift control dict in MSMM attributes
        self._stShifts()
        #Biomass in each metabolic state
        Bg = Blist[0] 
        Bm = Blist[1]
        #import metabolic shift controls and rates
        theta = self.MSctrls.copy()
        eta = self.specMSrate
        #compute biomass transfers
        Rm_g = Bm * eta * theta[metabolism]['GxM']
        Rg_m = Bg * eta * (1 - theta[metabolism]['GxM'])
        Rm_rip = Bm * eta * (1 - theta[metabolism]['M-RIP'])
        Rlist = [Rm_g, Rg_m, Rm_rip]
        return Rlist
    
    def _stShifts(self):
        """
        Function to compute shift controls between metabolic states.
                   
        """
        CSPdict = self.CSP.copy()
        st = self.st
        #Initialize thetaDict before shift controls (theta) calculations
        thetaDict = {}
        for metabolism in self.metabolisms:
            #Compute cell specific powers
            Pcat = CSPdict[metabolism]['Pcat']
            Ps = CSPdict[metabolism]['Ps']
            Pcell = CSPdict[metabolism]['Pcell']
            thetaDict[metabolism] = {}
            thetaDict[metabolism]['GxM'] = 1 / (np.exp((-Pcat + Pcell)/(st[0] * Pcell)) +1)
            thetaDict[metabolism]['M-RIP'] = 1 / (np.exp((-Pcat + Ps)/(st[1] * Ps)) +1)
        self.MSctrls = thetaDict

    def solveODE(self, Bini, tSpan, dt = 1, solExport = False):
        """
        Function to solve the MSMM ODE system and export the results as .xlsx document.
        
        Parameters
        ----------
        
        Bini : LIST of INT or FLOAT
            Initial biomass in each state (Growth, Maintenance, Death)
        tSpan : LIST or np.array
            Time range over which the microbial dynamic is computed, in hours
        dt : INT or FLOAT, optional (default : 1h)
            Time step for the integration.
        solExport : BOOL, optional (default : False)
            Command to export integrated biomass values as Excel document if set to True.
        
        Returns
        -------
        
        None 
        ODE solutions (numpy.ndarray of shape [3*len(metabolisms), tSpan+1]) are saved as MSMM attribute ('Bsol')
        If solExport is set to True, creates an Excel document of the results.
        
        """
        # check Bini
        if not isinstance(Bini, np.ndarray): Bini = np.array(Bini)
        if len(Bini) != 3*len(self.metabolisms):
            raise ValueError(f'Bini must contain {3*len(self.metabolisms)} elements, current length: {len(Bini)}.')
        # create time array for later plotting
        self.t_plot = np.linspace(0, tSpan, int(tSpan/dt)+1)
        #Initialize Bint matrix
        Bint = np.empty(3*len(self.metabolisms))
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
        Bstates = ['Growth', 'Maintenance', 'Death']
        nameSheet_B = 'MSMM biomass'
        metabolismnames = ', '.join([self.metabolisms[i] for i in range(len(self.metabolisms))])
        # import solutions of the ODE and time array from MSMM attributes
        time = pd.DataFrame(self.t_plot, columns = ['time (h)| states :'])

        # adapt header to environment model
        if self.envModel == 'ISA':
            alt = self.coord[0]
            introRowB = pd.DataFrame(np.array([f'States biomass [cell/m³ air] | Metabolisms: {metabolismnames} | Altitude: {alt}m | Environment: {self.envModel}']))
        elif self.envModel in ['ISAMERRA2', 'CAMSMERRA2']:
            alt = self.coord[0]
            lon = self.coord[1]
            lat = self.coord[2]
            introRowB = pd.DataFrame(np.array([f'States biomass [cell/m³ air] | Metabolisms: {metabolismnames} | Coordinates: {lon}LON;{lat}LAT | Altitude: {alt}m | Environment: {self.envModel}']))
        # write excel document 
        if not os.path.isfile(self.fullPathSave):
            with pd.ExcelWriter(self.fullPathSave) as writer:
                introRowB.to_excel(writer, sheet_name = nameSheet_B, index = False, header = False)
        with pd.ExcelWriter(self.fullPathSave, engine='openpyxl', mode = 'a', if_sheet_exists='overlay') as writer:
            time.to_excel(writer, sheet_name = nameSheet_B, startrow = 2, startcol = 1, index = False, header = True)
            #        Bdf = pd.DataFrame()
            for i in range(len(self.metabolisms)):    
                 Bdf_i = pd.DataFrame({self.metabolisms[i]+' '+Bstates[k] : self.Bsol[3*i+k] for k in range(3)})
             #           Bdf = pd.concat([Bdf,Bdf_i])   
                 Bdf_i.to_excel(writer, sheet_name = nameSheet_B, startrow = 2, startcol = 2+3*i, index = False, header = True)    
        
    def plotMSMM(self):
        
        """
        Function to plot MSMM microbial dynamic of a single point in the environment space.
        
        """
        # import solutions of the ODE and time array from MSMM attributes
        Bplot = getattr(self, 'Bsol', None)
        if Bplot is None:
            raise AttributeError('MSMM attribute "Bsol" could not be found. Please first use MSMM.solveODE().')
        # plotting of metabolic state curves
        if len(self.metabolisms)!=1:
            fig, axs =plt.subplots(nrows=len(self.metabolisms),sharex=True,layout='tight')#,squeeze=False)
            fig.set_size_inches(8.3, 12)
        
            for i in range(len(self.metabolisms)):
        
                axs[i].plot(self.t_plot, Bplot[0+3*i,:],'g-', linewidth=2.0, label='Growth')    #growth state curve
                axs[i].plot(self.t_plot, Bplot[1+3*i,:],'k-', linewidth=2.0, label='Maintenance')    #maintenance state curve
                axs[i].plot(self.t_plot, Bplot[2+3*i,:],'r--', linewidth=2.0 ,label='Dead cells')  #death state curve
                axs[i].set_xlim(0)
                #axs[i].set_ylim(0)
                theta1 = np.round(np.squeeze(self.MSctrls[self.metabolisms[i]]['GxM']),8)
                theta2 = np.round(np.squeeze(self.MSctrls[self.metabolisms[i]]['M-RIP']),8)
                axs[i].set_title(f'{self.communityNames[self.metabolisms[i]]}, θ(GxM)={theta1}, θ(M-RIP)={theta2}')
                axs[i].grid()
            axs[max(len(self.metabolisms)-2,0)].set_ylabel(self.plotYlabel)
            axs[len(self.metabolisms)-1].legend(bbox_to_anchor = (1.35, 0.8), title = 'Metabolic states:', title_fontproperties = {'size': 'large', 'weight': 'bold'})
            axs[len(self.metabolisms)-1].set_xlabel('time (hours)')
        else: 
            fig, ax =plt.subplots(layout='tight')#,squeeze=False)
            fig.set_size_inches(8.3, 4)
            ax.plot(self.t_plot, Bplot[0,:],'g-', linewidth=2.0, label='Growth')    #growth state curve
            ax.plot(self.t_plot, Bplot[1,:],'k-', linewidth=2.0, label='Maintenance')    #maintenance state curve
            ax.plot(self.t_plot, Bplot[2,:],'r--', linewidth=2.0 ,label='Dead cells')  #death state curve
            ax.set_xlim(0)
            #ax.set_ylim(0)
            theta1 = np.round(np.squeeze(self.MSctrls[self.metabolisms[0]]['GxM']),8)
            theta2 = np.round(np.squeeze(self.MSctrls[self.metabolisms[0]]['M-RIP']),8)            
            ax.set_title(f'{self.communityNames[self.metabolisms[i]]}, θ(GxM)={theta1}, θ(M-RIP)={theta2}')
            ax.grid()
            ax.set_ylabel(self.plotYlabel)
            ax.legend(bbox_to_anchor = (1.35, 0.8), title = 'Metabolic states:', title_fontproperties = {'size': 'large', 'weight': 'bold'})
            ax.set_xlabel('time (hours)')
        plt.suptitle(self.plotTitle)
        plt.show()
        return