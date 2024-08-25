# -*- coding: utf-8 -*-
"""
envdef - define and create environment(s)
=========================================

Summary
-------

Lorem ipsum...

Classes
-------

    :py:class:`Environment` - Common ancestor. Abstract class.
    
    :py:class:`ISA` - International Standard Atmosphere (ISA).

Functions
--------

    :py:fun:`lorem_ipsum` - Lorem ipsum.

------------------------------------------------------------------------------
"""
from thermodynamics import ThEq as eQ

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

import warnings
warnings.simplefilter(action = 'ignore')

class Environment:
    """
    Abstract base class for environment definition file.
    
    """
    def __init__(self, temperature, pressure, pH, compounds, composition):
        self.temperature = temperature
        self.pressure = pressure
        self.pH = pH
        self.compounds = compounds
        self.compositions = composition
        
    def setT(self, new_T):
        self.temperature = new_T
        
    def setP(self, new_P):
        self.pressure = new_P
    
    def setpH(self, new_pH):
        self.pH = new_pH
    
    def setComposition(self, compound, composition):
        """
        Modify composition of existing compounds or add new components with
        their respective compositions.

        Parameters
        ----------
        compound : STR or LIST
            Set of new and/or existing compounds.
        composition : FLOAT or LIST
            Composition(s) of new and/or existing compound(s).
            
        """
        if not isinstance(compound, list): compound = [compound]
        if not isinstance(composition, list): composition = [composition]
        pre_comp = self.compositions
        new_comp = dict(zip(compound, composition))
        pre_comp.update(new_comp)
        self.compositions = pre_comp
        ## Sorting doesn't work when setComposition() entry is a list(s).
        # comp_sorted = dict(sorted(pre_comp.items(), key=lambda x: x[1], reverse = True))
        # self.compositions = comp_sorted
        
class ISA(Environment):
    """
    International Standard Atmosphere (ISA). Static atmospheric model of how
    pressure, temperature, density and viscosity of the Earth's atmosphere
    change over a wide range of altitudes.
    --------------------------------------------------------------------------
    References: - ISO 2533:1975
                - National Oceanic and Atmospheric Administration (NOAA)
                - Goody & Yung (1989), doi: 10.1093/oso/9780195051346.001.0001
    
    """
    def __init__(self, layers, H2O, pH):
        # Layers of the Earth's atmosphere (ISA)
        ISAproperties = {
                        'Layer': [0, 1, 2, 3, 4, 5, 6, 7],
                        'Layer name': ['Troposphere', 'Tropopuse', 'Stratosphere', 'Stratosphere', 'Stratopause', 'Mesosphere', 'Mesosphere', 'Mesopause'],
                        'Start altitude': [0, 11000, 20000, 32000, 47000, 51000, 71000, 84852],         # [m (above MSL)]
                        'End altitude': [11000, 20000, 32000, 47000, 51000, 71000, 84852, 85852],       # [m (above MSL)]
                        'Lapse rate': [-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0],                      # [C/km)]
                        'Base temperature': [15, -56.5, -56.5, -44.5, -2.5, -2.5, -58.5, -86.204],      # [C]
                        'Base pressure': [101325, 22632, 5474.9, 868.02, 110.91, 66.939, 3.9564, 0],    # [Pa]
                        }
        dISA = pd.DataFrame(data = ISAproperties)
        # Dry composition
        DryComposition = {
                         'Compounds': ['N2', 'O2', 'Ar', 'CO2', 'Ne', 
                                       'He', 'CH4', 'Kr', 'H2', 'N2O', 
                                       'CO', 'Xe','O3', 'NO2', 'I2', 
                                       'NH3','NO', 'SO2', 'H2S'], # Formula
                         'Compositions': [7.8084e-1, 2.0946e-1, 9.34e-3, 4.2e-4, 1.8182e-5, 
                                          5.24e-6, 1.92e-6, 1.14e-6, 5.5e-7, 3.3e-7, 
                                          1e-7, 9e-8, 7e-8, 2e-8, 1e-8, 
                                          4e-9, 1e-9, 1e-9, 5e-11] # [%vol]
                          }
        dDC = pd.DataFrame(data = DryComposition)
        # super().__init__(temperature=None, pressure=None, pH=self, compounds=None, composition=None) # !!!
        self.layers = layers
        self.pH = pH
        self.computeTandP(layers, dISA)
        self.compounds = dDC['Compounds']
        self.compositions = pd.Series(dDC.Compositions.values, index = dDC.Compounds).to_dict()
        self.computeWaterContent(H2O, dDC)
    
    def computeTandP(self, layers, dISA):
        """
        Compute the change of temperaature and pressure of the Earth's
        atmosphere over the range of altitudes. Based on ISA (ISO 2533:1975).
        
        """
        # Constants
        R = 8.3144598   # Universal gas constant [J/mol/K]
        g0 = 9.80665    # Gravitational acceleration [m/s^2]
        M0 = 0.0289644  # Molar mass of Earth's air
        # Data from ISA to compute temperature
        if layers == 'All':
            layers = range(8)
        elif type(layers) == int:
            layers = range(layers, layers+1)
        lapse_rate = (dISA.loc[layers]['Lapse rate']).to_numpy()
        base_T = (dISA.loc[layers]['Base temperature']).to_numpy()
        base_P = (dISA.loc[layers]['Base pressure']).to_numpy()
        start_alt = (dISA.loc[layers]['Start altitude']).to_numpy()
        end_alt = (dISA.loc[layers]['End altitude']).to_numpy()
        # Temperature calculation (ISA)
        alt = t = p = []
        for i in range(lapse_rate.size):
            alt_aux = np.array(range(start_alt[i], end_alt[i], 1000))
            # Temperature
            t_aux = base_T[i] + lapse_rate[i] * ((alt_aux - start_alt[i]) / 1000)
            # Pressure
            if lapse_rate[i] != 0:
                p_aux = base_P[i] * (1 + ((lapse_rate[i]/1000) / (base_T[i]+273.15)) * (alt_aux - start_alt[i])) ** (-(g0 * M0) / (R * lapse_rate[i]/1000))
            else:
                p_aux = base_P[i] * np.exp((-(g0 * M0) / (R * (base_T[i]+273.15))) * (alt_aux - start_alt[i]))
            # Appending...
            alt = np.append(alt, alt_aux)
            t = np.append(t, t_aux)
            p = np.append(p, p_aux)
        # Last value of arrays (Altitude, Temperature and Pressure)
        alt_end = alt_aux[-1] + 1000
        alt = np.append(alt, alt_end)
        t = np.append(t, base_T[-1] + lapse_rate[-1] * ((alt_end - start_alt[-1]) / 1000)) # Last T value
        if lapse_rate[-1] != 0:
            p = np.append(p,  base_P[-1] * (1 + ((lapse_rate[-1]/1000) / (base_T[-1]+273.15)) * (alt_end - start_alt[-1])) ** (-(g0 * M0) / (R * lapse_rate[-1]/1000)))
        else:
            p = np.append(p, base_P[-1] * np.exp((-(g0 * M0) / (R * (base_T[-1]+273.15))) * (alt_end - start_alt[-1])))
        # Safe Altitude, Temperature and Pressure
        self.altitude = alt # [m]
        self.temperature = t # [C]
        self.pressure = p # [Pa]
        
    def computeWaterContent(self, H2O, dDC):
        """
        It is assumed that atmospheric compositon of minor components (<0.9%, 
        such as CO2) does not change significantly when including water vapor 
        (i.e., wet composition). Therefore, compositions of N2, O2 and Ar are
        recalculated (%vol).
        
        """
        self.H2O = H2O
        newComp = np.array(dDC['Compositions'])
        # self.composition = dryComposition
        if H2O > 0:
            newComp[0] -= 0.78100 * H2O # N2_wet
            newComp[1] -= 0.20925 * H2O # O2_wet
            newComp[2] -= 0.01100 * H2O # Ar_wet
            self.compositions = pd.Series(newComp, index = dDC.Compounds).to_dict()
    
    def getVerticalProfiles(self, phase, compound = None):
        """
        Computation of vertical profiles of compounds (parcial pressure, Pi;
        gas concentration, Ci_G; liquid concentration in fresh water, Ci_L-FW;
        and liquid concentration in sea water, Ci_L-SW).
        Gas concentrations (Ci_G) are calculated using Dalton's law and the 
        ideal gas law, and liquid concentration (Ci_LFW and Ci_LSW) with 
        Henry's law.
        
        Parameters
        ----------
        phase : STR ('G', 'L-FW', 'L-SW', 'L' or 'All')
            DESCRIPTION. Selection of phase of vertical profile.
                        'G' - Gas.
                        'L-FW' - Liquid fresh water.
                        'L-SW' - Liquid sea water.
                        'L' - Both liquid phases (L-FW, L-SW).
                        'All' - All phaes (G, L-FW, L-SW).
        compound : STR or LIST, optional
            DESCRIPTION. Interested compounds. Default: None -> All compounds.

        Returns
        -------
        Vertical profile/s of compounds as NumPy array (altitude)x(compounds).
        
        """
        # Data
        p = np.c_[self.pressure]        # [Pa] `np.c_[]` -> Column array
        t = self.temperature + 273.15   # [K]
        compounds = self.compounds
        compositions = np.array(list(self.compositions.values()))
        if compound:
            if type(compound) is str: compound = [compound]
            findC = compounds.reset_index().set_index('Compounds').loc[compound].reset_index().set_index('index').index
            compositions = compositions[findC]
            compounds = compounds[findC]
        # Constants
        R_g = 8314.46261815324  # Universal gas constant [(L·Pa)/(K·mol)]
        Hs_FW, notNaN_HsFW = eQ.solubilityHenry(compounds, 'FW', t)
        Hs_SW, notNaN_HsSW = eQ.solubilityHenry(compounds, 'SW', t)
        # Gas phase - Partial pressure (Pi)
        Pi = p * compositions
        # Gas concentration (Ci_G)
        Ci_G = (Pi.T / (R_g * t)).T
        # Re-selection of compounds (parcial pressures, Pi)
        Pi_FW = Pi[:, notNaN_HsFW]
        Pi_SW = Pi[:, notNaN_HsSW]
        # Liquid concentrations fresh water (Ci_LFW)
        Ci_LFW = Pi_FW * Hs_FW * (1/1000) # [mol/L]
        # Liquid concentrations sea water (Ci_LSW)
        Ci_LSW = Pi_SW * Hs_SW * (1/1000) # [mol/L]
        if phase == 'G':
            return Pi, Ci_G
        elif phase == 'L-FW':
            return Ci_LFW
        elif phase == 'L-SW':
            return Ci_LSW
        elif phase == 'L':
            return Ci_LFW, Ci_LSW
        elif phase == 'All':
            return Pi, Ci_G, Ci_LFW, Ci_LSW
        else:
            sys.exit('!EcosysEM.Error: No phase selected. Use one of the following string:\n'+
                     '                             \'G\'       - Gas.\n'+
                     '                             \'L-FW\'    - Liquid fresh water.\n'+
                     '                             \'L-SW\'    - Liquid sea water.\n'+
                     '                             \'L\'       - Both liquid phases (L-FW, L-SW).\n'+
                     '                             \'All\'     - All phases (G, L-FW, L-SW).')
    
    ## Plotting functions 
    def plotTandP(self):
        """
        Plotting of pressure (in atm) and temperature (in Kelvin) along the
        atmosphere (altitude in km).
        
        """
        # Variables
        alt = self.altitude / 1000      # [km]
        t = self.temperature + 273.15   # [K]
        p = self.pressure / 101325      # [atm]
        # Temperature
        fig, ax1 = plt.subplots(figsize = (3, 6))
        ax1.set_ylabel('Altitude (km)')
        ax1.set_xlabel('Temperature (K)', color = 'tab:red')
        ax1.set_xlim([0, 300])
        ax1.set_ylim([0, alt[-1]])
        ax1.plot(t, alt, color = 'tab:red')
        # Pressure
        ax2 = ax1.twiny()
        ax2.set_xlabel('Pressure (atm)', color = 'tab:blue')
        ax2.set_xlim([0, 1])
        ax2.plot(p, alt, color = 'tab:blue')
        # Plot properties and showing
        fig.tight_layout()
        plt.show()
    
    def plotCompsProfiles(self, C, xLabel, logCLabel = False, compounds = None):
        """
        Plotting of composition profiles along the atmosphere (altitude in km).
        
        """
        # Variables
        alt = self.altitude / 1000     # [km]
        # Plotting
        fig, ax = plt.subplots(figsize = (5,10))
        ax.set_ylabel('Altitude (km)')
        ax.set_xlabel(xLabel)
        ax.plot(C, alt)
        if logCLabel == True:
            ax.set_xscale('log')
        if not compounds:
            compoundsName = self.compounds
        else:
            compoundsName = compounds
        plt.legend(compoundsName, loc = 'center left', bbox_to_anchor = (1, 0.5))
        plt.show()
        
# Define Earth's atmosphere
envISA = ISA(0, 0.00 , 5.0) # ISA(Layer/s, H2O (%), pH)

#- DEBUGGING -#
# pG, cG, L_FWp, L_SWp = envISA.getVerticalProfiles('All', ['O2','N2'])
pG, cG, L_FWp, L_SWp = envISA.getVerticalProfiles('All')
# print('Parcial pressure (Pa)')
# print(pG)
# print('')
# print('Gas concentration (M)')
# print(cG)
# print('')
# print('Fresh water concentration (M)')
# print(L_FWp)
# print('')
# print('Sea water concentration (M)')
# print(L_SWp)
envISA.plotCompsProfiles(pG / 101325, 'Pressure (atm)', False)
envISA.plotCompsProfiles(cG, 'Concentration in gas (M)', True)
envISA.plotCompsProfiles(L_FWp * 10**(6), 'Concentration in FW (µM)', True)
envISA.plotCompsProfiles(L_SWp * 10**(6), 'Concentration in SW (µM)', True, 
                          ['N2','O2','Ar','CO2','CH4','H2','N2O','CO','NH3','NO','SO2','H2S'])

# print(envISA.altitude)
# print(envISA.temperature)
# print(envISA.pressure)
# print(envISA.compounds)
# print('----')
# print(envISA.compositions)
# print('')
# print(envISA.H2O)
# print(envISA.compositions)
# print('')
#-------------#

#- Info of functions and examples -#
### Modify composition of existing compounds or add new components (with respective compositions)
#> setComposition(compounds, compositions) /lists or floats/
# envISA.setComposition(['Test', 'N2'], [[0.5, 0.3], 0.49]) 
# print(envISA.compositions)
### Get Vertical Profiles of compounds in ISA
#> setComposition(phase, compounds = None), where `phase`: 'G', 'L-FW', 'L-SW', 'L', 'All'.
# envISA.getVerticalProfiles('L', ['N2', 'O2', 'Ar'], plot = None)  
### Plotting Temperature and Pressure
#> envISA.plotTandP()
### Plotting profile of components
#> plotCompsProfiles(verticalProfiles, yLabel {string}, logScale = True/False, compounds = None), logScale = True -> logarithmic scale in concentrations; compouds = list of strings
#! If all profiles are obtained, it is not necessary to define the name of compounds
# pG, cG,  L_FWp, L_SWp = envISA.getVerticalProfiles('All')
# envISA.plotCompsProfiles(pG / 101325, 'Pressure (atm)', False)
# envISA.plotCompsProfiles(cG, 'Concentration in gas (M)', True)
# envISA.plotCompsProfiles(L_FWp * 10**(6), 'Concentration in FW (µM)', True)
# envISA.plotCompsProfiles(L_SWp * 10**(6), 'Concentration in SW (µM)', True, 
#                          ['N2','O2','Ar','CO2','CH4','H2','N2O','CO','NH3','NO','SO2','H2S'])
#----------------------------------#

#- Old function versions -#
#-------------------#