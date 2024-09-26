# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 09:13:26 2024

@author: 2424069M
"""
from envdef import ISA
from thermodynamics import ThSA

import numpy as np

print('> Running test.py')
print('')

#- DEBUGGING -#
## Plot DeltaGr (Thermodynamic state analysis)
envISA = ISA(0, 0.00, 7.0, 500) # ISA(Layer/s, H2O (%), pH, resolution = 1000 {size of altitude nodes per layer, [m]})
altitude = envISA.altitude
temperature = envISA.temperature + 273.15
minpH = 3.0
maxpH = 8.0
pH = np.linspace(minpH, maxpH, len(altitude))
conc = envISA.getDictConc('L-FW', ['CO', 'O2', 'CO2', 'NH3'])
concx1 = {'CO': [1.08e-10],
          'O2': [3.40e-4],
          'CO2': [2.78e-5],
          'NH3': [1.33e-6]}
# ThSA.plotDeltaGr(temperature, pH, 'metabolisms', ['COOB', 'AOM', 'CMX'], ['CO', 'NH3', 'NH3'], concx1)
# ThSA.plotDeltaGr(temperature, pH, 'metabolisms', ['COOB', 'AOM', 'CMX'], ['CO', 'NH3', 'NH3'], conc, 'x')
# ThSA.plotDeltaGr(temperature, pH, 'metabolisms', ['COOB', 'AOM', 'CMX'], ['CO', 'NH3', 'NH3'], conc, 'y')
ThSA.plotDeltaGr(temperature, pH, 'metabolisms', ['COOB', 'AOM', 'CMX'], ['CO', 'NH3', 'NH3'], conc, 'xy')

## Get DeltaGr (Thermodynamic state analysis)
# envISA = ISA(0, 0.00, 7.0, 5000) # ISA(Layer/s, H2O (%), pH, resolution = 1000 {size of altitude nodes per layer, [m]})
# altitude = envISA.altitude
# print(f'·Altitude: {altitude}')
# print('')
# temperature = [216.65, 242.65, 255.65, 288.15]
# pH = [3.0, 5.0, 6.0, 8.0]
# nT = len(temperature)
# npH = len(pH)
# vT = list(range(nT))
# vpH = list(range(npH))
# conc = envISA.getDictConc('L-FW', ['CO', 'O2', 'CO2', 'NH3'])
# DGr, rInfoRxn = ThSA.getDeltaGr('metabolisms', ['COOB', 'AOM', 'CMX'], ['CO', 'NH3', 'NH3'], temperature, conc, pH)
# print(f'·DeltaGr ({rInfoRxn}):')
# print(DGr)
# print('---------------------------------------------------')
# print('')
# for idRxn, iRxn in enumerate(rInfoRxn):
#     DGr_T = DGr[idRxn, vT, :, vT]
#     DGr_pH = DGr[idRxn, :, vpH, vpH]
#     if iRxn == 'CMX':
#         print(f'·Concentration associated with temperature ({iRxn}):')
#         print(DGr_T)
#         print('')
#         print(f'·Concentration associated with pH, Transpose ({iRxn}):')
#         print(DGr_pH.T)
#         print('')
#-------------#

#- Info of functions and examples -#

# envdef.py
### Modify composition of existing compounds or add new components (with respective compositions)
#> setComposition(compounds, compositions) /lists or floats/
# envISA.setComposition(['Test', 'N2'], [[0.5, 0.3], 0.49]) 
# print(envISA.compositions)
### Get Vertical Profiles of compounds in ISA
#> getVerticalProfiles(phase, compounds = None), where `phase`: 'G', 'L-FW', 'L-SW', 'L', 'All'.
# envISA.getVerticalProfiles('L', ['N2', 'O2', 'Ar'])
#> getDictConc(phase, compounds = None), where `phase`: 'G', 'L-FW', 'L-SW', 'L', 'All'.
# envISA.getDictConc('L', ['N2', 'O2', 'Ar'])
### Plotting Temperature and Pressure
#> envISA.plotTandP()
### Plotting profile of components
#> plotCompsProfiles(verticalProfiles, yLabel {string}, logScale = True/False, compounds = None), logScale = True -> logarithmic scale in concentrations; compouds = list of strings
#! If all profiles are obtained, it is not necessary to define the name of compounds
# pG, cG,  L_FWp, L_SWp,nameCompounds_G, nameCompounds_FW, nameCompounds_SW = envISA.getVerticalProfiles('All')
# envISA.plotCompsProfiles(pG / 101325, 'Pressure (atm)', False)
# envISA.plotCompsProfiles(cG, 'Concentration in gas (M)', True)
# envISA.plotCompsProfiles(L_FWp * 10**(6), 'Concentration in FW (µM)', True)
# envISA.plotCompsProfiles(L_SWp * 10**(6), 'Concentration in SW (µM)', True, nameCompounds_SW)

# thermodynamics.py
### Get Henry's law solubility constants (Hs)
#> ThEq.solubilityHenry(compounds, wType, temperature), where wType: 'FW', 'SW'
# Hs, notNaN = ThEq.solubilityHenry(['O2', 'Ne', 'N2', 'Kr'], 'SW', [293.15, 298.15, 303.15])
# print(Hs)
### Calculate pH (or ion) speciation of selected compounds.
#> pHSpeciation(compounds, pH, temperature, Ct, rAllConc = False)
#> return.shape: (species)x(pH)x(total concentration)x(temperature)x(compounds)
# t = ThEq.pHSpeciation(['HCO3-', 'NH3', 'HNO2', 'HNO3'],     # Compounds: 4
#                       [6.0, 6.5, 7.0, 7.5, 8.0, 8.5],       # pH: 6
#                       [293.15, 298.15],                     # T: 2
#                       [0.0, 0.25, 0.5, 0.75, 1.0],          # Ct: 5
#                       True)                                 # Species: 4 (if True)
# print(t)
### Plot pH (or ion) speciation
#> plotpHSpeciation(compounds, pH, temperature), where temperature must be a FLOAT
# pH = np.arange(0, 14, 0.5)
# ThEq.plotpHSpeciation(['NH3', 'HNO2', 'HNO3', 'H2SO4', 'H2S'], pH, 298.15)
#> getDeltaGr(typeRxn, input_, specComp, temperature, concentrations, pH, assumptions for products, warnings), where input_ can be compounds or reactions
# conc = {'CO': [1.08e-10],
#         'O2': [3.40e-4],
#         'CO2': [2.78e-5],
#         'NH3': [1.33e-6]}
# temperature = 255.65
# pH = 7.0
# DGr, rInfoRxn = ThSA.getDeltaGr('metabolisms', ['CO', 'NH3'], True, temperature, conc, pH)
# DGr, rInfoRxn = ThSA.getDeltaGr('metabolisms', ['COOB', 'AOM', 'CMX'], ['CO', 'NH3', 'NH3'], temperature, conc, pH)

# reactions.py
### Function to get reactions by their name or involing component (calling getRxnByComp() and getRxnByName() functions)
#> Reactions.getRxn(typeRxn, input_, warnings = False), where input_ can be a string (one compound/reaction) or list (multiple compounds/reactions)
# rComp, mRxn, infoRxn = Reactions.getRxn('metabolisms', ['AOM', 'CMX'])
### Get reactions involving one or more compounds
#> Reactions.getRxnByComp(typeRxn, compounds, warnings = False), where compounds can be a string (one compound) or list (multiple compounds)
# rComp, mRxn, infoRxn = Reactions.getRxnByComp('metabolisms', ['NH3', 'CH4'])
### Get reactions from their name (parenthesis in Excel/csv files)
#> Reactions.getRxnByName(typeRxn, nameRxn, warnings = False), where nameRxn can be a string (one reaction) or list (multiple reactions)
# rComp, mRxn, infoRxn = Reactions.getRxnByName('metabolisms', ['AOM', 'CMX'])
#----------------------------------#

#- DEBUGGING -#

# envdef.py
# pG, cG, L_FWp, L_SWp = envISA.getVerticalProfiles('All', ['O2','N2'])
# pG, cG, L_FWp, L_SWp = envISA.getVerticalProfiles('All')
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
# envISA.plotCompsProfiles(pG / 101325, 'Pressure (atm)', False)
# envISA.plotCompsProfiles(cG, 'Concentration in gas (M)', True)
# envISA.plotCompsProfiles(L_FWp * 10**(6), 'Concentration in FW (µM)', True)
# envISA.plotCompsProfiles(L_SWp * 10**(6), 'Concentration in SW (µM)', True, 
#                           ['N2','O2','Ar','CO2','CH4','H2','N2O','CO','NH3','NO','SO2','H2S'])

# thermodynamic.py
## Henry's solubility
# Hs, notNaN = ThEq.solubilityHenry(['O2', 'Ne', 'N2', 'Kr'], 'SW', [293.15, 298.15, 303.15]) # T (K)
## Get thermodynamic paramether
# ThP.getThP('Hs', ['O2', 'Ne', 'N2', 'Kr'], 'SW')
## pH (ion) speciation
# t = ThEq.pHSpeciation(['HCO3-', 'NH3', 'HNO2', 'HNO3'],     # Compounds: 4
#                       [6.0, 6.5, 7.0, 7.5, 8.0, 8.5],       # pH: 6
#                       [293.15, 298.15],                     # T: 2
#                       [0.0, 0.25, 0.5, 0.75, 1.0],          # Ct: 5
#                       True)                                 # Species: 4 (if True)
# print(t.shape)
# print(t)
# print('')
# t_2 = t[:, 2, :, -1, :]
# print(t_2.shape)
# print(t_2)
# a = ThEq.pHSpeciation('HCO3-',  # Compounds: 1
#                       7.0,                # pH: 1
#                       298.15,             # T: 1
#                       1.0,                # Ct: 1
#                       False)              # Species: 4 (if True)
# print(a.shape)
# print(a)
## Plot pH Speciation
# pH = np.arange(0, 14, 0.25)
# ThEq.plotpHSpeciation(['NH4+', 'NO2-', 'HNO3', 'H2SO4', 
#                         'H2S', 'H2SO3', 'H2CO3'], 
#                         pH, 298.15)
#-------------#