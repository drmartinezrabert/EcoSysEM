# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 18:09:14 2024

@author: 2424069M
"""
from reactions import Reactions as Rxn

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
# import os
# import matplotlib.pyplot as plt

class ThP:
    """
    Class for thermodynamic parameters.
    
    """
    # Directory of databases
    path = 'db\\'
    
    def checkThP(typeParam, db, compounds, phase):
        """
        Function to check not available parameters.

        Parameters
        ----------
        typeParam : STR
            What parameter will be checked. {Only for warning message}.
        db : pandas
            Database of thermodynamic parameters.
        compounds : np.array, LIST or STR
            Name of compounds.
        phase : STR
            Fluid phase. Depending on typeParam.

        Returns
        -------
        notNaN : np.array or LIST
            Indices of compounds for which requested parameters are available.

        """
        if type(compounds) != 'numpy.ndarray':
            compounds = np.array(compounds)
        db_check = list(db.loc[db['Phase'] == phase, 'Formula'])
        findNaN = np.isin(compounds, list(set(compounds) - set(db_check)))
        compNaN = compounds[findNaN]
        notNaN = np.logical_not(findNaN)
        compnotNaN = compounds[notNaN]
        if compNaN.size != 0:
            print(f'!EcoSysEM.Warning: {typeParam} for {phase} not found for: {compNaN}.')
            print(f'>> Returned compounds: {compnotNaN}.\n')
        return notNaN
    
    def getThP(typeParam, compounds, phase):
        """
        Function to get thermodynamic parameters from csv file. 

        Parameters
        ----------
        typeParam : STR
            What parameters are requested, matching with csv name. E.g.:
                - 'Hs': Henry's solubility constant.
                - 'B': Temperature dependecy of Henry's solubility constant.
                - '...': ...
        compounds : np.array, LIST or STR
            Name of requested compounds.
        phase : STR
            Fluid phase. Depending on parameter type.
                - 'Hs' and 'B': FW - Fresh Water; SW - Sea Water.
                - '...': ...

        Returns
        -------
        Param : FLOAT or np.array
            Requested parameters. Shape: (temperature)x(compounds).
        notNaN : np.array
            Indices of parameters that are available.

        """
        dParam = pd.read_csv(ThP.path + typeParam + '.csv')
        dParam = dParam.set_index('Formula').loc[compounds].reset_index()
        Param = np.array(dParam.loc[dParam['Phase'] == phase, 'Value'])
        notNaN = ThP.checkThP(typeParam, dParam, compounds, phase)
        return Param, notNaN
    
    def getDeltaG0r(deltaG0f, mRxn):
        """
        Function to get gibbs free energy of reaction from DeltaG0f.

        Parameters
        ----------
        deltaG0f : LIST or np.array
            Standard Gibbs free energy of formation.
        mRxn : np.array
            Reaction matrix. (compounds)x(reactions)

        Returns
        -------
        deltaG0r : np.array
            Gibbs free energy of reaction.
    
        """
        deltaG0r = deltaG0f @ mRxn 
        return deltaG0r
    
    def getDeltaH0r(deltaH0f, mRxn):
        """
        Function to get enthalpy of reaction from DeltaH0f.

        Parameters
        ----------
        deltaH0f : LIST or np.array
            Standard enthalpy of formation.
        mRxn : np.array
            Reaction matrix. (compounds)x(reactions)

        Returns
        -------
        deltaHr : np.array
            Enthalpy of reaction.
    
        """
        deltaH0r = deltaH0f @ mRxn
        return deltaH0r
    
    def getKeq(compounds, mRxn, temperature):
        """
        Function to get equilibrium constants from DeltaG0f.

        Parameters
        ----------
        compounds : STR or np.array
            Requested compound(s).
        mRxn : np.array
            Reaction matrix. (compounds)x(reactions)
        temperature: FLOAT
            Absolute temperature [K]

        Returns
        -------
        Keq : np.array
            Equilibrium constants.
    
        """
        # Constants
        R = 0.0083144598   # Universal gas constant [kJ/mol/K]
        deltaG0f, notNaN = ThP.getThP('deltaG0f', compounds, 'L')
        deltaG0r = ThP.getDeltaG0r(deltaG0f, mRxn)
        Keq = np.exp(deltaG0r / (-R * temperature))
        return Keq

class ThEq:
    """
    Class for calulation of chemical, ion and interphase (G-L) equilibriums.
    
    """
    
    def solubilityHenry(compounds, wType = 'FW', temperature = []):
        """
        Get Henry's law solubility constants (Hs) in function of temperature.

        Parameters
        ----------
        compounds : STR or LIST
            Requested compound(s).
        wType : STR ('FW' or 'SW')
            Water type (Phase). FW - Fresh Water; SW - Sea Water.
        temperature : FLOAT or LIST
            Set temperature for Henry's law solubility constant(s).

        Returns
        -------
        Hs : FLOAT or np.array
            Henry's law solubility constant(s). Shape: (temperature)x(compounds).
            If no temperature is given, Hs is (1)x(compounds).
        notNaN : np.array
            Indices of parameters that are available.

        """
        if wType == 'FW' or wType == 'SW':
            # Henry's solubilities (Hs and B)
            Hs, notNaN_Hs = ThP.getThP('Hs', compounds, wType)
            B, notNaN_B = ThP.getThP('B', compounds, wType)
            if notNaN_Hs.sum() > notNaN_B.sum():
                notNaN = notNaN_B
            else:
                notNaN = notNaN_Hs
            if any(temperature):
                # Temperatures
                t = np.c_[temperature]                                  # Absolute temperature [K]
                Ts = 298.15                                             # Standard temperature [K]
                Hs = Hs * np.exp(B.astype(float) * ((1 / t) - (1/Ts)))  # [mol/m3]
            return Hs, notNaN
        else:
            print('!EcosysEM.Error: No water type selected. Use one of the following wType:\n'+
              '                 \'FW\'      - Fresh Water.\n'+
              '                 \'SW\'      - Sew Water.')
            return None, None
    
    def pHSpeciation(compounds, pH, temperature, Ct, rAllConc = False):
        """
        Calculate pH (or ion) speciation of selected compounds.

        Parameters
        ----------
        compounds : STR or LIST
            Requested compound(s).
        pH : FlOAT, LIST or np.array
                Set of pH.
        temperature : FLOAT, LIST or np.array
            Set of temperature for pH speciation [K].
        Ct : LIST or np.array
            Total concentrations of compounds.
        rAllConc : BOOL
            Option to select if return all compound species or only preferred.

        Returns
        -------
        rSpec : np.array
            Concentrations of selected compound species.
            If rAllConc = False: (pH)x(total concentration)x(temperature)x(compounds).
            If rAllConc = True: (species)x(pH)x(total concentration)x(temperature)x(compounds).
            Where species: [B], [B-], [B-2], [B-3].
    
        """
        if isinstance(compounds, str): compounds = [compounds]
        if isinstance(temperature, float or list): temperature = [temperature]
        if isinstance(pH, float or int): pH = [pH]
        if isinstance(Ct, float or int): Ct = [Ct]
        nCompounds = len(compounds)
        nTemperature = len(temperature)
        npH = len(pH)
        nCt = len(Ct)
        H = np.power(10, -abs(np.array(pH)))
        if rAllConc:
            rSpec = np.empty([4, npH, nCt, nTemperature, nCompounds])
        else:
            rSpec = np.empty([npH, nCt, nTemperature, nCompounds])
        for idCompound, iCompound in enumerate(compounds):
            # Reactions.getRxn -> only returns 1 reaction
            rComp, mRxn, infoRxn = Rxn.getRxn('pHSpeciation', iCompound)
            c_rComp = iCompound in rComp
            if not c_rComp:
                # print('!EcoSysEM.Warning: Check `pHSpeciation` file(s) (`reactions/` folder).')
                return Ct
            reqSp = rComp.index(iCompound) - 1 # Requested chemical species
            for idTemperature, iTemperature in enumerate(temperature):
                # Acid dissociation constants (Ka)
                transDict = {'First deP': 0, 'Second deP': 1, 'Third deP': 2}
                intRxn = [transDict[letter] for letter in infoRxn]
                Kai = ThP.getKeq(rComp, mRxn, iTemperature)
                Ka = np.zeros(3)
                Ka[intRxn] = Kai
                # Speciation
                theta = H**3 + (Ka[0] * H**2) + (Ka[0] * Ka[1] * H) + Ka[0] * Ka[1] * Ka[2]
                for idC, iC in enumerate(Ct):
                    mSpec_aux = np.array([iC * H**3 / theta,                    # [B]
                                          Ka[0] * iC * H**2 / theta,            # [B-]
                                          Ka[0] * Ka[1] * iC * H / theta,       # [B-2]
                                          Ka[0] * Ka[1] * Ka[2] * iC / theta])  # [B-3]
                    # Speciation matrix
                    if rAllConc:
                        rSpec[:,:,idC,idTemperature,idCompound] = mSpec_aux
                    else:
                        rSpec[:,idC,idTemperature,idCompound] = mSpec_aux[reqSp]
        rSpec = np.squeeze(rSpec)
        return rSpec

    def plotpHSpeciation(compounds, pH, temperature):
        """
        Plotting of pH (or ion) speciation of requested compound(s).

        Parameters
        ----------
        compounds : STR or LIST
            Requested compound(s).
        temperature : FLOAT
            Temperature for pH speciation [K].
        pH : LIST
            Set of pH.

        Returns
        -------
        None.

        """
        Ct = 1.0 # Total concentration [M]
        if not isinstance(temperature, float):
            print('!EcoSysEM.Warning: Temperature must be a FLOAT.')
            sys.exit()
        if isinstance(compounds, str): compounds = [compounds]
        for iCompound in compounds:
            Spec = ThEq.pHSpeciation(iCompound, pH, temperature, Ct, True)
            nFrac = (Spec.T / Ct) * 100 # Molar fraction [%]
            # Selection of involved chemical species
            c_nFrac = np.nonzero(np.sum(nFrac, axis = 0))
            nFrac = nFrac[:, c_nFrac[0]]
            # Get name of chemical species
            nCompounds = Rxn.getRxn('pHSpeciation', iCompound)[0][1:]
            # Plotting
            fig, ax = plt.subplots()
            ax.plot(pH, nFrac)
            ax.set_ylabel('Molar fraction (%)')
            ax.set_xlabel('pH')
            ax.set_xticks(np.arange(pH[0], pH[-1], 1))
            ax.set_yticks(np.arange(0, 110, 10))
            ax.margins(x = 0)
            plt.legend(nCompounds, loc = 'center left', bbox_to_anchor = (1, 0.5))
            plt.show()

class ThSA:
    """
    Class for thermodynamic state analysis of environment.
    
    """
    
    def getDeltaGr(compounds, T = 298.15, Ct = [], pH = [], asm = 'stoich'):
        """
        Calculate DeltaGr in function of pH, temperature and compound
        concentrations.

        Parameters
        ----------
        compounds : STR or LIST
            Requested compound(s).
        T : FLOAT
            Set of temperature [K].
        Ct : DICT
            Total concentrations of compounds {'compounds': [concentrations]}.
            All compounds of a reaction with the same number of concentrations.
        pH : LIST
            Set of pH.
        asm : STRING
            Assumption when products are not present in the environment.
            By default: 'stoich' - stoichiometric concentrations.

        Returns
        -------
        DGr : np.array
            Non-standard Gibbs free energy values.
            E.g., (species)x(pH)x(total concentration)x(temperature)x(compounds).
        infoRxn : LIST
            Name of reactions given by the user (see reaction/metabolisms.csv).

        """
        if isinstance(compounds, str): compounds = [compounds]
        if isinstance(T, float or int): T = [T]
        if isinstance(pH, float or int): pH = [pH]
        nT = len(T)
        npH = len(pH)
        # Check compound concentrations
        npCt = np.array(list(Ct.items()), dtype = object).T
        lenCt = [len(i) for i in npCt[1, :]]
        nCt = list(set(lenCt))
        cLen = len(nCt) == 1
        if not cLen:
            print('!EcoSysEM.Error: All compounds must have same number of concentrations.')
            sys.exit()
        nCt = nCt[0]
        # Initialise variables
        DGr = np.empty(0)
        rInfoRxn = np.empty(0)
        # Reactions.getRxn -> only 1 compound in arguments
        for idCompound, iCompound in enumerate(compounds):
            rComp, mRxn, infoRxn = Rxn.getRxn('metabolisms', iCompound)
            c_rComp = iCompound in rComp
            if not c_rComp:
                print(f'!EcoSysEM.Warning: {iCompound} not found. Check `metabolisms` file(s) (`reactions/` folder).')
                sys.exit()
            # Separate reactions that share the compound
            for iC in range(mRxn.shape[1]):
                mRxn_aux = np.c_[mRxn[:, iC]]
                cNonZ = np.nonzero(mRxn_aux)[:][0]
                rComp_aux = np.array(rComp)[cNonZ]
                mRxn_aux = mRxn_aux[cNonZ]
                infoRxn_aux = infoRxn[iC]
                # Save of reaction name
                rInfoRxn = np.concatenate((rInfoRxn, infoRxn_aux), axis = None)
                nRxn = len(rInfoRxn)
                #- DEBUGGING -#
                print(infoRxn_aux)
                print('-------------------------------------------')
                #-------------#
                # Reactions with iCompound as substrate
                cSubs = np.squeeze(mRxn_aux[np.where(rComp_aux == iCompound)])
                if cSubs < 0:
                    vSelected = abs(mRxn_aux[cSubs.astype(int)]) # Stoichiometric parameter of selectec compound
                    # print(vSelected)
                    # Calculate DeltaG0r
                    deltaG0f = ThP.getThP('deltaG0f', rComp_aux, 'L')[0]
                    deltaG0r = ThP.getDeltaG0r(deltaG0f, mRxn_aux)              # kJ
                    deltaG0r = deltaG0r / vSelected                             # kJ/mol
                    # Calculate DeltaH0r
                    deltaH0f = ThP.getThP('deltaH0f', rComp_aux, 'L')[0]
                    deltaH0r = ThP.getDeltaH0r(deltaH0f, mRxn_aux)          # kJ
                    deltaH0r = deltaH0r / vSelected  # kJ/mol
                    
                    rDGr = np.empty([nT, npH, nCt])
                    for idT, iT in enumerate(T):
                        for idpH, ipH in enumerate(pH):
                            #- DEBUGGING -#
                            print(iT)
                            print(ipH)
                            #-------------#
                            # Temperature influence (Gibbs–Helmholtz relationship)
                            Ts = 298.15                                             # Standard temperature [K]
                            deltaGr = deltaG0r * (iT / Ts) + deltaH0r * ((Ts - iT) / Ts)
                            #- DEBUGGING -#
                            print('-------------------------------------------')
                            print(f'DeltaG0r: {deltaG0r}')
                            print(f'DeltaGr (T={iT}): {deltaGr}')
                            #-------------#
                            if Ct:
                                # Calculate reaction quotient (Qr)
                                Qr = 0
                                R = 0.0083144598   # Universal gas constant [kJ/mol/K]
                                uComp = np.array(list(Ct.keys())) # Compounds given by user (their concentrations)
                                findVComp = np.argwhere(uComp == iCompound).squeeze()
                                if findVComp.size == 0:
                                    print(f'!EcoSysEM.Error: Compound concentration(s) for {iCompound} not found.')
                                    sys.exit()
                                for idComp, iComp in enumerate(rComp_aux):
                                    findComp = np.argwhere(uComp == iComp)
                                    vi = mRxn_aux[idComp] / vSelected # Stoichiometric parameter (per unit of selected compound)
                                    if findComp.size == 0:
                                        if vi < 0:
                                            print(f'!EcoSysEM.Error: Substrate concentration(s) for {iComp} not found.')
                                            sys.exit()
                                        else:
                                            if iComp == 'H+':
                                                if not pH:
                                                    # [H+] will be calculated assuming neutral pH (pH:7).
                                                    iConc = 10**(-7) * np.ones((1, nCt))
                                                else:
                                                    iConc = 10**(-ipH) * np.ones((1, nCt))
                                            elif iComp == 'H2O':
                                                # It is assumed a water activity of 1.0, as a pure liquid (it can be less).
                                                iConc = 1.0 * np.ones((1, nCt))
                                            else:
                                                # [P] is calculated based on stoichiometry.
                                                if asm == 'stoich':
                                                    iConc = (vi) * Ct[uComp[findVComp]]
                                    else:
                                        iConc = Ct[iComp]
                                    # pH speciation
                                    if pH and iComp != 'H+' and iComp != 'H2O' and iComp != 'OH-':
                                        iConc = ThEq.pHSpeciation(iComp, ipH, iT, iConc)
                                    # Calculation of reaction quotient (Qr)
                                    Qr += vi * np.log(iConc)
                                # Concentration influence (Nernst relationship)
                                deltaGr = deltaGr + R * iT * Qr
                    
                            #- DEBBUGING -#
                            print(f'DeltaGr (Nernst): {deltaGr}')
                            print('===========================================')
                            print('')
                            #-------------#
                            print(deltaGr)
                            rDGr[idT, idpH, :] = deltaGr
                    #- DEBUGGING -#
                    print(rDGr)
                    #-------------#
                    if DGr.size == 0:
                        DGr = rDGr
                    else:
                        DGr = np.concatenate((DGr, rDGr), axis = 0)
        DGr = DGr.reshape((nRxn, nT, npH, nCt))
        DGr = np.squeeze(DGr)
        #- DEBUGGING -#
        print('--------------')
        
        print('')
        print('List DGr:')
        print(DGr)
        print(f'List Info: {rInfoRxn}')
        #-------------#
        return DGr, rInfoRxn
    
#- DEBUGGING -#
## Henry's solubility
# Hs, notNaN = ThEq.solubilityHenry(['O2', 'Ne', 'N2', 'Kr'], 'SW', [293.15, 298.15, 303.15]) # T (K)
# print(Hs)
# print(notNaN)
# print(Hs)
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
# a = ThEq.pHSpeciation('HCO3-',  # Compounds: 2
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
## Get DeltaGr (Thermodynamic state analysis)
conc = {'CO': [1.08e-10, 1],
        'O2': [3.40e-4, 1],
        'CO2': [2.78e-5, 1],
        'NH3': [1.33e-6, 1]}
temperature = 255.65
pH = 7.0
# ---------------------------------------
# conc = {'CO': [1.08e-10, 1],
#         'O2': [3.40e-4, 1],
#         'CO2': [2.78e-5, 1],
#         'NH3': [1.33e-6, 1]}
# temperature = [216.65, 255.65, 288.15]
# pH = [2.0, 5.0, 8.0]
DGr, rInfoRxn = ThSA.getDeltaGr(['CO', 'NH3'], temperature, conc, pH)
# #-------------#

#- Info of functions and examples -#
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
#----------------------------------#

#- PREV. CODE -#
# if Ct:
#     # Calculate rection quotient (Qr)
#     npCt = np.array(list(Ct.items()), dtype = object).T
#     lenCt = [len(i) for i in npCt[1, :]]
#     cLen = len(set(lenCt)) == 1
#     if not cLen:
#         print('!EcoSysEM.Error: All compounds must have same number of concentrations.')
#         sys.exit()
#     else:
#         uComp = np.array(npCt[0, :], ndmin=2)               # Compounds given by user
#         print(uComp)
#         iCt = np.array([i for i in npCt[1, :]])
#         print(iCt)
#         print('-----------')
#         #- TEST -#
#         stoich = np.empty(rComp_aux.shape)
#         for idComp, iComp in enumerate(rComp_aux):
#             print(iComp)
#             findComp = np.argwhere(np.squeeze(uComp) == iComp)
#             if findComp.size == 0:
#                 if mRxn_aux[idComp] < 0:
#                     print(f'!EcoSysEM.Error: Substrate concentrations for {iComp} not found.')
#                     sys.exit()
#                 elif mRxn_aux[idComp] > 0 and iComp == 'H+':
#                     if pH:
#                         print('pH!')
#                     else:
#                         pass
#                         # If no pH, neutral pH is assumed
#                         cH = [10**(-7) * np.ones(lenCt[0])]
#                         iCt = np.concatenate((iCt, cH), axis = 0)
#                 elif mRxn_aux[idComp] > 0 and iComp != 'H+':
#                     if asm == 'stoich':
#                         print('Product no H+!')
#                 uComp = np.concatenate((uComp, iComp) , axis = None)
#                 print(iCt)
#                 print(uComp)
#             iF = (mRxn_aux[np.argwhere(np.squeeze(uComp) == iComp)] / vSelected ).item() # .item() to avoid Decrepated NumPy 1.25 
#             stoich[idComp] = iF
#             print(stoich)
#         #--------#
#         print(iCt)
#         npCt = np.concatenate((uComp, stoich, iCt.T), axis=0)
#         print(npCt)
#     #- DEBUGGING -#
#     print('-------------')
#     print('')
#     #-------------#
#--------------#