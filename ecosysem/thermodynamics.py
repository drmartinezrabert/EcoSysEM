# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 18:09:14 2024

@author: 2424069M
"""
from reactions import Reactions as Rxn

import pandas as pd
import numpy as np
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
    
    def getDeltaGr(deltaGf, mRxn):
        """
        Function to get gibbs free energy of reaction from gibbs formation.

        Parameters
        ----------
        deltaGf : LIST or np.array
            Gibbs free energy of formation.
        mRxn : np.array
            Reaction matrix. (compounds)x(reactions)

        Returns
        -------
        deltaGr : np.array
            Gibbs free energy of reaction.
    
        """
        deltaGr = deltaGf @ mRxn
        
        return deltaGr
    
    def getDeltaHr(deltaHf, mRxn):
        """
        Function to get enthalpy of reaction from enthalpy of formation.

        Parameters
        ----------
        deltaHf : LIST or np.array
            Enthalpy of formation.
        mRxn : np.array
            Reaction matrix. (compounds)x(reactions)

        Returns
        -------
        deltaHr : np.array
            Enthalpy of reaction.
    
        """
        deltaHr = deltaHf @ mRxn
        
        return deltaHr
    
    def getKeq(compounds, mRxn, temperature):
        """
        Function to get equilibrium constants from DeltaG.

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
        deltaGf, notNaN = ThP.getThP('deltaG', compounds, 'L')
        deltaGr = ThP.getDeltaGr(deltaGf, mRxn)
        Keq = np.exp(deltaGr / (-R * temperature))
        
        return Keq

class ThEq:
    """
    Class for calulation of chemical/ion and interphase (G-L) equilibriums.
    
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
                t = np.c_[temperature]  # Absolute temperature [K]
                Ts = 298.15             # Standard temperature [K]
                Hs = Hs * np.exp(B.astype(float) * ((1 / t) - (1/Ts)))
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
        if not isinstance(temperature, list): temperature = [temperature]
        if not isinstance(pH, list): pH = [pH]
        if not isinstance(Ct, list): Ct = [Ct]
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
            if not rComp:
                print('!EcoSysEM.Warning: Check `pHSpeciation` file(s) (`reactions/` folder).')
                sys.exit()
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
            
class ThSA:
    """
    Class for thermodynamic state analysis of environment.
    
    """
    
#- DEBUGGING -#
# Hs, notNaN = ThEq.solubilityHenry(['O2', 'Ne', 'N2', 'Kr'], 'SW', [293.15, 298.15, 303.15]) # T (K)
# print(Hs)
# print(notNaN)
# print(Hs)
# ThEq.pHSpeciation()
# ThP.getThP('Hs', ['O2', 'Ne', 'N2', 'Kr'], 'SW')

# (species)x(pH)x(total concentration)x(temperature)x(compounds)
t = ThEq.pHSpeciation(['HCO3-', 'NH3', 'HNO2', 'HNO3'],     # Compounds: 4
                      [6.0, 6.5, 7.0, 7.5, 8.0, 8.5],       # pH: 6
                      [293.15, 298.15],                     # T: 2
                      [0.0, 0.25, 0.5, 0.75, 1.0],          # Ct: 5
                      True)                                 # Species: 4 (if True)
# print(t.shape)
# print(t)
# print('')
t_2 = t[:, 2, :, -1, :]
print(t_2.shape)
print(t_2)
# a = ThEq.pHSpeciation(['HCO3-', 'NH4+'],  # Compounds: 2
#                       7.0,                # pH: 1
#                       298.15,             # T: 1
#                       1.0,                # Ct: 1
#                       False)              # Species: 4 (if True)
# print(a.shape)
# print(a)

# , 'NH5'
#-------------#

#- Info of functions and examples -#
### Get Henry's law solubility constants (Hs)
#> ThEq.solubilityHenry(compounds, wType, temperature), where wType: 'FW', 'SW'
# Hs, notNaN = ThEq.solubilityHenry(['O2', 'Ne', 'N2', 'Kr'], 'SW', [293.15, 298.15, 303.15])
# print(Hs)
#----------------------------------#
