# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 18:09:14 2024

@author: 2424069M
"""

import pandas as pd
import numpy as np
# import os
# import matplotlib.pyplot as plt
# import sys

class ThP:
    """
    Class for thermodynamic parameters.
    
    """
    # Database directory
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
        Function to get thermodynamic parameters. 

        Parameters
        ----------
        typeParam : STR
            What parameters are requested. 
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
        
    def setThP():
        pass

class Equilibrium:
    """
    Class for chemical/ion and interphase (G-L) equilibrium.
    
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
    
    def pHSpeciation():
        """
        Lorem ipsum...
    
        Parameters
        ----------
        PARAM : TYPE
            DESCRIPTION.
        PARAM : TYPE
            DESCRIPTION.
    
        Returns
        -------
        PARAM : TYPE
            DESCRIPTION.
    
        """
        pass
    
#- DEBUGGING -#
# Hs, notNaN = Equilibrium.solubilityHenry(['O2', 'Ne', 'N2', 'Kr'], 'SW', [293.15, 298.15, 303.15]) # T (K)
# print(Hs)
# print(notNaN)
# print(Hs)
# Equilibrium.pHSpeciation()
# ThP.getThP('HenrySolub', ['O2', 'Ne', 'N2', 'Kr'], 'SW')
#-------------#

#- Info of functions and examples -#
### Get Henry's law solubility constants (Hs)
#> Equilibrium.solubilityHenry(compounds, wType, temperature), where wType: 'FW', 'SW'
# Hs, notNaN = Equilibrium.solubilityHenry(['O2', 'Ne', 'N2', 'Kr'], 'SW', [293.15, 298.15, 303.15])
# print(Hs)
#----------------------------------#
