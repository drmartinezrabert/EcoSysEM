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
import os.path

def _rhoWater(T, S):
    """
    Density of water.
    --------------------------------------------------------------------------
    References: 
        - Millero & Chen (1973). doi: 10.1016/0198-0149(80)90016-3

    Parameters
    ----------
    T : FLOAT, LIST or np.ndarray
        Temperature [K].
    S : FLOAT, LIST or np.ndarray
        Salinity [ppt or g/L].

    Returns
    -------
    rho : FLOAT, LIST or np.ndarray
        Density of water [kg/L].

    """
    if not isinstance(T, np.ndarray): T = np.array(T)
    if not isinstance(S, np.ndarray): S = np.array(S)
    if not np.shape(T) == np.shape(S):
        print('!EcoSysEM.Warning: T and S arguments must have the same shape.')
        sys.exit()
    T = T - 273.15 # [°C]
    rho = (0.999841594 + 6.793952e-5 * T - 9.095290e-6 * T**2 + 1.001685e-7 * T**3 - 1.120083e-9 * T**4 + \
           6.536332e-12 * T**5) + (8.25917e-4 - 4.4490e-6 * T + 1.0485e-7 * T**2 - 1.2580e-9 * T**3 + \
           3.315e-12 * T**4) * S + (-6.33761e-6 + 2.8441e-7 * T - 1.6871e-8 * T**2 + 2.83258e-10 * T **3) * S**(3/2) + \
          (5.4705e-7 - 1.97975e-8 * T + 1.6641e-9 * T**2 - 3.1203e-11 * T**3) * S**2
    return rho

def _rhoSCWater(T, S):
    """
    Density of supercooled water.
    --------------------------------------------------------------------------
    References: 
        - Hare & Sorensen (1987). doi: 10.1063/1.453710

    Parameters
    ----------
    T : FLOAT, LIST or np.ndarray
        Temperature [K].
    S : FLOAT, LIST or np.ndarray
        Salinity [ppt or g/L].

    Returns
    -------
    rho : FLOAT, LIST or np.ndarray
        Density of water [kg/L].

    """
    if not isinstance(T, np.ndarray): T = np.array(T)
    if not isinstance(S, np.ndarray): S = np.array(S)
    if not np.shape(T) == np.shape(S):
        print('!EcoSysEM.Warning: T and S arguments must have the same shape.')
        sys.exit()
    rho = 0.99986 + 6.69e-5 * T - 8.486e-6 * T**2 + 1.518e-7 * T**3 - 6.9484e-9 * T**4 - 3.6449e-10 * T**5 - 7.497e-12 * T**6 + \
         (8.25917e-4 - 4.449e-6 * T + 1.0485e-7 * T**2 - 1.258e-9 * T**3 + 3.315e-12 * T**4) * S + (-6.33761e-6 + 2.8441e-7 * T - \
          1.6871e-8 * T**2 + 2.83258e-10 * T**3) * S**(3/2) + (5.4705e-7 - 1.97975e-8 * T + 1.6641e-9 * T**2 - 3.1203e-11 * T**3) * S**2
    return rho

def density(T, S, compound = 'H2O'):
    if compound == 'H2O':
        rho = np.where(T >= 0, _rhoWater(T, S), _rhoSCWater(T, S))
    return rho

class ThP:
    """
    Class for thermodynamic parameters.
    
    """
    # Directory of databases
    path = 'db\\'
    
    def checkThP(typeParam, db, compounds, phase, warnings = False):
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
        warnings : BOOL, optional
            Display function warnings. The default is False.

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
            if warnings:
                print(f'!EcoSysEM.Warning: {typeParam} for {phase} not found for: {compNaN}.')
                print(f'>> Returned compounds: {compnotNaN}.\n')
            if typeParam == 'deltaG0f' or typeParam == 'deltaH0f':
                print(f'!EcoSysEM.Error: {typeParam} for {phase} phase not found for: {compNaN}.')
                sys.exit()
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
            Requested parameters.
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
        Function to get Gibbs free energy of reaction from DeltaG0f.

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
    
    def getKeq(compounds, mRxn, t, phase):
        """
        Function to get equilibrium constants from DeltaG0f.

        Parameters
        ----------
        compounds : STR or np.array
            Requested compound(s).
        mRxn : np.array
            Reaction matrix. (compounds)x(reactions)
        temperature: FLOAT, LIST or np.ndarray
            Absolute temperature [K]
        phase : STR
            Fluid phase. Depending on typeParam.

        Returns
        -------
        Keq : np.array
            Equilibrium constants. Shape: (Z)x(Y)x(X)x(compounds).
    
        """
        R = 0.0083144598   # Universal gas constant [kJ/mol/K]
        deltaG0f, notNaN = ThP.getThP('deltaG0f', compounds, phase)
        deltaG0r = ThP.getDeltaG0r(deltaG0f, mRxn)
        dG0r = [np.nan, np.nan, np.nan]
        dG0r[0:len(deltaG0r)] = deltaG0r
        Keq_ = [np.exp(DG0r*np.ones(t.shape) / (-R * t)) for DG0r in dG0r]
        Keq = np.stack(Keq_, axis = -1)
        Keq = np.nan_to_num(Keq, nan = 0.0)
        return np.squeeze(Keq)

    def _sumI(composition):
        """
        Function to compute the sum concentration * (charge)^2 of all compounds.

        Parameters
        ----------
        composition : DICT
            Composition of solution.
            {'compound': concentration}

        Returns
        -------
        I : FLOAT
            Sum of concentration * (charge)^2.

        """
        if not isinstance(composition, dict):
            print('!EcosysEM.Error: Argument `composition` must be a dictionary.')
            sys.exit()
        # Initialize I
        I = 0
        for compound in composition:
            Ct = composition[compound]
            if not isinstance(Ct, np.ndarray): Ct = np.array(Ct)
            if '-' in compound:
                if len(compound) == compound.index('-')+1:
                    z = -1 * np.ones(Ct.shape)
                else:
                    z = int(compound[compound.index('-'):]) * np.ones(Ct.shape)
                c = Ct
                I_ = c * z**2
            elif '+' in compound:
                if len(compound) == compound.index('+')+1:
                    z = 1 * np.ones(Ct.shape)
                else:
                    z = int(compound[compound.index('+'):]) * np.ones(Ct.shape)
                c = Ct
                I_ = c * z**2
            else:
                pd_lyte = pd.read_csv('reactions/electrolytes.csv')
                try:
                    stoic = pd_lyte[compound].dropna()
                except:
                    continue
                else:
                    comp = pd_lyte['Compounds'][stoic.index].values
                    stoic = [stoic.values[i] * np.ones(Ct.shape) for i, c in enumerate(comp) if '-' in c or '+' in c]
                    stoic = np.squeeze(stoic)
                    comp = [c for c in comp if '-' in c or '+' in c]
                    # Composition dictionary with ions of salt
                    composition_ = dict(zip(comp, stoic * Ct))
                    I_ = 0
                    for c in composition_:
                        I_s = ThP._sumI({c: composition_[c]})
                        I_ += I_s
            I += I_
        return I
    
    def ionicStrength(composition):
        """
        Function to compute ionic strength of solution.

        Parameters
        ----------
        composition : DICT
            Composition of solution.
            {'compound': concentration}

        Returns
        -------
        I : FLOAT
            Ionic strength.

        """
        I = ThP._sumI(composition)
        return 0.5 * I
    
    def _debyeHuckel(composition, T, salinity = None, molality = True, solvent = 'H2O', selComp = None):
        """
        Function to estimate activity coefficients with Debye-Hückel theory.
        --------------------------------------------------------------------------
        References: 
            - Speight, J. (2005). Lange's Handbook of Chemistry, 16th ed
            - Meissner & Peppas (1973). doi: 10.1002/aic.690190419 
            
        Parameters
        ----------
        composition : DICT
            Composition of environment {'compound': concentration} [mol/L].
        T : FLOAT, LIST or np.ndarray
            Temperature [K].
        salinity : FLOAT, LIST or np.ndarray, optional
            Salinity of solvent [ppt or g/L]. The default is 0.0.
        molality : BOOL, optional
            Select if activity calculation in molality (True) or molarity (False). The default is True.
        solvent : STRING, optional
            Solvent name. The default is 'H2O' (water).
        selComp: STR, optional
            Select component from composition dictionary. The default is None.
            
        Returns
        -------
        actCoeff : DICT
            Activity coefficient of compounds.

        """
        # Solvent density
        if not isinstance(T, np.ndarray): T = np.array(T)
        if salinity is None:
            salinity = 0.0 * np.ones(T.shape)
        rho_solv = density(T, salinity, solvent)
        a = 4.6 # Effective ionic radius (Å)
        if molality:
            # Constants of solvents (unit weight of solvent)
            A = 4.495e-6 * T**2 - 1.848e-3 * T + 6.618e-1
            B = 5.779e-7 * T**2 - 1.966e-4 * T + 3.357e-1
            # B-dot parameter (unit weight of solvent)
            by = -8.825e-13 * T**4 - 2.121e-9 * T**3 + 2.733e-6 * T**2 - 9.075e-4 * T + 1.317e-1
            for comp in composition:
                composition[comp] = composition[comp] * (1/rho_solv)
        else:
            # Constants of solvents (unit volum of solvent)
            A = 5.135e-6 * T**2 - 2.154e-3 * T + 6.973e-1
            B = 8.800e-7 * T**2 - 3.292e-4 * T + 3.491e-1
            # B-dot parameter (unit weight of solvent)
            by = (-8.825e-13 * T**4 - 2.121e-9 * T**3 + 2.733e-6 * T**2 - 9.075e-4 * T + 1.317e-1) * (1/rho_solv)
        # Ionic strength
        I = ThP.ionicStrength(composition)
        # Select compound (if `selComp != None`)
        if selComp is not None:
            composition = {selComp: composition[selComp]}
        # Initialize dictionary of activity coefficients
        actCoeff = {}
        for compound in composition:
            if '-' in compound:
                if len(compound) == compound.index('-')+1:
                    z = -1 * np.ones(I.shape)
                else:
                    z = int(compound[compound.index('-'):]) * np.ones(I.shape)
                # Activity coefficient
                Y = 10**(-(A * z**2 * I**0.5) / (1 + B * a * I**0.5) + by * I)
            elif '+' in compound:
                if len(compound) == compound.index('+')+1:
                    z = 1 * np.ones(I.shape)
                else:
                    z = int(compound[compound.index('+'):]) * np.ones(I.shape)
                # Activity coefficient
                Y = 10**(-(A * z**2 * I**0.5) / (1 + B * a * I**0.5) + by * I)
            else:
                pd_lyte = pd.read_csv('reactions/electrolytes.csv')
                try:
                    stoic = pd_lyte[compound].dropna()
                except:
                    # Electrolyte (pH speciation)
                    # Not defined in `electrolytes.csv`. Check if it is defined in `pHSpeciation.csv`
                    comp, stoic, infoRxn = Rxn.getRxnpH(compound)
                    if not comp:    
                        # If compound is not an electrolyte (pH speciation), assume ideal behaviour (Y = 1.0)
                        Y = 1.0 * np.ones(T.shape)
                    else:
                        stoic = np.squeeze(stoic)
                        if np.ndim(stoic) > 1:
                            indexRxn = np.squeeze(np.argwhere(np.char.find(infoRxn, 'First deP') == 0))
                            stoic = stoic[:, indexRxn]
                else:
                    # Electrolyte (no pH speciation)
                    comp = pd_lyte['Compounds'][stoic.index].values
                    stoic = stoic.values
                dEly = {}
                if comp is not None:
                    for iComp in comp:
                        if '-' in iComp:
                            if len(iComp) == iComp.index('-')+1:
                                z = -1 * np.ones(I.shape)
                            else:
                                z = int(iComp[iComp.index('-'):]) * np.ones(I.shape)
                            dEly['s_neg'] = [stoic[i] for i, c in enumerate(comp) if '-' in c][0]
                            # Activity coefficient
                            dEly['Y_neg'] = 10**(-(A * z**2 * I**0.5) / (1 + B * a * I**0.5) + by * I)
                        elif '+' in iComp:
                            if len(iComp) == iComp.index('+')+1:
                                z = 1 * np.ones(I.shape)
                            else:
                                z = int(iComp[iComp.index('+'):]) * np.ones(I.shape)
                            dEly['s_pos'] = [stoic[i] for i, c in enumerate(comp) if '+' in c][0]
                            # Activity coefficient
                            dEly['Y_pos'] = 10**(-(A * z**2 * I**0.5) / (1 + B * a * I**0.5) + by * I)
                    try:
                        dEly['Y_neg'], dEly['Y_pos']
                    except:
                        Y = 1.0 * np.ones(T.shape)
                    else:
                        Y = (dEly['Y_pos']**dEly['s_pos'] * dEly['Y_neg']**dEly['s_neg'])**(1 / (dEly['s_pos'] + dEly['s_neg']))
                else:
                    Y = 1.0 * np.ones(T.shape)
            actCoeff[compound] = Y
        return actCoeff
    
    def _setschenowShumpe(composition, T, salinity = None, molality = True, solvent = 'H2O', selComp = None):
        """
        Function to estimate activity coefficients with Setschenow-Shumpe equation.
        --------------------------------------------------------------------------
        References: 
            - Weisenberger & Schumpe (1996), doi: 10.1002/aic.690420130
            - F. Millero (2000). doi: 10.1016/S0304-4203(00)00011-6
            
        Parameters
        ----------
        composition : DICT
            Composition of environment {'compound': concentration} [mol/L].
        T : FLOAT, LIST or np.ndarray
            Temperature [K].
        salinity : FLOAT, LIST or np.ndarray, optional
            Salinity of solvent [ppt or g/L]. The default is 0.0.
        molality : BOOL, optional
            Select if activity calculation in molality (True) or molarity (False). The default is True.
        solvent : STR, optional
            Solvent name. The default is 'H2O' (water).
        selComp: STR, optional
            Select component from composition dictionary. The default is None.    
        
        Returns
        -------
        actCoeff : DICT
            Activity coefficient of compounds.

        """
        # Model parameters
        hi = {'H+': 0, 'Li+': 0.0754, 'Na+': 0.1143, 'K+': 0.0992,'Rb+': 0.0839,
              'Cs+': 0.0759, 'NH4+': 0.0556, 'Mg+2': 0.1694, 'Ca+2': 0.1762,
              'Sr+2': 0.1881, 'Ba+2': 0.2168, 'Mn+2': 0.1463, 'Fe+2': 0.1523,
              'Co+2': 0.1680, 'Ni+2': 0.1654, 'Cu+2': 0.1675, 'Zn+2': 0.1537,
              'Cd+2': 0.1869, 'Al+3': 0.2174, 'Cr+3': 0.0648, 'Fe+3': 0.1161,
              'La+3': 0.2297, 'Ce+3': 0.2406, 'Th+4': 0.2709, 'OH-': 0.0839,
              'HS-': 0.0851, 'F-': 0.0920, 'Cl-': 0.0318, 'Br-': 0.0269,
              'I-': 0.0039, 'NO2-': 0.0795, 'NO3-': 0.0128, 'ClO3-': 0.1348,
              'BrO3-': 0.1116, 'IO3-': 0.0913, 'ClO4-': 0.0492, 'IO4-': 0.1464,
              'CN-': 0.0679, 'SCN-': 0.0627, 'HCrO4-': 0.0401, 'HCO3-': 0.0967,
              'H2PO4-': 0.0906, 'HSO3-': 0.0549, 'CO3-2': 0.1423, 'HPO4-2': 0.1499,
              'SO3-2': 0.1270, 'SO4-2': 0.1117, 'S2O3-2': 0.1149, 'PO4-3': 0.2119,
              '[Fe(CN)6]-4': 0.3574}
        h0g = {'H2': -0.0218, 'He': -0.0353, 'Ne': -0.0080, 'Ar': 0.0057, 
               'Kr': -0.0071, 'Xe': 0.0133, 'Rn': 0.0447, 'N2': -0.0010, 
               'O2': 0.0000, 'NO': 0.0060, 'N2O': -0.0085, 'NH3': -0.0481, 
               'CO2': -0.0172, 'CH4': 0.0022, 'C2H2': -0.0159, 'C2H4': 0.0037,
               'C2H6': 0.0120, 'C3H8': 0.0240, 'C4H10': 0.0297, 'H2S': -0.0333,
               'SO2': -0.0817, 'SF6': 0.0100}
        hT = {'H2': -2.99e-4, 'He': 4.64e-4, 'Ne': -9.13e-4, 'Ar': -4.85e-4, 
              'Kr': 0.0000, 'Xe': -3.29e-4, 'Rn': -1.38e-4, 'N2': -6.05e-4, 
              'O2': -3.34e-4, 'NO': 0.0000, 'N2O': -4.79e-4, 'NH3': 0.0000, 
              'CO2': -3.38e-4, 'CH4': -5.24e-4, 'C2H2': 0.0000, 'C2H4': 0.0000,
              'C2H6': -6.01e-4, 'C3H8': -7.02e-4, 'C4H10': -7.26e-4, 'H2S': 0.0000,
              'SO2': 2.75e-4, 'SF6': 0.0000}
        hi_comp = list(hi.keys())
        h0g_comp = list(h0g.keys())
        # Solvent density
        if not isinstance(T, np.ndarray): T = np.array(T)
        if salinity is None:
            salinity = 0.0 * np.ones(T.shape)
        rho_solv = density(T, salinity, solvent)
        # Electrolyte contributions
        dict_hi = {}
        dict_ci = {}
        for compEly in composition:                    
            if compEly in h0g_comp:
                continue
            elif compEly in hi_comp:
                hi_ = hi[compEly] * np.ones(T.shape)
                ci_ = composition[compEly]
                try:
                    dict_ci[compEly]
                except:
                    dict_ci[compEly] = ci_
                    dict_hi[compEly] = hi_
                else:
                    dict_ci[compEly] += ci_
            else:
                # Check `electrolytes.csv`
                pd_lyte = pd.read_csv('reactions/electrolytes.csv')
                try:
                    stoic = pd_lyte[compEly].dropna()
                except:
                    continue
                else:
                    # Electrolyte (no pH speciation)
                    compSalt = pd_lyte['Compounds'][stoic.index].values
                    stoic = stoic.values
                    c_ = composition[compEly]
                    for idComp, iComp in enumerate(compSalt):
                        if ('-' in iComp) or ('+' in iComp):
                            ci_ = c_ * abs(stoic[idComp])
                            if iComp in hi_comp:
                                hi_ = hi[iComp] * np.ones(T.shape)
                            try:
                                dict_ci[iComp]
                            except:
                                dict_ci[iComp] = ci_
                                dict_hi[iComp] = hi_
                            else:
                                dict_ci[iComp] += ci_
        # Non-electrolyte gas contribution
        actCoeff = {}
        # Select compound (if `selComp != None`)
        if selComp is not None:
            composition = {selComp: composition[selComp]}
        for comp in composition:
            # Gas compound (non-electrolyte)
            if comp in h0g_comp:
                h0g_ = h0g[comp] * np.ones(T.shape)
                hT_ = hT[comp] * np.ones(T.shape)
                hg_ = h0g_ + hT_ * (T - 298.15)
                if molality:
                    hg_ = hg_ * (rho_solv)
                logYn = 0
                for ely in dict_hi:
                    hi_ = dict_hi[ely]
                    ci_ = dict_ci[ely]
                    if molality:
                        hi_ = hi_ * (rho_solv)
                    logYn += (hi_ + hg_) * ci_
                Yn = 10**logYn
            else:
                Yn = 1.0 * np.ones(T.shape)
            actCoeff[comp] = Yn
        return actCoeff
        
    def activity(methods, composition, T = None, pH = None, salinity = None, molality = True, solvent = 'H2O', selComp = None):
        """
        Function to compute activities of compounds.

        Parameters
        ----------
        methods : DICT
            Method for coefficient activity estimation.
                'DH-ext'    - Debye-Hückel equation extended version.
                'SS'        - Setschenow-Shumpe equation.
        composition : DICT [mol/L; molarity]
            Composition of environment {'compound': concentration} [mol/L].
        T : FLOAT, LIST or np.ndarray, optional
            Temperature [K]. The default is None.
        pH : FLOAT, optional
            pH [-]. The default is None.
        salinity : FLOAT, LIST or np.ndarray, optional
            Salinity of solvent [ppt or g/L]. The default is None.
        molality : BOOL, optional
            Select if activity calculation in molality (True) or molarity (False). The default is True.
        solvent : STRING, optional
            Solvent name. The default is 'H2O' (water).
        selComp : STRING, optional
            Selected compound (must be in `methods` and `composition`). The default is None.
            
        Returns
        -------
        act : TYPE
            Activity of compounds.

        """
        # Temperature definition
        if T is None:
            compounds = list(composition.keys())
            Ct = composition[compounds[0]]
            T = 298.15 * np.ones(Ct.shape)
        # Solvent density
        if not isinstance(T, np.ndarray): T = np.array(T)
        if salinity is None:
            salinity = 0.0 * np.ones(T.shape)
        rho_solv = density(T, salinity, solvent)
        if pH is not None:
            composition_aux = {}
            methods_aux = {}
            for iComp in composition:
                if iComp == 'CO2': 
                    iComp_ = 'H2CO3'
                else:
                    iComp_ = iComp
                # Check if {iComp} is involved in ion speciation
                rComp, _, _ = Rxn.getRxnpH(iComp_)
                if rComp is not None:
                    rComp = rComp[1:]
                    C = composition[iComp]
                    if iComp != 'H+':
                        cSpec = ThEq.pHSpeciation(iComp, pH, T, C, True)
                        for idC, C in enumerate(rComp):
                            composition_aux[C] = cSpec[..., idC]
                            if C == 'H2CO3':
                                composition_aux['CO2'] = cSpec[..., idC]
                            try:
                                methods[C]
                            except:
                                # By default, Debye-Hückel theory to estimate activity coefficient
                                methods_aux[C] = 'DH-ext'
                            else:
                                continue
            composition = {**composition, **composition_aux}
            methods = {**methods, **methods_aux}
        act = {}
        if selComp is not None:
            try:
                composition[selComp]
            except:
                print(f'!EcosysEM.Error: Selected compound ({selComp}) was not found in `composition` argument, or pH is not defined (if {selComp} is a pH-related chemical species).')
                sys.exit()
            else:
                selCompounds = [selComp]
        else:
            selCompounds = list(composition.keys())
        for comp in selCompounds:
            if molality:
                c = composition[comp] * (1/rho_solv)
            else:
                c = composition[comp]
            try:
                methods[comp]
            except:
                # By default, Debye Huckel theory is used to estimate coefficient activity
                actCoeff = ThP._debyeHuckel(composition, T, salinity, molality, solvent, selComp = comp)
                act[comp] = c * actCoeff[comp]
            else:
                if methods[comp] == 'DH-ext':
                    actCoeff = ThP._debyeHuckel(composition, T, salinity, molality, solvent, selComp = comp)
                elif methods[comp] == 'SS':
                    actCoeff = ThP._setschenowShumpe(composition, T, salinity, molality, solvent, selComp = comp)
                elif methods[comp] == 'ideal':
                    actCoeff = {comp: 1.0 * np.ones(T.shape)}
                else:
                    print(f'!EcosysEM.Error: Method to estimate activity coefficient of {comp} not defined.')
                    sys.exit()
                act[comp] = c * actCoeff[comp]
        return act
    
class ThEq:
    """
    Class for calulation of chemical, ion and interphase (G-L) equilibriums.
    
    """
    
    def solubilityHenry(compounds, wType = 'FW', t = None):
        """
        Get Henry's law solubility constants (Hs) in function of temperature.

        Parameters
        ----------
        compounds : STR or LIST
            Requested compound(s).
        wType : STR ('FW' or 'SW')
            Water type (Phase). FW - Fresh Water; SW - Sea Water.
        temperature : FLOAT or LIST, optional
            Set temperature for Henry's law solubility constant(s). The default is None.

        Returns
        -------
        Hs : np.array
            Henry's law solubility constant(s). Shape: (Z)x(Y)x(X)x(compounds).
            If no temperature is given, Hs is (1)x(compounds).
        notNaN : np.array
            Indices of parameters that are available.

        """
        if isinstance(compounds, str): compounds = [compounds]
        # Standard temperature [K]
        Ts = 298.15
        if wType == 'FW' or wType == 'SW':
            # Henry's solubilities (Hs and B)
            Hs, notNaN_Hs = ThP.getThP('Hs', compounds, wType)
            B, notNaN_B = ThP.getThP('B', compounds, wType)
            if notNaN_Hs.sum() > notNaN_B.sum():
                notNaN = notNaN_B
            else:
                notNaN = notNaN_Hs
            if t.all():
                # Initialize Hs result
                Hs_r = np.empty(t.shape)
                Hs_r = Hs_r[..., np.newaxis]
                Hs_r = np.repeat(Hs_r, len(compounds), axis = -1)
                Hs_ = [H_*np.ones(t.shape) * np.exp(B[idH]*np.ones(t.shape) * ((1 / t) - (1/Ts))) for idH, H_ in enumerate(Hs)]
                Hs = np.stack(Hs_, axis = -1)
                c = 0
                for idComp in range(len(compounds)):
                    if notNaN[idComp] == True:
                        Hs_r[..., idComp] = Hs[..., c]
                        c += 1
                    else:
                        Hs_r[..., idComp] = np.zeros(Hs[..., 0].shape)
            return Hs_r, notNaN
        else:
            print('!EcosysEM.Error: No water type selected. Use one of the following wType:\n'+
              '                 \'FW\'      - Fresh Water.\n'+
              '                 \'SW\'      - Sew Water.')
            return None, None
    
    def pHSpeciation(iCompound, pH, t, Ct, rAllConc = False):
        """
        Calculate pH (or ion) speciation of selected compounds.

        Parameters
        ----------
        compounds : STR
            Requested compound.
        pH : FlOAT
                Set of pH.
        t : FLOAT, LIST or np.array
            Set of temperature for pH speciation [K].
        Ct : LIST or np.array
            Total concentrations of compounds.
        rAllConc : BOOL, optional
            Option to select if return all compound species or only preferred. The default is False.

        Returns
        -------
        rSpec : np.array
            Concentrations of selected compound species.
            If rAllConc = False: (Z)x(Y)x(X).
            If rAllConc = True: (Z)x(Y)x(X)x(species).
            Where species: [B], [B-], [B-2], [B-3].
    
        """
        if not isinstance(t, np.ndarray): t = np.array(t)
        if not isinstance(Ct, np.ndarray): Ct = np.array(Ct)
        # Simplification hydration/dehydration equil.: [(CO2)aq] >>>> [H2CO3]
        if iCompound == 'CO2': iCompound = 'H2CO3'
        H = np.power(10, -abs(np.array(pH)))
        rComp, mRxn, infoRxn = Rxn.getRxnpH(iCompound)
        # Return same concentration(s) if compound doesn't have acid-base equilibrium
        if not rComp:
            return Ct
        c_rComp = iCompound in rComp
        if not c_rComp:
            return Ct
        reqSp = rComp.index(iCompound) - 1 # Requested chemical species
        Ka = ThP.getKeq(rComp, mRxn, t, 'L')
        # Speciation
        theta = H**3 + (Ka[..., 0] * H**2) + (Ka[..., 0] * Ka[..., 1] * H) + Ka[..., 0] * Ka[..., 1] * Ka[..., 2]
        mSpec__ = np.array([Ct * H**3 / theta,                                    # [B]
                            Ka[..., 0] * Ct * H**2 / theta,                       # [B-]
                            Ka[..., 0] * Ka[..., 1] * Ct * H / theta,             # [B-2]
                            Ka[..., 0] * Ka[..., 1] * Ka[..., 2] * Ct / theta])   # [B-3]
        # Reshape matrix of chemical species
        mSpec_ = [mSpec__[i, ...] for i in range(4)]
        mSpec_aux = np.stack(mSpec_, axis = -1)
        if rAllConc == False:
            rSpec = mSpec_aux[..., reqSp]
        else:
            nzeros = np.nonzero(mSpec_aux)
            shapes = t.shape + (len(rComp) - 1,)
            rSpec = np.reshape(mSpec_aux[nzeros], shapes)
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
            Spec_ = [ThEq.pHSpeciation(iCompound, pH_, temperature, Ct, True) for pH_ in pH]
            Spec = np.stack(Spec_, axis = 0)
            nFrac = (Spec / Ct) * 100 # Molar fraction [%]
            # Selection of involved chemical species
            c_nFrac = np.nonzero(np.sum(nFrac, axis = 0))
            nFrac = nFrac[:, c_nFrac[0]]
            # Get name of chemical species
            nCompounds = Rxn.getRxnpH(iCompound)[0][1:]
            # Simplification hydration/dehydration equil.: [(CO2)aq] >>>> [H2CO3]
            nCompounds = [nC.replace('H2CO3', 'H2CO3 (CO2)') for nC in nCompounds]
            text_ = f'Chemical (or ion) speciation at {float(f'{temperature:.2f}')}K.'
            # Plotting
            fig, ax = plt.subplots()
            ax.plot(pH, nFrac)
            ax.set_ylabel('Molar fraction (%)')
            ax.set_xlabel('pH')
            ax.set_yticks(np.arange(0, 110, 10))        
            fig.text(0.5, 0.025, text_, horizontalalignment = 'center', wrap = True)
            fig.tight_layout(rect=(0,.05,1,1)) 
            ax.margins(x = 0)
            plt.legend(nCompounds, loc = 'center left', bbox_to_anchor = (1, 0.5))
            plt.show()

class ThSA:
    """
    Class for thermodynamic state analysis of environment.
    
    """

    def getDeltaGr(typeRxn, input_, phase, specComp = False, T = 298.15, pH = 7.0, S = None, Ct = 1.0,
                   fluidType = 'ideal', molality = True, methods = None, solvent = 'H2O', asm = 'stoich', 
                   warnings = False, printDG0r = False, printDH0r = False):
        """
        Calculate DeltaGr in function of pH, temperature and compound
        concentrations.

        Parameters
        ----------
        typeRxn : STR
            What reaction(s) type are requested, matching with csv name. E.g.:
                - 'metabolisms': metabolic activities.
        input_ : STR or LIST
            Name(s) of requested compound(s) or reaction(s).
        phase: STR
            Phase in which reaction(s) ocurr. 'G' - Gas, 'L' - Liquid.
        specComp : (if input_ is reactions; STR or LIST) or (if input_ is compounds; BOOL - True), optional
            Name(s) of compound(s) to calculate specific deltaGr (kJ/mol-compound). The default is False.
        T : FLOAT or LIST, optional
            Set of temperature [K]. The default is 298.15 K (standard temperature).
        pH : INT or FLOAT, optional
            Set of pH. The default is 7.0 (neutral pH).
        S : FLOAT, LIST or np.array
            Salinity [ppt]. The default is None.
        Ct : DICT
            Total concentrations of compounds {'compounds': [concentrations]}.
            All compounds of a reaction with the same number of concentrations.
        fluidType : STR, optional
            Type of fluid (ideal or non-ideal). The default is ideal.
        molality : BOOL, optional
            Select if activity units are in molality (True) or molarity (False). The default is True.
        methods : DICT, optional
            Method for coefficient activity estimation. The default is None.
                'DH-ext'    - Debye-Hückel equation extended version.
                'SS'        - Setschenow-Shumpe equation.
        solvent : STRING, optional
            Solvent name. The default is 'H2O' (water).
        asm : STRING, optional
            Assumption when products are not present in the environment.
            The default is 'stoich' (stoichiometric concentrations).
        warnings : BOOL, optional
            Display function warnings. The default is False.
        printDG0r : BOOL, optional
            Print in console the values of standard Gibbs free energy of reactions. The default is False.
        printDH0r : BOOL, optional
            Print in console the values of standard enthalpy of reactions. The default is False.

        Returns
        -------
        DGr : np.array
            Non-standard Gibbs free energy values.
            Shape: (Z)x(Y)x(X)x(reactions).
        infoRxn : LIST
            Name of reactions given by the user (see reaction/{typeRxn}.csv).

        """
        if not isinstance(typeRxn, str): typeRxn = str(typeRxn)
        if not isinstance(input_, list): input_ = [input_]
        if phase != 'G' and phase != 'L':
            print('!EcoSysEM.Error: `phase` argument must be "G" (gas) or "L" (liquid)')
            sys.exit()
        if isinstance(T, int): T = float(T)
        if isinstance(T, float): T = [T]
        if not isinstance(T, np.ndarray): T = np.array(T)
        if isinstance(pH, int): pH = float(pH)
        # Check shape of temperature (if `Ct` is a dictionary)
        if isinstance(Ct, dict):
            shapeT = T.shape
            shapeC = Ct[list(Ct.keys())[0]].shape
            if shapeT[0] == 1:
                T = T * np.ones(shapeC)
            else:
                print(f'!EcoSysEM.Error: `T` shape {shapeT} and `Ct` keys shape {shapeC} doesn''t match.')
                sys.exit()
        # Get reactions
        rComp, mRxn, infoRxn = Rxn.getRxn(typeRxn, input_, warnings)
        nRxn = infoRxn.size
        # Initialise variables
        Ts = 298.15                                                             # Standard temperature [K]
        R = 0.0083144598                                                        # Universal gas constant [kJ/mol/K]
        # Initializate DGr matrix
        DGr = np.empty(T.shape)
        DGr = DGr[..., np.newaxis]
        DGr = np.repeat(DGr, len(input_), axis = -1)
        if specComp:
            if specComp == True:
                tSpecComp = 'compounds'
                specComp = input_
                for iSpecComp in specComp:
                    c_specComp = np.squeeze(np.where(np.array(rComp) == iSpecComp))
                    if c_specComp.size == 0:
                        print(f'!EcoSysEM.Error: {iSpecComp} was not found as a compound in {typeRxn}.csv file. ' 
                              'Use the `specComp` argument to specify compounds or set it to False.')
                        sys.exit()
            else:
                tSpecComp = 'reactions'
                if isinstance(specComp, str): specComp = [specComp]
                n_specComp = len(specComp)
                if nRxn != n_specComp:
                    print(f'!EcoSysEM.Error: Different number of reactions and specific compounds was found. # reactions: {nRxn}; # specific compounds: {n_specComp}')
                    sys.exit()
        # Select reactions w/ requested compounds as substrates (if input_ = compounds) (How ?)
        for idRxn, iRxn in enumerate(infoRxn):
            # iVariables definition
            rDGr = np.empty(T.shape)            # Initialise auxiliar result matrix (rDGr)
            i_mRxn = np.c_[mRxn[:, idRxn]]      # mRxn_aux
            cNonZ = np.nonzero(i_mRxn)[:][0]
            i_rComp = np.array(rComp)[cNonZ]    # rComp_aux
            i_mRxn = i_mRxn[cNonZ]              # mRxn_aux
            if specComp:
                vSelected = 0
                if tSpecComp == 'compounds':
                    for specComp_aux in specComp:
                        c_specComp = np.squeeze(np.where(i_rComp == specComp_aux))
                        if c_specComp.size > 0:
                            if vSelected == 0:
                                i_specComp = specComp_aux
                            else:
                                print('!EcoSysEM.Error: More than one specific compound has been used. Only one by reaction.')
                                sys.exit()
                else:
                    i_specComp = specComp[idRxn]
                    # Check if i_specComp are in reaction
                    c_specComp = np.squeeze(np.where(i_rComp == i_specComp))
                    if c_specComp.size == 0:
                        print(f'!EcoSysEM.Error: {i_specComp} was not found in reaction {iRxn}.')
                        sys.exit()
                # Stoichiometric parameter of selected compound
                id_specComp = np.squeeze(np.where(i_rComp == i_specComp))
                vSelected = abs(np.squeeze(i_mRxn[id_specComp]))
            else:
                vSelected = 1.0
            # Calculate DeltaG0r
            deltaG0f = ThP.getThP('deltaG0f', i_rComp, phase)[0]
            deltaG0r = ThP.getDeltaG0r(deltaG0f, i_mRxn)                        # kJ
            deltaG0r = deltaG0r / vSelected                                     # kJ/mol-i (if specDGr)
            if printDG0r:
                print(f'· DG0r of {iRxn}: {deltaG0r}.')
            # Calculate DeltaH0r
            deltaH0f = ThP.getThP('deltaH0f', i_rComp, phase)[0]
            deltaH0r = ThP.getDeltaH0r(deltaH0f, i_mRxn)                        # kJ
            deltaH0r = deltaH0r / vSelected                                     # kJ/mol-i (if specDGr)
            if printDH0r:
                print(f'· DH0r of {iRxn}: {deltaH0r}.')
            if printDG0r or printDH0r:
                print('')
            deltaGTr = deltaG0r * (T / Ts) + deltaH0r * ((Ts - T) / Ts)
            if isinstance(Ct, dict):
                # Calculate reaction quotient (Qr)
                Qr = 0
                uComp = np.array(list(Ct.keys()))                       # Compounds given by user (Ct)
                findSpecComp = np.argwhere(uComp == i_specComp).squeeze()
                for idComp, iComp in enumerate(i_rComp):
                    findComp = np.argwhere(uComp == iComp)
                    vi = i_mRxn[idComp] / vSelected
                    if fluidType == 'ideal':
                        if findComp.size == 0:
                            if iComp == 'H+':
                                iConc = 10**(-pH) * np.ones(T.shape)
                            elif iComp == 'H2O':                # It is assumed a water activity of 1.0, as a pure liquid (it can be less).
                                iConc = 1.0 * np.ones(T.shape)
                            else:
                                # Check if product is a species of iComp
                                rComp_pH, _, _ = Rxn.getRxnpH(iComp)
                                uIComp = np.isin(uComp, rComp_pH)
                                if (rComp is None) or (any(uIComp) == False):
                                    if asm == 'stoich': # [P] is calculated based on stoichiometry.
                                        iConc = (vi) * Ct[uComp[findSpecComp]]
                                else:
                                    rxn_iComp =  uComp[uIComp][0]
                                    # uComp = np.char.replace(uComp, rxn_iComp, iComp) # ???
                                    iConc = Ct[rxn_iComp]
                        else:
                            iConc = Ct[iComp]
                        # pH speciation
                        if iComp != 'H+' and iComp != 'H2O' and iComp != 'OH-':
                            # Only pH speciation in liquid
                            if phase == 'L':
                                if iComp == 'CO2': iComp = 'H2CO3'      # Simplification hydration/dehydration equil.: [(CO2)aq] >>>> [H2CO3]
                                iConc = ThEq.pHSpeciation(iComp, pH, T, iConc)
                        iAct = iConc
                    elif fluidType == 'non-ideal':
                        if methods is None:
                            print("!EcoSysEM.Error: `methods` must be defined to calculate activities of species.")
                            sys.exit()
                        if iComp == 'H+':
                            Ct[iComp] = 10**(-pH) * np.ones(T.shape)
                        if (iComp == 'H2O' and solvent == 'H2O') or ('(s)' in iComp):
                            iAct = 1.0 * np.ones(T.shape)
                        else:
                            if phase == 'G':
                                iAct = iConc
                                print(f"!EcoSysEM.Warning: Estimation of fugacities not included. Ideal behaviour of {iComp} is assumed.")
                            elif phase == 'L':
                                activity = ThP.activity(methods, Ct, T, pH, S, molality, solvent, iComp)
                                iAct = activity[iComp]
                    else:
                        print("!EcoSysEM.Error: `fluidType` argument must be 'ideal' {Qr = f(concentrations)} or 'non-ideal' {Qr = f(activities)}.")
                        sys.exit()
                    # Calculation of reaction quotient (Qr)
                    Qr += vi * np.log(iAct)
                # Activity/Concentration influence (Nernst equation)
                deltaGr = deltaGTr + R * T * Qr
                rDGr = deltaGr
            else:
                rDGr = deltaGTr
            DGr[..., idRxn] = rDGr
        return np.squeeze(DGr), infoRxn
    
    def exportDeltaGr(modeExport, typeRxn, input_, phase, T, pH = 7.0, S = None, Ct = 1.0,
                      specComp = False, altitude = False, fluidType = 'ideal', molality = True, 
                      methods = None, solvent = 'H2O', asm = 'stoich', warnings = False,
                      printDG0r = False, printDH0r = False):
        """
        Export DeltaGr in function of pH and temperature or concentrations.

        Parameters
        ----------
        modeExport : STR
            How DeltaGr is exported: 'plot': plot in Spyder; 'Excel': write in Excel document.
        typeRxn : STR
            What reaction(s) type are requested, matching with csv name. E.g.:
                - 'metabolisms': metabolic activities.
        input_ : STR or LIST
            Name(s) of requested compound(s) or reaction(s).
        phase: STR
            Phase in which reaction(s) ocurr. 'G' - Gas, 'L' - Liquid.
        T : LIST or np.array
            Set of temperatures [K].
        pH : LIST or np.array, optional
            Set of pH values. The default is 7.0.
        S : FLOAT, LIST or np.array
            Salinity [ppt]. The default is None.
        Ct : DICT
            Total concentrations of compounds {'compounds': [concentrations]}.
            All compounds of a reaction with the same number of concentrations.
        specComp : (if input_ is reactions; STR or LIST) or (if input_ is compounds; BOOL - True)
           Name(s) of compound(s) to calculate specific deltaGr (kJ/mol-compound). Default: False.
           The default is 1.0.
        altitude : BOOL, optional
            Altitude in coordenate y of plots. The default is False.
        fluidType : STR, optional
            Type of fluid (ideal or non-ideal). The default is ideal.
        molality : BOOL, optional
            Select if activity calculation in molality (True) or molarity (False). The default is True.
        methods : DICT, optional
            Method for coefficient activity estimation. The default is None.
                'DH-ext'    - Debye-Hückel equation extended version.
                'SS'        - Setschenow-Shumpe equation.
        solvent : STRING, optional
            Solvent name. The default is 'H2O' (water).
        asm : STRING, optional
            Assumption when products are not present in the environment.
            The default is 'stoich' (stoichiometric concentrations).
        warnings : BOOL, optional
            Display function warnings. The default is False.
        printDG0r : BOOL, optional
            Print in console the values of standard Gibbs free energy of reactions. The default is False.
        printDH0r : BOOL, optional
            Print in console the values of standard enthalpy of reactions. The default is False.

        Returns
        -------
        Plot(s).

        """
        # Checking arguments and variable definition
        if not isinstance(modeExport, str): modeExport = str(modeExport)
        ### Calculations here
        if isinstance(pH, int): pH = float(pH)
        if isinstance(pH, float): pH = [pH]
        if isinstance(T, np.ndarray):
            if np.ndim(T) > 1:
                print('!EcoSysEM.Error: temperature (`T`) must be a list or 1D np.ndarray.')
                sys.exit()
        ### If modeExport == 'plot', T must be a list or np.ndarray with ndim = 1
        # Initializate DGr matrix
        DGr = np.empty(T.shape)
        DGr = DGr[..., np.newaxis]
        DGr = np.repeat(DGr, len(input_), axis = -1) # Compounds
        DGr = DGr[..., np.newaxis]
        DGr = np.repeat(DGr, len(pH), axis = -1) # pH
        for idpH, pH_ in enumerate(pH):
            DGr_, infoRxn  = ThSA.getDeltaGr(typeRxn, input_, phase, specComp, T, pH_, S, Ct,
                                             fluidType, molality, methods, solvent, asm, 
                                             warnings, printDG0r, printDH0r)
            DGr[..., idpH] = DGr_
        if modeExport == 'Excel':
            path = 'results/'
            nameDocument = input(' > Name of result document: ')
            fullPathSave = path + nameDocument + '.xlsx'
            if os.path.isfile(fullPathSave):
                val = input(' > Do you want to overwrite `' + nameDocument + '.xlsx`? [Y/N]: ')
                if val == 'Y' or val == 'y':       
                    os.remove(fullPathSave)
            if not isinstance(altitude, bool):
                ThSA._writeExcel(DGr, infoRxn, fullPathSave, Ct, pH, altitude, True)
            else:
                ThSA._writeExcel(DGr, infoRxn, fullPathSave, Ct, pH, T)
        elif modeExport == 'plot':
            text_ = 'Concentrations and temperature associated with altitude.'
            for idRxn, rxn in enumerate(infoRxn):
                if not isinstance(altitude, bool):
                    ThSA._contourfAlt(pH, altitude, DGr[:, idRxn, :], rxn, text_)
                else:
                    ThSA._contourfT(pH, T, DGr[:, idRxn, :], rxn, text_)
        else:
            print("!EcoSysEM.Error: `modeExport` argument must be 'plot' (for plotting in Spyder) or 'Excel' (for writing in Excel document).")
            sys.exit()
    
    def _writeExcel(DGr, infoRxn, fullPathSave, Ct, pH, y, altitude = False):
        """
        Write calculated DeltaGr in Excel document.

        """
        nameSheet_DGr = 'DGr'
        introRowDGr = pd.DataFrame(np.array([f'∆Gr in Excel | Reactions: {infoRxn}']))
        introRowCt = pd.DataFrame(np.array(['Aerosol concentrations (mol/L)']))
        if not os.path.isfile(fullPathSave):
            with pd.ExcelWriter(fullPathSave) as writer:
                introRowDGr.to_excel(writer, sheet_name = nameSheet_DGr, index = False, header = False)
        pH = pd.DataFrame(np.round(np.array([pH]), 3))
        if altitude:
            yCt = pd.DataFrame({'Alt. (km)': np.round(y / 1000, 3)})
            yDGr = pd.DataFrame({'Alt. (km)| pH': np.round(y / 1000, 3)})
        else:
            yDGr = pd.DataFrame({'T (K)|pH': np.round(y, 3)})
            yCt = pd.DataFrame({'T (K)': np.round(y, 3)})
        if isinstance(Ct, dict):
            nameSheet_Ct = 'Ct'
            dFrame_Ct = pd.DataFrame(Ct)
            with pd.ExcelWriter(fullPathSave, engine='openpyxl', mode = 'a', if_sheet_exists='overlay') as writer: 
                introRowCt.to_excel(writer, sheet_name = nameSheet_Ct, index = False, header = False)
                yCt.to_excel(writer, sheet_name = nameSheet_Ct, startrow = 2, startcol = 1, index = False, header = True)
                dFrame_Ct.to_excel(writer, sheet_name = nameSheet_Ct, startrow = 2, startcol = 2, index = False, header = True)
        for idRxn, rxn in enumerate(infoRxn):
            infoR = pd.DataFrame(np.array([rxn]))
            df = pd.DataFrame(DGr[:, idRxn, :])
            with pd.ExcelWriter(fullPathSave, engine='openpyxl', mode = 'a', if_sheet_exists='overlay') as writer:
                sRow = writer.sheets[nameSheet_DGr].max_row
                infoR.to_excel(writer, sheet_name = nameSheet_DGr, startrow = sRow+1, startcol = 1, index = False, header = False)
                pH.to_excel(writer, sheet_name = nameSheet_DGr, startrow = sRow+2, startcol = 2, index = False, header = False)
                yDGr.to_excel(writer, sheet_name = nameSheet_DGr, startrow = sRow+2, startcol = 1, index = False, header = True)
                df.to_excel(writer, sheet_name = nameSheet_DGr, startrow = sRow+3, startcol = 2, index = False, header = False)
        
    def _contourfT(pH, T, DGr_plot, iRxn, text_):
        """
        Specific `contourf()` function (from matplotlib.pyplot) for plotting DGr.
        
        """
        fig, ax = plt.subplots()
        F = ax.contourf(pH, T, DGr_plot, cmap = 'hot')
        clb = fig.colorbar(F)
        clb.set_label('∆Gr (kJ/mol)')
        ax.set_xlabel('pH')
        ax.set_ylabel('Temperature (K)')
        fig.text(0.5, 0.025, text_, horizontalalignment = 'center', wrap = True)
        fig.tight_layout(rect=(0,.05,1,1)) 
        plt.title(iRxn)
        plt.gca().invert_yaxis()
        plt.show()
    
    def _contourfAlt(pH, Alt, DGr_plot, iRxn, text_):
        """
        Specific `contourf()` function (from matplotlib.pyplot) for plotting DGr.
        
        """
        fig, ax = plt.subplots()
        F = ax.contourf(pH, Alt/1000, DGr_plot, cmap = 'hot')
        clb = fig.colorbar(F)
        clb.set_label('∆Gr (kJ/mol)')
        ax.set_xlabel('pH')
        ax.set_ylabel('Altitude (km)')
        fig.text(0.5, 0.025, text_, horizontalalignment = 'center', wrap = True)
        fig.tight_layout(rect=(0,.05,1,1)) 
        plt.title(iRxn)
        plt.show()