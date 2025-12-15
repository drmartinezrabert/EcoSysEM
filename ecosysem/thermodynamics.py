# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 18:09:14 2024

@author: 2424069M
"""
from reactions import Reactions as Rxn

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import Locator
from matplotlib.patches import Rectangle
from matplotlib import colors as clr
import matplotlib.ticker as tkr
import os.path
from itertools import combinations
from scipy import stats
from scipy.stats import qmc
import warnings
#-DEBUGGING-#
import time
#-----------#

class MinorSymLogLocator(Locator):
    """
    Dynamically find minor tick positions based on the positions of
    major ticks for a symlog scaling.
    """
    def __init__(self, linthresh):
        """
        Ticks will be placed between the major ticks.
        The placement is linear for x between -linthresh and linthresh,
        otherwise its logarithmically
        """
        self.linthresh = linthresh

    def __call__(self):
        'Return the locations of the ticks'
        majorlocs = self.axis.get_majorticklocs()

        # iterate through minor locs
        minorlocs = []

        # handle the lowest part
        for i in range(1, len(majorlocs)):
            majorstep = majorlocs[i] - majorlocs[i-1]
            if abs(majorlocs[i-1] + majorstep/2) < self.linthresh:
                ndivs = 10
            else:
                ndivs = 9
            minorstep = majorstep / ndivs
            locs = np.arange(majorlocs[i-1], majorlocs[i], minorstep)[1:]
            minorlocs.extend(locs)

        return self.raise_if_exceeds(np.array(minorlocs))

    def tick_values(self, vmin, vmax):
        raise NotImplementedError('Cannot get tick locations for a '
                                  '%s type.' % type(self))

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
    if not np.shape(T) == np.shape(S): raise ValueError(f' Argument T ({T.shape}) and argument S ({S.shape}) must have the same shape.')
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
    if not np.shape(T) == np.shape(S): raise ValueError(f' Argument T ({T.shape}) and argument S ({S.shape}) must have the same shape.')
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
    path = 'db/'
    
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
            if typeParam == 'deltaG0f' or typeParam == 'deltaH0f': raise ValueError(f'{typeParam} for {phase} phase not found for: {compNaN}.')
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
        try:
            dParam = dParam.set_index('Formula').loc[compounds].reset_index()
        except:
            Param = np.empty(0)
        else:
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
    
    def getDeltaCp(Cpi, mRxn):
        """
        Function to get heat capacity of reaction from Cpi.

        Parameters
        ----------
        Cpi : LIST or np.array
            Specific heat capacity.
        mRxn : np.array
            Reaction matrix. (compounds)x(reactions)

        Returns
        -------
        deltaCpi : np.array
            Heat capacity of reaction.
    
        """
        deltaCpi = Cpi @ mRxn
        return deltaCpi
    
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
        if not isinstance(composition, dict): raise TypeError('Argument `composition` must be a dictionary.')
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
                    if iComp != 'H+' and iComp != 'H2O' and iComp != 'OH-':
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
            except: raise ValueError(f'Selected compound ({selComp}) was not found in `composition` argument, or pH is not defined (if {selComp} is a pH-related chemical species).')
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
                else: raise ValueError(f'Method to estimate activity coefficient of {comp} not defined.')
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
        nnzeros_Ka = np.nonzero(Ka)
        shape_nnz_Ka = Ka[nnzeros_Ka].shape
        shape_nnz_Ka = (shape_nnz_Ka[0]+1,)
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
            if np.isnan(mSpec_aux).any():
                rSpec = np.nan * np.ones(shape_nnz_Ka)
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
        if not isinstance(temperature, float): raise ValueError('Argument `temperature` must be a float.')
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
    def _check_arguments(varNames, varValues, varTypes):
        for idVar, varName in enumerate(varNames):
            varValue = varValues[idVar]
            varType = varTypes[idVar]
            varType_str = str(varType).replace('<class \'','').replace('\'>', '').replace(',', ' or').replace('(', '').replace(')','')
            if not isinstance(varValue, varType): raise TypeError(f'Argument `{varName}` must be {varType_str}.')
    
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
        if phase != 'G' and phase != 'L': raise ValueError('Argument `phase` must be \'G\' (gas) or \'L\' (liquid).')
        if isinstance(T, int): T = float(T)
        if isinstance(T, float): T = [T]
        if not isinstance(T, np.ndarray): T = np.array(T)
        if isinstance(pH, int): pH = float(pH)
        # Check shape of temperature (if `Ct` is a dictionary)
        if isinstance(Ct, dict):
            shapeT = T.shape
            shapeC = Ct[list(Ct.keys())[0]].shape
            if shapeT != shapeC:
                if shapeT[0] == 1:
                    T = T * np.ones(shapeC)
                else: raise ValueError(f' Argument T ({T.shape}) and argument `Ct` keys ({shapeC}) must have the same shape.')
        # Get reactions
        rComp, mRxn, infoRxn = Rxn.getRxn(typeRxn, input_, warnings)
        nRxn = infoRxn.size
        # Initialize variables
        Ts = 298.15                                                             # Standard temperature [K]
        R = 0.0083144598                                                        # Universal gas constant [kJ/mol/K]
        # Initialize DGr matrix
        DGr = np.empty(T.shape)
        DGr = DGr[..., np.newaxis]
        DGr = np.repeat(DGr, len(input_), axis = -1)
        if specComp:
            if specComp == True:
                tSpecComp = 'compounds'
                specComp = input_.copy()
                for iSpecComp in specComp:
                    c_specComp = np.squeeze(np.where(np.array(rComp) == iSpecComp))
                    if c_specComp.size == 0: 
                        raise ValueError(f'{iSpecComp} was not found as a compound in {typeRxn}.csv file. Use the `specComp` argument to specify compounds or set it to False.')
            else:
                tSpecComp = 'reactions'
                if isinstance(specComp, str): specComp = [specComp]
                n_specComp = len(specComp)
                if nRxn != n_specComp: raise ValueError(f'Different number of reactions and specific compounds was found. # reactions: {nRxn}; # specific compounds: {n_specComp}.')
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
                            else: raise ValueError('More than one specific compound has been used. Only one compound per reaction.')
                else:
                    i_specComp = specComp[idRxn]
                    # Check if i_specComp are in reaction
                    c_specComp = np.squeeze(np.where(i_rComp == i_specComp))
                    if c_specComp.size == 0: raise ValueError(f'{i_specComp} was not found in reaction {iRxn}.')
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
            # Calculate DeltaHr (if necessary)
            Cpi, notNan_Cpi = ThP.getThP('Cpi', i_rComp, phase)                 # J/mol-i/K
            check_Cpi = notNan_Cpi.all()
            if check_Cpi:
               deltaCp = ThP.getDeltaCp(Cpi, i_mRxn) / 1000                     # kJ/mol-i/K
               deltaHr = deltaH0r + deltaCp * (T - Ts)                          # kJ/mol-i (if specDGr)
            else:
               deltaHr = deltaH0r
            deltaGTr = deltaG0r * (T / Ts) + deltaHr * ((Ts - T) / Ts)
            if isinstance(Ct, dict):
                # Calculate reaction quotient (Qr)
                Qr = 0
                uComp = np.array(list(Ct.keys()))                       # Compounds given by user (Ct)
                if specComp:
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
                                if 'H2CO3' in rComp_pH:
                                    rComp_pH += ['CO2']
                                uIComp = np.isin(uComp, rComp_pH)
                                if not any(uIComp):
                                    if asm == 'stoich': 
                                        # [P] is calculated based on stoichiometry.
                                        if specComp:
                                            iConc = (vi) * Ct[uComp[findSpecComp]]
                                        else: raise ValueError('`specComp` must be given to calculate the stoichiometric concentration of {iComp}.')
                                else:
                                    rxn_iComp =  uComp[uIComp][0]
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
                        if methods is None: raise ValueError('Argument `methods` must be defined to calculate activities of species.')
                        if iComp == 'H+':
                            Ct[iComp] = 10**(-pH) * np.ones(T.shape)
                        # Activity estimation
                        if (iComp == 'H2O' and solvent == 'H2O') or ('(s)' in iComp):
                            iAct = 1.0 * np.ones(T.shape)
                        else:
                            if phase == 'G':
                                iAct = iConc
                                print(f"!EcoSysEM.Warning: Estimation of fugacities not included. Ideal behaviour of {iComp} is assumed.")
                            elif phase == 'L':
                                activity = ThP.activity(methods, Ct, T, pH, S, molality, solvent, iComp)
                                iAct = activity[iComp]
                    else: raise ValueError(f'Unknown fluidType ({fluidType}). Existing fluidType: \'ideal\' or \'non-ideal\'.')
                    # Calculation of reaction quotient (Qr)
                    Qr += vi * np.log(iAct)
                # Activity/Concentration influence (Nernst equation)
                deltaGr = deltaGTr + R * T * Qr
                rDGr = deltaGr
            else:
                rDGr = deltaGTr
            DGr[..., idRxn] = rDGr
        return DGr, infoRxn
    
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
            if np.ndim(T) > 1: raise TypeError('Temperature argument (`T`) must be a list or 1D np.ndarray.')
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
        else: raise ValueError(f'Unknown modeExport ({modeExport}). Existing modeExport: \'plot\' (for plotting in Spyder) or \'Excel\' (for writing in Excel document).')
    
    def smmryDeltaGr(typeRxn, input_, specComp, phase, T, pH, S, Ct, fluidType = 'ideal', molality = True, 
                     methods = None, renameRxn = None, write_results_csv = False, logScaleX = True, vmin = None,
                     vmax = None, printDG0r = False, printDH0r = False, showMessage = True):
        """
        Create a summary plot with the range of Gibbs free energy of a set of reactions at a specific
        range of T and pH. 

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
        T : LIST
            Set of temperatures [K].
        pH : LIST, optional
            Set of pHs.
        S : LIST or np.array
            Salinity [ppt].
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
        renameRxn : None or DICT, optional
            If it's a DICT, change de name of reactions of .csv file in the plot. {'originalName': 'NewName'}
            The default is None.
        write_results_csv : BOOL, optional
            Write DGr values in a .csv file. The default is False.
        logScaleX : BOOL, optional
            If True, DGr is plotted using symmetrical log coordinate. The default is True.
        vmin : None or FLOAT, optional
            Set minimum value (left) of coordinate-X. The default is None.
        vmax : None or FLOAT, optional
                Set maximum value (right) of coordinate-X. The default is None.
        printDG0r : BOOL, optional
            Print in console the values of standard Gibbs free energy of reactions. The default is False.
        printDH0r : BOOL, optional
            Print in console the values of standard enthalpy of reactions. The default is False.
        showMessage : BOOL, optional
             Boolean to set whether informative messages are displayed in Console. The default is True.

        Returns
        -------
        None.

        """
        if showMessage:
            print('  > Creating summary of thermodynamic state analysis...')
        if renameRxn:
            if not isinstance(renameRxn, dict):
                raise TypeError(f'Argument \'renameRxn\' must be a dictionary ({type(renameRxn)}).')
        # Variable initialization    
        DGr = np.empty((len(pH), len(T), 1))
        DGr_m = np.empty((len(input_), len(pH), len(T)))    
        iC = 0
        for id_, i_ in enumerate(input_):
            iC += 1
            specCp_ = specComp[id_]
            for idpH, pH_ in enumerate(pH):
                if idpH == 0:
                    if printDG0r: 
                        pDG0r = True
                    else:
                        pDG0r = False
                    if printDH0r: 
                        pDH0r = True
                    else:
                        pDH0r = False
                else:
                    pDG0r = False
                    pDH0r = False
                try:
                    DGr_, infoRxn = ThSA.getDeltaGr(typeRxn = typeRxn,
                                                    input_ = i_, 
                                                    phase = phase, 
                                                    specComp = specCp_,
                                                    T = T,
                                                    pH = pH_,
                                                    S = S, 
                                                    Ct = Ct,
                                                    fluidType = fluidType,
                                                    molality = molality,
                                                    methods = methods,
                                                    printDG0r = pDG0r,
                                                    printDH0r = pDH0r)
                except:
                    DGr[idpH, ...] = np.nan
                    infoRxn = i_
                    if idpH == 0:
                        print(f'Error calculating DGr of {i_}.')
                else:
                    if np.ndim(DGr_) == 1: DGr_ = np.expand_dims(DGr_, axis = 1)
                    DGr[idpH, ...] = DGr_
            DGr_m[id_, ...] = np.squeeze(DGr)
        if not vmin:
            vmin = np.sign(np.nanmin(DGr_m)) * 1.1 * abs(np.nanmin(DGr_m))
        if not vmax:
            vmax = np.sign(np.nanmax(DGr_m)) * 0.95 * abs(np.nanmax(DGr_m))
        # Numpy to csv
        if write_results_csv:
            header = pd.DataFrame(np.array([infoRxn]))
            df = pd.DataFrame(np.squeeze(DGr), index = pH_)
            header.to_csv('DGr.csv', mode='a', header=False, index=False)
            df.to_csv('DGr.csv', mode='a', header=False)
        # Plotting Summary of DGr
        fig, ax = plt.subplots(figsize = (5.0, 11.0))
        for idRxn, rxn in enumerate(input_):
            DGr_rxn = DGr_m[idRxn, ...]
            if np.all(np.isnan(DGr_rxn)):
                ax.axhline(y = idRxn+1, color = 'darkgrey', linewidth = 0.5, linestyle = '--')
                continue
            min_DGr = np.min(DGr_rxn)
            max_DGr = np.max(DGr_rxn)
            if min_DGr > 0 and max_DGr > 0:
                c_ = 'dodgerblue'
            elif min_DGr < 0 and max_DGr > 0:
                c_ = ['tomato', 'dodgerblue']
            else:
                c_ = 'tomato'
            ax.axhline(y = idRxn+1, color = 'darkgrey', linewidth = 0.5, linestyle = '--')
            h = 0.5
            if isinstance(c_, str):
                ax.add_patch(Rectangle((min_DGr, idRxn+1-(h/2)),
                                       width = max_DGr - min_DGr,
                                       height = h,
                                       fc = c_,
                                       ec = c_,
                                       lw = 1.0,
                                       zorder = 2.5))
            else:
                ax.add_patch(Rectangle((min_DGr, idRxn+1-(h/2)),
                                       width = abs(min_DGr),
                                       height = h,
                                       fc = c_[0],
                                       ec = c_[0],
                                       lw = 1.0,
                                       zorder = 2.5))
                ax.add_patch(Rectangle((0.0, idRxn+1-(h/2)),
                                       width = max_DGr,
                                       height = h,
                                       fc = c_[1],
                                       ec = c_[1],
                                       lw = 1.0,
                                       zorder = 2.5))
        ax.vlines(0.0, ymin = 0.0, ymax = len(input_), color = 'k', linewidth = 0.8, zorder = 2.6)
        # x-label
        ax.set_xlabel('∆Gr (kJ/mol)', labelpad = 7.0)
        ax.set_xlim(left = vmin, right = vmax)
        ax.xaxis.set_ticks_position('top')
        ax.xaxis.set_label_position('top')
        # y-label
        positionsY = np.arange(1, len(input_)+1)
        if renameRxn:
            for rxn_rename in renameRxn:
                try:
                    ind = input_.index(rxn_rename)
                except:
                    continue
                else:
                    input_[ind] = renameRxn[rxn_rename] 
        labelsY = input_
        ax.yaxis.set_major_locator(tkr.FixedLocator(positionsY))
        ax.yaxis.set_major_formatter(tkr.FixedFormatter(labelsY))
        ax.set_ylim(top = len(input_)+2, bottom = 0)
        ax.invert_yaxis()
        ax.margins(y = 0)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        if logScaleX:
            plt.xscale('symlog')
            xaxis = plt.gca().xaxis
            xaxis.set_minor_locator(MinorSymLogLocator(1e-1))
        else:
            # plt.minorticks_on()
            ax.xaxis.set_minor_locator(tkr.AutoMinorLocator())
        plt.show()
        if showMessage:
            print('  > Done.')
    
    def conc_var_DeltaGr(typeRxn, input_, specComp, Ct, range_val, T = 298.15, pH = 7.0, S = 0.0, num = 50, 
                         phase = 'L', fluidType = 'ideal', molality = True, methods = None, marker = 'o', 
                         mec = 'k', mew = 1, mfc = 'w', ms = 8, figsize = (9.0, 6.0), fontsize_label = 12,
                         savePlot = False, printDG0r = False, printDH0r = False, showMessage = True):
        """
        Show the variation of Gibbs free energy for a set of reactions at a specific range of substrate
        and product concentrations. If `savePlot=True`, the plots are saved in `results/#. rxnName` folder.
        
        Parameters
        ----------
        typeRxn : STR
            What reaction(s) type are requested, matching with csv name. E.g.:
                - 'metabolisms': metabolic activities.
        input_ : STR or LIST
            Name(s) of requested compound(s) or reaction(s).
        specComp : (if input_ is reactions; STR or LIST) or (if input_ is compounds; BOOL - True), optional
            Name(s) of compound(s) to calculate specific deltaGr (kJ/mol-compound). The default is False.
        Ct : DICT
            Total concentrations of compounds {'compounds': [concentrations]}.
            All compounds of a reaction with the same number of concentrations.
        range_val : (min_value, max_value)
            Set minimum and maximum concentration values.
        T : FLOAT, optional
            Set of temperature [K]. The default is 298.15 K (standard temperature).
        pH : FLOAT, optional
            Set of pH. The default is 7.0 (neutral pH).
        S : FLOAT
            Salinity [ppt]. The default is None.
        num : INT, optional
            Number of concentration to generate between min_value and max_value. The default is 50.
        phase : STR, optional
            Phase of fluid. The default is 'L'.
        fluidType : STR, optional
            Type of fluid (ideal or non-ideal). The default is ideal.
        molality : BOOL, optional
            Select if activity units are in molality (True) or molarity (False). The default is True.
        methods : DICT, optional
            Method for coefficient activity estimation. The default is None.
                'DH-ext'    - Debye-Hückel equation extended version.
                'SS'        - Setschenow-Shumpe equation.
        marker : STR, optional
            Set the line marker. The default is 'o'.
        mec : STR, optional
            Set the marker edge color. The default is 'k'.
        mew : FLOAT, optional
            Set the marker edge width in points. The default is 1.0.
        mfc : STR, optional
            Set the marker face color. The default is 'w'.
        ms : FLOAT, optional
            Set the marker size in points. The default is 8.0.
        figsize : (FLOAT, FLOAT), optional
            Figure size. (Width, Height) in inches. The default is (9.0, 6.0).
        fontsize_label : FLOAT, optional
            Font size label. The default is 12.
        savePlot : BOOL, optional
            Save resultant plot in `results/` folder. The default is False.
        printDG0r : BOOL, optional
            Print in console the values of standard Gibbs free energy of reactions. The default is False.
        printDH0r : BOOL, optional
            Print in console the values of standard enthalpy of reactions. The default is False.
        showMessage : BOOL, optional
             Boolean to set whether informative messages are displayed in Console. The default is True.

        Returns
        -------
        Plot in Spyder or in `results/` folder.

        """
        if showMessage:
            print('  > Running concentration sensitivity analysis of Gibbs free energy...')
        savePath = 'results/'  
        if range_val[0] > range_val[1]:
            raise ValueError(f'Invalid \'range_val\' ({range_val}). It must be (min_val, max_val).')
        min_val = range_val[0]
        max_val = range_val[1]
        log_min_val = int(np.ceil(np.log10(min_val)))
        log_max_val = int(np.ceil(np.log10(max_val)))
        gShape = np.ones((num, num))
        list_val = {}
        for comp in Ct:
            list_val[comp] = np.logspace(log_min_val, log_max_val, num = num)
        if not isinstance(T, (int, float)) and (isinstance(T, (list, np.ndarray)) and len(T) != 1):
            raise TypeError(f'Temperature (argument `T`) must be an integer or float ({type(T)})')
        if not isinstance(S, (int, float)) and (isinstance(S, (list, np.ndarray)) and len(T) != 1):
            raise TypeError(f'Salinity (argument `S`) must be an integer or float ({type(S)})')
        T = T * gShape
        S = S * gShape
        for c in Ct:
            Ct[c] = Ct[c] * gShape
        remove_comp = ['H2O', 'H+']
        iC = 0
        for idRxn, rxn in enumerate(input_):
            iC += 1
            specCp_ = specComp[idRxn]
            sa_comp, _, _ = Rxn.getRxn(typeRxn, rxn)
            for r in remove_comp:
                try:
                    sa_comp.remove(r)
                except:
                    continue
            comb_comp = combinations(sa_comp, 2)
            for comb in comb_comp:
                fix_comp = sa_comp[:]
                C = Ct.copy()
                # Analyzed compounds
                for idc, c in enumerate(comb):
                    c_prev = c
                    try:
                        list_val[c]
                    except:
                        # pH speciation
                        pH_comp, _, _ = Rxn.getRxnpH(c)
                        if 'H2CO3' in pH_comp:
                            pH_comp = [w.replace('H2CO3', 'CO2') for w in pH_comp]
                        u, counts = np.unique(np.hstack((pH_comp, list(Ct.keys()))), return_counts=True)
                        c = u[counts > 1][0]
                    if idc == 0:
                        C[c] = np.transpose(np.tile(list_val[c], (num, 1)))
                    elif idc == 1:
                        C[c] = np.tile(list_val[c], (num, 1))
                    fix_comp.remove(c_prev)
                # Fixed compounds
                for f in fix_comp:
                    try:
                        Ct[f]
                    except:
                        # pH speciation
                        pH_comp, _, _ = Rxn.getRxnpH(f)
                        if 'H2CO3' in pH_comp:
                            pH_comp = [w.replace('H2CO3', 'CO2') for w in pH_comp]
                        u, counts = np.unique(np.hstack((pH_comp, list(Ct.keys()))), return_counts=True)
                        f = u[counts > 1][0]
                    C[f] = np.tile(Ct[f][0, 0], (num, num))
                # DGr calculation
                DGr, _ = ThSA.getDeltaGr(typeRxn = 'microprony',
                                         input_ = rxn, 
                                         phase = phase, 
                                         specComp = specCp_,
                                         T = T,
                                         pH = pH,
                                         S = S, 
                                         Ct = C,
                                         fluidType = fluidType,
                                         molality = molality,
                                         methods = methods,
                                         printDG0r = printDG0r,
                                         printDH0r = printDH0r)
                # Plotting
                try:
                    list_val[comb[1]]
                except:
                    # pH speciation
                    pH_comp, _, _ = Rxn.getRxnpH(comb[1])
                    if 'H2CO3' in pH_comp:
                        pH_comp = [w.replace('H2CO3', 'CO2') for w in pH_comp]
                    u, counts = np.unique(np.hstack((pH_comp, list(Ct.keys()))), return_counts=True)
                    c = u[counts > 1][0]
                    x = list_val[c]
                else:
                    x = list_val[comb[1]]
                try:
                    list_val[comb[0]]
                except:
                    # pH speciation
                    pH_comp, _, _ = Rxn.getRxnpH(comb[0])
                    if 'H2CO3' in pH_comp:
                        pH_comp = [w.replace('H2CO3', 'CO2') for w in pH_comp]
                    u, counts = np.unique(np.hstack((pH_comp, list(Ct.keys()))), return_counts=True)
                    c = u[counts > 1][0]
                    y = list_val[c]
                else:
                    y = list_val[comb[0]]
                fig, ax = plt.subplots(figsize = figsize)
                F = ax.contourf(x, y, DGr, num, norm=clr.CenteredNorm(), cmap = 'coolwarm_r')
                clb = fig.colorbar(F, format=tkr.FormatStrFormatter('%.2f'))
                eq_line = (np.max(DGr)-np.min(DGr)) * 1e-3
                levels = [-eq_line, eq_line]
                plt.contourf(x, y, DGr, levels=levels, colors='k')
                # yLabel
                try:
                    Ct[comb[0]]
                except:
                    # pH speciation
                    pH_comp, _, _ = Rxn.getRxnpH(comb[0])
                    if 'H2CO3' in pH_comp:
                        pH_comp = [w.replace('H2CO3', 'CO2') for w in pH_comp]
                    u, counts = np.unique(np.hstack((pH_comp, list(Ct.keys()))), return_counts=True)
                    c = u[counts > 1][0]
                else:
                    c = comb[0]
                try:
                    Ct[c]
                except:
                    print(f'Warning: ∆Gr of {rxn} cannot be estimated. Missing concentration of {c}.')
                    continue
                else:
                    y = Ct[c]
                    yLabel = c
                try:
                    Ct[comb[1]]
                except:
                    # pH speciation
                    pH_comp, _, _ = Rxn.getRxnpH(comb[1])
                    if 'H2CO3' in pH_comp:
                        pH_comp = [w.replace('H2CO3', 'CO2') for w in pH_comp]
                    u, counts = np.unique(np.hstack((pH_comp, list(Ct.keys()))), return_counts=True)
                    c = u[counts > 1][0]
                else:
                    c = comb[1]
                try:
                    Ct[c]
                except:
                    print(f'Warning: ∆Gr of {rxn} cannot be estimated. Missing concentration of {c}.')
                    continue
                else:
                    x = Ct[c]
                    xLabel = c
                plt.plot(x, y, marker = marker, mec = mec, mew = mew,
                         mfc = mfc, ms = ms, alpha = 0.8)
                clb.set_label('∆Gr (kJ/mol)', fontsize = fontsize_label)
                clb.ax.axhline(0.0, c='k')
                ax.set_xscale('log')
                ax.set_yscale('log')
                ax.set_xlabel(f'{xLabel} (M)')
                ax.set_ylabel(f'{yLabel} (M)')
                ticks = np.logspace(log_min_val, log_max_val, int(abs(log_max_val-log_min_val))+1)
                ax.set_xticks(ticks)
                ax.set_yticks(ticks)
                xaxis = plt.gca().xaxis
                xaxis.set_minor_locator(MinorSymLogLocator(1e-1))
                yaxis = plt.gca().yaxis
                yaxis.set_minor_locator(MinorSymLogLocator(1e-1))
                plt.title(f'{rxn} | x: {xLabel}; y: {yLabel}')
                if savePlot:
                    path = savePath + str(iC) + '.' + rxn + '/'
                    if not os.path.isdir(path):
                        os.makedirs(path)
                    plt.savefig(f'{path}{rxn}_{xLabel}_{yLabel}', bbox_inches='tight')
                else:
                    plt.show()
                plt.close()
            print(f'    > {iC}. {rxn} done.')
        if showMessage:
            print('  > Done.')
    
    def local_sa_DeltaGr(typeRxn, input_, specComp, list_var, Ct, T_sv = 298.15, pH_sv = 7.0, S_sv = 0.0, fontsize = 11, 
                         rangeType = 'VR', range_ = None, num = 50, num_pH = 10, sensitivity_method = 'sigma-norm',
                         concLog10 = False, phase = 'L', fluidType = 'ideal', molality = True, methods = None, 
                         renameRxn = None, figsize = (12.0, 8.0), cb_limit = False, vmin = None, vmax = None, 
                         cb_fontsize = 12, cb_orientation = 'horizontal', marker = '*', mec = 'k', mew = 0.75, 
                         mfc = 'gold', ms = 8, printDG0r = False, printDH0r = False, showMessage = True):
        """
        Perform the local sensitivity analysis of Gibbs free energy for a set of reactions at a specific
        range of temperature, pH and concentrations of substrates and products.

        Parameters
        ----------
        typeRxn : STR
            What reaction(s) type are requested, matching with csv name. E.g.:
                - 'metabolisms': metabolic activities.
        input_ : STR or LIST
            Name(s) of requested compound(s) or reaction(s).
        specComp : (if input_ is reactions; STR or LIST) or (if input_ is compounds; BOOL - True), optional
            Name(s) of compound(s) to calculate specific deltaGr (kJ/mol-compound). The default is False.
        list_var : LIST of STR
            List of variables. Temperature as 'T', pH as 'pH' and concentrations as 'conc_compoundSymbol'
            (e.g., 'conc_H2').
        Ct : DICT
            Total concentrations of compounds {'compounds': [concentrations]}.
            All compounds of a reaction with the same number of concentrations.
        T_sv : FLOAT, optional
            Set of temperature [K]. The default is 298.15 K (standard temperature).
        pH_sv : FLOAT, optional
            Set of pH. The default is 7.0 (neutral pH).
        S_sv : FLOAT
            Salinity [ppt]. The default is 0.0.
        fontsize : FLOAT, optional
            Set font size. The default is 11.
        rangeType : STR ('DR' or 'VR'), optional
            Type of range. 'DR' means 'defined range' and the user gives the maximum and minimum values of
            each variable. 'VR' measn 'value range' and the range are defined based on original values. 
            The default is 'VR'.
        range_ : DICT, optional
            Range of values or upper and lower order of magnitudes/difference. The default is None.
            If rangeType = 'DR': {'var': [min_val, max_val]}
            If rangeType = 'VR': {'T' or 'pH': [lower_diff, upper_diff]}; {'conc_compound': [lower_oom, upper_oom]}
        num : INT, optional
            Number of temperature and concentration to generate between min_value and max_value. The default is 50.
        num_pH : INT, optional
            Number of pH values to generate between min_value and max_value. The default is 10.
        sensitivity_method : STR ('local', 'sigma-norm'), optional
            Set sensitivity method. The default is 'sigma-norm'.
        concLog10 : BOOL, optional
            Establish whether concentration is analysed using log10(Ci). The default is True.
        phase : STR, optional
            Phase of fluid. The default is 'L'.
        fluidType : STR, optional
            Type of fluid (ideal or non-ideal). The default is ideal.
        molality : BOOL, optional
            Select if activity units are in molality (True) or molarity (False). The default is True.
        methods : DICT, optional
            Method for coefficient activity estimation. The default is None.
                'DH-ext'    - Debye-Hückel equation extended version.
                'SS'        - Setschenow-Shumpe equation.
        renameRxn : None or DICT, optional
            If it's a DICT, change de name of reactions of .csv file in the plot. {'originalName': 'newName'}
            The default is None.
        figsize : (FLOAT, FLOAT), optional
            Figure size. (Width, Height) in inches. The default is (12.0, 8.0).
        cb_limit : BOOL, optional
            Active/inactive limits of colorbar. The default is False.
        cb_vmin : FLOAT or None, optional
            Set minimum value of colorbar. The default is None.
        cb_vmax : FLOAT or None, optional
            Set maximum value of colorbar. The default is None.
        cb_fontsize : FLOAT, optional
            Set size of colorbar font. The default is 12.
        cb_orientation : STR ('vertical', 'horizontal'), optional
            Set orientation of colorbar. The default is 'horizontal'.
        marker : STR, optional
            Set the line marker. The default is '*'.
        mec : STR, optional
            Set the marker edge color. The default is 'k'.
        mew : FLOAT, optional
            Set the marker edge width in points. The default is 0.75.
        mfc : STR, optional
            Set the marker face color. The default is 'gold'.
        ms : FLOAT, optional
            Set the marker size in points. The default is 8.
        printDG0r : BOOL, optional
            Print in console the values of standard Gibbs free energy of reactions. The default is False.
        printDH0r : BOOL, optional
            Print in console the values of standard enthalpy of reactions. The default is False.
        showMessage : BOOL, optional
            Boolean to set whether informative messages are displayed in Console. The default is True.

        Returns
        -------
        Plot in Spyder.

        """
        font = {'size': fontsize}
        plt.rc('font', **font)
        sensitivity_methods = ('local', 'sigma-norm', 'ref-norm', 'var-based', 'pearson')
        if showMessage:
            print('  > Running local sensitivity analysis of Gibbs free energy...')
        if not isinstance(T_sv, (int,float)) or (isinstance(T_sv, (list, np.ndarray)) and len(T_sv) > 1):
            raise ValueError('Temperature (argument \'T_sv\') must be an integer or float.')
        if not isinstance(pH_sv, (int,float)) or (isinstance(pH_sv, (list, np.ndarray)) and len(pH_sv) > 1):
            raise ValueError('pH (argument \'pH_sv\') must be an integer or float.')
        if not isinstance(S_sv, (int,float)) or (isinstance(S_sv, (list, np.ndarray)) and len(S_sv) > 1):
            raise ValueError('Salinity (argument \'S_sv\') must be an integer or float.')
        if vmin and vmax and cb_limit == False:
            cb_limit = True
        if rangeType == 'VR':
            if not range_:
                range_val = {}
                for var in list_var:
                    if var == 'T':
                        range_val[var] = [float(T_sv)-10, float(T_sv)+10]
                    elif var == 'pH':
                        range_val[var] = [max(float(pH_sv)-5, 0.0), min(float(pH_sv)+5, 14.0)]
                    elif 'conc_' in var:
                        indx = var.find('_') + 1
                        comp = var[indx:]
                        range_val[var] = [float(Ct[comp])/100, float(Ct[comp])*100]
            else:
                range_val = {}
                for var in list_var:
                    lower_oom = range_[var][0]
                    upper_oom = range_[var][1]
                    if var == 'T':
                        range_val[var] = [float(T_sv) - lower_oom, float(T_sv) + upper_oom]
                    elif var == 'pH':
                        range_val[var] = [float(pH_sv) - lower_oom, float(pH_sv) + upper_oom]
                    elif 'conc_' in var:
                        indx = var.find('_') + 1
                        comp = var[indx:]
                        range_val[var] = [float(Ct[comp])/lower_oom, float(Ct[comp])*upper_oom]
        elif rangeType == 'DR':
            if not range_:
                raise ValueError('Range values must be defined. {\'var\': [lower_value, upper_value]}')
            else:
                range_val = range_
        res = num
        res_pH = num_pH
        # Matrix initialization
        m_eqch = np.empty((len(input_), len(list_var))) 
        m_diff = np.empty((len(input_), len(list_var)))
        for idRxn, rxn in enumerate(input_):
            print(f'  {idRxn+1}.{rxn}')
            specCp_ = specComp[idRxn]
            C = Ct.copy()
            for c in Ct:
                C[c] = np.array([C[c]])
            DGr_ref, _ = ThSA.getDeltaGr(typeRxn = typeRxn,
                                         input_ = rxn, 
                                         phase = phase, 
                                         specComp = specCp_,
                                         T = [T_sv],
                                         pH = pH_sv,
                                         S = [S_sv], 
                                         Ct = C,
                                         fluidType = fluidType,
                                         methods = methods,
                                         printDG0r = printDG0r,
                                         printDH0r = printDH0r)
            xLabels = []
            for idVar, var in enumerate(list_var):
                if 'conc_' in var:
                    indx = var.find('_') + 1
                    comp = var[indx:]
                    xLabels += [f'[{comp}]']
                else:
                    xLabels += [var]
                if var == 'T':
                    ref_val = T_sv
                    T_min = range_val['T'][0]
                    T_max = range_val['T'][1]
                    T = np.linspace(T_min, T_max, res)
                    list_val = T
                    S = S_sv * np.ones(np.array(list_val).shape)
                elif var == 'pH':
                    ref_val = pH_sv
                    pH_min = range_val['pH'][0]
                    pH_max = range_val['pH'][1]
                    pH = np.linspace(pH_min, pH_max, res_pH)
                    list_val = pH
                    T = T_sv * np.ones(np.array(list_val).shape)
                    S = S_sv * np.ones(np.array(list_val).shape)
                elif 'conc_' in var:
                    C_min = range_val[var][0]
                    C_max = range_val[var][1]
                    if concLog10:
                        list_val = np.logspace(np.log10(C_min), np.log10(C_max), num = res)
                    else:
                        list_val = np.linspace(C_min, C_max, num = res)
                    T = T_sv * np.ones(np.array(list_val).shape)
                    S = S_sv * np.ones(np.array(list_val).shape)
                    indx = var.find('_') + 1
                    comp = var[indx:]
                    ref_val = Ct[comp]
                    # pH speciation of 'ref_val'?
                    rxn_comp, _, _ = Rxn.getRxn(typeRxn, rxn)
                    all_comp = rxn_comp.copy()
                    for c in all_comp:
                        if c == 'CO2': c = 'H2CO3'
                        pH_comp, _, _ = Rxn.getRxnpH(c)
                        if pH_comp and c != 'H+' and c != 'H2O':
                            if 'CO2' in pH_comp:
                                pH_comp = [w.replace('H2CO3', 'CO2') for w in pH_comp]
                            all_comp = np.unique(np.hstack((all_comp, pH_comp)))
                    if not comp in all_comp:
                        m_diff[idRxn, idVar] = np.NaN
                        m_eqch[idRxn, idVar] = np.NaN
                        continue
                if var != 'pH':
                    C = Ct.copy()
                    for c in Ct:
                        if 'conc_' in var:
                            if c == comp:
                                C[c] = list_val
                                continue
                        C[c] = Ct[c] * np.ones(np.array(list_val).shape)
                    # DGr calculation
                    DGr, _ = ThSA.getDeltaGr(typeRxn = typeRxn,
                                             input_ = rxn, 
                                             phase = phase, 
                                             specComp = specCp_,
                                             T = T,
                                             pH = pH_sv,
                                             S = S, 
                                             Ct = C,
                                             fluidType = fluidType,
                                             methods = methods,
                                             printDG0r = printDG0r,
                                             printDH0r = printDH0r)
                else:
                    DGr = []
                    C = Ct.copy()
                    for c in C:
                        C[c] = np.array([C[c]])
                    for idpH, pH_ in enumerate(pH):
                        DGr_pH, _ = ThSA.getDeltaGr(typeRxn = typeRxn,
                                                    input_ = rxn, 
                                                    phase = phase, 
                                                    specComp = specCp_,
                                                    T = [T_sv],
                                                    pH = pH_,
                                                    S = [S_sv], 
                                                    Ct = C,
                                                    fluidType = fluidType,
                                                    methods = methods,
                                                    printDG0r = printDG0r,
                                                    printDH0r = printDH0r)
                        DGr += [DGr_pH]
                if sensitivity_method == 'local':
                    sa_method = 'Local sensitivity analysis'
                    m_diff[idRxn, idVar] = np.squeeze((DGr[-1] - DGr[0]) / (list_val[-1] - list_val[0]))
                elif sensitivity_method == 'sigma-norm':
                    sa_method = 'Sigma-normalized Derivative'
                    sigma_x = np.nanstd(list_val)
                    sigma_y = np.nanstd(DGr)
                    m_diff[idRxn, idVar] = np.squeeze((sigma_x / sigma_y) * ((DGr[-1] - DGr[0]) / (list_val[-1] - list_val[0])))
                elif sensitivity_method == 'ref-norm':
                    sa_method = 'Reference-normalized Derivative'
                    m_diff[idRxn, idVar] = np.squeeze((ref_val / DGr_ref) * ((DGr[-1] - DGr[0]) / (list_val[-1] - list_val[0])))
                elif sensitivity_method == 'var-based':
                    sa_method = 'Variance-normalized Derivative'
                    N = len(DGr)
                    varEDGr_X = np.sum(abs(DGr - DGr_ref)**2) / N
                    varEX = np.sum(abs(list_val - ref_val)**2) / N
                    m_diff[idRxn, idVar] = np.squeeze((varEX / varEDGr_X) * ((DGr[-1] - DGr[0]) / (list_val[-1] - list_val[0])))
                elif sensitivity_method == 'pearson':
                    sa_method = 'Pearson\'s correlation'
                    cov_xy = np.cov(list_val, DGr) # np.cov(a,b) = [[cov(a,a), cov(a,b)][cov(a,b), cov(b,b)]]
                    Sxy = cov_xy[0][1]
                    sigma_x = np.nanstd(list_val)
                    sigma_y = np.nanstd(DGr)
                    m_diff[idRxn, idVar] = np.squeeze((Sxy / (sigma_x * sigma_y)))
                else:
                    raise ValueError(f'Unknown sensitivity method ({sensitivity_method}). Existing sensitivity methods: {sensitivity_methods}')
                if min(DGr) < 0 and max(DGr) > 0:
                    m_eqch[idRxn, idVar] = 1
                else:
                    m_eqch[idRxn, idVar] = 0
                if 'conc_' in var:
                    print(f'    ·[{comp}] done.')
                else:
                    print(f'    ·{var} done.')
        # #-Plotting local_sa_DGr
        ThSA._plot_mesh_sa(input_, list_var, m_diff, m_eqch, None, sa_method , cb_limit, cb_orientation, cb_fontsize, vmin, vmax, figsize, 
                           renameRxn, marker, mec, mew, mfc, ms)
        if showMessage:
            print('  > Done.')
    
    def sobol_indices_DeltaGr(typeRxn, input_, specComp, list_var, Ct = 1.0, T_sv = 298.15, pH_sv = 7.0, S_sv = 0.0, 
                              phase = 'L', fluidType = 'ideal', methods = None, molality = True, rangeType = 'VR', 
                              range_ = None, num = 64, rng = None, sobol_simpl = 'saltelli', check_negatives = True, 
                              stat_sign_method = 'kendall', plotMode = False, renameRxn = None, cb_limit = False, 
                              cb_orientation = 'horizontal', cb_fontsize = 12, vmin = None, vmax = None, figsize = (12.0, 8.0),
                               marker = '*', mec = 'k', mew = 0.75, mfc = 'gold', ms = 8, showMessage = True):
        """
        Perform the global sensitivity analysis (variance-based sensitivity analysis or Sobol' indices) of Gibbs free 
        energy for a set of reactions at a specific range of temperature, pH and concentrations of substrates and products.

        Parameters
        ----------
        typeRxn : STR
            What reaction(s) type are requested, matching with csv name. E.g.:
                - 'metabolisms': metabolic activities.
        input_ : STR or LIST
            Name(s) of requested compound(s) or reaction(s).
        specComp : (if input_ is reactions; STR or LIST) or (if input_ is compounds; BOOL - True), optional
            Name(s) of compound(s) to calculate specific deltaGr (kJ/mol-compound).
        list_var : LIST of STR
            List of variables. Temperature as 'T', pH as 'pH' and concentrations as 'conc_compoundSymbol'
            (e.g., 'conc_H2').
        Ct : DICT
            Total concentrations of compounds {'compounds': [concentrations]}.
            All compounds of a reaction with the same number of concentrations. The default is 1.0.
        T_sv : FLOAT, optional
            Set of temperature [K]. The default is 298.15 K (standard temperature).
        pH_sv : FLOAT, optional
            Set of pH. The default is 7.0 (neutral pH).
        S_sv : FLOAT
            Salinity [ppt]. The default is 0.0.
        phase : STR, optional
            Phase of fluid. The default is 'L'.
        fluidType : STR, optional
            Type of fluid (ideal or non-ideal). The default is ideal.
        methods : DICT, optional
            Method for coefficient activity estimation. The default is None.
                'DH-ext'    - Debye-Hückel equation extended version.
                'SS'        - Setschenow-Shumpe equation.
        molality : BOOL, optional
            Select if activity units are in molality (True) or molarity (False). The default is True.
        rangeType : STR ('DR' or 'VR'), optional
            Type of range. 'DR' means 'defined range' and the user gives the maximum and minimum values of
            each variable. 'VR' measn 'value range' and the range are defined based on original values. 
            The default is 'VR'. The default is 'VR'.
        range_ : DICT, optional
            Range of values or upper and lower order of magnitudes/difference. The default is None.
            If rangeType = 'DR': {'var': [min_val, max_val]}
            If rangeType = 'VR': {'T' or 'pH': [lower_diff, upper_diff]}; {'conc_compound': [lower_oom, upper_oom]}
        num : INT, optional
            Number of temperature and concentration to generate between min_value and max_value. The default is 64.
        rng : INT, optional
            Pseudorandom number generator state (rng seed). The default is None.
            When rng is None, a new numpy.random.Generator is created using entropy from the operating system.
        sobol_simpl : STR, optional
            Estimator method to calculate Sobol' indices. The default is 'saltelli'.
                'saltelli'    - Estimators from Saltelli et al. (2010), doi: 10.1016/j.cpc.2009.09.018
        check_negatives : BOOL, optional
            Boolean to set whether negative Sobol' indices are displayed in Console (if any). The default is True.
        stat_sign_method : STR, optional
            Set statistical method to compute correlation between variable and DGr and to obtain statistical sign. The default is 'kendall'.
        plotMode : BOOL, optional
            Boolean to set whether color mesh plot is created. The default is False.
        showMessage : BOOL, optional
            Boolean to set whether informative messages are displayed in Console. The default is True.

        Returns
        -------
        si : np.ndarray
            Matrix of first order indices. Shape: (reactions, variables).
        st : np.ndarray
            Matrix of total indices. Shape: (reactions, variables).
        stat_sign : np.ndarray
            Matrix of statistical sign. Shape: (reactions, variables).
                If positive (+1), DGr and variable are positively correlated (both increase or decrease together).
                If negative (-1), DGr and variable are negatively correlated (one increase and the other decrease or vice versa).
                If zero (0), DGr and variable are not correlated.
        eqch : np.ndarray
            Matrix of equilbirum change (endergonic <-> exergonic). Shape: (reactions, variables).
                If one (1), an equilibrium change is observed in the range of values of variables.
                If zero (0), no equilibrium change is observed in the range of values of variables.
        Plot in Spyder if 'plotMode = True'.
        
        """
        #-Private functions
        def _N_power_of_two(num, maxN = 15):
            power_of_two = np.array([int(2**n) for n in list(np.linspace(0, maxN, maxN+1))])
            indx = np.argmin(abs(power_of_two - num))
            return power_of_two[indx]
            
        def _u_l_bounds(list_var, rangeType, range_, sv):
            Ct = sv['Ct'].copy()
            T_sv = sv['T']
            pH_sv = sv['pH']
            S_sv = sv['S']
            if rangeType:
                l_bounds = np.nan * np.ones(len(list_var))
                u_bounds = np.nan * np.ones(len(list_var))
                if rangeType == 'VR':
                    if not range_:
                        for idVar, var in enumerate(list_var):
                            if var == 'T':
                                l_bounds[idVar] = float(T_sv) - 10
                                u_bounds[idVar] = float(T_sv) + 10
                            elif var == 'pH':
                                l_bounds[idVar] = max(float(pH_sv)-5, 0.0)
                                u_bounds[idVar] = min(float(pH_sv)+5, 14.0)
                            elif var == 'S':
                                l_bounds[idVar] = max(float(S_sv)-10, 0.0)
                                u_bounds[idVar] = float(S_sv) + 10
                            elif 'conc_' in var:
                                indx = var.find('_') + 1
                                comp = var[indx:]
                                try:
                                    l_bounds[idVar] = float(Ct[comp]) / 100
                                    u_bounds[idVar] = float(Ct[comp]) * 100
                                except:
                                    raise ValueError(f'Compound {comp} was not included in Ct dictionary:'+' Ct = {\'compound\': concentration}.')
                    else:
                        for idVar, var in enumerate(list_var):
                            lower_oom = range_[var][0]
                            upper_oom = range_[var][1]
                            if var == 'T':
                                l_bounds[idVar] = float(T_sv) - lower_oom
                                u_bounds[idVar] = float(T_sv) + upper_oom
                            elif var == 'pH':
                                l_bounds[idVar] = max(float(pH_sv)-lower_oom, 0.0)
                                u_bounds[idVar] = min(float(pH_sv)+upper_oom, 14.0)
                            elif 'conc_' in var:
                                indx = var.find('_') + 1
                                comp = var[indx:]
                                try:
                                    l_bounds[idVar] = float(Ct[comp]) / lower_oom
                                    u_bounds[idVar] = float(Ct[comp]) * upper_oom
                                except:
                                    raise ValueError(f'Compound {comp} was not included in Ct dictionary:'+' Ct = {\'compound\': concentration}.')
                elif rangeType == 'DR':
                    if not range_:
                        raise ValueError('Range values must be defined with \'range_\' argument. {\'variable\': [lower_value, upper_value]}.')
                    else:
                        for idVar, var in enumerate(list_var):
                            l_bounds[idVar] = range_[var][0]
                            u_bounds[idVar] = range_[var][1]
                else: raise ValueError(f'Unknown rangeType ({rangeType}). Existing rangeType: \'VR\' (Value Range) and \'DR\' (Defined Range).')
            else:
                l_bounds = None
                u_bounds = None
            return l_bounds, u_bounds
        
        def _sample_A_B_AB(N, k, l_bounds = None, u_bounds = None, rng = None):
            A_B = qmc.Sobol(d = 2*k, seed = rng, bits = 64).random(N)
            if (l_bounds is not None) and (u_bounds is not None):
                if ~np.isnan(l_bounds).any() and ~np.isnan(u_bounds).any():
                    l_bounds = np.log10(l_bounds)
                    u_bounds = np.log10(u_bounds)
                    l_bounds = np.tile(l_bounds, 2)
                    u_bounds = np.tile(u_bounds, 2)
                    A_B = qmc.scale(A_B, l_bounds, u_bounds)
                    A_B = 10 ** A_B
                else:
                    warnings.warn(f'NaN value(s) present in \'l_bounds\' ({l_bounds}) and/or \'u_bounds\' ({u_bounds}).')
            A = A_B[:, 0:k]
            B = A_B[:, k:]
            AB = np.tile(A, (k, 1, 1))
            for i in np.arange(k):
                AB[i, :, i] = B[:, i]
            return A, B, AB
        
        def _sobol_DGr(typeRxn, input_, specComp, phase, fluidType, methods, molality, sv, list_var, sample, AB = False):
            if AB:
                d, k, d = sample.shape
                sample = sample.reshape(-1, sample.shape[-1])
            warnMessage = False
            non_DGr_parameters = []
            N, _ = sample.shape
            C = sv['Ct'].copy()
            T_ = sv['T']
            pH_ = sv['pH']
            S_ = sv['S']
            f = []
            for n in range(N):
                for idVar, var in enumerate(list_var):
                    if var == 'T':
                        T_ = sample[n, idVar]
                    elif var == 'pH':
                        pH_ = sample[n, idVar]
                    elif var == 'S':
                        S_ = sample[n, idVar]
                    elif var.startswith('conc_'):
                        indx = var.find('_') + 1
                        comp = var[indx:]
                        try:
                            C[comp] = sample[n, idVar]
                        except:
                            raise ValueError(f'Compound {comp} was not included in Ct dictionary:'+' Ct = {\'compound\': concentration}.')
                    else:
                        warnMessage = True
                        if not var in non_DGr_parameters:
                            non_DGr_parameters += var
                DGr, _ = ThSA.getDeltaGr(typeRxn = typeRxn,
                                         input_ = rxn,
                                         specComp = specComp, 
                                         phase = phase,
                                         T = [T_],
                                         pH = pH_,
                                         S = [S_],
                                         Ct = C,
                                         fluidType = fluidType,
                                         methods = methods,
                                         molality = molality)
                f += [np.squeeze(DGr)]
            f = np.array(f).T
            if warnMessage:
                warnings.warn(f'Some variables given in \'list_var\' are not used for DGr calculations: {non_DGr_parameters}.')
            if AB:
                f = f.reshape((-1, k))
            else:
                f = np.squeeze(f)
            return f
        
        def _indices(f_A, f_B, f_AB, method):
            valid_methods = {'saltelli'}
            I, _ = f_AB.shape
            first_order = []
            total_order = []
            var = np.nanvar([f_A, f_B], axis = (0, -1))
            for i in range(I):
                fi_AB = f_AB[i, :]
                method = method.lower()
                if not method in valid_methods:
                    raise NameError(f'Invalid method ({method}) to compute sobol\' indices. Valid methods: {valid_methods}.')
                if method == 'saltelli':
                    first_order += [np.nanmean(f_B * (fi_AB - f_A)) / var]
                    total_order += [0.5 * np.nanmean((f_A - fi_AB) ** 2) / var]
            return np.array(first_order), np.array(total_order)
        
        def _stat_sign(typeRxn, input_, specComp, phase, fluidType, methods, molality, sv, sample, list_var, stat_sign_method):
            valid_methods = {'kendall', 'spearman', 'pearson'}
            if not stat_sign_method in valid_methods:
                raise NameError(f'Invalid method ({stat_sign_method}) to compute sobol\' indices. Valid methods: {valid_methods}.')
            N, _ = sample.shape
            Ct = sv['Ct'].copy()
            T_sv = sv['T']
            pH_sv = sv['pH']
            S_sv = sv['S']
            # Matrix initialization
            eqch = np.nan * np.ones((len(input_), len(list_var))) 
            stat_sign = np.nan * np.ones((len(input_), len(list_var)))
            for idRxn, rxn in enumerate(input_):
                specC_ = specComp[idRxn]
                for idVar, var in enumerate(list_var):
                    var_values = sample[:, idVar]
                    shape_var = var_values.shape
                    if var == 'T':
                        T = var_values
                        pH = pH_sv
                        S = S_sv * np.ones(shape_var)
                        C = Ct.copy()
                        for c in C:
                            C[c] = C[c] * np.ones(shape_var)
                    elif var == 'pH':
                        T = T_sv
                        pH = var_values
                        S = S_sv
                        C = Ct.copy()
                        for c in C:
                            C[c] = np.array(C[c])
                    elif var == 'S':
                        T = T_sv * np.ones(shape_var)
                        pH = pH_sv
                        S = sample[:, idVar]
                        C = Ct.copy()
                        for c in C:
                            C[c] = C[c] * np.ones(shape_var)
                    elif var.startswith('conc_'):
                        indx = var.find('_') + 1
                        comp = var[indx:]
                        rxn_comp, _, _ = Rxn.getRxn(typeRxn, rxn)
                        all_comp = rxn_comp.copy()
                        for c in all_comp:
                            if c == 'CO2': c = 'H2CO3'
                            pH_comp, _, _ = Rxn.getRxnpH(c)
                            if pH_comp and c != 'H+' and c != 'H2O':
                                if 'CO2' in pH_comp:
                                    pH_comp = [w.replace('H2CO3', 'CO2') for w in pH_comp]
                                all_comp = np.unique(np.hstack((all_comp, pH_comp)))
                        if not comp in all_comp:
                            stat_sign[idRxn, idVar] = np.nan
                            eqch[idRxn, idVar] = np.nan
                            continue
                        C = Ct.copy()
                        for c in C:
                            if c == comp:
                                C[comp] = var_values
                            else:
                                C[c] = C[c] * np.ones(shape_var)
                    if var != 'pH':
                        DGr, _ = ThSA.getDeltaGr(typeRxn = typeRxn,
                                                 input_ = rxn, 
                                                 phase = phase, 
                                                 specComp = specC_,
                                                 T = T,
                                                 pH = pH,
                                                 S = S, 
                                                 Ct = C,
                                                 fluidType = fluidType,
                                                 methods = methods)
                    else:
                        DGr = np.empty(shape_var)
                        for idpH, pH_ in enumerate(pH):
                            DGr_pH, _ = ThSA.getDeltaGr(typeRxn = typeRxn,
                                                        input_ = rxn, 
                                                        phase = phase, 
                                                        specComp = specC_,
                                                        T = [T],
                                                        pH = pH_,
                                                        S = [S], 
                                                        Ct = C,
                                                        fluidType = fluidType,
                                                        methods = methods)
                            DGr[idpH] = np.squeeze(DGr_pH)
                    # Thermodynamic change analysis (endergonic <-> exergonic)
                    if np.nanmin(DGr) < 0 and np.nanmax(DGr) > 0:
                        eqch[idRxn, idVar] = 1
                    else:
                        eqch[idRxn, idVar] = 0
                    # Statistical sign from DGr = f(var)
                    if stat_sign_method == 'kendall':
                        stat, _ = stats.kendalltau(var_values, DGr)
                    elif stat_sign_method == 'spearman':
                        stat, _ = stats.spearmanr(var_values, DGr)
                    elif stat_sign_method == 'pearson':
                        stat, _ = stats.pearsonr(var_values, DGr)
                    stat_sign[idRxn, idVar] = np.squeeze(np.sign(stat))
            return stat_sign, eqch
        
        if showMessage:
            print('  > Running variance-based sensitivity analysis (Sobol\' method) of Gibbs free energy...')
        #-Checking arguments
        ThSA._check_arguments(('list_var', 'T_sv', 'pH_sv', 'S_sv'), 
                              (list_var, T_sv, pH_sv, S_sv),
                              (list, (int, float), (int, float), (int, float)))
        sv = {'Ct': Ct.copy(),
              'T': T_sv,
              'pH': pH_sv,
              'S': S_sv}
        #-Definiion of `l_bounds` and `u_bounds` based on rangeType ('VR' or 'DR') and range_. Shape: (k,)
        k = len(list_var)
        l_bounds, u_bounds = _u_l_bounds(list_var, rangeType, range_, sv)
        #-The balance properties of Sobol' points require N to be a power of 2.
        N = _N_power_of_two(num, maxN = 15)
        # Creation of sampling matrices (N, k)
        A, B, AB = _sample_A_B_AB(N, k, l_bounds, u_bounds, rng)
        # Variable initialization
        si = np.zeros((len(input_), len(list_var)))
        st = np.zeros((len(input_), len(list_var)))
        if check_negatives:
            rxns_negative_first_order = []
            negatives_first_order = []
            rxns_negative_total_order = []
            negatives_total_order = []
        for idRxn, rxn in enumerate(input_):
            print(f'  {idRxn+1}.{rxn}')
            specC_ = specComp[idRxn]
            start_time = time.time()
            #-Calculation of DGr with samples A, B and AB: f(A), f(B), f(AB)
            f_A = _sobol_DGr(typeRxn, input_, specC_, phase, fluidType, methods, molality, sv, list_var, A, AB = False)
            print('    ·Sample matrix A done (%s seconds).' % (time.time() - start_time))
            start_time = time.time()
            f_B = _sobol_DGr(typeRxn, input_, specC_, phase, fluidType, methods, molality, sv, list_var, B, AB = False)
            print('    ·Sample matrix B done (%s seconds).' % (time.time() - start_time))
            start_time = time.time()
            f_AB = _sobol_DGr(typeRxn, input_, specC_, phase, fluidType, methods, molality, sv, list_var, AB, AB = True)
            print('    ·Sample matrix AB done (%s seconds).' % (time.time() - start_time))
            #-Normalization by mean - empirically centered function (for sake of stability)
            mean = np.nanmean([f_A, f_B], axis = (0, -1))
            f_A -= mean
            f_B -= mean
            f_AB -= mean
            #-Compute indices
            s_i, s_t = _indices(f_A = f_A, f_B = f_B, f_AB = f_AB, method = sobol_simpl)
            if check_negatives: 
                if (s_i < 0).any():
                    rxns_negative_first_order += [rxn]
                    negatives_first_order += [np.nanmin(s_i)]
                if (s_t < 0).any():
                    rxns_negative_total_order += [rxn]
                    negatives_total_order += [np.nanmin(s_t)]
            si[idRxn, :] = s_i
            st[idRxn, :] = s_t
        print('')
        if check_negatives:
            if rxns_negative_first_order:
                zip_negative_first_order = zip(rxns_negative_first_order, np.round(negatives_first_order, decimals = 5))
                print(f'[!] Reactions with negative first order indices [\'reaction\', min(index)]: {list(zip_negative_first_order)}.')
                print('')
            if rxns_negative_total_order:
                zip_negative_total_order = zip(rxns_negative_total_order, np.round(negatives_total_order, decimals = 5))
                print(f'[!] Reactions with negative total order indices [\'reaction\', min(index)]: {list(zip_negative_total_order)}.')
                print('')
        # Thermodynamic change analysis (endergonic <-> exergonic) & Statistical sign from DGr = f(var)
        print('  Executing thermodynamic change analysis & correlation analysis between DGr and variables...')
        stat_sign, eqch = _stat_sign(typeRxn, input_, specComp, phase, fluidType, methods, molality, sv, A, list_var, stat_sign_method)
        print('  Done.')
        if plotMode:
            ThSA._plot_mesh_sa(input_, list_var, st, eqch, stat_sign, "[±] Sobol' indices" , cb_limit, cb_orientation, cb_fontsize, vmin, vmax, figsize, 
                               renameRxn, marker, mec, mew, mfc, ms)
        else:
            return si, st, stat_sign, eqch
    
    def _plot_mesh_sa(input_, list_var, si, m_eqch, stat_sign, sa_method, cb_limit, cb_orientation, cb_fontsize, vmin, vmax, figsize, 
                      renameRxn, marker, mec, mew, mfc, ms, fontFamily = 'Arial'):
        plt.rcParams["font.family"] = fontFamily
        fig, ax = plt.subplots(figsize = figsize)
        xLabels = []
        for idVar, var in enumerate(list_var):
            if 'conc_' in var:
                indx = var.find('_') + 1
                comp = var[indx:]
                xLabels += [f'[{comp}]']
            else:
                xLabels += [var]
        yLabels = input_.copy()
        if renameRxn:
            for rxn_rename in renameRxn:
                try:
                    ind = yLabels.index(rxn_rename)
                except:
                    continue
                else:
                    yLabels[ind] = renameRxn[rxn_rename] 
        x = np.arange(-0.5, len(list_var), 1)
        y = np.arange(-0.5, len(yLabels), 1)
        if isinstance(stat_sign, np.ndarray):
            z = np.where(stat_sign != 0, stat_sign * si, si)
        else:
            z = si
        #-DEBUGGING-#
        print(f'min(z): {np.nanmin(z)}')
        print(f'max(z): {np.nanmax(z)}')
        #-----------#
        if cb_limit:
            pc = ax.pcolormesh(x, y, z, edgecolor = 'k', snap = True, vmin = vmin, vmax = vmax,
                               cmap = 'coolwarm_r')
            if (np.nanmin(z) < vmin) and (np.nanmax(z) > vmax):
                extend_ = 'both'
            elif (np.nanmin(z) < vmin) and not (np.nanmax(z) > vmax):
                extend_ = 'min'
            elif (np.nanmax(z) > vmax) and not (np.nanmin(z) < vmin):
                extend_ = 'max'
            else:
                extend_ = 'neither'
            clb = fig.colorbar(pc, extend=extend_, orientation = cb_orientation)
        else:
            pc = ax.pcolormesh(x, y, z, edgecolor = 'k', snap = True,
                               norm = clr.CenteredNorm(), cmap = 'coolwarm_r')
            clb = fig.colorbar(pc, orientation = cb_orientation)
        a = f'Sᵢ ({sa_method})'
        if cb_orientation == 'vertical':
            clb.set_label(a, fontsize = cb_fontsize) # , y = 0.53, labelpad = 9
        else:
            clb.set_label(a, fontsize = cb_fontsize)
        # Equilibrium change (endergonic <-> exergonic)
        eqch = np.argwhere(m_eqch == 1)
        for idx_eqch in eqch:
            plt.plot(idx_eqch[1], idx_eqch[0], marker = marker, mec = mec, mew = mew,
                     mfc = mfc, ms = ms)
        ax.xaxis.set_ticks_position('top')
        ax.xaxis.set_label_position('top')
        ax.set_xticks(np.arange(0, len(list_var), 1))
        ax.set_xticklabels(xLabels, rotation = 90)
        ax.set_yticks(np.arange(0, len(yLabels), 1))
        ax.set_yticklabels(yLabels)
        ax.set_facecolor('dimgrey') # dimgrey
        ax.invert_yaxis()
        fig.tight_layout()
        plt.gca().set_aspect('equal', adjustable='box')  # Ensures equal aspect ratio
        plt.show()
    
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