# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 10:00:16 2025

@author: zLemaire
"""

from thermodynamics import ThSA
from reactions import KinRates as KR

import numpy as np

class CSP:
    
    Rj = 8.314      #universal gas constant [J/mol.K]
    T0 = 298.15     #standard temperature [K]
    
    def getPcat(typeMetabo, reaction, pH = 7., S = None, Ct = 1., T = 298.15, phase = 'L'):
        """
        Function to compute the catabolic cell-specific power.
        
        Parameters
        ----------
        
        typeMetabo : STR
            Requested metabolism type, matching with csv name. E.g.:
                - 'metabolisms' : aerobic metabolisms
                - 'AnMetabolisms' : anaerobic metabolisms
        reaction : STR
            Requested reaction name. E.g.:
                -'COOB' : carbon monoxide oxidation
                -'HOB' : hydrogen oxidation
        pH : INT or FLOAT, optional
            Set of pH. The default is 7.0 (neutral pH).
        S : FLOAT, LIST or np.array
            Salinity [ppt]. The default is None.
        Ct : DICT
            Total concentrations of compounds {'compounds': [concentrations]}.
            All compounds of a reaction with the same number of concentrations.
        T : FLOAT or LIST
            Set of temperature [K]. The default is 298.15 K (standard temperature).
        phase : STR
            Phase in which reaction(s) ocurr. 'G' - Gas, 'L' - Liquid. The default is 'L'.
                   
        Returns
        -------
        Pcat : FLOAT or LIST
        	Catabolic cell-specific power: energy flux produced by the cell, using environmental resources
            or internal reservoirs.
            
        """
        
        # Get non-standard Gibbs free energy and cell-specific uptake rate     
        DGr = ThSA.getDeltaGr(typeMetabo, reaction, phase, T = T)
        Rs = KR.getRs('MM-Arrhenius', 'ArrhCor_AtmMicr', reaction, Ct)
        # compute Pcat
        Pcat = abs(Rs * DGr) * 1e18
        return Pcat #[fJ/cell.h]
    
    def getPana(typeMetabo, reaction, pH = 7., S = None, Ct = 1., T = 298.15, phase = 'L', DGsynth = 9.54E-11):
        """
        Function to compute the anabolic cell-specific power.
        
        Parameters
        ----------
        
        typeMetabo : STR
            Requested metabolism type, matching with csv name. E.g.:
                - 'metabolisms' : aerobic metabolisms
                - 'AnMetabolisms' : anaerobic metabolisms
        reaction : STR
            Requested reaction name. E.g.:
                -'COOB' : carbon monoxide oxidation
                -'HOB' : hydrogen oxidation
        pH : INT or FLOAT, optional
            Set of pH. The default is 7.0 (neutral pH).
        S : FLOAT, LIST or np.array
            Salinity [ppt]. The default is None.
        Ct : DICT
            Total concentrations of compounds {'compounds': [concentrations]}.
            All compounds of a reaction with the same number of concentrations.
        T : FLOAT or LIST
            Set of temperature [K]. The default is 298.15 K (standard temperature).
        phase : STR
            Phase in which reaction(s) ocurr. 'G' - Gas, 'L' - Liquid. The default is 'L'.
        DGsynth : FLOAT
            Energy necessary to synthesize a cell [J/cell], by default: 9.54E-11
        
        Returns
        -------
        Pana : FLOAT or LIST
        	Anabolic cell-specific power: energy flux associated with the synthesis of
            cellular components.
            
        """
        
        # Get non-standard Gibbs free energy and cell-specific uptake rate
        DGr = ThSA.getDeltaGr(typeMetabo, reaction, phase, T = T)
        Rs = KR.getRs('MM-Arrhenius', 'ArrhCor_AtmMicr', reaction, Ct)
        # compute Pana
        Pana = DGr * Rs * DGsynth * 1e15
        return Pana #[fJ/cell.h]
    
    def getPmg(T = 298.15):
        """
        Function to compute the growth-based maintenance power.
        
        Parameters
        ----------
        
        T : FLOAT or LIST
            Set of temperature [K]. The default is 298.15 K (standard temperature).
        
        Returns
        -------
        Pmg : FLOAT or LIST
        	Growth-based maintenance power: energy flux that microbes use that does not
            result in growth while they are growing (Pirt et al., 1965).
            
        """
        # initialise variables
        R = CSP.Rj
        T0 = CSP.T0
        # compute Pmg
        Pmg = (8.96 * np.exp((-6.40e4 / R) * ((1 / T) - (1 / T0)))) * 1e15
        return Pmg #[fJ/cell.h]
    
    def getPm0(T = 298.15):
        """
        Function to compute the basal maintenance power.
        
        Parameters
        ----------
        
        T : FLOAT or LIST
            Set of temperature [K]. The default is 298.15 K (standard temperature).
        
        Returns
        -------
        Pm0 : FLOAT or LIST
        	Basal maintenance power: energy flux associated with the minimal set of functions
            required to sustain a basal functional state (Hoehler et al., 2013).
            
        """
        # initialise variables
        R = CSP.Rj
        T0 = CSP.T0
        # compute Pm0
        Pm0 = (1.70e-3 * np.exp((-9.07e4/R) * ((1/T) - (1/T0)))) * 1e15
        return Pm0 #[fJ/cell.h]
    
    def getPs(T = 298.15):
        """
        Function to compute the survival power.
        
        Parameters
        ----------
        
        T : FLOAT or LIST
            Set of temperature [K]. The default is 298.15 K (standard temperature).
        
        Returns
        -------
        Ps : FLOAT or LIST
        	Survival power: minimal energy flux for preservation of membrane integrity and
            key macromolecules (e.g., enzymes), as well as other maintenance costs, such
            as maintaining energized membranes or the conservation of catabolic energy. 
        """
        
        # initialise variables
        R = CSP.Rj
        T0 = CSP.T0
        # compute Ps
        Ps = (2.67e-6 * np.exp((-5.68e4 / R) * ((1 / T) - (1 / T0)))) * 1e15
        return Ps #[fJ/cell.h]
    
        
    def getAllCSP(typeMetabo, reaction, Ct, T = 298.15, phase = 'L', DGsynth = 9.54E-11):
        """
        Function to compute every cell-specific powers.
        
        Parameters
        ----------
        
        typeMetabo : STR
            Requested metabolism type, matching with csv name. E.g.:
                - 'metabolisms' : aerobic metabolisms
                - 'AnMetabolisms' : anaerobic metabolisms
        reaction : STR
            Requested reaction name. E.g.:
                -'COOB' : carbon monoxide oxidation
                -'HOB' : hydrogen oxidation
        Ct : DICT
            Total concentrations of compounds {'compounds': [concentrations]}.
            All compounds of a reaction with the same number of concentrations.
        T : FLOAT or LIST
            Set of temperature [K]. The default is 298.15 K (standard temperature).
        phase : STR
            Phase in which reaction(s) occur. 'G' - Gas, 'L' - Liquid. The default is 'L'.
        DGsynth : FLOAT
            Energy necessary to synthesize a cell [J/cell], by default: 9.54E-11
        
        Returns
        -------
        Pcat : FLOAT or LIST
        	Catabolic cell-specific power: energy flux produced by the cell, using environmental resources
            or internal reservoirs.
        Pana : FLOAT or LIST
        	Anabolic cell-specific power: energy flux associated with the synthesis of
            cellular components.
        Pmg : FLOAT or LIST
        	Growth-based maintenance power: energy flux that microbes use that does not
            result in growth while they are growing (Pirt et al., 1965).
        Pm0 : FLOAT or LIST
        	Basal maintenance power: energy flux associated with the minimal set of functions
            required to sustain a basal functional state (Hoehler et al., 2013).
        Ps : FLOAT or LIST
        	Survival power: minimal energy flux for preservation of membrane integrity and
            key macromolecules (e.g., enzymes), as well as other maintenance costs, such
            as maintaining energized membranes or the conservation of catabolic energy.
            
        """
        
        Pcat = CSP.getPcat(typeMetabo, reaction, T, Ct, phase)
        Pana = CSP.getPana(typeMetabo, reaction, T, Ct, phase, DGsynth)
        Pmg = CSP.getPmg(T)
        Pm0 = CSP.getPm0(T)
        Ps = CSP.getPs(T)
        Pcell = Pana + Pmg
        
        
        return Pcat, Pana, Pmg, Pm0, Ps, Pcell  #[fJ/cell.h]