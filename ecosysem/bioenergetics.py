# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 10:00:16 2025

@author: zLemaire
"""

from thermodynamics import ThSA
from reactions import KinRates as KR

import numpy as np
import pandas as pd

class CSP:
    
    Rj = 8.314      #universal gas constant [J/mol.K]
    T0 = 298.15     #standard temperature [K]
    
    def getPcat(paramDB, typeKin, typeMetabo, reaction, specComp, Ct, 
                T = 298.15, pH = 7., S = None, phase = 'L', sample = 'All',
                fluidType = 'ideal', molality = 'True', methods = 'None',
                solvent = 'H2O', asm = 'stoich', DGsynth = 9.54E-11):
        """   
        Function to compute the catabolic cell-specific power.
        
        Parameters
        ----------
        paramDB : STR or LIST
            Name of parameter database, matching with csv name, E.g. 'ArrhCor_AtmMicr'
        typeKin : STR
            Type of kinetic equations 
                MM - 'Michaelis-Menten equation'.
                MM-Arrhenius - 'Michaelis-Menten-Arrhenius equation'
        typeMetabo : STR
            Requested metabolism type, matching with csv name. E.g.:
                - 'metabolisms' : aerobic metabolisms
                - 'AnMetabolisms' : anaerobic metabolisms
        reaction : STR
            Requested reaction name. E.g.:
                -'COOB' : carbon monoxide oxidation
                -'HOB' : hydrogen oxidation
        specComp : (if input_ is reactions; STR or LIST) or (if input_ is compounds; BOOL - True), optional
            Name(s) of compound(s) to calculate specific deltaGr (kJ/mol-compound). The default is False.
        Ct : DICT
            Total concentrations of compounds {'compounds': [concentrations]}.
            All compounds of a reaction with the same number of concentrations.
        T : FLOAT or LIST
            Set of temperature [K]. The default is 298.15 K (standard temperature).
        pH : INT or FLOAT, optional
            Set of pH. The default is 7.0 (neutral pH).
        S : FLOAT, LIST or np.array, optional
            Salinity [ppt]. The default is None.
        phase : STR, optional
            Phase in which reaction(s) ocurr. 'G' - Gas, 'L' - Liquid. The default is 'L'.
        sample : STR or LIST, optional
            Requested samples (rows of `paramDB.csv`). The default is 'All'.
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
        asm : STR, optional
            Assumption when products are not present in the environment.
            The default is 'stoich' (stoichiometric concentrations).
        DGsynth : FLOAT
            Energy necessary to synthesize a cell [J/cell], by default: 9.54E-11.         
                   
        Returns
        -------
        Pcat : pandas.core.series.Series
        	Catabolic cell-specific power: energy flux produced by the cell,
            using environmental resources or internal reservoirs.
            
        """
        # Set requested arguments for ThSA.getDeltaGr & KR.getRs
        DGr_args = {'typeRxn' : typeMetabo,
                    'input_' : reaction,
                    'phase' : phase,
                    'specComp' : specComp,
                    'T' : T,
                    'pH' : pH,
                    'S' : S,
                    'Ct' : Ct,
                    'fluidType' : fluidType,
                    'molality' : molality,
                    'methods' : methods,
                    'solvent' : solvent,
                    'asm' : asm, 
                    'warnings' : False,
                    'printDG0r' : False,
                    'printDH0r' : False
                    }       
        Rs_args = {'typeKin' : typeKin,
                   'paramDB' : paramDB,
                   'reactions' : reaction,
                   'T' : T,
                   'pH' : pH,
                   'Ct' : Ct, 
                   'sample' : sample
                   }
        # Get non-standard Gibbs free energy and cell-specific uptake rate
        DGr = ThSA.getDeltaGr(**DGr_args)[0] * 1000  #[J/moleD]
        Rs = pd.DataFrame((KR.getRs(**Rs_args)[0])[reaction]).mean(axis=1) / 3600   #[moleD/(cell.s)]
        # Compute Pcat
        Pcat = -(Rs * DGr * 1e15)
        return Pcat      #[fW/cell]

    def getPana(paramDB, typeKin, typeMetabo, reaction, specComp, Ct, 
                T = 298.15, pH = 7., S = None, phase = 'L', sample = 'All',
                fluidType = 'ideal', molality = 'True', methods = 'None',
                solvent = 'H2O', asm = 'stoich', DGsynth = 9.54E-11):
        """
        Function to compute the anabolic cell-specific power.
        
        Parameters
        ----------
        paramDB : STR or LIST
            Name of parameter database, matching with csv name, E.g. 'ArrhCor_AtmMicr'
        typeKin : STR
            Type of kinetic equations 
                MM - 'Michaelis-Menten equation'.
                MM-Arrhenius - 'Michaelis-Menten-Arrhenius equation'
        typeMetabo : STR
            Requested metabolism type, matching with csv name. E.g.:
                - 'metabolisms' : aerobic metabolisms
                - 'AnMetabolisms' : anaerobic metabolisms
        reaction : STR
            Requested reaction name. E.g.:
                -'COOB' : carbon monoxide oxidation
                -'HOB' : hydrogen oxidation
        specComp : (if input_ is reactions; STR or LIST) or (if input_ is compounds; BOOL - True), optional
            Name(s) of compound(s) to calculate specific deltaGr (kJ/mol-compound). The default is False.
        Ct : DICT
            Total concentrations of compounds {'compounds': [concentrations]}.
            All compounds of a reaction with the same number of concentrations.
        T : FLOAT or LIST
            Set of temperature [K]. The default is 298.15 K (standard temperature).
        pH : INT or FLOAT, optional
            Set of pH. The default is 7.0 (neutral pH).
        S : FLOAT, LIST or np.array, optional
            Salinity [ppt]. The default is None.
        phase : STR, optional
            Phase in which reaction(s) ocurr. 'G' - Gas, 'L' - Liquid. The default is 'L'.
        sample : STR or LIST, optional
            Requested samples (rows of `paramDB.csv`). The default is 'All'.
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
        asm : STR, optional
            Assumption when products are not present in the environment.
            The default is 'stoich' (stoichiometric concentrations).
        DGsynth : FLOAT
            Energy necessary to synthesize a cell [J/cell], by default: 9.54E-11.
        
        Returns
        -------
        Pana : pandas.core.series.Series
        	Anabolic cell-specific power: energy flux associated with the synthesis of
            cellular components.
        """
        # Set requested arguments for ThSA.getDeltaGr & KR.getRs
        DGr_args = {'typeRxn' : typeMetabo,
                    'input_' : reaction,
                    'phase' : phase,
                    'specComp' : specComp,
                    'T' : T,
                    'pH' : pH,
                    'S' : S,
                    'Ct' : Ct,
                    'fluidType' : fluidType,
                    'molality' : molality, 
                    'methods' : methods,
                    'solvent' : solvent,
                    'asm' : asm, 
                    'warnings' : False,
                    'printDG0r' : False,
                    'printDH0r' : False
                    }       
        Rs_args = {'typeKin' : typeKin,
                   'paramDB' : paramDB,     
                   'reactions' : reaction,
                   'T' : T,
                   'pH' : pH,
                   'Ct' : Ct, 
                   'sample' : sample
                   }
        # Get non-standard Gibbs free energy and cell-specific uptake rate
        DGr = ThSA.getDeltaGr(**DGr_args)[0] * 1000 #[J/moleD]
        Rs = pd.DataFrame((KR.getRs(**Rs_args)[0])[reaction]).mean(axis=1) / 3600   #[moleD/cell/s]      
        # compute cell-growth yield and Pana
        Yx = -(DGr * (0.5/1.04e-10))   # [cell/moleD]
        Pana = Yx * Rs * DGsynth * 1e15
        return Pana     #[fW/cell]
    
    def getPmg(T = 298.15):
        """
        Function to compute the growth-based maintenance power.
        
        Parameters
        ----------
        T : FLOAT or LIST
            Set of temperature [K]. The default is 298.15 K (standard temperature).
        
        Returns
        -------
        Pmg : pandas.core.series.Series
        	Growth-based maintenance power: energy flux that microbes use that does not
            result in growth while they are growing (Pirt et al., 1965).
        """
        # initialise variables
        R = CSP.Rj
        T0 = CSP.T0
        # compute Pmg
        Pmg = pd.Series(8.96 * np.exp((-6.40e4 / R) * ((1 / T) - (1 / T0))))
        return Pmg      #[fW/cell]
    
    def getPm0(T = 298.15):
        """
        Function to compute the basal maintenance power.
        
        Parameters
        ----------
        T : FLOAT or LIST
            Set of temperature [K]. The default is 298.15 K (standard temperature).
        
        Returns
        -------
        Pm0 : pandas.core.series.Series
        	Basal maintenance power: energy flux associated with the minimal set of functions
            required to sustain a basal functional state (Hoehler et al., 2013).
        """
        # initialise variables
        R = CSP.Rj
        T0 = CSP.T0
        # compute Pm0
        Pm0 = pd.Series(1.70e-3 * np.exp((-9.07e4/R) * ((1/T) - (1/T0))))
        return Pm0      #[fW/cell]
    
    def getPs(T = 298.15):
        """
        Function to compute the survival power.
        
        Parameters
        ----------
        T : FLOAT or LIST
            Set of temperature [K]. The default is 298.15 K (standard temperature).
        
        Returns
        -------
        Ps : pandas.core.series.Series
        	Survival power: minimal energy flux for preservation of membrane integrity and
            key macromolecules (e.g., enzymes), as well as other maintenance costs, such
            as maintaining energized membranes or the conservation of catabolic energy. 
        """
        # initialise variables
        R = CSP.Rj
        T0 = CSP.T0
        # compute Ps
        Ps = pd.Series(2.67e-6 * np.exp((-5.68e4 / R) * ((1 / T) - (1 / T0))))
        return Ps   #[fW/cell]
    
        
    def getAllCSP(paramDB, typeKin, typeMetabo, reaction, specComp, Ct, 
                T = 298.15, pH = 7., S = None, phase = 'L', sample = 'All',
                fluidType = 'ideal', molality = 'True', methods = 'None',
                solvent = 'H2O', asm = 'stoich', DGsynth = 9.54E-11):
        """
        Function to compute every cell-specific powers.
        
        Parameters
        ----------
        paramDB : STR or LIST
            Name of parameter database, matching with csv name, E.g. 'ArrhCor_AtmMicr'
        typeKin : STR
            Type of kinetic equations 
                MM - 'Michaelis-Menten equation'.
                MM-Arrhenius - 'Michaelis-Menten-Arrhenius equation'
        typeMetabo : STR
            Requested metabolism type, matching with csv name. E.g.:
                - 'metabolisms' : aerobic metabolisms
                - 'AnMetabolisms' : anaerobic metabolisms
        reaction : STR
            Requested reaction name. E.g.:
                -'COOB' : carbon monoxide oxidation
                -'HOB' : hydrogen oxidation
        specComp : (if input_ is reactions; STR or LIST) or (if input_ is compounds; BOOL - True), optional
            Name(s) of compound(s) to calculate specific deltaGr (kJ/mol-compound). The default is False.
        Ct : DICT
            Total concentrations of compounds {'compounds': [concentrations]}.
            All compounds of a reaction with the same number of concentrations.
        T : FLOAT or LIST
            Set of temperature [K]. The default is 298.15 K (standard temperature).
        pH : INT or FLOAT, optional
            Set of pH. The default is 7.0 (neutral pH).
        S : FLOAT, LIST or np.array, optional
            Salinity [ppt]. The default is None.
        phase : STR, optional
            Phase in which reaction(s) ocurr. 'G' - Gas, 'L' - Liquid. The default is 'L'.
        sample : STR or LIST, optional
            Requested samples (rows of `paramDB.csv`). The default is 'All'.
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
        asm : STR, optional
            Assumption when products are not present in the environment.
            The default is 'stoich' (stoichiometric concentrations).
        DGsynth : FLOAT
            Energy necessary to synthesize a cell [J/cell], by default: 9.54E-11.
        
        Returns
        -------
        dfCSP : pandas.core.frame.DataFrame
            Contains all cell specific power series :
                - 'Pcat' : Catabolic cell-specific power: energy flux produced by the cell, using environmental resources or internal reservoirs.
                - 'Pana' : Anabolic cell-specific power: energy flux associated with the synthesis of cellular components.
                - 'Pmg' : Growth-based maintenance power: energy flux that microbes use that does not result in growth while they are growing (Pirt et al., 1965).
                - 'Pm0' : Basal maintenance power: energy flux associated with the minimal set of functions required to sustain a basal functional state (Hoehler et al., 2013).
                - 'Ps' : Survival power: minimal energy flux for preservation of membrane integrity and key macromolecules (e.g., enzymes), as well as other maintenance costs, such as maintaining energized membranes or the conservation of catabolic energy.
        """
        
        Pcat = CSP.getPcat(paramDB, typeKin, typeMetabo, reaction, specComp, Ct, 
                    T, pH, S, phase, sample, fluidType, molality, methods,
                    solvent, asm, DGsynth)
        Pana = CSP.getPana(paramDB, typeKin, typeMetabo, reaction, specComp, Ct, 
                    T, pH, S, phase, sample, fluidType, molality, methods,
                    solvent, asm, DGsynth)
        Pmg = CSP.getPmg(T)
        Pm0 = CSP.getPm0(T)
        Ps = CSP.getPs(T)
        Pcell = Pana + Pmg
        dfCSP = pd.concat([Pcat, Pana, Pmg, Pm0, Ps, Pcell], axis = 1)
        dfCSP.columns = ['Pcat', 'Pana', 'Pmg', 'Pm0', 'Ps', 'Pcell']
        return dfCSP  # [fW/cell]



 
    