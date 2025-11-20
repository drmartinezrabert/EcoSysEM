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
    
    def getPcat(paramDB, typeKin, typeMetabo, reaction, specComp, Ct, 
                T = 298.15, pH = 7., S = None, phase = 'L', sample = 'All',
                fluidType = 'ideal', molality = 'True', methods = 'None',
                solvent = 'H2O', asm = 'stoich', DGsynth = 9.54E-11, Rs = None, DGr = None):
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
        specComp : (if input_ is reaction; STR) or (if input_ is compound; BOOL - True), optional
            Name of compound to calculate specific deltaGr (kJ/mol-compound). The default is False.
        Ct : DICT
            Total concentrations of compound {'compound': [concentrations]}.
        T : INT or FLOAT or LIST or numpy.ndarray
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
        Rs : np.ndarray, LIST, INT, FLOAT, optional
            Cell-specific uptake rate(s) in [moleD/cell.h]. It can be either given or kept to default (None) to be computed inside of the function.
                NB : Input units must be respected.
                     If given, Rs must have the same shape as DGr.
        DGr : np.ndarray, LIST, INT, FLOAT, optional
            Non-standard Gibbs free energy in [kJ/moleD]. It can be either given or kept to default (None) to be computed inside of the function.
                NB : Input units must be respected.
                     If given, DGr must have the same shape as Rs.
                   
        Returns
        -------
        Pcat : numpy.ndarray
        	Catabolic cell-specific power [fW/cell]: energy flux produced by the cell,
            using environmental resources or internal reservoirs.
            
        """
        # Get cell-specific uptake rate (Rs)
        if Rs is None or np.all(Rs == 0):
            # Set requested arguments for KR.getRs (no Rs given in arguments)
            Rs_args = {'typeKin' : typeKin,
                       'paramDB' : paramDB,
                       'reactions' : reaction,
                       'T' : T,
                       'pH' : pH,
                       'Ct' : Ct, 
                       'sample' : sample
                       }
            #compute Rs dict
            Rs = (KR.getRs(**Rs_args)[0])[reaction]
            #extract Rs values and compute mean of sample combinations
            _rs = np.array(list(val for val in Rs.values()))
            _Rs = np.nanmean(_rs, axis = 0) / 3600      #[moleD/(cell.s)]
        else:
            #use Rs argument
            if not isinstance(Rs,np.ndarray): 
                if isinstance(Rs, (int, float, list)): Rs = np.array(Rs)
                else: raise TypeError(f'Rs must be a non-zero float, int, list or array. current type: {type(Rs)})')
            _Rs = Rs / 3600  #[moleD/(cell.s)]
        # Get non-standard Gibbs free energy (DGr)
        if DGr is None or np.all(DGr == 0):
            # Set requested arguments for ThSA.getDeltaGr (no DGr given in arguments)
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
            _DGr = np.squeeze(ThSA.getDeltaGr(**DGr_args)[0]) * 1000  #[J/moleD]
        else: 
            #use DGr argument
            if not isinstance(DGr, np.ndarray):
                if isinstance(DGr,(float, int, list)): DGr = np.array(DGr)
                else: raise TypeError(f'DGr must be a non-zero float, int, list or array. current type: {type(DGr)}')
            _DGr = DGr * 1000  #[J/moleD]
        #check shape of _DGr and _Rs arrays
        if not _DGr.shape == _Rs.shape :
            raise IndexError(f'Arrays of different shape cannot be used as operands (_Rs:{_Rs.shape} ; _DGr:{_DGr.shape}).')
        # Compute Pcat
        Pcat = -(_Rs * _DGr * 1e15)
        return Pcat      #[fW/cell]

    def getPana(paramDB, typeKin, typeMetabo, reaction, specComp, Ct, 
                T = 298.15, pH = 7., S = None, phase = 'L', sample = 'All',
                fluidType = 'ideal', molality = 'True', methods = 'None',
                solvent = 'H2O', asm = 'stoich', DGsynth = 9.54E-11, Rs = None, DGr = None):
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
        specComp : (if input_ is reaction; STR) or (if input_ is compound; BOOL - True), optional
            Name of compound to calculate specific deltaGr (kJ/mol-compound). The default is False.
        Ct : DICT
            Total concentrations of compound {'compound': [concentrations]}.
        T : INT or FLOAT or LIST or numpy.ndarray
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
        Rs : np.ndarray, LIST, INT, FLOAT, optional
            Cell-specific uptake rate(s) in [moleD/cell.h]. It can be either given or kept to default (None) to be computed inside of the function.
                NB : Input units must be respected.
                     If given, Rs must have the same shape as DGr.
        DGr : np.ndarray, LIST, INT, FLOAT, optional
            Non-standard Gibbs free energy in [kJ/moleD]. It can be either given or kept to default (None) to be computed inside of the function.
                NB : Input units must be respected.
                     If given, DGr must have the same shape as Rs.
        
        Returns
        -------
        Pana : numpy.ndarray
        	Anabolic cell-specific power [fW/cell]: energy flux associated with the synthesis of
            cellular components.
        """
        # Get cell-specific uptake rate (Rs)
        if Rs is None or np.all(Rs == 0):
            # Set requested arguments for KR.getRs (no Rs given in arguments)
            Rs_args = {'typeKin' : typeKin,
                       'paramDB' : paramDB,
                       'reactions' : reaction,
                       'T' : T,
                       'pH' : pH,
                       'Ct' : Ct, 
                       'sample' : sample
                       }
            #compute Rs dict
            Rs = (KR.getRs(**Rs_args)[0])[reaction]
            #extract Rs values and compute mean of sample combinations
            _rs = np.array(list(val for val in Rs.values()))
            _Rs = np.nanmean(_rs, axis = 0) / 3600      #[moleD/(cell.s)]
        else:
            #use Rs argument
            if not isinstance(Rs,np.ndarray): 
                if isinstance(Rs, (int, float, list)): Rs = np.array(Rs)
                else: raise TypeError(f'Rs must be a non-zero float, int, list or array. current type: {type(Rs)})')
            _Rs = Rs / 3600  #[moleD/(cell.s)]
        # Get non-standard Gibbs free energy (DGr)
        if DGr is None or np.all(DGr == 0):
            # Set requested arguments for ThSA.getDeltaGr (no DGr given in arguments)
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
            _DGr = np.squeeze(ThSA.getDeltaGr(**DGr_args)[0]) * 1000  #[J/moleD]
        else: 
            #use DGr argument
            if not isinstance(DGr, np.ndarray):
                if isinstance(DGr,(float, int, list)): DGr = np.array(DGr)
                else: raise TypeError(f'DGr must be a non-zero float, int, list or array. current type: {type(DGr)}')
            _DGr = DGr * 1000  #[J/moleD]
        #check shape of _DGr and _Rs arrays
        if not _DGr.shape == _Rs.shape :
            raise IndexError(f'Arrays of different shape cannot be used as operands (_Rs:{_Rs.shape} ; _DGr:{_DGr.shape}).')
        # compute cell-growth yield and Pana
        Yx = -(_DGr * (0.5/1.04e-10))   # [cell/moleD]
        Pana = Yx * _Rs * DGsynth * 1e15
        return Pana     #[fW/cell]
    
    def getPmg(T = 298.15):
        """
        Function to compute the growth-based maintenance power.
        
        Parameters
        ----------
        T : INT or FLOAT or LIST or numpy.ndarray
            Set of temperature [K]. The default is 298.15 K (standard temperature).
        
        Returns
        -------
        Pmg : numpy.ndarray
        	Growth-based maintenance power [fW/cell]: energy flux that microbes use that does not
            result in growth while they are growing (Pirt et al., 1965).
        """
        # initialise variables
        R = CSP.Rj
        T0 = CSP.T0
        # compute Pmg
        Pmg = np.array(8.96 * np.exp((-6.40e4 / R) * ((1 / T) - (1 / T0))))
        return Pmg      #[fW/cell]
    
    def getPm0(T = 298.15):
        """
        Function to compute the basal maintenance power.
        
        Parameters
        ----------
        T : INT or FLOAT or LIST or numpy.ndarray
            Set of temperature [K]. The default is 298.15 K (standard temperature).
        
        Returns
        -------
        Pm0 : numpy.ndarray
        	Basal maintenance power [fW/cell]: energy flux associated with the minimal set of functions
            required to sustain a basal functional state (Hoehler et al., 2013).
        """
        # initialise variables
        R = CSP.Rj
        T0 = CSP.T0
        # compute Pm0
        Pm0 = np.array(1.70e-3 * np.exp((-9.07e4/R) * ((1/T) - (1/T0))))
        return Pm0      #[fW/cell]
    
    def getPs(T = 298.15):
        """
        Function to compute the survival power.
        
        Parameters
        ----------
        T : INT or FLOAT or LIST or numpy.ndarray
            Set of temperature [K]. The default is 298.15 K (standard temperature).
        
        Returns
        -------
        Ps : numpy.ndarray
        	Survival power [fW/cell]: minimal energy flux for preservation of membrane integrity and
            key macromolecules (e.g., enzymes), as well as other maintenance costs, such
            as maintaining energized membranes or the conservation of catabolic energy. 
        """
        # initialise variables
        R = CSP.Rj
        T0 = CSP.T0
        # compute Ps
        Ps = np.array(2.67e-6 * np.exp((-5.68e4 / R) * ((1 / T) - (1 / T0))))
        return Ps   #[fW/cell]
    
        
    def getAllCSP(paramDB, typeKin, typeMetabo, reaction, specComp, Ct,
                  T = 298.15, pH = 7., S = None, phase = 'L', sample = 'All',
                  fluidType  = 'ideal', molality = 'True', methods = 'None',
                  solvent = 'H2O', asm = 'stoich', DGsynth = 9.54E-11,
                  Rs = None, DGr = None):
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
        specComp : (if input_ is reaction; STR) or (if input_ is compound; BOOL - True), optional
            Name of compound to calculate specific deltaGr (kJ/mol-compound). The default is False.
        Ct : DICT
            Total concentrations of compound {'compound': [concentrations]}.
        T : INT or FLOAT or LIST or numpy.ndarray
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
        Rs : np.ndarray, LIST, INT, FLOAT, optional
            Cell-specific uptake rate(s) in [moleD/cell.h]. It can be either given or kept to default (None) to be computed inside of the function.
                NB : Input units must be respected.
                     If given, Rs must have the same shape as DGr.
        DGr : np.ndarray, LIST, INT, FLOAT, optional
            Non-standard Gibbs free energy in [kJ/moleD]. It can be either given or kept to default (None) to be computed inside of the function.
                NB : Input units must be respected.
                     If given, DGr must have the same shape as Rs.
        
        Returns
        -------
        CSP_dict : DICT of numpy.ndarray
            Contains all cell specific power series + Pcell in fW/cell:
                - 'Pcat' : Catabolic cell-specific power: energy flux produced by the cell, using environmental resources or internal reservoirs.
                - 'Pana' : Anabolic cell-specific power: energy flux associated with the synthesis of cellular components.
                - 'Pmg' : Growth-based maintenance power: energy flux that microbes use that does not result in growth while they are growing (Pirt et al., 1965).
                - 'Pm0' : Basal maintenance power: energy flux associated with the minimal set of functions required to sustain a basal functional state (Hoehler et al., 2013).
                - 'Ps' : Survival power: minimal energy flux for preservation of membrane integrity and key macromolecules (e.g., enzymes), as well as other maintenance costs, such as maintaining energized membranes or the conservation of catabolic energy.
                - 'Pcell' : Growth power: energy flux of a growing cell (sum of Pana & Pmg).
        """
        # Compute CSP through existing methods                    
        Pcat = CSP.getPcat(paramDB, typeKin, typeMetabo, reaction, specComp, Ct, 
                    T, pH, S, phase, sample, fluidType, molality, methods,
                    solvent, asm, DGsynth, Rs, DGr)
        Pana = CSP.getPana(paramDB, typeKin, typeMetabo, reaction, specComp, Ct, 
                    T, pH, S, phase, sample, fluidType, molality, methods,
                    solvent, asm, DGsynth, Rs, DGr)
        Pmg = CSP.getPmg(T)
        Pm0 = CSP.getPm0(T)
        Ps = CSP.getPs(T)
        Pcell = Pana + Pmg
        # Create DataFrame of CSP results
        CSP_dict = {'Pcat': Pcat, 'Pana': Pana, 'Pmg': Pmg, 'Pm0': Pm0, 'Ps': Ps, 'Pcell': Pcell}
        return CSP_dict  # [fW/cell]