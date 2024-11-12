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
import re
from scipy.stats import kendalltau as kTau
import os.path

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
    
    def getKeq(compounds, mRxn, temperature, phase):
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
        phase : STR
            Fluid phase. Depending on typeParam.

        Returns
        -------
        Keq : np.array
            Equilibrium constants.
    
        """
        # Constants
        R = 0.0083144598   # Universal gas constant [kJ/mol/K]
        deltaG0f, notNaN = ThP.getThP('deltaG0f', compounds, phase)
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
        if isinstance(temperature, int): temperature = float(temperature)
        if isinstance(temperature, float): temperature = [temperature]
        if isinstance(pH, int): pH = float(pH)
        if isinstance(pH, float): pH = [pH]
        if isinstance(Ct, int): Ct = float(Ct)
        if isinstance(Ct, float): Ct = [Ct]
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
            rComp, mRxn, infoRxn = Rxn.getRxnpH(iCompound)
            if not rComp:
                return Ct
            c_rComp = iCompound in rComp
            if not c_rComp:
                return Ct
            reqSp = rComp.index(iCompound) - 1 # Requested chemical species
            for idTemperature, iTemperature in enumerate(temperature):
                # Acid dissociation constants (Ka)
                transDict = {'First deP': 0, 'Second deP': 1, 'Third deP': 2}
                intRxn = [transDict[letter] for letter in infoRxn]
                Kai = ThP.getKeq(rComp, mRxn, iTemperature, 'L')
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
        textT = float(f'{temperature:.2f}')
        for iCompound in compounds:
            Spec = ThEq.pHSpeciation(iCompound, pH, temperature, Ct, True)
            nFrac = (Spec.T / Ct) * 100 # Molar fraction [%]
            # Selection of involved chemical species
            c_nFrac = np.nonzero(np.sum(nFrac, axis = 0))
            nFrac = nFrac[:, c_nFrac[0]]
            # Get name of chemical species
            nCompounds = Rxn.getRxnpH(iCompound)[0][1:]
            # Simplification hydration/dehydration equil.: [(CO2)aq] >>>> [H2CO3]
            nCompounds = [nC.replace('H2CO3', 'H2CO3 (CO2)') for nC in nCompounds]
            text_ = f'Chemical (or ion) speciation at {textT}K.'
            # Plotting
            fig, ax = plt.subplots()
            ax.plot(pH, nFrac)
            ax.set_ylabel('Molar fraction (%)')
            ax.set_xlabel('pH')
            ax.set_xticks(np.arange(pH[0], pH[-1]+1, 1))
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
    
    def getDeltaGr(typeRxn, input_, phase, specComp = False, T = 298.15, Ct = 1.0, pH = 7.0, asm = 'stoich', warnings = False):
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
        specComp : (if input_ is reactions; STR or LIST) or (if input_ is compounds; BOOL - True)
            Name(s) of compound(s) to calculate specific deltaGr (kJ/mol-compound). Default: False.
        T : FLOAT or LIST
            Set of temperature [K]. By default: standard temperature (298.15 K).
        Ct : DICT
            Total concentrations of compounds {'compounds': [concentrations]}.
            All compounds of a reaction with the same number of concentrations.
        pH : FLOAT or LIST
            Set of pH. By default: neutral pH (pH: 7.0)
        asm : STRING
            Assumption when products are not present in the environment.
            By default: 'stoich' - stoichiometric concentrations.
        warnings : BOOL
            Display function warnings. Default: False.

        Returns
        -------
        DGr : np.array
            Non-standard Gibbs free energy values.
            Shape: (reactions)x(temperature)x(pH)x(total concentration)..
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
        if isinstance(pH, int): pH = float(pH)
        if isinstance(pH, float): pH = [pH]
        # Variable lenghts
        nT = len(T)
        npH = len(pH)
        # Check compound concentrations
        if isinstance(Ct, float):
            nCt = 1
        else:
            npCt = np.array(list(Ct.items()), dtype = object).T
            lenCt = [len(i) for i in npCt[1, :]]
            nCt = list(set(lenCt))
            cLen = len(nCt) == 1
            if not cLen:
                print('!EcoSysEM.Error: All compounds must have same number of concentrations.')
                sys.exit()
            if not specComp and asm == 'stoich':
                print('!EcoSysEM.Error: You have to give the main compound (with `specComp` argument) to estimate product concentrations.')
                sys.exit()
            nCt = nCt[0]
        # Get reactions
        rComp, mRxn, infoRxn = Rxn.getRxn(typeRxn, input_, warnings)
        nRxn = infoRxn.size
        # Initialise variables
        Ts = 298.15                                                             # Standard temperature [K]
        R = 0.0083144598                                                        # Universal gas constant [kJ/mol/K]
        DGr = np.empty([nRxn, nT, npH, nCt])
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
                    print('!EcoSysEM.Error: Different number of reactions and specific compounds was found.')
                    sys.exit()
        # Select reactions w/ requested compounds as substrates (if input_ = compounds) (How ?)
        for idRxn, iRxn in enumerate(infoRxn):
            # iVariables definition
            rDGr = np.empty([nT, npH, nCt])                                     # Initialise auxiliar result matrix (rDGr)
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
            # Calculate DeltaH0r
            deltaH0f = ThP.getThP('deltaH0f', i_rComp, phase)[0]
            deltaH0r = ThP.getDeltaH0r(deltaH0f, i_mRxn)                        # kJ
            deltaH0r = deltaH0r / vSelected                                     # kJ/mol-i (if specDGr)
            for idT, iT in enumerate(T):
                # Temperature influence (Gibbs–Helmholtz relationship)
                deltaGTr = deltaG0r * (iT / Ts) + deltaH0r * ((Ts - iT) / Ts)
                for idpH, ipH in enumerate(pH):
                    if isinstance(Ct, dict):
                        # Calculate reaction quotient (Qr)
                        Qr = 0
                        uComp = np.array(list(Ct.keys()))                       # Compounds given by user (Ct)
                        findSpecComp = np.argwhere(uComp == i_specComp).squeeze()
                        for idComp, iComp in enumerate(i_rComp):
                            findComp = np.argwhere(uComp == iComp)
                            vi = i_mRxn[idComp] / vSelected
                            if findComp.size == 0:
                                if vi < 0:
                                    # iComp is a substrate and its concentration was not given.
                                    print(f'!EcoSysEM.Error: Concentration for {iComp} not found.')
                                else:
                                    if iComp == 'H+':
                                        iConc = 10**(-ipH) * np.ones((1, nCt))
                                    elif iComp == 'H2O':                # It is assumed a water activity of 1.0, as a pure liquid (it can be less).
                                        iConc = 1.0 * np.ones((1, nCt))
                                    else: # [P] is calculated based on stoichiometry.
                                        if asm == 'stoich':
                                            iConc = (vi) * Ct[uComp[findSpecComp]]
                            else:
                                iConc = Ct[iComp]
                            # pH speciation
                            if iComp != 'H+' and iComp != 'H2O' and iComp != 'OH-':
                                # Only pH speciation in liquid
                                if phase == 'L':
                                    if iComp == 'CO2': iComp = 'H2CO3'      # Simplification hydration/dehydration equil.: [(CO2)aq] >>>> [H2CO3]
                                    iConc = ThEq.pHSpeciation(iComp, ipH, iT, iConc)
                            # Calculation of reaction quotient (Qr)
                            Qr += vi * np.log(iConc)
                        # Concentration influence (Nernst relationship)
                        deltaGr = deltaGTr + R * iT * Qr
                        rDGr[idT, idpH, :] = deltaGr
                    else:
                        rDGr[idT, idpH, :] = deltaGTr
            DGr[idRxn, :, :, :] = rDGr
        DGr = np.squeeze(DGr)
        return DGr, infoRxn
    
    def exportDeltaGr(modeExport, T, pH, phase, typeRxn, input_, specComp = False, Ct = 1.0, Ct_associated = None, asm = 'stoich', warnings = False):
        """
        Export DeltaGr in function of pH, temperature and/or compound
        concentrations.

        Parameters
        ----------
        modeExport : STR
            How DeltaGr is exported: 'plot': plot in Spyder; 'Excel': write in Excel document.
        T : LIST or np.array
            Set of temperatures [K].
        pH : LIST or np.array
            Set of pH values.
        phase: STR
            Phase in which reaction(s) ocurr. 'G' - Gas, 'L' - Liquid.
        typeRxn : STR
            What reaction(s) type are requested, matching with csv name. E.g.:
                - 'metabolisms': metabolic activities.
        input_ : STR or LIST
            Name(s) of requested compound(s) or reaction(s).
        specComp : (if input_ is reactions; STR or LIST) or (if input_ is compounds; BOOL - True)
           Name(s) of compound(s) to calculate specific deltaGr (kJ/mol-compound). Default: False.
        Ct : DICT
            Total concentrations of compounds {'compounds': [concentrations]}.
            All compounds of a reaction with the same number of concentrations.
            The default is 1.0.
        Ct_associated : STR ('x', 'y' or 'xy')
            Set if concentration is associated with coordinate x, y or both ('xy'). The default is None.
            If None and Ct != 1.0, Ct is z automatically.
        asm : STRING
            Assumption when products are not present in the environment.
            By default: 'stoich' - stoichiometric concentrations.
        warnings : BOOL
            Display function warnings. Default: False.

        Returns
        -------
        Plot(s).

        """
        # Default values
        oneT = False
        onepH = False
        # Checking arguments and variable definition
        if not isinstance(modeExport, str): modeExport = str(modeExport)
        if modeExport == 'Excel':
            path = 'results/'
            nameDocument = input(' > Name of result document: ')
            fullPathSave = path + nameDocument + '.xlsx'
            if os.path.isfile(fullPathSave):
                val = input(' > Do you want to overwrite `' + nameDocument + '.xlsx`? [Y/N]: ')
                if val == 'Y' or val == 'y':       
                    os.remove(fullPathSave)
        elif modeExport != 'plot' and modeExport != 'Excel':
            print("`modeExport` argument must be 'plot' (for plotting in Spyder) or 'Excel' (for writing in Excel document).")
            sys.exit()
        if not isinstance(typeRxn, str): typeRxn = str(typeRxn)
        if not isinstance(input_, list): input_ = [input_]
        if T is None: oneT = True; T = [298.15] # Standard conditions (K)
        if isinstance(T, int): T = float(T)
        if isinstance(T, float): oneT = True; T = [T]
        if pH is None: onepH = True; pH = [7.0]  # Standard conditions
        if isinstance(pH, int): pH = float(pH)
        if isinstance(pH, float): onepH = True; pH = [pH]
        # Variable lenghts
        nT = len(T)
        npH = len(pH)
        if phase != 'G' and phase != 'L':
            print('!EcoSysEM.Error: `phase` argument must be "G" (gas) or "L" (liquid)')
            sys.exit()
        if isinstance(Ct, float):
            nCt = 1
        else:
            npCt = np.array(list(Ct.items()), dtype = object).T
            lenCt = [len(i) for i in npCt[1, :]]
            nCt = list(set(lenCt))
            cLen = len(nCt) == 1
            if not cLen:
                print('!EcoSysEM.Error: All compounds must have same number of concentrations.')
                sys.exit()
            nCt = nCt[0]
        # Exporting DGr
        if Ct_associated:
            if Ct_associated == 'x' and nCt != npH:
                print('!EcoSysEM.Error: Total concentrations and pH must have the same size.')
                sys.exit()
            elif Ct_associated == 'y' and nCt != nT:
                print('!EcoSysEM.Error: Total concentrations and temperature must have the same size.')
                sys.exit()
            elif Ct_associated == 'xy' and len(list({nT, npH, nCt})) != 1:
                print('!EcoSysEM.Error: Total concentrations, temperature and pH must have the same size.')
                sys.exit()
            elif Ct_associated != 'x' and Ct_associated != 'y' and Ct_associated != 'xy':
                print('!EcoSysEM.Error: Ct_associated must be "x", "y" or "xy".')
                sys.exit()
            DGr, infoRxn = ThSA.getDeltaGr(typeRxn, input_, phase, specComp, T, Ct, pH, asm, warnings)
            nDimDGr = DGr.ndim
            for idRxn, iRxn in enumerate(infoRxn):
                # Selection of DGr for plotting
                if Ct_associated == 'xy':
                    index_Ct = list(range(nCt))
                    if nDimDGr == 4:
                        DGr_export = DGr[idRxn, index_Ct, index_Ct, index_Ct]
                    else:
                        DGr_export = DGr[index_Ct, index_Ct, index_Ct]
                    text_ = 'Concentrations (in mol/L) associated with pH and T.'
                elif Ct_associated == 'x':
                    index_pH = list(range(npH))
                    if not oneT:
                        if nDimDGr == 4:
                            DGr_export = DGr[idRxn, :, index_pH, index_pH].T
                        elif nDimDGr == 3:
                            DGr_export = DGr[:, index_pH, index_pH]
                        text_ = 'Concentrations (in mol/L) associated with pH.'
                    else:
                        if nDimDGr == 3:
                            DGr_export = DGr[idRxn, index_pH, index_pH]
                        elif nDimDGr == 2:
                            DGr_export = DGr[index_pH, index_pH]
                        Ct_associated = 'OnlypH'
                        text_ = f'Concentrations (in mol/L) associated with pH (T = {str(T[0])}K).'
                elif Ct_associated == 'y':
                    index_T = list(range(nT))
                    if not onepH:
                        if nDimDGr == 4:
                            DGr_export = DGr[idRxn, index_T, :, index_T]
                        elif nDimDGr == 3:
                            DGr_export = DGr[index_T, :, index_T]
                        text_ = 'Concentrations (in mol/L) associated with T.'
                    else:
                        if nDimDGr == 3:
                            DGr_export = DGr[idRxn, index_T, index_T]
                        elif nDimDGr == 2:
                            DGr_export = DGr[index_T, index_T]
                        Ct_associated = 'OnlyT'
                        text_ = f'Concentrations (in mol/L) associated with T (pH = {str(pH[0])}).'
                else:
                    print("`Ct_associated` argument must be 'x' (pH, only one pH or pH = None), 'y' (temperature, only one temperature or temperature = None) or 'xy' (pH and temperature).")
                    sys.exit()
                # Exporting
                if Ct_associated == 'xy':
                    if modeExport == 'plot':
                        ThSA.plot_(pH, T, DGr_export, iRxn, text_)
                    elif modeExport == 'Excel':
                        ThSA.writeExcel_(Ct, pH, T, DGr_export, iRxn, text_, Ct_associated, fullPathSave)
                elif Ct_associated == 'OnlyT':
                    if modeExport == 'plot':
                        ThSA.plotOnlyT(T, DGr_export, iRxn, text_)
                    elif modeExport == 'Excel':
                        ThSA.writeExcel_(Ct, pH, T, DGr_export, iRxn, text_, Ct_associated, fullPathSave)
                    Ct_associated = 'y'
                elif Ct_associated == 'OnlypH':
                    if modeExport == 'plot':
                        ThSA.plotOnlypH(pH, DGr_export, iRxn, text_)
                    elif modeExport == 'Excel':
                        ThSA.writeExcel_(Ct, pH, T, DGr_export, iRxn, text_, Ct_associated, fullPathSave)
                    Ct_associated = 'x'
                else:
                    if modeExport == 'plot':
                        ThSA.contourf_(pH, T, DGr_export, iRxn, text_)
                    elif modeExport == 'Excel':
                        ThSA.writeExcel_(Ct, pH, T, DGr_export, iRxn, text_, Ct_associated, fullPathSave)
        else:
            if isinstance(Ct, dict):
                nameComp = list(Ct)
                conc = np.array(list(Ct.values()))
                toRemove = ["{", "}", "[", "]"]
                pattern = '[' + ''.join(toRemove) + ']'
            for idCt in range(nCt):
                if isinstance(Ct, dict):
                    iCt = {nameComp [i]: [float(conc[i][idCt])] for i in range(len(nameComp))}
                    text_ = re.sub(pattern, '', str(iCt)).replace("'", "[").replace("[:", "]:") + ' (in mol/L)'
                else:
                    iCt = 1.0
                    text_ = 'Compound concentrations were not given.'
                DGr, infoRxn = ThSA.getDeltaGr(typeRxn, input_, phase, specComp, T, iCt, pH, asm, warnings)
                nDimDGr = DGr.ndim
                for idRxn, iRxn in enumerate(infoRxn):
                    # Selection of DGr for exporting
                    if not onepH and not oneT:
                        if nDimDGr == 3:
                            DGr_export = DGr[idRxn, :, :]
                        else:
                            DGr_export = DGr
                        if modeExport == 'plot':
                            ThSA.contourf_(pH, T, DGr_export, iRxn, text_)
                        elif modeExport == 'Excel':
                            ThSA.writeExcel_(Ct, pH, T, DGr_export, iRxn, text_, 'T-pH-Ct_noAssociated', fullPathSave)
                    elif onepH:
                        if nDimDGr == 2:
                            DGr_export = DGr[idRxn, :]
                        else:
                            DGr_export = DGr
                        if modeExport == 'plot':
                            ThSA.plotOnlyT(T, DGr_export, iRxn, text_)
                        elif modeExport == 'Excel':
                            ThSA.writeExcel_(Ct, pH, T, DGr_export, iRxn, text_, 'OnlyT-Ct_noAssociated', fullPathSave)
                    elif oneT:
                        if nDimDGr == 2:
                            DGr_export = DGr[idRxn, :]
                        else:
                            DGr_export = DGr
                        if modeExport == 'plot':
                            ThSA.plotOnlypH(pH, DGr_export, iRxn, text_)
                        elif modeExport == 'Excel':
                            ThSA.writeExcel_(Ct, pH, T, DGr_export, iRxn, text_, 'OnlypH-Ct_noAssociated', fullPathSave)
    
    def writeExcel_(Ct, pH, T, DGr_export, iRxn, text_, typeData, fullPathSave):
        """
        Write calculated DeltaGr in Excel document.

        """
        # Excel properties
        nameSheet = '∆Gr'
        startCol = 1
        # Compounds and concentrations
        if isinstance(Ct, dict):
            nComp = np.array(list(Ct.keys()))
            Concs = np.array(list(Ct.values()))
        # Info reaction
        infoR = pd.DataFrame(np.array([f'·Reaction: {iRxn}.']))
        if typeData == 'OnlyT' and not isinstance(Ct, dict):
            introRow = pd.DataFrame(np.array(['∆Gr in Excel | Only temperature dependency. ' + text_]))
            colNames = ['Temperature (K)', '∆Gr (kJ/mol)']
            dFram = np.array([T, DGr_export]).T
        elif typeData == 'OnlypH' and not isinstance(Ct, dict):
            introRow = pd.DataFrame(np.array(['∆Gr in Excel | Only pH dependency. ' + text_]))
            colNames = ['pH', '∆Gr (kJ/mol)']
            dFram = np.array([pH, DGr_export]).T
        elif typeData == 'OnlyT' and isinstance(Ct, dict):
            introRow = pd.DataFrame(np.array(['∆Gr in Excel | Only temperature dependency. ' + text_]))
            colNames = np.append(nComp, ['Temperature (K)', '∆Gr (kJ/mol)'])
            dFram = np.append(Concs, [T, DGr_export], axis = 0).T
        elif typeData == 'OnlypH' and isinstance(Ct, dict):
            introRow = pd.DataFrame(np.array(['∆Gr in Excel | Only pH dependency. ' + text_]))
            colNames = np.append(nComp, ['pH', '∆Gr (kJ/mol)'])
            dFram = np.append(Concs, [pH, DGr_export], axis = 0).T
        elif typeData == 'y' and isinstance(Ct, dict):
            introRow = pd.DataFrame(np.array(['∆Gr in Excel | Temperature & pH dependency. ' + text_]))
            colNames = np.append(nComp, ['Temperature (K)/pH'])
            colNames = np.append(colNames, np.float16(pH))
            dFram = np.append(Concs, [T], axis = 0).T
            dFram = np.append(dFram, DGr_export, axis = 1)
        elif typeData == 'x' and isinstance(Ct, dict):
            introRow = pd.DataFrame(np.array(['∆Gr in Excel | Temperature & pH dependency. ' + text_]))
            colNames = np.append(nComp, ['pH/Temperature (K)'])
            colNames = np.append(colNames, np.float16(T))
            dFram = np.append(Concs, [np.float16(pH)], axis = 0).T
            dFram = np.append(dFram, DGr_export, axis = 1)
        elif typeData == 'xy' and isinstance(Ct, dict):
            introRow = pd.DataFrame(np.array(['∆Gr in Excel | Temperature & pH dependency. ' + text_]))
            colNames = np.append(nComp, ['Temperature (K)', 'pH', '∆Gr (kJ/mol)'])
            dFram = np.append(Concs, [T, pH, DGr_export], axis = 0).T
        elif typeData == 'T-pH-Ct_noAssociated':
            introRow = pd.DataFrame(np.array(['∆Gr in Excel']))
            infoR = pd.DataFrame(np.array([f'·Reaction: {iRxn}. | {text_}']))
            colNames = np.append('Temperature (K)/pH', np.float16(pH))
            dFram = np.append(np.array([T]), DGr_export, axis = 0).T
        elif typeData == 'OnlyT-Ct_noAssociated':
            introRow = pd.DataFrame(np.array([f'∆Gr in Excel | Only temperature dependency (pH = {str(pH[0])}).']))
            infoR = pd.DataFrame(np.array([f'·Reaction: {iRxn}. | {text_}']))
            colNames = ['Temperature (K)', '∆Gr (kJ/mol)']
            dFram = np.array([T, DGr_export]).T
        elif typeData == 'OnlypH-Ct_noAssociated':
            introRow = pd.DataFrame(np.array([f'∆Gr in Excel | Only pH dependency (T = {str(T[0])}K).']))
            infoR = pd.DataFrame(np.array([f'·Reaction: {iRxn}. | {text_}']))
            colNames = ['pH', '∆Gr (kJ/mol)']
            dFram = np.array([pH, DGr_export]).T
        df = pd.DataFrame(dFram, columns = colNames)
        # Time to write
        if not os.path.isfile(fullPathSave):
            with pd.ExcelWriter(fullPathSave) as writer:
                introRow.to_excel(writer, sheet_name = nameSheet, index = False, header = False)
        with pd.ExcelWriter(fullPathSave, mode = 'a', if_sheet_exists='overlay') as writer:
            sRow = writer.sheets[nameSheet].max_row
            infoR.to_excel(writer, sheet_name = nameSheet, startrow = sRow+1, startcol = 0, index = False, header = False)
            df.to_excel(writer, sheet_name = nameSheet, startrow = sRow+2, startcol = startCol, index = False)

    def contourf_(pH, T, DGr_plot, iRxn, text_):
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
    
    def plot_(pH, T, DGr_plot, iRxn, text_):
        """
        Specific `plot()` function (from matplotlib.pyplot) for plotting DGr.
        
        """
        # Check if temperature axis must be inverted using Kendall's Tau.
        tau_pH, p_pH = kTau(pH, DGr_plot)
        tauSign_pH = tau_pH > 0
        tau_T, p_T = kTau(T, DGr_plot)
        tauSign_T = tau_T > 0
        if p_pH < 0.05 and p_T < 0.05:
            cTaus = tauSign_pH == tauSign_T
        else:
            cTaus = True
        fig = plt.figure()
        ax1 = fig.add_subplot()
        ax2 = fig.add_subplot()
        ax1.plot(pH, DGr_plot, 'k')
        ax1.set_xlabel('pH')
        ax1.set_ylabel('∆Gr (kJ/mol)')
        ax2.plot(T, DGr_plot, 'k')
        ax2.set_xlabel('Temperature (K)')
        ax2.xaxis.set_ticks_position('bottom')
        ax2.xaxis.set_label_position('bottom')
        ax2.axes.get_yaxis().set_visible(False)
        ax2.spines['bottom'].set_position(('axes', -0.20))
        if not cTaus:
            ax2.invert_xaxis()
        fig.text(0.5, 0.025, text_, horizontalalignment = 'center', wrap = True)
        fig.tight_layout(rect=(0,0.025,1,1)) 
        plt.title(iRxn)
        plt.show()
        
    def plotOnlyT(T, DGr_plot, iRxn, text_):
        """
        Specific `plot()` function (from matplotlib.pyplot) for plotting DGr.
        
        """
        fig = plt.figure()
        ax1 = fig.add_subplot()
        ax1.plot(T, DGr_plot, 'k')
        ax1.set_xlabel('Temperature (K)')
        ax1.set_ylabel('∆Gr (kJ/mol)')
        fig.text(0.5, 0.025, text_, horizontalalignment = 'center', wrap = True)
        fig.tight_layout(rect=(0,0.025,1,1)) 
        plt.title(iRxn)
        plt.gca().invert_xaxis()
        ax1.margins(x=0)
        plt.show()
    
    def plotOnlypH(pH, DGr_plot, iRxn, text_):
        """
        Specific `plot()` function (from matplotlib.pyplot) for plotting DGr.
        
        """
        fig = plt.figure()
        ax1 = fig.add_subplot()
        ax1.plot(pH, DGr_plot, 'k')
        ax1.set_xlabel('pH')
        ax1.set_ylabel('∆Gr (kJ/mol)')
        fig.text(0.5, 0.025, text_, horizontalalignment = 'center', wrap = True)
        fig.tight_layout(rect=(0,0.025,1,1)) 
        plt.title(iRxn)
        ax1.margins(x=0)
        plt.show()