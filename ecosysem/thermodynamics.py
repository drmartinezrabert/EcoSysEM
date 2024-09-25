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
# import matplotlib.pyplot as plt

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
            rComp, mRxn, infoRxn = Rxn.getRxnByComp('pHSpeciation', iCompound)
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
            nCompounds = Rxn.getRxnByComp('pHSpeciation', iCompound)[0][1:]
            # Simplification hydration/dehydration equil.: [(CO2)aq] >>>> [H2CO3]
            nCompounds = [nC.replace('H2CO3', 'H2CO3 (CO2)') for nC in nCompounds]
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
    
    def getDeltaGr(typeRxn, input_, specComp = False, T = 298.15, Ct = 1.0, pH = 7.0, asm = 'stoich', warnings = False):
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
        if isinstance(T, float or int): T = [T]
        nT = len(T)
        if isinstance(pH, float or int): pH = [pH]
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
            deltaG0f = ThP.getThP('deltaG0f', i_rComp, 'L')[0]
            deltaG0r = ThP.getDeltaG0r(deltaG0f, i_mRxn)                        # kJ
            deltaG0r = deltaG0r / vSelected                                     # kJ/mol-i (if specDGr)
            # Calculate DeltaH0r
            deltaH0f = ThP.getThP('deltaH0f', i_rComp, 'L')[0]
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
                                    # iComp is a substrate
                                    print(f'!EcoSysEM.Error: Substrate concentration(s) for {iComp} not found.')
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
    
    def plotDeltaGr(T, pH, typeRxn, input_, specComp = False, Ct = 1.0, Ct_associated = None, asm = 'stoich', warnings = False):
        """
        Plot DeltaGr in function of pH, temperature and/or compound
        concentrations.

        Parameters
        ----------
        T : LIST or np.array
            Set of temperatures [K].
        pH : LIST or np.array
            Set of pHs.
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
        Ct_associated : STR ('x' or 'y')
            Set if concentration is associated to coordinate x or y. The default is None.
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
        # Checking arguments and variable definition
        if not isinstance(typeRxn, str): typeRxn = str(typeRxn)
        if not isinstance(input_, list): input_ = [input_]
        if isinstance(T, float or int): T = [T]
        if isinstance(pH, float or int): pH = [pH]
        # Variable lenghts
        nT = len(T)
        npH = len(pH)
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
            # Concentration values
            npCt = np.array(list(Ct.items()), dtype = object)
        if nT != npH:
            print('!EcoSysEM.Error: Temperature and pH must have the same size.')
            sys.exit()
        # Plotting DGr
        if Ct_associated:
            if Ct_associated == 'x' and nCt != npH:
                print('!EcoSysEM.Error: Total concentrations and pH must have the same size.')
                sys.exit()
            elif Ct_associated == 'y' and nCt != nT:
                print('!EcoSysEM.Error: Total concentrations and temperature must have the same size.')
                sys.exit()
            elif Ct_associated != 'x' and Ct_associated != 'y':
                print('!EcoSysEM.Error: Ct_associated must be "x" or "y".')
                sys.exit()
            DGr, infoRxn = ThSA.getDeltaGr(typeRxn, input_, specComp, T, Ct, pH, asm, warnings)
            for idRxn, iRxn in enumerate(infoRxn):
                if Ct_associated == 'x':
                    index_pH = list(range(npH))
                    DGr_plot = DGr[idRxn, :, index_pH, index_pH]
                    text_ = 'Concentrations (in mol/L) associated to pH.'
                elif Ct_associated == 'y':
                    index_T = list(range(nT))
                    DGr_plot = DGr[idRxn, index_T, :, index_T]
                    text_ = 'Concentrations (in mol/L) associated to T.'
                else:
                    print("`Ct_associated` argument must be 'x' or 'y'.")
                    sys.exit()
                ThSA.contourf_(pH, T, DGr_plot, iRxn, text_)
        else:
            nameComp = list(Ct)
            conc = np.array(list(Ct.values()))
            toRemove = ["\{", "\}", "\[", "\]"]
            pattern = '[' + ''.join(toRemove) + ']'
            for idCt in range(nCt):
                iCt = {nameComp [i]: [float(conc[i][idCt])] for i in range(len(nameComp))}
                text_ = re.sub(pattern, '', str(iCt)).replace("'", "[").replace("[:", "]:") + ' (in mol/L)'
                DGr, infoRxn = ThSA.getDeltaGr(typeRxn, input_, specComp, T, iCt, pH, asm, warnings)
                for idRxn, iRxn in enumerate(infoRxn):
                    DGr_plot = DGr[idRxn, :, :]
                    ThSA.contourf_(pH, T, DGr_plot, iRxn, text_)            
    
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
## Get DeltaGr (Thermodynamic state analysis)
# Concentrations in mol/L
# conc = {'CO': [1.08e-10],
#         'O2': [3.40e-4],
#         'CO2': [2.78e-5],
#         'NH3': [1.33e-6]}
# temperature = 255.65
# pH = 7.0
# ---------------------------------------
# conc = {'CO': [1.08e-10, 1],
#         'O2': [3.40e-4, 1],
#         'CO2': [2.78e-5, 1],
#         'NH3': [1.33e-6, 1]}
# temperature = [216.65, 255.65, 288.15]
# pH = [2.0, 5.0, 8.0]
# DGr, rInfoRxn = ThSA.getDeltaGr('metabolisms', ['CO', 'NH3'], True, temperature, conc, pH)
# print(DGr)
# print(rInfoRxn)
# print('')
# DGr, rInfoRxn = ThSA.getDeltaGr('metabolisms', ['COOB', 'AOM', 'CMX'], ['CO', 'NH3', 'NH3'], temperature, conc, pH)
# print(DGr)
# print(rInfoRxn)

## Plot DeltaGr (Thermodynamic state analysis)
# Concentrations in mol/L
# temperature = [288.15, 281.65, 275.15, 268.65, 262.15, 255.65, 249.15, 242.65, 236.15, 229.65, 223.15, 216.65]
# pH = []
#-------------#