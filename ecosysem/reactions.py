# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 07:44:51 2024

@author: 2424069M
"""
# from thermodynamics import ThEq
import importlib
import pandas as pd
import numpy as np

class KinP:
    """
    Class for kinetic parameters.
    
    """
    # Direction of kinetic parameters
    path = 'kinetics\\'
    # Parameters associated to Compound (like 'Km', 'Ks' or 'Ki')
    assocParamComp = ['Km', 'Ks', 'Ki']
    
    def checkKinP(typeParam, dParam, params):
        """
        Function to check if parameters are in .csv file.

        Parameters
        ----------
        typeParam : STR
            Name of parameter database, matching with csv name.
        dParam : DICT
            Dictionary from `typeParam.csv` file.
        params : TYPE
            Requested parameters.

        Returns
        -------
        None.

        """
        col = dParam.columns.values
        for iParam in params:
            if not iParam in KinP.assocParamComp:
                if not iParam in col: raise ValueError(f'Parameter ({iParam}(`) not found in {typeParam}.csv. Please, add a column with {iParam} values.')
    
    def getKinP(paramDB, params, reaction, sample = 'All', compounds = None):
        """
        Function to get kinetic parameters from csv file. 

        Parameters
        ----------
        paramDB : STR
            Name of parameter database, matching with csv name.
        params : STR or LIST
            Requested parameters.
        reaction : STR
            Requested reaction.
        sample : STR or LIST, optional
            Requested samples (rows of `typeParam.csv`). The default is 'All'.
        comp : STR or LIST, optional
            Compounds of parameters associated to compounds (e.g., Km). The default is None.

        Returns
        -------
        dictR : DICT
            Dictionary with requested parameters.
        sampleNames : STR or LIST
            Names of samples (rows of `typeParam.csv`.
        """
        if not isinstance(params, list): params = [params]
        dParam = pd.read_csv(KinP.path + paramDB + '.csv', encoding_errors='ignore')
        KinP.checkKinP(paramDB, dParam, params)
        dParam = dParam.loc[dParam['Reaction'] == reaction]
        if sample != 'All':
            dParam = dParam.loc[dParam['Sample'].isin(sample)]
        sampleNames = dParam['Sample'].values
        dictR = {}
        for iParam in params:
            if not iParam in KinP.assocParamComp:
                dictR[iParam] = dParam[iParam].values
            else:
                if compounds:
                    for comp in compounds:
                        try:
                            dParam[comp]
                        except:
                            pass
                        else:
                            dictR[f'{iParam}_{comp}'] = dParam[comp].fillna(0).values
                            cVal = np.all(f'{iParam}_{comp}' == 0.0)
                            if cVal: raise ValueError(f'Limiting substrate ({comp}) not found. Check `Ct` argument in `getRs()`, and Km values in csv file.')
                else: raise ValueError(f'You must to define the compounds for {iParam} with `comp` parameter.')
        return dictR, sampleNames
    
class KinRates:
    """
    Class for kinetic rate calculations (biotic and abiotic transformations).
    
    """
    # Direction of kinetic parameters
    path = 'kinetics\\'
    
    def getRs(typeKin, paramDB, reactions, Ct, sample = 'All', pH = None, T = None):
        """
        Function to compute reaction rates.

        Parameters
        ----------
        paramDB : STR
            Type of kinetic equations 
                MM - 'Michaelis-Menten equation'.
                MM-Arrhenihus - 'Michaelis-Menten-Arrhenius equation'
        paramDB : STR
            Name of parameter database, matching with csv name.
        reactions : STR
            Requested reaction.
        Ct : DICT
            Concentration of substrates, products and/or inhibitors.
        sample : STR or LIST, optional
            Requested samples (rows of `paramDB.csv`). The default is 'All'.
        pH : FLOAT, optional
            pH values. The default is None.
        T : FLOAT, LIST or np.ndarray, optional
            Temperatures. The default is None.
        
        Returns
        -------
        Rs : DICT
            Resultant rates. Shape: (Z)x(Y)x(X)x(reactions).
        combNames : DICT
            Combination of samples (rows of `typeParam.csv`).
        orderComb : STR
            Order of samples in `combNames`.

        """
        # Dynamic import of ThEq
        MyModule = importlib.import_module('thermodynamics')
        ThEq = MyModule.ThEq
        # Check variables
        if not isinstance(T, np.ndarray): T = np.ndarray(T)
        if not isinstance(reactions, list): reactions = [reactions]
        if not isinstance(Ct, dict): raise TypeError('`Ct` argument must be a dictionary.')
        else:
            # Check shapes of Ct dictionary
            checkArray_compounds = []
            compounds = list(Ct.keys())
            for c in Ct:
                if not isinstance(Ct[c], np.ndarray): raise TypeError('All concentrations in `Ct` argument must be np.ndarray.')
                else:
                    checkArray_compounds += [Ct[c].shape]
            uShpConc = list(set(checkArray_compounds))
            if not len(uShpConc) == 1: raise ValueError('All compounds must have same shape.')
        if sample != 'All':
            if not isinstance(sample, list): sample = [sample]
        ## Concentrations (pH speciation)
        if pH:
            if T is None: raise ValueError('Temperature not defined (argument `T`).')
            for comp in compounds:
                Ct[comp] = ThEq.pHSpeciation(comp, pH, T, Ct[comp])
        if typeKin == 'MM':
            params = ['qmax', 'Km']
            # Initialize results
            Rs = {}
            combNames = {}
            for idRxn, Rxn in enumerate(reactions):
                # Initialize combNames[Rxn]
                combNames[Rxn] = []
                p, sampleNames = KinP.getKinP(paramDB, params, Rxn, sample, compounds)
                qmax = p['qmax']
                # Select Km
                Km = {k : p[k] for k in p if 'Km' in k}
                rs = {}
                for idMM, sample_MM in enumerate(sampleNames):
                    # Rates
                    qmax_ = qmax[idMM]
                    Km_ = {k : Km[k][idMM] for k in Km}
                    rs[f'comb_{idMM+1}'] = KinRates._rMM(qmax_, Km_, Ct)
                    # CombNames
                    combNames[Rxn].append(sample_MM)
                Rs[Rxn] = rs
            orderComb = 'MM'
            return Rs, combNames, orderComb
        elif typeKin == 'MM-Arrhenius':
            if len(paramDB) != 2: raise ValueError('For MM-Arrhenius equation, 2 databases must be given: paramDB = ["Michaelis-Menten DB", "Arrhenius DB"]')
            paramDB_MM = paramDB[0]
            paramDB_Arrh = paramDB[1]
            # Initialize results
            Rs = {}
            combNames = {}
            for idRxn, Rxn in enumerate(reactions):
                # Initialize combNames[Rxn]
                combNames[Rxn] = []
                ## Michaelis-Menten parameters
                p, sampleNames_MM = KinP.getKinP(paramDB_MM, ['qmax', 'Km'], Rxn, sample, compounds)
                qmax = p['qmax']
                # Select Km
                Km = {k : p[k] for k in p if 'Km' in k}
                ## Arrhenius parameters
                p, sampleNames_Arrh = KinP.getKinP(paramDB_Arrh, 'Coeff', Rxn)
                O = p['Coeff']
                p, _ = KinP.getKinP(paramDB_MM, 'Texp', Rxn, sample)
                Texp = p['Texp']
                ## Rate calculation (combinations of MM and Arrhenius)
                c = 1
                rs = {}
                for idMM, sample_MM in enumerate(sampleNames_MM):
                    for idArrh, sample_Arrh in enumerate(sampleNames_Arrh):
                        #-DEBUGGING-#
                        sampleComb = f'{sample_MM} - {sample_Arrh}'
                        #-----------#
                        ## Michaelis-Menten
                        qmax_ = qmax[idMM]
                        Km_ = {k : Km[k][idMM] for k in Km}
                        rs_ = KinRates._rMM(qmax_, Km_, Ct)
                        ## Arrhenius
                        O_ = O[idArrh]
                        Texp_ = Texp[idMM]
                        rs[f'comb_{c}'] = KinRates._arrhCorr(rs_, O_, Texp_, T)
                        c += 1
                        # CombNames
                        combNames[Rxn].append(sampleComb)
                Rs[Rxn] = rs
            orderComb = 'MM - Arrhenius'
            return Rs, combNames, orderComb
        else: raise ValueError(f'Unknown typeKin ({typeKin}). Existing typeKin: "MM" (Michaelis-Menten eq.); "MM-Arrhenius" (Michaelis-Menten-Arrhenius eq.).')
    
    def _rMM(qmax, Km, C):
        """
        Function to compute Michaelis-Menten equation.

        Parameters
        ----------
        qmax : FLOAT
            Value of maximum uptake rate.
        Km : DICT
            Values of half-saturation constants (or Michaelis constants).
        C : DICT
            Concentrations of limiting substrates.

        Returns
        -------
        r : DICT
            Resultant substrate uptake rates.

        """
        if not isinstance(qmax, np.ndarray): qmax = np.array(qmax)
        if qmax.ndim == 0: qmax = np.array([qmax])
        if not isinstance(Km, dict): raise TypeError('`Km` argument must be a dictionary. Km = {"Km_compound": [values]}')
        limCompds = [str(np.char.replace(k, 'Km_', '')) for k in Km if 'Km' in k]
        M = 1.0
        r = {}
        for limComp in limCompds:
            c = C[limComp]
            k = Km[f'Km_{limComp}']
            if not isinstance(k, (list, np.ndarray)): k = np.array(k)
            if isinstance(k, np.ndarray):
                if k.ndim != 1:
                    k = np.array([k])
            K = k * np.ones(c.shape)
            M *= (c) / (K + c + 1e-25)
        r = qmax * M
        return r
    
    def _arrhCorr(rateBase, O, tempBase, temp):
        """
        Function to compute Arrhenius correlation.

        Parameters
        ----------
        rateBase : FLOAT, LIST or np.ndarray
            Reaction rate at base (measured) temparature (tempBase).
        O : FLOAT
            Arrhenius coefficient.
        tempBase : FLOAT
            Temperature base, that is, the original temperature of substrate rate.
        temp : FLOAT, LIST or np.ndarray
            Set of temperatures.

        Returns
        -------
        rT : DICT
            Resultant substrate uptake rates as function of temperature.

        """
        if not isinstance(rateBase, np.ndarray): rateBase = np.array(rateBase)
        if not isinstance(temp, np.ndarray): rateBase = np.array(temp)
        if rateBase.shape != temp.shape: raise ValueError(f'Argument `rateBase` ({rateBase.shape}) and argument `temp` ({temp.shape}) must have the same shape.')
        rT = rateBase * O ** (temp - tempBase)
        return rT

class Reactions:
    # Directory of reactions
    path = 'reactions\\'
    
    def getRxn(typeRxn, input_ = 'All', warnings = False):
        """
        Function to get reaction(s) involving the requested compound or with the reaction name.
    
        Parameters
        ----------
        typeRxn : STR
            What reaction(s) type are requested, matching with csv name. E.g.:
                - 'pHSpeciation': pH equilibrium
                - 'metabolisms': metabolic activities
        input_ : STR or LIST
            Name(s) of requested compound(s) or reaction(s).
        warnings : BOOL, optional
            Display function warnings. The default is False.
        
        Returns
        -------
        rComp : LIST
            Involving compounds of reaction(s).
        mRxn : np.array
            Reaction matrix (compounds)x(reactions).        
        infoRxn : LIST
            Information of reaction if any (in parenthesis in csv/Excel file)
        
        """
        if not isinstance(typeRxn, list): 
            if not isinstance(typeRxn, str): typeRxn = str(typeRxn)
        if input_ == 'All':
            return Reactions.getAllRxn(typeRxn, warnings)
        else:
            if not isinstance(input_, list): input_ = [input_]
            gRbC_rComp, gRbC_mRxn, gRbC_infoRxn = Reactions.getRxnByComp(typeRxn, input_, warnings)
            gRbN_rComp, gRbN_mRxn, gRbN_infoRxn = Reactions.getRxnByName(typeRxn, input_, warnings)
            if (gRbC_rComp is None) and (gRbN_rComp is None):
                print(f'!EcoSysEM.Warning: Requested reaction(s) or compound(s) not found in {typeRxn}.csv file.')
                return None, None, None
            elif gRbC_rComp is None:
                rComp = gRbN_rComp
                mRxn = gRbN_mRxn
                infoRxn = gRbN_infoRxn
            elif gRbN_rComp is None:
                rComp = gRbC_rComp
                mRxn = gRbC_mRxn
                infoRxn = gRbC_infoRxn
            elif (gRbC_rComp is not None) and (gRbN_rComp is not None):
                print(f'!EcoSysEM.Warning: Requested input(s) found as a reaction and compound in {typeRxn}.csv file. '
                      'All inputs must be the same type: compound or reaction. If a same name is used for a compound and a reaction, '
                      'use `getRxnByComp()` or `getRxnByName()` instead.')
                return None, None, None
            return rComp, mRxn, infoRxn
    
    def getRxnpH(compounds, typeRxn='pHSpeciation', warnings = False):
        """
        Specific function to get reaction(s) involving the requested compound 
        on acid-base equilibrium (defined in `pHSpeciation.csv`)
    
        Parameters
        ----------
        compounds : STR or LIST
            Name(s) of requested compound(s).
            STR - one compound; LIST - multiple compounds.
        typeRxn : STR, optional
            What reaction(s) type are requested, matching with csv name. The default is 'pHSpeciation'.
        warnings : BOOL, optional
            Display function warnings. The default is False.
            
        Returns
        -------
        rComp : LIST
            Involving compounds of reaction(s).
        mRxn : np.array
            Reaction matrix (compounds)x(reactions).        
        infoRxn : LIST
            Information of reaction if any (in parenthesis in csv/Excel file)
        
        """
        if not isinstance(typeRxn, list): 
            if not isinstance(typeRxn, str): typeRxn = str(typeRxn)
            typeRxn = [typeRxn]
        if not isinstance(compounds, list): compounds = [compounds]
        dRxn = pd.read_csv(Reactions.path + typeRxn[0] + '.csv')
        rComp = list(dRxn['Compounds'])
        diffComp = []
        if len(typeRxn) > 1:
            typeRxn = typeRxn[1:]
            for iDB in typeRxn:
                dRxn_aux = pd.read_csv(Reactions.path + iDB + '.csv')
                rComp_aux = list(dRxn_aux['Compounds'])
                diffComp_aux = list(set(rComp_aux)-set(rComp))
                diffComp += diffComp_aux
            dRxn = pd.concat([dRxn, pd.DataFrame({'Compounds': diffComp})], ignore_index=True)
            for iDB in typeRxn:
                dRxn_aux = pd.read_csv(Reactions.path + iDB + '.csv')
                dRxn = pd.merge(dRxn, dRxn_aux, on='Compounds', how='left')
        headers = np.empty(0)
        infoRxn = np.empty(0)
        for iCompound in compounds:
            dRxnAux1 = dRxn.filter(regex = f'^{iCompound}/')
            dRxnAux2 = dRxn.filter(regex = f'/{iCompound}/')
            dRxnAux3 = dRxn.filter(regex = f'/{iCompound} *')
            dRxnSizes = [dRxnAux1.shape[1], dRxnAux2.shape[1], dRxnAux3.shape[1]]
            ind_dRxnMax = dRxnSizes.index(max(dRxnSizes))
            if ind_dRxnMax == 0:
                dRxnAux = dRxnAux1
            elif ind_dRxnMax == 1:
                dRxnAux = dRxnAux2
            elif ind_dRxnMax == 2:
                dRxnAux = dRxnAux3
            if dRxnAux.empty:
                if warnings:
                    print(f'!EcoSysEM.Warning: Reactions with {iCompound} as a compound not found.')
            else:
                iHeader = dRxnAux.columns.values[0:]
                for iHead in iHeader:
                    # Check if reaction has been already selected
                    c_Head = np.where(headers == iHead)[0].size
                    if c_Head == 0:
                        iRxn = iHead[iHead.find('(')+1:iHead.find(')')]
                        infoRxn = np.append(infoRxn, iRxn)
                        headers = np.append(headers, iHead)
        if infoRxn.size == 0:
            return None, None, None
        else:
            dRxnF = dRxn[headers].dropna(how = 'all')
            dRxnF = pd.concat([dRxn.iloc[dRxnF.index, [0]], dRxnF], axis=1).fillna(0)
            rComp = list(dRxnF['Compounds'])
            mRxn = np.array(dRxnF.loc[:, dRxnF.columns != 'Compounds'])
            return rComp, mRxn, infoRxn
        
    def getRxnByComp(typeRxn, compounds, warnings = False):
        """
        Function to get reaction(s) involving the requested compound.
    
        Parameters
        ----------
        typeRxn : STR
            What reaction(s) type are requested, matching with csv name. E.g.:
                - 'pHSpeciation': pH equilibrium
                - 'metabolisms': metabolic activities
        compounds : STR or LIST
            Name(s) of requested compound(s).
            STR - one compound; LIST - multiple compounds.
        warnings : BOOL, optional
            Display function warnings. The default is False.
            
        Returns
        -------
        rComp : LIST
            Involving compounds of reaction(s).
        mRxn : np.array
            Reaction matrix (compounds)x(reactions).        
        infoRxn : LIST
            Information of reaction if any (in parenthesis in csv/Excel file)
        
        """
        if not isinstance(typeRxn, list): 
            if not isinstance(typeRxn, str): typeRxn = str(typeRxn)
            typeRxn = [typeRxn]
        if not isinstance(compounds, list): compounds = [compounds]
        dRxn = pd.read_csv(Reactions.path + typeRxn[0] + '.csv')
        rComp = list(dRxn['Compounds'])
        diffComp = []
        if len(typeRxn) > 1:
            typeRxn = typeRxn[1:]
            for iDB in typeRxn:
                dRxn_aux = pd.read_csv(Reactions.path + iDB + '.csv')
                rComp_aux = list(dRxn_aux['Compounds'])
                diffComp_aux = list(set(rComp_aux)-set(rComp))
                diffComp += diffComp_aux
            dRxn = pd.concat([dRxn, pd.DataFrame({'Compounds': diffComp})], ignore_index=True)
            for iDB in typeRxn:
                dRxn_aux = pd.read_csv(Reactions.path + iDB + '.csv')
                dRxn = pd.merge(dRxn, dRxn_aux, on='Compounds', how='left')
        headers = np.empty(0)
        infoRxn = np.empty(0)
        for iCompound in compounds:
            dRxnAux = dRxn[dRxn['Compounds'].isin([iCompound])].dropna(axis='columns')
            if dRxnAux.empty:
                if warnings:
                    print(f'!EcoSysEM.Warning: Reactions with {iCompound} as a compound not found.')
            else:
                iHeader = dRxnAux.columns.values[1:]
                for iHead in iHeader:
                    # Check if reaction has been already selected
                    c_Head = np.where(headers == iHead)[0].size
                    if c_Head == 0:
                        iRxn = iHead[iHead.find('(')+1:iHead.find(')')]
                        infoRxn = np.append(infoRxn, iRxn)
                        headers = np.append(headers, iHead)
        if infoRxn.size == 0:
            return None, None, None
        else:
            dRxnF = dRxn[headers].dropna(how = 'all')
            dRxnF = pd.concat([dRxn.iloc[dRxnF.index, [0]], dRxnF], axis=1).fillna(0)
            rComp = list(dRxnF['Compounds'])
            mRxn = np.array(dRxnF.loc[:, dRxnF.columns != 'Compounds'])
            return rComp, mRxn, infoRxn
    
    def getRxnByName(typeRxn, nameRxn, warnings = False):
        """
        Function to get reaction(s) with the info/name (parenthesis in Excel).
        Do not use this function for pHSpeciation(). Instead use getRxnByComp().
        
        Parameters
        ----------
        typeRxn : STR
            What reaction(s) type are requested, matching with csv name. E.g.:
                - 'pHSpeciation': pH equilibrium
                - 'metabolisms': metabolic activities
        compounds : STR or LIST
            Name(s) of requested compound(s).
            STR - one compound; LIST - multiple compounds.
        warnings : BOOL, optional
            Display function warnings. The default is False.
        
        Returns
        -------
        rComp : LIST
            Involving compounds of reaction(s).
        mRxn : np.array
            Reaction matrix (compounds)x(reactions).        
        infoRxn : LIST
            Information of reaction if any (in parenthesis in csv/Excel file)
        """
        if not isinstance(typeRxn, list): 
            if not isinstance(typeRxn, str): typeRxn = str(typeRxn)
            typeRxn = [typeRxn]
        if not isinstance(nameRxn, list): nameRxn = [nameRxn]
        dRxn = pd.read_csv(Reactions.path + typeRxn[0] + '.csv')
        rComp = list(dRxn['Compounds'])
        diffComp = []
        if len(typeRxn) > 1:
            typeRxn = typeRxn[1:]
            for iDB in typeRxn:
                dRxn_aux = pd.read_csv(Reactions.path + iDB + '.csv')
                rComp_aux = list(dRxn_aux['Compounds'])
                diffComp_aux = list(set(rComp_aux)-set(rComp))
                diffComp += diffComp_aux
            dRxn = pd.concat([dRxn, pd.DataFrame({'Compounds': diffComp})], ignore_index=True)
            for iDB in typeRxn:
                dRxn_aux = pd.read_csv(Reactions.path + iDB + '.csv')
                dRxn = pd.merge(dRxn, dRxn_aux, on='Compounds', how='left')
        headers = np.empty(0)
        infoRxn = np.empty(0)
        for iRxn in nameRxn:
            dRxnAux1 = dRxn.filter(regex = f'^{iRxn} *').dropna(how='all')
            dRxnAux2 = dRxn.filter(like = f'({iRxn})').dropna(how = 'all')
            dRxnSizes = [dRxnAux1.shape[1], dRxnAux2.shape[1]]
            ind_dRxnMax = dRxnSizes.index(max(dRxnSizes))
            if ind_dRxnMax == 0:
                dRxnAux = dRxnAux1
                iHead = str(dRxnAux.columns.values)
                iRxn = iHead[iHead.find('(')+1:iHead.find(')')]
            elif ind_dRxnMax == 1:
                dRxnAux = dRxnAux2
            iHeader = dRxnAux.columns.values[0:]
            headers = np.append(headers, iHeader)
            if not dRxnAux.empty:
                infoRxn = np.append(infoRxn, iRxn)
            else:
                if warnings:
                    print(f'!EcoSysEM.Warning: Reaction {iRxn} not found.')
        if infoRxn.size == 0:
            return None, None, None
        else:
            dRxnF = dRxn[headers].dropna(how = 'all')
            dRxnF = pd.concat([dRxn.iloc[dRxnF.index, [0]], dRxnF], axis=1).fillna(0)
            rComp = list(dRxnF['Compounds'])
            mRxn = np.array(dRxnF.loc[:, dRxnF.columns != 'Compounds'])
            return rComp, mRxn, infoRxn

    def getAllRxn(typeRxn, warnings = False):
        """
        Function to get all reaction from reaction database (`typeRxn.csv`).

        Parameters
        ----------
        typeRxn : STR
            What reaction(s) type are requested, matching with csv name. E.g.:
                - 'pHSpeciation': pH equilibrium
                - 'metabolisms': metabolic activities
        warnings : BOOL, optional
            Display function warnings. The default is False.

        Returns
        -------
        rComp : LIST
            Involving compounds of reaction(s).
        mRxn : np.array
            Reaction matrix (compounds)x(reactions).        
        infoRxn : LIST
            Information of reaction if any (in parenthesis in csv/Excel file)

        """
        if not isinstance(typeRxn, list): typeRxn = [typeRxn]
        dRxn = pd.read_csv(Reactions.path + typeRxn[0] + '.csv')
        rComp = list(dRxn['Compounds'])
        diffComp = []
        if len(typeRxn) > 1:
            typeRxn = typeRxn[1:]
            for iDB in typeRxn:
                dRxn_aux = pd.read_csv(Reactions.path + iDB + '.csv')
                rComp_aux = list(dRxn_aux['Compounds'])
                diffComp_aux = list(set(rComp_aux)-set(rComp))
                diffComp += diffComp_aux
            dRxn = pd.concat([dRxn, pd.DataFrame({'Compounds': diffComp})], ignore_index=True)
            for iDB in typeRxn:
                dRxn_aux = pd.read_csv(Reactions.path + iDB + '.csv')
                dRxn = pd.merge(dRxn, dRxn_aux, on='Compounds', how='left')
        # Get compounds list
        rComp = list(dRxn['Compounds'])
        # Get info reactions (or names)
        infoRxn = np.empty(0)
        headers = dRxn.columns.values[1:]
        for iHead in headers:
            iRxn = iHead[iHead.find('(')+1:iHead.find(')')]
            infoRxn = np.append(infoRxn, iRxn)
        # Get Stoichiometric matrix
        dRxnF = dRxn.dropna(how = 'all')
        dRxnF = pd.concat([dRxn.iloc[dRxnF.index, [0]], dRxnF], axis=1).fillna(0)
        mRxn = np.array(dRxnF.loc[:, dRxnF.columns != 'Compounds'])
        if infoRxn.size == 0:
            return None, None, None
        else:
            return rComp, mRxn, infoRxn