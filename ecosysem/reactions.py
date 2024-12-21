# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 07:44:51 2024

@author: 2424069M
"""

import pandas as pd
import numpy as np
import sys

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
                if not iParam in col:
                    print(f'!EcoSysEM.Error: Parameter `{iParam}` not found in `{typeParam}.csv`. Please, add a column with `{iParam}` values.')
                    sys.exit()
    
    def getKinP(typeParam, params, reaction, sample = 'All', comp = None):
        """
        Function to get kinetic parameters from csv file. 

        Parameters
        ----------
        typeParam : STR
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
        dParam = pd.read_csv(KinP.path + typeParam + '.csv')
        KinP.checkKinP(typeParam, dParam, params)
        dParam = dParam.loc[dParam['Reaction'] == reaction]
        if sample != 'All':
            dParam = dParam.loc[dParam['Sample'].isin(sample)]
        sampleNames = dParam['Sample'].values
        dictR = {}
        for iParam in params:
            if not iParam in KinP.assocParamComp:
                dictR[iParam] = dParam[iParam].values
            else:
                if comp:
                    dictR[iParam] = np.squeeze(dParam[comp].fillna(0).values)
                    cVal = np.all(dictR[iParam] == 0.0)
                    if cVal:
                        print('!EcoSysEM.Error: Limiting substrate not found. Check `Ct` parameter in `getRs()` or `exportRs()`, and Km values in csv file.')
                        sys.exit()
                else:
                    print(f'!EcoSysEM.Error: You must to define the compounds for {iParam} with `comp` parameter.')
                    sys.exit()
        return dictR, sampleNames
    
class KinRates:
    """
    Class for kinetic rate calculations (biotic and abiotic transformations).
    
    """
    # Direction of kinetic parameters
    path = 'kinetics\\'
    
    def getRs(typeKin, typeParam, Rxn, Ct, sample = 'All', T = None, Tcorr = None):
        """
        Lorem ipsum...

        Parameters
        ----------
        typeKin : STR
            Type of kinetic equations (e.g., rMM - 'Michaelis-Menten equation').
        typeParam : STR
            Name of parameter database, matching with csv name.
        Rxn : STR
            Requested reaction.
        Ct : FLOAT, LIST or DICT
            Concentration of substrates, products and/or inhibitors.
        sample : STR or LIST, optional
            Requested samples (rows of `typeParam.csv`). The default is 'All'.
        T : FLOAT or LIST, optional
            Temperatures. The default is None.
        Tcorr : STR, optional
            Type of temperature correlation (e.g., Arrhenius - 'Arrhenius correlation'). The default is None.

        Returns
        -------
        rs : np.ndarray
            Resultant rates.
        sampleNames : STR or LIST
            Names of samples (rows of `typeParam.csv`).

        """
        if not isinstance(typeKin, str): typeKin = str(typeKin)
        if not isinstance(typeParam, str): typeParam = str(typeParam)
        if not isinstance(Rxn, str): Rxn = str(Rxn)
        if not isinstance(Ct, dict): 
            print('!EcoSysEM.Error: `Ct` parameter must be a dictionary.')
            sys.exit()
        else:
            npCt = np.array(list(Ct.items()), dtype = object).T
            lenCt = [len(i) for i in npCt[1, :]]
            nCt = list(set(lenCt))
            cLen = len(nCt) == 1
            if not cLen:
                print('!EcoSysEM.Error: All compounds must have same number of concentrations.')
                sys.exit()
            nCt = nCt[0]
        if sample != 'All':
            if not isinstance(sample, list): sample = [sample]
        if T:
            if isinstance(T, int): T = float(T)
            if isinstance(T, float): T = [T]
            if not Tcorr:
                print('!EcoSysEM.Error: You must to define the temperature correlation with `Tcorr` parameter.')
                sys.exit()
            else:
                if not isinstance(Tcorr, str): Tcorr = str(Tcorr)
        # Type Kinetics (Equation to calculate rs)   
        if typeKin == 'rMM':
            params = ['qmax', 'Km']
            comp = list(Ct.keys())
            Cs = pd.DataFrame(Ct).values
            p, sampleNames = KinP.getKinP(typeParam, params, Rxn, sample, comp)
            qmax = p['qmax']
            Km = p['Km']
            rs = KinRates.rMM(qmax, Km, Cs)
            if T:
                if Tcorr == 'Arrhenius':
                    p, _ = KinP.getKinP('ArrhCor', 'Coeff', Rxn)
                    O = p['Coeff']
                    # If more than one Arrhenius parameter, average value is taken.
                    if len(O) > 1:
                        O = np.mean(O)
                    p, _ = KinP.getKinP(typeParam, 'Texp', Rxn, sample)
                    Texp = p['Texp']
                    rs = KinRates.arrhCorr(rs, O, Texp, T)
                else:
                    print(f'!EcoSysEM.Error: {Tcorr} correlation not found. Available temperature correlations: "Arrhenius".')
                    sys.exit()
        else:
            print(f'!EcoSysEM.Error: {typeKin} equation not found. Available reaction equations: "rMM" (Michaelis-Menten eq.).')
            sys.exit()
        return np.squeeze(rs), sampleNames
    
    def rMM(qmax, Km, C):
        """
        Function to compute Michaelis-Menten equation.

        Parameters
        ----------
        qmax : FLOAT, LIST, np.ndarray
            Values of maximum uptake rate.
        Km : FLOAT, LIST, np.ndarray
            Values of half-saturation constants (or Michaelis constants).
        C : FLOAT, LIST, np.ndarray, DICT
            Concentrations of limiting substrates.

        Returns
        -------
        r : np.ndarray
            Resultant substrate uptake rates.

        """
        if not isinstance(qmax, np.ndarray): qmax = np.array(qmax)
        if qmax.ndim == 0: qmax = np.array([qmax])
        if not isinstance(Km, np.ndarray): Km = np.array(Km)
        if Km.ndim == 0: Km = np.array([Km])
        if isinstance(C, dict): 
            C = np.array(pd.DataFrame(C))
        else:
            if not isinstance(C, np.ndarray): C = np.array(C)
            if C.ndim == 0: C = np.array([C])
        # qmax
        nqmax = len(qmax)
        # Km
        if Km.ndim == 1:
            Km = np.array([Km])
            if Km.shape[0] == 1 and nqmax > 1:
                Km = Km.T
        nKm = Km.shape[0]
        # Checking number of qmax and Km values
        cN = nqmax == nKm
        if not cN:
            print('!EcoSysEM.Error: Same number of qmax and Km must be given.')
            sys.exit()
        else:
            nParam = nKm
        # C
        if C.ndim == 1:
            C = np.array([C])
            if C.shape[0] == 1 and Km.shape[1] == 1:
                C = C.T
        nConc = C.shape[0]
        # Checking number of limiting substrates
        nLS_Km = Km.shape[1]
        nLS_C = C.shape[1]
        cN = nLS_Km == nLS_C
        if not cN:
            print('!EcoSysEM.Error: Number of limiting substrates is not consistent. Km and C must have the same number of columns.')
            sys.exit()
        # Initializing result variable
        r = np.empty([nParam, nConc])
        # Calculation of rate(s)
        for iP in range(nParam):
            for iC in range(nConc):
                K = Km[iP,:]
                c = C[iC,:]
                M = np.prod((c) / (K + c + 1e-25)) # + 1e-25 to prevent NaN when conc == 0 and Ks == 0
                r[iP, iC] = qmax[iP] * M
        return np.squeeze(r) # (samples)x(concs), (samples) = (qmax) and (Km)
    
    def arrhCorr(rateBase, O, tempBase, temp):
        """
        Function to compute Arrhenius correlation.

        Parameters
        ----------
        rateBase : FLOAT or LIST or np.ndarray
            Reaction rate at base (measured) temparature (tempBase).
        O : FLOAT
            Arrhenius coefficient.
        tempBase : FLOAT or LIST
            Temperature base, that is, the original temperature of substrate rate.
        temp : FLOAT or LIST
            Set of temperatures.

        Returns
        -------
        rT : np.ndarray
            Resultant substrate uptake rates as function of temperature.

        """
        rBaslist = False
        if not isinstance(rateBase, np.ndarray): rateBase = np.array(rateBase)
        if rateBase.ndim == 0: rateBase = np.array([rateBase])
        if rateBase.ndim == 1: nRateBase = len(rateBase); rateBase = np.array([rateBase]); rBaslist = True
        if not isinstance(tempBase, np.ndarray): tempBase = np.array(tempBase)
        if tempBase.ndim == 0: tempBase = np.array([tempBase])
        if not isinstance(temp, list): temp = [temp]
        # Check number of rates and base temperatures
        if not rBaslist:
            nRateBase = rateBase.shape[0]
        nTempBase = len(tempBase)
        cN = nRateBase == nTempBase
        if not cN:
            print('!EcoSysEM.Error: Same number of rateBase and tempBase must be given.')
            sys.exit()
        else:
            nSamples = rateBase.shape[0]
        # Number of concentrations
        nConc = rateBase.shape[1]
        # Number of temperatures
        nT = len(temp)
        # Initializing result variable
        rT = np.empty([nSamples, nConc, nT])
        for iS in range(nSamples):
            for iC in range(nConc):
                for iT in range(nT):
                    rT[iS, iC, iT] = rateBase[iS, iC] * O ** (temp[iT] - tempBase[iS])
        return np.squeeze(rT) # (samples)x(concs)x(temp)

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
        warnings : BOOL
            Display function warnings. Default: False.
        
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
                print(f'!EcoSysEM.Error: Requested reaction(s) or compound(s) not found in {typeRxn}.csv file.')
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
                print(f'!EcoSysEM.Error: Requested input(s) found as a reaction and compound in {typeRxn}.csv file. '
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
        typeRxn : STR
            What reaction(s) type are requested, matching with csv name.
            Default: 'pHSpeciation'.
        warnings : BOOL
            Display function warnings. Default: False.
            
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
        warnings : BOOL
            Display function warnings. Default: False.
            
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
        warnings : BOOL
            Display function warnings. Default: False.
        
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
        warnings : BOOL
            Display function warnings. Default: False.

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