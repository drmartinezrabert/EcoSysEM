# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 07:44:51 2024

@author: 2424069M
"""

import pandas as pd
import numpy as np
import sys

class Reactions:
    # Directory of reactions
    path = 'reactions\\'
    
    def getRxn(typeRxn, input_, warnings = False):
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
        if not isinstance(typeRxn, str): typeRxn = str(typeRxn)
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
        if not isinstance(typeRxn, str): typeRxn = str(typeRxn)
        if not isinstance(compounds, list): compounds = [compounds]
        dRxn = pd.read_csv(Reactions.path + typeRxn + '.csv')
        headers = np.empty(0)
        infoRxn = np.empty(0)
        for iCompound in compounds:
            dRxnAux1 = dRxn.filter(regex = f'^{iCompound}/')
            dRxnAux2 = dRxn.filter(regex = f'/{iCompound}/')
            dRxnAux3 = dRxn.filter(like = f'/{iCompound} ')
            if not dRxnAux1.empty:
                dRxnAux = dRxnAux1
            elif not dRxnAux2.empty:
                dRxnAux = dRxnAux2
            elif not dRxnAux3.empty:
                dRxnAux = dRxnAux3
            else:
                dRxnAux = dRxnAux1
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
        
        """
        if not isinstance(typeRxn, str): typeRxn = str(typeRxn)
        if not isinstance(nameRxn, list): nameRxn = [nameRxn]
        dRxn = pd.read_csv(Reactions.path + typeRxn + '.csv')
        headers = np.empty(0)
        infoRxn = np.empty(0)
        for iRxn in nameRxn:
            dRxnAux = dRxn.filter(like = f'({iRxn})').dropna(how = 'all')
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
