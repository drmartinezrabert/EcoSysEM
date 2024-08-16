# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 07:44:51 2024

@author: 2424069M
"""

import pandas as pd
import numpy as np

class Reactions:
    # Directory of reactions
    path = 'reactions\\'
    
    def getRxn(typeRxn, compounds):
        """
        Function to get reaction(s) involving the compound from csv file.
    
        Parameters
        ----------
        typeRxn : STR
            What reaction(s) type are requested, matching with csv name. E.g.:
                - 'pHSpeciation': pH equilibrium
        compounds : STR or LIST
            Name(s) of requested compound(s).
            STR - one compound; LIST - multiple compounds.
        
        Returns
        -------
        rComp : LIST
            Involving compounds of reaction(s).
        mRxn : np.array
            Reaction matrix (compounds)x(reactions).        
        infoRxn : LIST
            Information of reaction if any (in parenthesis in csv/Excel file)
        
        """
        if not isinstance(compounds, list): compounds = [compounds]
        dRxn = pd.read_csv(Reactions.path + typeRxn + '.csv')
        dRxnAux = dRxn
        for iCompound in compounds:
            dRxnAux = dRxnAux.filter(like = iCompound).dropna(how = 'all')
        dRxn = pd.concat([dRxn.iloc[dRxnAux.index, [0]], dRxnAux], axis=1).fillna(0)
        
        if dRxn.empty:
            print(f'!EcoSysEM.Warning: Reaction involving {compounds} not found.')
            return None, None, None
        else:
            rComp = list(dRxn['Compounds'])
            mRxn = np.array(dRxn.loc[:, dRxn.columns != 'Compounds'])
            headers = dRxn.columns.values[1:]
            infoRxn = np.empty(0)
            for iHead in headers:
                if iHead.find('(') != -1:
                    newInfo = iHead[iHead.find('(')+1:iHead.find(')')]
                    infoRxn = np.append(infoRxn, newInfo)
            if not infoRxn.any():
                infoRxn = None
            return rComp, mRxn, infoRxn

#- DEBUGGING -#
# dR = Reactions.getRxn('pHSpeciation', 'H2CO3')
# print(dR)
#-------------#

#- Info of functions and examples -#
### Get reactions involving one or more compounds
#> Reactions.getRxn(typeRxn, compound, path), where compound can be a string (one compound) or list (multiple compounds)
# dR = Reactions.getRxn('pHSpeciation', 'NH3')
# dR = Reactions.getRxn('pHSpeciation', ['HCO3-', 'CO3-2'])
#----------------------------------#