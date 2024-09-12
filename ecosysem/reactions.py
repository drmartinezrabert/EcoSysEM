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
    
    def getRxnByComp(typeRxn, compound):
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
        
        Returns
        -------
        rComp : LIST
            Involving compounds of reaction(s).
        mRxn : np.array
            Reaction matrix (compounds)x(reactions).        
        infoRxn : LIST
            Information of reaction if any (in parenthesis in csv/Excel file)
        
        """
        if isinstance(compound, list):
            if len(compound) > 1:
                print('!EcoSysEM.Error: Only one compound can be requested.')
                sys.exit()
            else:
                compound = compound[0]
        if not isinstance(compound, str): compound = str(compound)
        dRxn = pd.read_csv(Reactions.path + typeRxn + '.csv')
        dRxnAux = dRxn
        infoRxn = np.empty(0)
        # for iCompound in compounds:
        dRxnAux1 = dRxn.filter(regex = f'^{compound}/').dropna(how = 'all')
        dRxnAux2 = dRxn.filter(regex = f'/{compound}/').dropna(how = 'all')
        dRxnAux3 = dRxn.filter(like = f'/{compound} ').dropna(how = 'all')
        if not dRxnAux1.empty:
            dRxnAux = dRxnAux1
        elif not dRxnAux2.empty:
            dRxnAux = dRxnAux2
        elif not dRxnAux3.empty:
            dRxnAux = dRxnAux3
        else:
            return None, None, None
        else:
            rComp = list(dRxn['Compounds'])
            mRxn = np.array(dRxn.loc[:, dRxn.columns != 'Compounds'])
            headers = dRxn.columns.values[1:]
            for iHead in headers:
                if iHead.find('(') != -1:
                    newInfo = iHead[iHead.find('(')+1:iHead.find(')')]
                    infoRxn = np.append(infoRxn, newInfo)
            if not infoRxn.any():
                infoRxn = None
            return rComp, mRxn, infoRxn

#- DEBUGGING -#

#-------------#

#- Info of functions and examples -#
### Get reactions involving one or more compounds
#> Reactions.getRxn(typeRxn, compounds), where compound can be a string (one compound) or list (multiple compounds)
# rComp, mRxn, infoRxn = Reactions.getRxn('pHSpeciation', 'NH4+')
# print(rComp)
# print(mRxn)
# print(infoRxn)
#----------------------------------#