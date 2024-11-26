# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 13:24:02 2024

@author: emartinez
"""

from reactions import Reactions as Rxn

import sys
import pandas as pd
from iteration_utilities import deepflatten

class EcoA():
    """
    Class for Ecosystem Analysis module (EcoA).
    
    """
    def rxnDirEdges(typeRxn, input_ = 'All', rxnDef = None, cycles = None, warnings = False):
        """
        Function to get directional hyperedges from stoichiometric matrix 
        (read from `typeRxn.csv`) and cycle hyperedges (if cycles are defined).

        Parameters
        ----------
        typeRxn : STR
            What reaction(s) type are requested, matching with csv name. E.g.:
                - 'pHSpeciation': pH equilibrium
                - 'metabolisms': metabolic activities
        input_ : STR or LIST
            Name(s) of requested compound(s) or reaction(s). The default is 'All'.
        rxnDef : DICT, optional
            A dictionary with abbreviations of reactions. The default is None.
            If None : rxnDef = {'Reaction 1': 'T1',..., 'Reaction N': 'TN'}.
            If 'Same' : Abbreviations from reaction names. rxnDef = {'Rxn1': 'Rxn1',..., 'RxnN': 'RxnN'}.
            If dict : {'Reaction 1': 'Abbr 1', ..., 'Reaction N': 'Abbr N'}
                The number of dict keys must coincide with the number of reactions.
        cycles : DICT, optional
            Classification of reactions by biogeochemical cycles or other grouping type. The default is None.
        warnings : BOOL
            Display function warnings. Default: False.

        Returns
        -------
        diedges : DICT
            A dictionary of directed hyperedges and their nodes.
            diedges = {'DirHyperedge': [['Arrowtail node(s)'], ['Arrowhead node(s)']]}
        edges : DICT
            A dictionary of hyperedges and their nodes.
            edges = {'Edge': ['Nodes']}
            If value of a cycle is nan `{'Cycle 1': [nan]}`, a wrong reaction name was given.

        """
        # Get reactions 
        compNames, StoykM, rxnNames = Rxn.getRxn(typeRxn, input_, warnings)
        # Reaction abbreviations
        if not rxnDef:
            rxnDef = {}
            for idRxn, iRxn in enumerate(rxnNames):
                abbrRxn = 'T' + str(idRxn+1)
                rxnDef[iRxn] = abbrRxn
            rxnAbbr = list(rxnDef.values()) # prev. reacSymb
        elif rxnDef == 'Same':
            rxnAbbr = rxnNames # prev. reacSymb
            rxnDef = dict(zip(rxnAbbr, rxnAbbr))
        else:
            if not isinstance(rxnDef, dict):
                print('!EcoSysEM.Warning: `rxnDef` must be a dictionary.')
                sys.exit()
            rxnAbbr = list(rxnDef.values()) # prev. reacSymb
            if len(rxnAbbr) != len(rxnNames):
                print(f'!EcoSysEM.Warning: The number of reactions ({len(rxnNames)}) and abbreviations ({len(rxnAbbr)}) mismatch.')
                print(f'> Reaction names: {rxnNames}.')
                sys.exit()
        # DirEdges definition from Stoichiometric matrix
        diedges = {}
        tSubs = []
        tProd = []
        for idColumn, column in enumerate(StoykM.T):
            diedges[rxnAbbr[idColumn]] = None
            for idRow, row in enumerate(column):
                if row == 0:
                    pass
                elif row > 0:
                    tProd.append(compNames[idRow])
                elif row < 0:
                    tSubs.append(compNames[idRow])
            diedges[rxnAbbr[idColumn]] = [tSubs, tProd]
            # Clean up variables
            tSubs = []
            tProd = []
        # Edge cycles definition
        if not cycles:
            return diedges
        else:
            if not isinstance(cycles, dict):
                print('!EcoSysEM.Warning: `cycles` must be a dictionary.')
                sys.exit()
            for cycle in cycles:
                iCycle = cycles[cycle]
                symbCycle = (pd.Series(iCycle)).map(rxnDef)
                combCycle = list((pd.Series(symbCycle)).map(diedges))  
                combCycleMerged = list(set(deepflatten(combCycle, types=list)))
                cycles[cycle] = combCycleMerged
            edges = cycles
            return diedges, edges

class EcoM():
    """
    Class for Ecosystem Modelling module (EcoM).
    
    """
    pass