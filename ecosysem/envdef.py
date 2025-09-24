# -*- coding: utf-8 -*-
"""
envdef - define and create environment(s)
=========================================

Summary
-------

Lorem ipsum...

Classes
-------
    
    :py:class:`ISA` - International Standard Atmosphere (ISA).

Functions
--------

    :py:fun:`lorem_ipsum` - Lorem ipsum.

------------------------------------------------------------------------------
"""
from thermodynamics import ThEq as eQ
from scipy.interpolate import RegularGridInterpolator
from molmass import Formula

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import earthaccess
import itertools
import xarray as xr
import datetime
import calendar
import os
import time
import cdsapi
import zipfile

import warnings
warnings.simplefilter(action = 'ignore')

class ISA:
    """
    International Standard Atmosphere (ISA). Static atmospheric model of how
    pressure, temperature, density and viscosity of the Earth's atmosphere
    change over a wide range of altitudes.
    --------------------------------------------------------------------------
    References: - ISO 2533:1975
                - National Oceanic and Atmospheric Administration (NOAA)
                - Goody & Yung (1989), doi: 10.1093/oso/9780195051346.001.0001
    
    """
    def __init__(self, layers, H2O, pH, resolution = 1000):
        # Layers of the Earth's atmosphere (ISA)
        _ISAproperties = {
                        'Layer': [0, 1, 2, 3, 4, 5, 6, 7],
                        'Layer name': ['Troposphere', 'Tropopuse', 'Stratosphere', 'Stratosphere', 'Stratopause', 'Mesosphere', 'Mesosphere', 'Mesopause'],
                        'Start altitude': [0, 11000, 20000, 32000, 47000, 51000, 71000, 84852],         # [m (above MSL)]
                        'End altitude': [11000, 20000, 32000, 47000, 51000, 71000, 84852, 85852],       # [m (above MSL)]
                        'Lapse rate': [-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0],                      # [C/km)]
                        'Base temperature': [15.0, -56.5, -56.5, -44.5, -2.5, -2.5, -58.5, -86.204],      # [C]
                        'Base pressure': [101325, 22632, 5474.9, 868.02, 110.91, 66.939, 3.9564, 0],    # [Pa]
                        }
        dISA = pd.DataFrame(data = _ISAproperties)
        # Dry composition
        dryComposition = {
                         'Compounds': ['N2', 'O2', 'Ar', 'CO2', 'Ne', 
                                       'He', 'CH4', 'Kr', 'H2', 'N2O', 
                                       'CO', 'Xe', 'O3', 'NO2', 'SO2', 
                                       'I2','NH3', 'HNO2', 'HNO3', 'H2S'], # Formula
                         'Compositions': [7.8084e-1, 2.0946e-1, 9.34e-3, 4.26e-4, 1.8182e-5, 
                                          5.24e-6, 1.92e-6, 1.14e-6, 5.5e-7, 3.3e-7, 
                                          1e-7, 9e-8, 7e-8, 2e-8, 1.5e-8, 
                                          1e-8, 6e-9, 1e-9, 1e-9, 3.3e-10] # [%vol]
                          }
        dDC = pd.DataFrame(data = dryComposition)
        self.layers = layers
        self.pH = pH
        self.resolution = resolution
        self._computeTandP_ISA(layers, dISA)
        self.compounds = dDC['Compounds']
        self.compositions = pd.Series(dDC.Compositions.values, index = dDC.Compounds).to_dict()
        self._computeWaterContent(H2O, dDC)
    
    def _computeTandP_ISA(self, layers, dISA):
        """
        Compute the change of temperature and pressure of the Earth's
        atmosphere over the range of altitudes. Based on ISA (ISO 2533:1975).
        
        """
        # Constants
        R = 8.3144598                   # Universal gas constant [J/mol/K]
        g0 = 9.80665                    # Gravitational acceleration [m/s^2]
        M0 = 0.0289644                  # Molar mass of Earth's air
        resolution = self.resolution    # Resolution (size of altitude nodes per layer) [m]
        # Data from ISA to compute temperature
        if layers == 'All':
            layers = range(8)
        elif type(layers) == int:
            layers = range(layers, layers+1)
        lapse_rate = (dISA.loc[layers]['Lapse rate']).to_numpy()
        base_T = (dISA.loc[layers]['Base temperature']).to_numpy()
        base_P = (dISA.loc[layers]['Base pressure']).to_numpy()
        start_alt = (dISA.loc[layers]['Start altitude']).to_numpy()
        end_alt = (dISA.loc[layers]['End altitude']).to_numpy()
        # Temperature calculation (ISA)
        alt = t = p = []
        for i in range(lapse_rate.size):
            alt_aux = np.array(range(int(start_alt[i]), int(end_alt[i]), resolution))
            # Temperature
            t_aux = base_T[i] + lapse_rate[i] * ((alt_aux - start_alt[i]) / 1000)
            # Pressure
            if lapse_rate[i] != 0:
                p_aux = base_P[i] * (1 + ((lapse_rate[i]/1000) / (base_T[i]+273.15)) * (alt_aux - start_alt[i])) ** (-(g0 * M0) / (R * lapse_rate[i]/1000))
            else:
                p_aux = base_P[i] * np.exp((-(g0 * M0) / (R * (base_T[i]+273.15))) * (alt_aux - start_alt[i]))
            # Appending...
            alt = np.append(alt, alt_aux)
            t = np.append(t, t_aux)
            p = np.append(p, p_aux)
        # Last value of arrays (Altitude, Temperature and Pressure)
        alt_end = min(alt_aux[-1] + resolution, end_alt[i])
        alt = np.append(alt, alt_end)
        t = np.append(t, base_T[-1] + lapse_rate[-1] * ((alt_end - start_alt[-1]) / 1000)) # Last T value
        if lapse_rate[-1] != 0:
            p = np.append(p,  base_P[-1] * (1 + ((lapse_rate[-1]/1000) / (base_T[-1]+273.15)) * (alt_end - start_alt[-1])) ** (-(g0 * M0) / (R * lapse_rate[-1]/1000)))
        else:
            p = np.append(p, base_P[-1] * np.exp((-(g0 * M0) / (R * (base_T[-1]+273.15))) * (alt_end - start_alt[-1])))
        # Safe Altitude, Temperature and Pressure
        self.ISAaltitude = alt # [m]
        self.ISAtemperature = t # [C]
        self.ISApressure = p # [Pa]
    
    def _computeWaterContent(self, H2O, dDC):
        """
        It is assumed that atmospheric compositon of minor components (<0.9%, 
        such as CO2) does not change significantly when including water vapor 
        (i.e., wet composition). Therefore, compositions of N2, O2 and Ar are
        recalculated (%vol).
        
        """
        self.H2O = H2O
        newComp = np.array(dDC['Compositions'])
        # self.composition = dryComposition
        if H2O > 0:
            newComp[0] -= 0.78100 * H2O # N2_wet
            newComp[1] -= 0.20925 * H2O # O2_wet
            newComp[2] -= 0.01100 * H2O # Ar_wet
            self.ISAcompositions = pd.Series(newComp, index = dDC.Compounds).to_dict()
    
    def getConcISA(self, phase, compound = None):
        """
        Computation of vertical profiles of compounds (parcial pressure, Pi;
        gas concentration, Ci_G; liquid concentration in fresh water, Ci_L-FW;
        and liquid concentration in sea water, Ci_L-SW).
        Gas concentrations (Ci_G) are calculated using Dalton's law and the 
        ideal gas law, and liquid concentration (Ci_LFW and Ci_LSW) with 
        Henry's law.
        
        Parameters
        ----------
        phase : STR ('G', 'L-FW', 'L-SW', 'L' or 'All')
            DESCRIPTION. Selection of phase of vertical profile.
                        'G' - Gas.
                        'L-FW' - Liquid fresh water.
                        'L-SW' - Liquid sea water.
                        'L' - Both liquid phases (L-FW, L-SW).
                        'All' - All phaes (G, L-FW, L-SW).
        compound : STR or LIST, optional
            DESCRIPTION. Interested compounds. Default: None -> All compounds.

        Returns
        -------
        Dictionaries with pressures/concentrations.
        dict_Pi, dict_Ci_G : DICT (if phase='G')
        dict_Ci_LFW : DICT (if phase='L-FW')
        dict_Ci_LSW : DICT (if phase='L-SW')
        dict_Ci_LFW, dictCi_LSW : DICT (if phase='L')
        dict_Pi, dict_Ci_G, dict_Ci_LFW, dict_Ci_LSW : DICT (if phase='All')
             dict_Pi : Parcial pressure of desired compounds.
             dict_Ci_G : Concentration in gas of desired compounds.
             dict_Ci_LFW : Concentration in liquid (freshwater) of desired compounds.
             dict_Ci_LSW : Concentration in liquid (seawater) of desired compounds.
        
        """
        # Data
        p = self.ISApressure
        t = self.ISAtemperature + 273.15   # [K]
        compounds = self.compounds
        compositions = np.array(list(self.compositions.values()))
        if compound:
            if type(compound) is str: compound = [compound]
            findC = compounds.reset_index().set_index('Compounds').loc[compound].reset_index().set_index('index').index
            compositions = compositions[findC]
            compounds = compounds[findC]
        # Constants
        R_g = 8314.46261815324  # Universal gas constant [(L·Pa)/(K·mol)]
        Hs_FW, notNaN_HsFW = eQ.solubilityHenry(compounds, 'FW', t)
        Hs_SW, notNaN_HsSW = eQ.solubilityHenry(compounds, 'SW', t)
        # Dictionaries initialization
        dict_Pi = {}
        dict_Ci_G = {}
        dict_Ci_LFW = {}
        dict_Ci_LSW = {}
        compounds = compounds.values
        for id_, composition in enumerate(compositions):
            # Gas phase - Partial pressure (Pi)
            Pi = p * composition # [Pa]
            # Gas phase - Gas concentration (Ci_G)
            Ci_G = (Pi / (R_g * t))
            # Liquid phase - Freshwater (Ci_LFW)
            if notNaN_HsFW[id_]:
                Ci_LFW = Pi * Hs_FW[..., id_] * (1/1000) # [mol/L]
            else:
                Ci_LFW = None
            # Liquid phase - Seawater (Ci_LSW)
            if notNaN_HsSW[id_]:
                Ci_LSW = Pi * Hs_SW[..., id_] * (1/1000) # [mol/L]
            else:
                Ci_LSW = None
            # Save data in dictionary
            dict_Pi[compounds[id_]] = Pi
            dict_Ci_G[compounds[id_]] = Ci_G
            dict_Ci_LFW[compounds[id_]] = Ci_LFW
            dict_Ci_LSW[compounds[id_]] = Ci_LSW
        if phase == 'G':
            return dict_Pi, dict_Ci_G
        elif phase == 'L-FW':
            return dict_Ci_LFW
        elif phase == 'L-SW':
            return dict_Ci_LSW
        elif phase == 'L':
            return dict_Ci_LFW, dict_Ci_LSW
        elif phase == 'All':
            return dict_Pi, dict_Ci_G, dict_Ci_LFW, dict_Ci_LSW
        else:
            print('!EcosysEM.Error: No phase selected. Use one of the following string:\n'+
                  '                             \'G\'       - Gas.\n'+
                  '                             \'L-FW\'    - Liquid fresh water.\n'+
                  '                             \'L-SW\'    - Liquid sea water.\n'+
                  '                             \'L\'       - Both liquid phases (L-FW, L-SW).\n'+
                  '                             \'All\'     - All phases (G, L-FW, L-SW).')
            sys.exit()
            
    def setComposition(self, compound, composition):
        """
        Modify composition of existing compounds or add new components with
        their respective compositions.

        Parameters
        ----------
        compound : STR or LIST
            Set of new and/or existing compounds.
        composition : FLOAT or LIST
            Composition(s) of new and/or existing compound(s).
        
        """
        if not isinstance(compound, list): compound = [compound]
        if not isinstance(composition, list): composition = [composition]
        pre_comp = self.compositions
        new_comp = dict(zip(compound, composition))
        pre_comp.update(new_comp)
        self.compositions = pre_comp
    
    def selectAltitude(self, selAlt):
        """
        Select a specific region of atmosphere based on a minimum and maximum 
        altitud value.
        
        """
        if isinstance(selAlt, list):
            if len(selAlt) > 2:
                print('!EcoSysEM.Error: `selAlt` argument must be a int, float or list [minAlt, maxAlt].')
                sys.exit()
            elif len(selAlt) == 1:
                selAlt = [0, selAlt[0]]
            minAlt = selAlt[0]
            maxAlt = selAlt[1]
        elif isinstance(selAlt, (int, float)):
            minAlt = 0
            maxAlt = selAlt
        # Previous altitude, temperature and pressure
        prevAlt = self.ISAaltitude
        prevT = self.ISAtemperature
        prevP = self.ISApressure
        # Correct min and max altitude, if out of previous range
        minAlt = max(minAlt, min(prevAlt))
        maxAlt = min(maxAlt, max(prevAlt))
        # Indexes of min and max altitudes
        imAlt = int(np.argwhere(prevAlt <= minAlt)[-1])
        iMAlt = int(np.argwhere(prevAlt >= maxAlt)[0]) + 1
        # New altitude, temperature and pressure
        self.ISAaltitude = prevAlt[imAlt:iMAlt]
        self.ISAtemperature = prevT[imAlt:iMAlt]
        self.ISApressure = prevP[imAlt:iMAlt]
    
    ## Plotting functions 
    def plotTandP_ISA(self):
        """
        Plotting of pressure (in atm) and temperature (in Kelvin) along the
        atmosphere (altitude in km).
        
        """
        # Variables
        alt = self.ISAaltitude / 1000      # [km]
        t = self.ISAtemperature + 273.15   # [K]
        p = self.ISApressure / 101325      # [atm]
        # Temperature
        fig, ax1 = plt.subplots(figsize = (4.2, 2))
        ax1.set_xlabel('Altitude (km)')
        ax1.set_ylabel('Temperature (K)', color = 'tab:red')
        ax1.set_ylim([150, 300])
        ax1.set_xlim([0, alt[-1]])
        ax1.plot(alt, t, color = 'tab:red')
        # Pressure
        ax2 = ax1.twinx()
        ax2.set_ylabel('Pressure (atm)', color = 'tab:blue')
        ax2.set_ylim([0, 1])
        ax2.plot(alt, p, color = 'tab:blue')
        # Plot properties and showing
        fig.tight_layout()
        plt.show()
    
    def plotCompsProfilesISA(self, C, xLabel, logCLabel = False, compounds = None):
        """
        Plotting of composition profiles along the atmosphere (altitude in km).
        
        """
        # Variables
        alt = self.ISAaltitude / 1000     # [km]
        # Plotting
        fig, ax = plt.subplots(figsize = (5,10))
        ax.set_ylabel('Altitude (km)')
        ax.set_xlabel(xLabel)
        ax.plot(C, alt)
        if logCLabel == True:
            ax.set_xscale('log')
        if not compounds:
            compoundsName = self.compounds
        else:
            compoundsName = compounds
        plt.legend(compoundsName, loc = 'center left', bbox_to_anchor = (1, 0.5))
        plt.show()

class MERRA2:
    """
    The Modern-Era Retrospective analysis for Research and Applications, 
    Version 2 (MERRA-2) provides data beginning in 1980. It was introduced 
    to replace the original MERRA dataset because of the advances made in 
    the assimilation system that enable assimilation of modern hyperspectral 
    radiance and microwave observations, along with GPS-Radio Occultation 
    datasets. It also uses NASA's ozone profile observations that began in 
    late 2004. Additional advances in both the GEOS model and the GSI 
    assimilation system are included in MERRA-2. Spatial resolution remains 
    about the same (about 50 km in the latitudinal direction) as in MERRA.
    --------------------------------------------------------------------------
    References: - Global Modeling and Assimilation Office (NASA)
                    https://gmao.gsfc.nasa.gov/reanalysis/merra-2/
                - Gelaro et al. (2017), doi: 10.1175/JCLI-D-16-0758.1
    """
    
    def _testMERRA2(self):
        print('MERRA2 class works.')

    def _closestCoord(self, lonR, latR, Coor):
        """
        Get closest coordinates from coordinates system of MERRA-2

        Parameters
        ----------
        lonR : LIST
            Reference longitude.
        latR : LIST
            Reference latitude.
        Coor : TUPLE
            Desired coordinates.
            (lower_left_lon, lower_left_lat, upper_right_lon, upper_right_lat)

        Returns
        -------
        bbox : TUPLE
            Closest coordinate from reference coordinates.
            (lower_left_lon, lower_left_lat, upper_right_lon, upper_right_lat)

        """
        bbox = (lonR[np.abs(np.array(lonR) - np.array(Coor[0])).argmin()],
                latR[np.abs(np.array(latR) - np.array(Coor[1])).argmin()],
                lonR[np.abs(np.array(lonR) - np.array(Coor[2])).argmin()],
                latR[np.abs(np.array(latR) - np.array(Coor[3])).argmin()])
        
        # bbox = (lonR[np.abs(lonR - Coor[0]).argmin()],
        #         latR[np.abs(latR - Coor[1]).argmin()],
        #         lonR[np.abs(lonR - Coor[2]).argmin()],
        #         latR[np.abs(latR - Coor[3]).argmin()])
        return bbox
    
    def _HfromP(self, P):
        """
        Get atmospheric altitude from pressure altitude (from NOAA).

        Parameters
        ----------
        P : INT, FLOAT or ndarray
            Pressure(s) values. [Pa]

        Returns
        -------
        INT, FLOAT or ndarray.
            Atmopspheric altitude value(s).

        """
        return np.maximum(145366.45 * (1 - (P / 101325) ** 0.190284) * (0.3048 / 1), np.zeros(P.shape))
    
    def _LR(self, T1, T2, H1, H2):
        """
        Estimate atmospheric lapse rate (K/km)

        Parameters
        ----------
        T1 : FLOAT, LIST or ndarray
            Temperature value(s) of minimum atmospheric altitude (H1). [K or °C]
        T2 : FLOAT, LIST or ndarray
            Temperature value(s) of maximum atmospheric altitude (H2). [K or °C]
        H1 : FLOAT, LIST or ndarray
            Minimum atmospheric altitude. [m]
        H2 : FLOAT, LIST or ndarray
            Maximum atmospheric altitude. [m]

        Returns
        -------
        FLOAT, LIST or ndarray
            Lapse rate value(s). (K/km)

        """
        return ((T2 - T1) / (H2 - H1) * 1000)
    
    def _saveNPZMERRA2(self, data, dataType, y, m, d = None):
        """
        Create .npz file with downladed data.

        Parameters
        ----------
        date : DICT
            Data in dictionary form.
        dataType : STR ('dly', 'mly' or 'cmly')
            Type of data.
        y : INT or LIST of INT
            Year(s) of data.
        m : INT or LIST of INT
            Month of data  
        
        Returns
        -------
        None.
        
        """
        # Path check and folder creation if this does not exist
        path = f'data/MERRA2/{dataType}/'
        if not os.path.isdir(path):
            os.makedirs(path)
        # File name (based on dataType)
        if dataType == 'dly':
            if not isinstance(y, int):
                print('\n!EcoSysEM.Error: argument \'y\' must be a integer')
                sys.exit()
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
                sys.exit()
            if not isinstance(d, int):
                print('\n!EcoSysEM.Error: argument \'d\' must be a integer')
                sys.exit()
            file = f'{y}_{m}_{d}_day.npz'
        elif dataType == 'mly':
            if not isinstance(y, int):
                print('\n!EcoSysEM.Error: argument \'y\' must be a integer')
                sys.exit()
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
                sys.exit()
            file = f'{y}_{m}_month.npz'
        elif dataType == 'cmly':
            if not isinstance(y, list):
                print('\n!EcoSysEM.Error: argument \'y\' must be a list: [start_year, end_year]')
                sys.exit()
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
                sys.exit()
            file = f'{y[0]}_{y[-1]}_{m}.npz'
        # Path generation
        pathfile = path + file
        # Save .npz file
        np.savez(pathfile, **data)
    
    def _openNPZMERRA2(self, dataType, y, m, d = None):
        """
        Open .npz file with downladed data.

        Parameters
        ----------
        date : DICT
            Data in dictionary form.
        dataType   : STR ('dly', 'mly' or 'cmly')
            Type of data.
        y : INT or LIST of INT
            Year(s) of data.
        m : INT or LIST of INT
            Month of data  
        
        Returns
        -------
        None.
        
        """
        path = f'data/MERRA2/{dataType}/'
        if dataType == 'dly':
            if not isinstance(y, int):
                print('\n!EcoSysEM.Error: argument \'y\' must be a integer')
                sys.exit()
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
                sys.exit()
            if not isinstance(d, int):
                print('\n!EcoSysEM.Error: argument \'d\' must be a integer')
                sys.exit()
            file = f'{y}_{m}_{d}_day.npz'
        if dataType == 'mly':
            if not isinstance(y, int):
                print('\n!EcoSysEM.Error: argument \'y\' must be a integer')
                sys.exit()
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
                sys.exit()
            file = f'{y}_{m}_month.npz'
        elif dataType == 'cmly':
            if not isinstance(y, list):
                print('\n!EcoSysEM.Error: argument \'y\' must be a list: [start_year, end_year]')
                sys.exit()
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
                sys.exit()
            file = f'{y[0]}_{y[-1]}_{m}.npz'
        return np.load(path + file)
    
    def _combDataMERRA2(self, dataType, year, month, days = None, keys = 'All', mlyDelete = False):
        """
        Get average and standard deviation from a group of data.
        
        Parameters
        ----------
        years : INT or LIST of INT
            Year(s) of data.
        month : INT
            Month of data
        dataType   : STR ('mly', 'dly')
            Type of data.
        keys  : LIST of STR
            List of requested variables. (Default: 'All')     
        
        Returns
        -------
        None
        
        """
        if not isinstance(month, int):
            print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
            sys.exit()
        if dataType == 'mly':
            if isinstance(year, int): year = [year]
            if isinstance(year, float): year = [year]
            if len(year) <= 1:
                print('\n!EcoSysEM.Error: Introduce at least 2 years to combine data.')
                sys.exit()
        # Get all files from data\npz
        folder = f'data/MERRA2/{dataType}/'
        allFiles = np.array(os.listdir(folder))
        # Test elements
        tEl = np.empty((0))
        if dataType == 'mly':
            for y in year:    
                el = f'{y}_{month}_month.npz'
                tEl = np.append(tEl, el)
        if dataType == 'dly':
            if not isinstance(year, int):
                print('\n!EcoSysEM.Error: `year` argument must be a integer.')
                sys.exit()
            if not isinstance(month, int):
                print('\n!EcoSysEM.Error: `month` argument must be a integer.')
                sys.exit()
            y = year
            m = month
            if days == 'All':
                last_day = calendar.monthrange(y, m)[1]
                first_day = 1
            else:
                if not isinstance(days, list): days = [days]
                days = [calendar.monthrange(y, m)[1] if d == -1 else d for d in days]
                last_day = max(days)
                if len(days) > 1:
                    first_day = min(days)
                else:
                    first_day = 1
            days = np.arange(first_day, last_day + 1, step = 1)  
            for d in days:
                el = f'{y}_{month}_{d}_day.npz'
                tEl = np.append(tEl, el)
        selFiles = allFiles[np.isin(allFiles, tEl)]
        # Stack matrices
        combData = {}
        monthData = {}
        for file in selFiles:
            path = folder + file
            f = np.load(path)
            # Select keys
            if keys == 'All':
                keys = f.files
            else:
                coor = []
                if not np.any(np.char.find('lat', keys) > -1):
                    coor += ['lat']
                if not np.any(np.char.find('lon', keys) > -1):
                    coor += ['lon']
                keys = coor + keys
            for key in keys:
                if file == selFiles[0]:
                    combData[key] = f[key]
                else:
                    combData[key] = np.dstack((combData[key], f[key]))
            f.close()
        # Monthly average and std
        for key in keys:
            if np.char.find(key, '_std') == -1:
                monthData[key] = np.nanmean(combData[key], axis = -1)
            else:
                monthData[key] = np.nanstd(combData[key], axis = -1)
            if key == 'lat' or key == 'lon':
                monthData[key] = np.squeeze(monthData[key])
        # Save numpy matrices in .npz format (v2)
        if dataType == 'mly':
            MERRA2._saveNPZMERRA2(self, data = monthData, dataType = 'cmly', y = [year[0], year[-1]], m = month)
        elif dataType == 'dly':
            MERRA2._saveNPZMERRA2(self, data = monthData, dataType = 'mly', y = y, m = m)
        # Delete monthly data (if necessary)
        if mlyDelete:
            for file in selFiles:
                path = folder + file
                os.remove(path)
            
    def getDataMERRA2(self, dataType, years, months,
                      days = 'All',
                      product = 'M2I1NXASM',
                      version = '5.12.4',
                      bbox = (-180, -90, 180, 90),
                      var = ['PS', 'TROPPB', 'T2M', 'TROPT', 'TROPH', 'LR']):
        """
        Download data from MERRA2 database.
        
        Parameters
        ----------
        dataType: STR or LIST of STR ('dly', 'mly', 'cmly' or 'All')
            Type of data
        years : INT or LIST of INT
            Year(s) of data.
        months : INT or LIST of INT
            Month(s) of data
        days : INT, LIST of INT, or STR ('All')
            Day(s) of data
        product  : STR
            Product of data (section of MERRA2 database).
        version : STR
            Version of data.
        bbox : TUPLE
            Earths region of data, the bounding box.
            (lower_left_longitude, lower_left_latitude, upper_right_longitude, upper_right_latitude)
        var : LIST of STR
            List of requested variables.   
        
        Returns
        -------
        None.
        
        """
        # Initialize `mlyDelete` for _combDataMERRA2()
        mlyDelete = False
        # Data checking
        if isinstance(dataType, str): dataType = [dataType]
        if not np.all(np.isin(dataType, ['dly', 'mly', 'cmly', 'All'])):
            print('\n!EcoSysEM.Error: dataType not found. Data type must be "dly", "mly", "cmly", list of data types or "All".')
            sys.exit()
        if np.any(np.isin(dataType, 'All')): dataType = ['dly', 'mly', 'cmly']
        if np.any(np.isin(dataType, 'cmly')) and not np.any(np.isin(dataType, 'mly')): 
            dataType += ['mly']
            mlyDelete = True
        if np.any(np.isin(dataType, 'All')): dataType = ['dly', 'mly', 'cmly']
        # Coordinates (bbox)
        if len(bbox) == 2:
            bbox = bbox + bbox
        elif len(bbox) == 4:
            bbox = bbox
        else:
            print('\n!EcoSysEM.Error: boundaries() requires 2 `(lon, lat)` or 4 `(lower_left_lon, lower_left_lat, upper_right_lon, upper_right_lat)` positional arguments.')
            sys.exit()
        # Check arguments
        if not isinstance(years, list): years = [years]
        if not isinstance(months, (list, np.ndarray)): months = [months]
        # Sort years and months
        years = sorted(years)
        months = sorted(months)
        if np.any(np.isin(dataType, 'cmly')) and len(years) <= 1:
            print('\n!EcoSysEM.Error: Introduce at least 2 years to combine monthly data.')
            sys.exit()
        start_time = time.time()
        # This will work if Earthdata prerequisite files have already been generated
        earthaccess.login()
        # Combinations [(year1, month1), (year1, month2), ...]
        dateList = list(itertools.product(years, months))
        # Initialize dictionaries
        monthData = {}
        for date in dateList:
            y = date[0]
            m = date[1]
            if days == 'All':
                last_day = calendar.monthrange(y, m)[1]
                first_day = 1
            else:
                if not isinstance(days, list): days = [days]
                days = [calendar.monthrange(y, m)[1] if d == -1 else d for d in days]
                last_day = max(days)
                if len(days) > 1:
                    first_day = min(days)
                else:
                    first_day = 1
            start_date = datetime.date(y, m, first_day)
            end_date = datetime.date(y, m, last_day)
            delta = datetime.timedelta(days=1)
            date = start_date
            # Initialize dictionaries
            hourData = {}
            av_dayData = {}
            dayData = {}
            ## Day matrices
            while(date <= end_date):
                print(f"> {date.strftime('%Y-%m-%d')}")
                # Open granules to local path
                results = earthaccess.search_data(
                    short_name = product, 
                    version = version,
                    temporal = (date, date),
                    bounding_box = bbox
                )
                # Open granules using Xarray
                print(f'>> Product: {product}')
                fs = earthaccess.open(results)
                ds = xr.open_mfdataset(fs)
                lonR = ds.lon.to_numpy()
                latR = ds.lat.to_numpy()
                # Get closest coordinates
                bbox =  MERRA2._closestCoord(self, lonR, latR, bbox)
                lat_slice = slice(bbox[1], bbox[3])
                lon_slice = slice(bbox[0], bbox[2])
                ds = ds.sel(lat = lat_slice, lon = lon_slice)
                # PHIS [Surface geopotential height] (local)
                d_PHIS = np.load('data/MERRA2/PHIS.npz')
                d_PHIS = MERRA2.selectRegion(self, d_PHIS, bbox)
                H = d_PHIS['PHIS']
                H = np.where(H < 0, 0, H)
                H = np.repeat(H[np.newaxis, ...], 24, axis = 0)
                hourData['H'] = H
                lon = ds.lon.to_numpy()
                lat = ds.lat.to_numpy()
                # User variables
                for key in var:
                    try:
                        ds[key]
                    except:
                        print(f'  {key} estimated.')
                    else:
                        hourData[key] = ds[key].values
                        print(f'  {key} done.')
                # Atmospheric altitudes
                varH = ['TROPH']
                varReqH = ['TROPPB', 'TROPPT', 'TROPPV']
                if np.any([np.char.find(var, iVar) == 0 for iVar in varH]):
                    c = 0
                    for iV in var:
                        if np.any(np.char.find(iV, varReqH) > -1):
                            varNames = {'TROPPB': 'TROPH',
                                        'TROPPT': 'TROPH',
                                        'TROPPV': 'TROPH'}
                            hourData[varNames[iV]] = MERRA2._HfromP(self, hourData[iV])
                            c += 1
                    if c == 0:
                        print('\n!EcoSysEM.Error: missing required variable missing for atmospheric altitude(s) - \'TROPPB\', \'TROPPT\', or \'TROPPV\'.')
                        sys.exit()
                # Lapse rate
                varReqLR = ['T2M', 'TROPT', 'TROPH']
                if np.any(np.char.find(var, 'LR') > -1):
                    if np.all([np.any(np.char.find(var, iVarReq) > -1) for iVarReq in varReqLR]):
                        hourData['LR'] = MERRA2._LR(self, hourData['T2M'], hourData['TROPT'], hourData['H'], hourData['TROPH'])
                    else:
                        print('\n!EcoSysEM.Error: missing required variable missing for lapse rate - \'T2M\', \'TROPT\', or \'TROPH\'.')
                        sys.exit()
                # PHIS [Surface geopotential height]
                av_dayData['H'] = np.round(np.nanmean(H, axis = 0))
                if date == start_date:
                    dayData['H'] = av_dayData['H']
                else:
                    dayData['H'] = np.dstack((dayData['H'], av_dayData['H'])) 
                # User variables
                for iV in var:
                    # Daily averages and std
                    av_dayData[iV] = np.nanmean(hourData[iV], axis = 0)
                    av_dayData[f'{iV}_std'] = np.nanstd(hourData[iV], axis = 0)
                    # Day matrices
                    if date == start_date:
                        dayData[iV] = av_dayData[iV]
                    else:
                        dayData[iV] = np.dstack((dayData[iV], av_dayData[iV]))
                # Save daily data
                if np.any(np.isin(dataType, 'dly')):
                    av_dayData['lat'] = lat
                    av_dayData['lon'] = lon
                    MERRA2._saveNPZMERRA2(self, data = av_dayData, dataType = 'dly', y = y, m = m, d = int(date.day))
                # Next date
                date += delta
            # Month matrices
            monthData['lat'] = lat
            monthData['lon'] = lon
            # PHIS [Surface geopotential height]
            if last_day != 1:
                monthData['H'] = np.round(np.nanmean(dayData['H'], axis = -1))
            else:
                monthData['H'] = dayData['H']
            # User variables
            for iV in var:
                if last_day != 1:
                    monthData[iV] = np.nanmean(dayData[iV], axis = -1)
                    monthData[f'{iV}_std'] = np.nanstd(dayData[iV], axis = -1)
                else:
                    monthData[iV] = dayData[iV]
                    monthData[f'{iV}_std'] = 0.0 * np.ones(dayData[iV].shape)
            # Save numpy matrices in .npz format
            if np.any(np.isin(dataType, 'mly')):
                MERRA2._saveNPZMERRA2(self, data = monthData, dataType = 'mly', y = y, m = m)
        # Combine monthly data if user required ('cmly')
        if np.any(np.isin(dataType, 'cmly')):
            var += ['H']
            for m in months:
                MERRA2._combDataMERRA2(self, dataType = 'mly',
                                       year = years,
                                       month = m,
                                       days = days,
                                       keys = var, 
                                       mlyDelete = mlyDelete)
        print("--- %s seconds ---" % (time.time() - start_time))
    
    def loadDataMERRA2(self, dataType, y, m, d = None, keys = 'All'):
        """
        Get data in dictionary form.

        Parameters
        ----------
        dataType : STR ('mly', 'cmly', 'dly')
            Type of data
        y : INT or LIST of INT
            Year(s) of data
        m : INT or LIST of INT
            Month of data
        d : INT or LIST of INT
            Day(s) of data
        keys : LIST of STR
            List of requested variables. (Default: 'All')
        
        Returns
        -------
        dictVar : DICT
            Dictionary with requested variables.

        """
        npz = MERRA2._openNPZMERRA2(self, dataType, y, m, d)
        if keys == 'All':
            keys = MERRA2.keysMERRA2(self, dataType, y, m, d)
            dictVar = {key: npz[key] for key in keys}
        else:
            coor = []
            if not np.any(np.char.find('lat', keys) > -1):
                coor += ['lat']
            if not np.any(np.char.find('lon', keys) > -1):
                coor += ['lon']
            keys = coor + keys
            dictVar = {key: npz[key] for key in keys}
        npz.close()
        return dictVar    
    
    def selectRegion(self, data, bbox): # !!! (redundant)
        """
        Select specific region of Earth of downloaded data

        Parameters
        ----------
        data : DICT
            Data
        bbox : Tuple
            Earths region of data, the bounding box.
            (lower_left_longitude, lower_left_latitude, upper_right_longitude, upper_right_latitude)

        Returns
        -------
        dataSel : DICT
            Data from specific region.

        """
        # Selected coordinates
        if len(bbox) == 2:
            bbox = bbox + bbox
            uniqueCoor = True
        elif len(bbox) == 4:
            bbox = bbox
            uniqueCoor = False
        else:
            print('\n!EcoSysEM.Error: boundaries() requires 2 `(lon, lat)` or 4 `(lower_left_lon, lower_left_lat, upper_right_lon, upper_right_lat)` positional arguments.')
            sys.exit()
        # Check user coordinates
        if bbox[0] > bbox[2]:
            print('\n!EcoSysEM.Error: `upper_right_longitude` (bbox[2]) must be higher than `lower_left_longitude` (bbox[0])')
            sys.exit()
        if bbox[1] > bbox[3]:
            print('\n!EcoSysEM.Error: `upper_right_latitude` (bbox[3]) must be higher than `lower_left_latitude` (bbox[1])')
            sys.exit()
        # BBox from data
        lonR = data['lon']
        latR = data['lat']
        bboxData = (lonR[0], latR[0], lonR[-1], latR[-1])
        # Get closest coordinates
        bbox =  MERRA2._closestCoord(self, lonR, latR, bbox)
        # Check if requested region is in data
        cBBOX = [(bboxData[i] <= bbox[i] and bboxData[i+2] >= bbox[i+2]) for i in np.arange(2)]
        if all(cBBOX) == True:
            # Get indices
            # (lower_left_lon, lower_left_lat, upper_right_lon, upper_right_lat)
            idx = (int(np.argwhere(lonR == bbox[0])),
                   int(np.argwhere(latR == bbox[1])),
                   int(np.argwhere(lonR == bbox[2])) + 1,
                   int(np.argwhere(latR == bbox[3])) + 1)
            # Initialize dictionary of selected data
            dataSel = {}
            for var in data:
                if var == 'lon':
                    if uniqueCoor:
                        dataSel[var] = data[var][idx[0]]
                    else:
                        if idx[0] == idx[2]:
                            dataSel[var] = data[var][idx[0]]
                        else:
                            dataSel[var] = data[var][idx[0]:idx[2]]
                elif var == 'lat':
                    if uniqueCoor:
                        dataSel[var] = data[var][idx[1]]
                    else:
                        if idx[1] == idx[3]:
                            dataSel[var] = data[var][idx[1]]
                        else:
                            dataSel[var] = data[var][idx[1]:idx[3]]
                else:
                    if uniqueCoor:
                        dataSel[var] = data[var][idx[1],idx[0]]
                    else:
                        if idx[0] == idx[2]:
                            dataSel[var] = data[var][idx[1]:idx[3],idx[0]]
                        elif idx[1] == idx[3]:
                            dataSel[var] = data[var][idx[1],idx[0]:idx[2]]
                        else:
                            dataSel[var] = data[var][idx[1]:idx[3],idx[0]:idx[2]]
                if not isinstance(dataSel[var], (list, np.ndarray)): dataSel[var] = [dataSel[var]]
            return dataSel
        else:
            print('\n!EcoSysEM.Error: selected region is outside the data boundaries.')
            sys.exit()
    
    def getTPAlt(self, dataType, year, month, day = None, bbox = (-180, -90, 180, 90), altArray = None, num = 50):
        """
        Compute the change of temperature and pressure of the Earth's
        atmosphere over the range of altitudes. Based on ISA (ISO 2533:1975).
        
        Parameters
        ----------
        dataType : STR ('mly', 'cmly', 'dly')
            Type of data
        year : INT or LIST of INT
            Year(s) of data
        month : INT or LIST of INT
            Month of data
        day : INT or LIST of INT
            Day(s) of data
        bbox : TUPLE, optional
            Earths region of data, the bounding box.
            (lower_left_longitude, lower_left_latitude, upper_right_longitude, upper_right_latitude)
        altArray : LIST or np.ndarray, optional
            List of altitudes
        num : INT, optional
            Number of altitude steps to generate.

        Returns
        -------
        T : FLOAT, LIST or ndarray
            Temperature in function of altitude. [K]
        P : FLOAT, LIST or ndarray
            Pressure in function of altitude. [Pa]
        H : FLOAT, LIST or ndarray
            Altitude. [m]
        """
        # Check argument format and dimension of data
        d = MERRA2.loadDataMERRA2(self, dataType, year, month, day)
        d = MERRA2.selectRegion(self, d, bbox)
        TS = np.array(d['T2M'])
        LR = np.array(d['LR'])
        PS = np.array(d['PS'])
        TROPH = np.array(d['TROPH'])
        HS = np.array(d['H'])
        HS = np.where(HS < 0, 0, HS)
        # Constants
        R = 8.3144598                   # Universal gas constant [J/mol/K]
        g0 = 9.80665                    # Gravitational acceleration [m/s^2]
        M0 = 0.0289644                  # Molar mass of Earth's air
        # Altitude. Shape: (alt, lat, lon)
        if not isinstance(altArray, (list, np.ndarray)):
            HS_min = 0.0 * np.ones(HS.shape)
            TROPH_max = np.nanmax(TROPH) * 0.99 * np.ones(TROPH.shape)
            H = np.linspace(start = HS_min, stop = TROPH_max, num = num)        # [m]
        else:
            max_TROPH = np.nanmax(TROPH) * 0.99
            altChk = altArray
            if not isinstance(altArray, (list, np.ndarray)): altArray = np.array(altArray)
            altArray = np.where(altArray < max_TROPH, altArray, np.NaN)
            altArray = altArray[~np.isnan(altArray)]
            if np.any(altChk > max_TROPH):
                altArray = np.append(altArray, max_TROPH)
            y = TS.shape[0]
            x = TS.shape[1]
            z = len(altArray)
            H = np.tile(altArray, y * x).reshape((z, y, x), order = 'F')
        # User variables
        # 3D matrix creation (TS, LR, HS)
        TS = np.repeat(TS[np.newaxis, ...], H.shape[0], axis = 0)              # [K]
        LR = np.repeat(LR[np.newaxis, ...], H.shape[0], axis = 0) / 1000       # [K/m]
        HS = np.repeat(HS[np.newaxis, ...], H.shape[0], axis = 0)              # [m]
        # Temperature profile
        T = TS + LR * (H - HS)
        # Pressure profile
        P = PS * (1 + ((LR) / (TS)) * (H - HS)) ** (-(g0 * M0) / (R * LR))
        # Temperature
        T = np.where(H < HS, np.NaN, T)
        T = np.where(H > TROPH, np.NaN, T)
        # Pressure
        P = np.where(H < HS, np.NaN, P)
        P = np.where(H > TROPH, np.NaN, P)
        return T, P, H
    
    def keysMERRA2(self, dataType, y, m, d = None):
        """
        Get variable list of data.

        Parameters
        ----------
        dataType : STR ('dly', 'mly' or 'cmly')
            Type of data.
        y : INT or LIST of INT
            Year(s) of data.
        m : INT or LIST of INT
            Month of data  
        
        Returns
        -------
        keys  : LIST of STR
            Variable list of data.

        """
        npz = MERRA2._openNPZMERRA2(self, dataType, y, m, d)
        keys = npz.files
        npz.close()
        return keys
    
    def deleteKeyMERRA2(self, keys, dataType, y, m, d = None):
        """
        Delete variable(s) from data.

        Parameters
        ----------
        keys : LIST of STR
            List of requested variables.
        dataType : STR ('dly', 'mly' or 'cmly')
            Type of data.
        y : INT or LIST of INT
            Year(s) of data.
        m : INT or LIST of INT
            Month of data  

        """
        dictVar = MERRA2.dictMERRA2(self, dataType, y, m, d)
        for key in keys:
            dictVar.pop(key, None)
        # Save new .npz file
        MERRA2._saveNPZMERRA2(self, dictVar, dataType, y, m, d)

class CAMS:
    """
    The Copernicus Atmosphere Monitoring Service (CAMS) provides continuous 
    data and information on atmospheric composition, supporting a wide range 
    of applications from air quality monitoring and forecasting to climate 
    change assessment. It combines observations from satellites and in-situ 
    measurements with sophisticated numerical models to provide consistent 
    and quality-controlled data records. CAMS data typically encompasses 
    global and regional analyses and forecasts of reactive gases, greenhouse 
    gases (GHGs), aerosols, and stratospheric ozone. The data records span 
    from approximately 2003 onwards for analyses and reanalyses, with 
    forecast data available on a rolling basis.
    --------------------------------------------------------------------------
    Reference: Copernicus CAMS Website - https://atmosphere.copernicus.eu/
    
    """

    def __init__(self):
        self.base_path = 'data/CAMS/'

    def _testCAMS(self):
        print('CAMS class works.')

    def getDataCAMS(self, dataType, years, months, days = 'All',
                    dataset = None,
                    pressure_levels = [50, 100, 200, 400, 
                                       600, 800, 900, 1000],
                    variables = None,
                    bbox = [90, -180, -90, 180],
                    mode = None,
                    method = 'linear'):
        """
        Download data from CAMS Global Greenhouse Gas Forecasts database.

        Parameters
        ----------
        dataType : STR or LIST of STR
            Type(s) of data ('dly', 'mly').
        years : INT or LIST of INT
            Year(s) of data.
        months : INT or LIST of INT
            Month(s) of data.
        days : INT, LIST of INT, or STR ('All')
            Day(s) of data.
        dataset : STR
            CAMS dataset name (e.g., "cams-global-greenhouse-gas-forecasts",
                               "cams-global-ghg-reanalysis-egg4", 
                               "cams-global-atmospheric-composition-forecasts").
        pressure_levels : INT or LIST of INT
            A list of pressure levels to download (e.g., [100, 200, ..., 1000]).
        variables : STR or LIST of STR or None
            A list of variables to download (Allowed: "co", "co2", "ch4").
            If None, uses the dataset defaults.
        bbox : LIST
            Earth's region of data, the bounding box (e.g., [90, -180, -90, 180]).
            (upper_right_latitude, lower_left_longitude, lower_left_latitude, upper_right_longitude)
        mode : STR or None
            Mode for the download of data (Allowed: "add").
            If "add", adds variable(s) to downloaded data. If None, downloads new data.
        method : STR, optional
            Method of interpolation (default: 'linear').
        """
        # Normalize inputs
        dataTypes = [dataType] if isinstance(dataType, str) else list(dataType)
        years = [years] if isinstance(years, int) else list(years)
        months = [months] if isinstance(months, int) else list(months)
    
        if days == 'All':
            days_list = None
        else:
            days_list = [days] if isinstance(days, int) else list(days)
    
        if isinstance(pressure_levels, (int, float)):
            pressure_levels = [pressure_levels]
        pressure_levels = [str(int(pl)) for pl in pressure_levels]
    
        if dataset == "cams-global-greenhouse-gas-forecasts":
            support  = {"co": True, "co2": True, "ch4": True}
            req_name = {"co": "carbon_monoxide", "co2": "carbon_dioxide", "ch4": "methane"}
        elif dataset == "cams-global-ghg-reanalysis-egg4":
            support  = {"co": False, "co2": True, "ch4": True}
            req_name = {"co2": "carbon_dioxide", "ch4": "methane"}
        elif dataset == "cams-global-atmospheric-composition-forecasts":
            support  = {"co": True, "co2": False, "ch4": True}
            req_name = {"co": "carbon_monoxide", "ch4": "methane"}
        elif dataset == None:
            raise ValueError("Dataset should be defined (e.g., 'cams-global-greenhouse-gas-forecasts', 'cams-global-ghg-reanalysis-egg4', 'cams-global-atmospheric-composition-forecasts')")
        else:
            raise ValueError(f"Unknown dataset: {dataset}")
    
        # Variables to request
        allowed_keys = ("co", "co2", "ch4")
        if variables is None:
            selected_mols = []
            for mol, ok in support.items():
                if ok:
                    selected_mols.append(mol)
        else:
            if isinstance(variables, str):
                variables = [variables]
        
            selected_mols = []
            for v in variables:
                mol = str(v).lower()
                if mol in allowed_keys and mol in support and support[mol]:
                    selected_mols.append(mol)

            if not selected_mols:
                raise ValueError("Selected variables are not available for this dataset. Use 'co', 'co2', 'ch4'.")
    
        req_vars = [req_name[m] for m in selected_mols]
    
        raw_folder = os.path.join(self.base_path, 'temporary')
        os.makedirs(raw_folder, exist_ok=True)
    
        client = cdsapi.Client()
    
        for y in years:
            for m in months:
                # Day ranges
                month_max_day = calendar.monthrange(y, m)[1]
                if days_list is None:
                    day_ranges = [(1, month_max_day)]
                else:
                    invalid_days = [d for d in days_list if d < 1 or d > month_max_day]
                    if invalid_days:
                        raise ValueError(f"Invalid day(s) {invalid_days} for year {y}, month {m}")
                        sys.exit()
                    day_ranges = [(d, d) for d in days_list]
    
                for d_start, d_end in day_ranges:
                    date_range = f"{datetime.date(y, m, d_start)}/{datetime.date(y, m, d_end)}"
    
                    # Unique filenames
                    if d_start == d_end:
                        stamp = f"{y}_{m}_{d_start}"
                    else:
                        stamp = f"{y}_{m}_{d_start}-{d_end}"
    
                    zip_path  = os.path.join(raw_folder, f"CAMS_{stamp}.zip")
                    target_nc = os.path.join(raw_folder, f"CAMS_{stamp}.nc")
    
                    if dataset == "cams-global-greenhouse-gas-forecasts":
                        req = {
                            'variable': req_vars,
                            'pressure_level': pressure_levels,
                            'date': [date_range],
                            'data_format': 'netcdf_zip',
                            'area': bbox,
                            "leadtime_hour": ["0"]
                        }
                    elif dataset == "cams-global-ghg-reanalysis-egg4":
                        req = {
                            'variable': req_vars,
                            'pressure_level': pressure_levels,
                            'date': [date_range],
                            'data_format': 'netcdf_zip',
                            'area': bbox,
                            "step": ["0"]
                        }
                    elif dataset == "cams-global-atmospheric-composition-forecasts":
                        req = {
                            'variable': req_vars,
                            'pressure_level': pressure_levels,
                            'date': [date_range],
                            'data_format': 'netcdf_zip',
                            'area': bbox,
                            "time": ["00:00"],
                            "leadtime_hour": ["0"],
                            "type": ["forecast"]
                        }
                    else:
                        raise ValueError(f"Unknown dataset: {dataset}")
    
                    print(f"[{dataset}] Downloading {date_range} for {req_vars}")
                    try:
                        client.retrieve(dataset, req, zip_path)
    
                        # Extract .nc
                        with zipfile.ZipFile(zip_path, 'r') as zf:
                            nc_file = zf.namelist()[0]
                            temp_nc = target_nc + ".tmp"
                            with zf.open(nc_file) as zipped_nc, open(temp_nc, 'wb') as local_nc:
                                local_nc.write(zipped_nc.read())

                        os.replace(temp_nc, target_nc)
                        os.remove(zip_path)
                        print("Extracted .nc")
    
                        # Process and clean up
                        for dt in dataTypes:
                            self._processDataCAMS('mly', dataset=dataset, molecules=selected_mols, month=m, mode=mode, method=method)
                        self._deleteTempDataCAMS(target_nc)
    
                    except Exception as e:
                        print(f"Failed for {date_range}: {e}")
                        sys.exit()
        
        print("\nAll downloads completed.")

    def _processDataCAMS(self, dataType, dataset, molecules, month = None, mode = None, method = 'linear'):
        """
        Process CAMS .nc files and save as .npz format.

        Parameters
        ----------
        dataType : STR 
            Type of data ('dly', 'mly', or 'cmly'). 
        dataset : STR
            CAMS dataset name (e.g., "cams-global-greenhouse-gas-forecasts",
                               "cams-global-ghg-reanalysis-egg4", or 
                               "cams-global-atmospheric-composition-forecasts").
        molecules : STR or LIST of STR
            Molecules to process: 'co', 'co2', 'ch4'. 
        month : INT
            Input for cmly dataType.
        mode : STR or None
            Mode for the download of data (Allowed: "add").
            If "add", adds variable(s) to downloaded data. If None, downloads new data.
        method : STR, optional
            Method of interpolation (default: 'linear').
        """
        input_folder = os.path.join(self.base_path, 'temporary')
        processed_folder = os.path.join(self.base_path, f"{dataType}")
        os.makedirs(processed_folder, exist_ok=True)
    
        try:
            nc_files = [f for f in os.listdir(input_folder) if f.endswith('.nc')]
        except FileNotFoundError:
            print(f"Input folder not found: '{input_folder}'")
            sys.exit()
    
        if not nc_files:
            print(f"No .nc files found in: '{input_folder}'")
            sys.exit()
    
        # Variable names in netCDF data
        if dataset == "cams-global-greenhouse-gas-forecasts":
            variable_names = {
                "co":  ["co"],
                "co2": ["co2"],
                "ch4": ["ch4"]
            }
        elif dataset == "cams-global-ghg-reanalysis-egg4":
            variable_names = {
                "co2": ["co2"],
                "ch4": ["ch4"]
                # 'co' not available
            }
        elif dataset == "cams-global-atmospheric-composition-forecasts":
            variable_names = {
                "co":  ["co"],
                "ch4": ["ch4_c"]
                # 'co2' not available
            }
        else:
            raise ValueError(f"Unknown dataset: {dataset}")
    
        mols_to_process = [molecules] if isinstance(molecules, str) else list(molecules)
        mols_to_process = [m for m in mols_to_process if m in variable_names]
    
        if dataType == "dly":
            for nc_file in nc_files:
                path_file = os.path.join(input_folder, nc_file)
                print(f"Processing: {path_file}")
                try:
                    ds = xr.open_dataset(path_file)
    
                    # Time dimension
                    if ("forecast_reference_time" in ds.dims) or ("forecast_reference_time" in ds.coords):
                        time_dim = "forecast_reference_time"
                    elif ("time" in ds.dims) or ("time" in ds.coords):
                        time_dim = "time"
                    elif ("valid_time" in ds.dims) or ("valid_time" in ds.coords):
                        time_dim = "valid_time"
                    else:
                        ds.close()
                        raise ValueError("Time dimension not found ('forecast_reference_time'/'valid_time').")
    
                    # Levels in Pa
                    if "pressure_level" in ds:
                        pressure_pa = ds["pressure_level"].values * 100
                    elif "level" in ds:
                        pressure_pa = ds["level"].values * 100
                    elif "pressure_level" in ds.coords:
                        pressure_pa = ds.coords["pressure_level"].values * 100
                    elif "level" in ds.coords:
                        pressure_pa = ds.coords["level"].values * 100
                    else:
                        ds.close()
                        raise ValueError("Pressure levels not found (pressure_level/level).")
    
                    # Coords
                    if 'latitude' in ds:
                        lat = ds['latitude'].values
                    elif 'lat' in ds:
                        lat = ds['lat'].values
                    else:
                        ds.close()
                        raise ValueError("Latitude coord not found (latitude/lat).")
                    
                    if 'longitude' in ds:
                        lon = ds['longitude'].values
                    elif 'lon' in ds:
                        lon = ds['lon'].values
                    else:
                        ds.close()
                        raise ValueError("Longitude coord not found (longitude/lon).")
    
                    # Daily means
                    daily_by_molecule = {}
                    first_key = None
    
                    for mol in mols_to_process:
                        names = variable_names[mol]
                        if isinstance(names, str):
                            names = [names]
                        vname = next((nm for nm in names if nm in ds.data_vars), None)
                        if vname is None:
                            print(f"  - {mol.upper()} variable not found in file; skip.")
                            continue
    
                        da = ds[vname]
                        if "step" in da.dims:
                            da = da.isel(step=0)
                        if "leadtime" in da.dims:
                            da = da.isel(leadtime=0)
    
                        da_daily = da.resample({time_dim: '1D'}).mean()
                        daily_by_molecule[mol] = da_daily
                        if first_key is None:
                            first_key = mol
    
                    if not daily_by_molecule:
                        ds.close()
                        print("No variables found to process for daily product.")
                        continue
    
                    days_coord = daily_by_molecule[first_key][time_dim].values
    
                    # Extract year/month from filename
                    parts = os.path.basename(nc_file).split('_')
                    year_int = int(parts[1]); month_int = int(parts[2].split('.')[0])
    
                    for i_day in range(len(days_coord)):
                        data_dict = {}
                        for mol, da_daily in daily_by_molecule.items():
                            data_dict[mol.upper()] = da_daily.isel({time_dim: i_day}).values
    
                        data_dict.update({
                            'alt': MERRA2._HfromP(self, pressure_pa), # !!! [redudant]
                            'lat': lat,
                            'lon': lon,
                            'P_level': pressure_pa
                        })
    
                        # Make sure arrays are 3D
                        for key_out, val_out in list(data_dict.items()):
                            if isinstance(val_out, np.ndarray):
                                arr = np.asarray(val_out).squeeze()
                                if arr.ndim == 4 and 1 not in arr.shape:
                                    arr = arr.mean(axis=0)
                                data_dict[key_out] = arr
                        
                        if mode == 'add':
                            NPZ = CAMS._openNPZCAMS(self, dataType, year_int, month_int, (i_day+1)) # !!! [redudant]
                            existing = {k: NPZ[k] for k in NPZ.files}
                            targ_lats = existing['lat']
                            targ_lons = existing['lon']
                            targ_alts = existing['alt']
                            NPZ.close()
                            
                            orig_lats = lat
                            orig_lons = lon
                            
                            mesh_lat, mesh_lon = np.meshgrid(targ_lats, targ_lons, indexing='ij')
                            points = np.stack([mesh_lat.ravel(), mesh_lon.ravel()], axis=-1)
                            
                            vars_to_add = []
                            for m in (mols_to_process if isinstance(mols_to_process, list) else [mols_to_process]):
                                up = m.upper()
                                if up in data_dict: 
                                    vars_to_add.append(up)
                                if f'{up}_std' in data_dict: 
                                    vars_to_add.append(f'{up}_std')
                            
                            for var in (vars_to_add):
                                data3d = np.asarray(data_dict[var])
                                
                                if data3d.shape[0] != targ_alts.size:
                                    raise ValueError(f"Vertical level mismatch: Source={data3d.shape[0]} vs Target={targ_alts.size}")
                                    
                                out = np.empty((targ_alts.size, targ_lats.size, targ_lons.size), dtype=data3d.dtype)
                                for k in range(targ_alts.size):
                                    interp = RegularGridInterpolator(
                                        (orig_lats, orig_lons),
                                        data3d[k, :, :],
                                        method=method,
                                        bounds_error=False,
                                        fill_value=np.nan
                                    )
                                    out[k] = interp(points).reshape(targ_lats.size, targ_lons.size)
                                existing[var] = out
                            CAMS._saveNPZCAMS(self, existing, "dly", year_int, month_int, (i_day+1)) # !!! [redudant]
                            print(f"{mols_to_process} is added to the file.")
                        else:
                            CAMS._saveNPZCAMS(self, data_dict, "dly", year_int, month_int, (i_day+1)) # !!! [redudant]
    
                    ds.close()
    
                except Exception as e:
                    print(f"Error occurred: {nc_file}, {e}")
                    try:
                        ds.close()
                    except:
                        pass
                    sys.exit()
    
            print(f"{dataType} files are processed and saved.")
    
        elif dataType == "mly":
            target_list = nc_files
            for nc_file in target_list:
                path_file = os.path.join(input_folder, nc_file)
                print(f"Processing: {path_file}")
                try:
                    ds = xr.open_dataset(path_file)
    
                    # Time dimension
                    if ("forecast_reference_time" in ds.dims) or ("forecast_reference_time" in ds.coords):
                        time_dim = "forecast_reference_time"
                    elif ("time" in ds.dims) or ("time" in ds.coords):
                        time_dim = "time"
                    elif ("valid_time" in ds.dims) or ("valid_time" in ds.coords):
                        time_dim = "valid_time"
                    else:
                        ds.close()
                        raise ValueError("Time dimension not found ('forecast_reference_time'/'valid_time').")
    
                    # Levels in Pa
                    if "pressure_level" in ds:
                        pressure_pa = ds["pressure_level"].values * 100
                    elif "level" in ds:
                        pressure_pa = ds["level"].values * 100
                    elif "pressure_level" in ds.coords:
                        pressure_pa = ds.coords["pressure_level"].values * 100
                    elif "level" in ds.coords:
                        pressure_pa = ds.coords["level"].values * 100
                    else:
                        ds.close()
                        raise ValueError("Pressure levels not found (pressure_level/level).")
    
                    # Coords
                    if 'latitude' in ds:
                        lat = ds['latitude'].values
                    elif 'lat' in ds:
                        lat = ds['lat'].values
                    else:
                        ds.close()
                        raise ValueError("Latitude coord not found (latitude/lat).")
                    
                    if 'longitude' in ds:
                        lon = ds['longitude'].values
                    elif 'lon' in ds:
                        lon = ds['lon'].values
                    else:
                        ds.close()
                        raise ValueError("Longitude coord not found (longitude/lon).")
    
                    data_dict = {}
                    for mol in mols_to_process:
                        var_name = variable_names[mol]
                        vname = None
                        for name in var_name:
                            if name in ds.data_vars:
                                vname = name
                                break
                        if vname is None:
                            print(f"  - {mol.upper()} variable not found in file; skip.")
                            continue
    
                        da = ds[vname]
                        if "step" in da.dims:
                            da = da.isel(step=0)
                        if "leadtime" in da.dims:
                            da = da.isel(leadtime=0)
    
                        data_dict[mol.upper()] = da.mean(dim=time_dim).values
                        data_dict[f"{mol.upper()}_std"] = da.std(dim=time_dim).values
    
                    data_dict.update({
                        'alt': MERRA2._HfromP(self, pressure_pa), # !!! [redudant]
                        'lat': lat,
                        'lon': lon,
                        'P_level': pressure_pa
                    })
    
                    # Make sure arrays are 3D
                    for key_out, val_out in list(data_dict.items()):
                        if isinstance(val_out, np.ndarray):
                            arr = np.asarray(val_out).squeeze()
                            if arr.ndim == 4 and 1 not in arr.shape:
                                arr = arr.mean(axis=0)
                            data_dict[key_out] = arr
    
                    parts = os.path.basename(nc_file).split('_')
                    year_int = int(parts[1]); month_int = int(parts[2].split('.')[0])
                    if mode == 'add':
                        NPZ = CAMS._openNPZCAMS(self, dataType, year_int, month_int) # !!! [redudant]
                        existing = {k: NPZ[k] for k in NPZ.files}
                        targ_lats = existing['lat']
                        targ_lons = existing['lon']
                        targ_alts = existing['alt']
                        NPZ.close()
                        
                        orig_lats = lat
                        orig_lons = lon
                        
                        mesh_lat, mesh_lon = np.meshgrid(targ_lats, targ_lons, indexing='ij')
                        points = np.stack([mesh_lat.ravel(), mesh_lon.ravel()], axis=-1)
                        
                        vars_to_add = []
                        for m in (mols_to_process if isinstance(mols_to_process, list) else [mols_to_process]):
                            up = m.upper()
                            if up in data_dict: 
                                vars_to_add.append(up)
                            if f'{up}_std' in data_dict: 
                                vars_to_add.append(f'{up}_std')
                        
                        for var in (vars_to_add):
                            data3d = np.asarray(data_dict[var])
                            
                            if data3d.shape[0] != targ_alts.size:
                                raise ValueError(f"Vertical level mismatch: Source={data3d.shape[0]} vs Target={targ_alts.size}")
                                
                            out = np.empty((targ_alts.size, targ_lats.size, targ_lons.size), dtype=data3d.dtype)
                            for k in range(targ_alts.size):
                                interp = RegularGridInterpolator(
                                    (orig_lats, orig_lons),
                                    data3d[k, :, :],
                                    method=method,
                                    bounds_error=False,
                                    fill_value=np.nan
                                )
                                out[k] = interp(points).reshape(targ_lats.size, targ_lons.size)
                            existing[var] = out
                        CAMS._saveNPZCAMS(self, existing, dataType, year_int, month_int) # !!! [redudant]
                        print(f"{mols_to_process} is added to the file.")
                    else:
                        CAMS._saveNPZCAMS(self, data_dict, dataType, year_int, month_int) # !!! [redudant]
    
                    ds.close()
    
                except Exception as e:
                    print(f"Error occurred: {nc_file}, {e}")
                    try:
                        ds.close()
                    except:
                        pass
                    sys.exit()
    
            print(f"{dataType} file is processed and saved.")
    
    def combDataCAMS(self, dataType, years, months, method = 'linear',
                     target_lats = None, target_lons = None):
        """
        Combine data as 'cmly', 'yly' or 'cyly'.

        Parameters
        ----------
        dataType : STR 
            Type of data ('cmly', 'yly', or 'cyly'). 
        years : INT or LIST of INT
            Year(s) of data.
        months : INT or LIST of INT
            Month(s) of data.
        method : STR, optional
            Method of interpolation (default: 'linear').
        target_lats : 1D array or None optional
            Desired latitudes for the CAMS grid (e.g.: np.arange(-90, 90.1, 0.5)).
            (default: None; uses last year/month grid).
        target_lons : 1D array, optional.
            Desired longitudes for the CAMS grid (e.g.: np.arange(-180, 179.375+0.001, 0.625)). 
            (default: None; uses last year/month grid).
        """ 
        
        out_dir = f"data/CAMS/{dataType}"
        os.makedirs(out_dir, exist_ok=True)
        
        years = [years] if isinstance(years, int) else list(years)
        months = [months] if isinstance(months, int) else list(months)
        source_type = 'mly'
        
        if dataType == 'cmly':
            opened = {}
            for m in months:
                opened[m] = {}
                for y in years:
                    try:
                        npz = self._openNPZCAMS(source_type, y, m) # !!! [redudant]
                        opened[m][y] = npz
                        print(f"[OK] opened: data/CAMS/{source_type}/{y}_{m}_month.npz")
                    except FileNotFoundError:
                        continue
            
            for m in months:
                missing_years = [y for y in years if y not in opened.get(m, {})]
                if missing_years:
                    raise RuntimeError(f"[ERROR] month {m}: missing source files in 'mly' for years {missing_years}.")
            
            print("Regridding for data combining...")
            regridded = {} 
    
            for m, yrs in opened.items():
                # Last year that requested for this month
                base_year = next((y for y in reversed(years) if y in yrs), None)
            
                base_npz  = yrs[base_year]
                if target_lats is not None:
                    targ_lats = np.flip(target_lats)
                else:
                    targ_lats = np.asarray(base_npz["lat"])
                if target_lons is not None:
                    targ_lons = target_lons
                else:
                    targ_lons = np.asarray(base_npz["lon"])
            
                targ_levels = np.asarray(base_npz["P_level"])
                level_key = "P_level"
            
                LAT, LON = np.meshgrid(targ_lats, targ_lons, indexing="ij")
                points = np.stack([LAT.ravel(), LON.ravel()], axis=-1)
            
                default_keys = {"lat", "lon", "P_level", "alt"}
                regridded[m] = {}
                for y, npz in yrs.items():
                    orig_lats = np.asarray(npz["lat"])
                    orig_lons = np.asarray(npz["lon"])
            
                    out_dict = {
                        "lat": targ_lats,
                        "lon": targ_lons,
                    }
                    out_dict[level_key] = targ_levels
                    out_dict["alt"] = np.asarray(npz["alt"])
            
                    for key in npz.files:
                        if key in default_keys:
                            continue
                        data = np.asarray(npz[key])
                        if data.ndim != 3:
                            raise ValueError(f"Expected 3D (alt,lat,lon) for '{key}', got shape={data.shape} in {y}-{m}")
                        if data.shape[0] != np.size(targ_levels):
                            raise ValueError(f"Vertical level mismatch for '{key}' in {y}-{m}: "f"source={data.shape[0]} vs target={np.size(targ_levels)}")
                        out = np.empty((np.size(targ_levels), targ_lats.size, targ_lons.size), dtype=data.dtype)
                        for klevel in range(np.size(targ_levels)):
                            interp = RegularGridInterpolator(
                                (orig_lats, orig_lons), data[klevel, :, :],
                                method=method, bounds_error=False, fill_value=np.nan
                            )
                            out[klevel] = interp(points).reshape(targ_lats.size, targ_lons.size)
                        out_dict[key] = out
            
                    regridded[m][y] = out_dict
                    
                    if hasattr(npz, "close"):
                        npz.close()
            
            print("Combining the data...")
            
            for m, yrs in regridded.items():        
                years_used = sorted(yrs.keys())
                sample = next(iter(yrs.values()))  # schema
                
                targ_lats = sample["lat"]
                targ_lons = sample["lon"]
            
                output_data = {"lat": targ_lats,
                               "lon": targ_lons}
                output_data["P_level"] = sample["P_level"]
                output_data["alt"] = sample["alt"]
            
                axis_keys = {"lat", "lon", "P_level", "alt"}
                all_keys  = set().union(*[d.keys() for d in yrs.values()])
                base_vars = [k for k in all_keys if (k not in axis_keys) and (not k.endswith("_std"))]
            
                for var in base_vars:
                    arr_list = [yrs[y][var] for y in years_used if var in yrs[y]]
                    if not arr_list:
                        continue
                    stack = np.stack(arr_list, axis=0).astype(float)  # (nyears, lev, lat, lon)
                    output_data[var] = np.nanmean(stack, axis=0)
                    output_data[f"{var}_std"] = np.nanstd(stack, axis=0)
                    
                    missing = [y for y in years_used if var not in yrs[y]]
                    if missing:
                        print(f"[INFO] month {m}: '{var}' averaged over {stack.shape[0]}/{len(years_used)} years; missing: {missing}")
                
                CAMS._saveNPZCAMS(self, data=output_data, dataType='cmly', y=years, m=m) # !!! [redudant]

        elif dataType == 'yly':
            opened = {}
            for y in years:
                opened[y] = {}
                for m in months:
                    try:
                        npz = self._openNPZCAMS(source_type, y, m) # !!! [redudant]
                        opened[y][m] = npz
                        print(f"[OK] opened: data/CAMS/{source_type}/{y}_{m}_month.npz")
                    except FileNotFoundError:
                        continue
            
            for y in years:
                missing_months = [m for m in months if m not in opened.get(y, {})]
                if missing_months:
                    raise RuntimeError(f"[ERROR] year {y}: missing source files in 'mly' for months {missing_months}.")
            
            print("Regridding for data combining...")
            regridded = {} 
    
            for y, mhs in opened.items():
                # Last year that requested for this month
                base_month = next((m for m in reversed(months) if m in mhs), None)
            
                base_npz  = mhs[base_month]
                if target_lats is not None:
                    targ_lats = np.flip(target_lats)
                else:
                    targ_lats = np.asarray(base_npz["lat"])
                if target_lons is not None:
                    targ_lons = target_lons
                else:
                    targ_lons = np.asarray(base_npz["lon"])
            
                targ_levels = np.asarray(base_npz["P_level"])
                level_key = "P_level"
            
                LAT, LON = np.meshgrid(targ_lats, targ_lons, indexing="ij")
                points = np.stack([LAT.ravel(), LON.ravel()], axis=-1)
            
                default_keys = {"lat", "lon", "P_level", "alt"}
                regridded[y] = {}
                for m, npz in mhs.items():
                    orig_lats = np.asarray(npz["lat"])
                    orig_lons = np.asarray(npz["lon"])
            
                    out_dict = {
                        "lat": targ_lats,
                        "lon": targ_lons,
                    }
                    out_dict[level_key] = targ_levels
                    out_dict["alt"] = np.asarray(npz["alt"])
            
                    for key in npz.files:
                        if key in default_keys:
                            continue
                        data = np.asarray(npz[key])
                        if data.ndim != 3:
                            raise ValueError(f"Expected 3D (alt,lat,lon) for '{key}', got shape={data.shape} in {y}-{m}")
                        if data.shape[0] != np.size(targ_levels):
                            raise ValueError(f"Vertical level mismatch for '{key}' in {y}-{m}: "f"source={data.shape[0]} vs target={np.size(targ_levels)}")
                        out = np.empty((np.size(targ_levels), targ_lats.size, targ_lons.size), dtype=data.dtype)
                        for klevel in range(np.size(targ_levels)):
                            interp = RegularGridInterpolator(
                                (orig_lats, orig_lons), data[klevel, :, :],
                                method=method, bounds_error=False, fill_value=np.nan
                            )
                            out[klevel] = interp(points).reshape(targ_lats.size, targ_lons.size)
                        out_dict[key] = out
            
                    regridded[y][m] = out_dict
                    
                    if hasattr(npz, "close"):
                        npz.close()
            
            print("Combining the data...")
            
            for y, mhs in regridded.items():        
                months_used = sorted(mhs.keys())
                sample = next(iter(mhs.values()))  # schema
                
                targ_lats = sample["lat"]
                targ_lons = sample["lon"]
            
                output_data = {"lat": targ_lats,
                               "lon": targ_lons}
                output_data["P_level"] = sample["P_level"]
                output_data["alt"] = sample["alt"]
            
                axis_keys = {"lat", "lon", "P_level", "alt"}
                all_keys  = set().union(*[d.keys() for d in mhs.values()])
                base_vars = [k for k in all_keys if (k not in axis_keys) and (not k.endswith("_std"))]
            
                for var in base_vars:
                    arr_list = [mhs[m][var] for m in months_used if var in mhs[m]]
                    if not arr_list:
                        continue
                    stack = np.stack(arr_list, axis=0).astype(float)  # (nmonths, lev, lat, lon)
                    output_data[var] = np.nanmean(stack, axis=0)
                    output_data[f"{var}_std"] = np.nanstd(stack, axis=0)
                    
                    missing = [mm for mm in months_used if var not in mhs[mm]]
                    if missing:
                        print(f"[INFO] year {y}: '{var}' averaged over {stack.shape[0]}/{len(months_used)} months; missing: {missing}")
                
                CAMS._saveNPZCAMS(self, data=output_data, dataType='yly', y=y) # !!! [redudant]
                
        elif dataType == 'cyly':
            opened = {}
            for y in years:
                opened[y] = {}
                for m in months:
                    try:
                        npz = self._openNPZCAMS(source_type, y, m) # !!! [redudant]
                        opened[y][m] = npz
                        print(f"[OK] opened: data/CAMS/{source_type}/{y}_{m}_month.npz")
                    except FileNotFoundError:
                        continue

            for y in years:
                missing_months = [m for m in months if m not in opened.get(y, {})]
                if missing_months:
                    raise RuntimeError(f"[ERROR] year {y}: missing source files in 'mly' for months {missing_months}.")

            print("Regridding for data combining...")
            regridded = {} 
            
            global_base_year = next((yy for yy in reversed(years) if opened.get(yy)), None)
            if global_base_year is None:
                raise RuntimeError("[ERROR] cyly: no source files found.")
            global_base_month = next((mm for mm in reversed(months) if mm in opened[global_base_year]), None)
            if global_base_month is None:
                raise RuntimeError(f"[ERROR] cyly: base year {global_base_year} has none of the requested months.")

            base_npz  = opened[global_base_year][global_base_month]
            if target_lats is not None:
                targ_lats = np.flip(target_lats)
            else:
                targ_lats = np.asarray(base_npz["lat"])
            if target_lons is not None:
                targ_lons = target_lons
            else:
                targ_lons = np.asarray(base_npz["lon"])
            targ_levels = np.asarray(base_npz["P_level"])
            level_key = "P_level"

            LAT, LON = np.meshgrid(targ_lats, targ_lons, indexing="ij")
            points = np.stack([LAT.ravel(), LON.ravel()], axis=-1)
            
            for y, mhs in opened.items():

                default_keys = {"lat", "lon", "P_level", "alt"}
                regridded[y] = {}
                for m, npz in mhs.items():
                    orig_lats = np.asarray(npz["lat"])
                    orig_lons = np.asarray(npz["lon"])

                    out_dict = {"lat": targ_lats, "lon": targ_lons}
                    out_dict[level_key] = targ_levels
                    out_dict["alt"] = np.asarray(npz["alt"])

                    for key in npz.files:
                        if key in default_keys:
                            continue
                        data = np.asarray(npz[key])
                        if data.ndim != 3:
                            raise ValueError(f"Expected 3D (alt,lat,lon) for '{key}', got shape={data.shape} in {y}-{m}")
                        if data.shape[0] != np.size(targ_levels):
                            raise ValueError(f"Vertical level mismatch for '{key}' in {y}-{m}: "f"source={data.shape[0]} vs target={np.size(targ_levels)}")
                        out = np.empty((np.size(targ_levels), targ_lats.size, targ_lons.size), dtype=data.dtype)
                        for klevel in range(np.size(targ_levels)):
                            interp = RegularGridInterpolator(
                                (orig_lats, orig_lons), data[klevel, :, :],
                                method=method, bounds_error=False, fill_value=np.nan
                            )
                            out[klevel] = interp(points).reshape(targ_lats.size, targ_lons.size)
                        out_dict[key] = out

                    regridded[y][m] = out_dict
                    if hasattr(npz, "close"):
                        npz.close()

            print("Combining the data...")

            annual_map = {}
            for y, mhs in regridded.items():        
                months_used = sorted(mhs.keys())
                if not months_used:
                    continue
                sample = next(iter(mhs.values()))
                targ_lats = sample["lat"]
                targ_lons = sample["lon"]

                axis_keys = {"lat", "lon", "P_level", "alt"}
                all_keys  = set().union(*[d.keys() for d in mhs.values()])
                base_vars = [k for k in all_keys if (k not in axis_keys) and (not k.endswith("_std"))]

                annual_map[y] = {"lat": targ_lats, "lon": targ_lons, 
                                 "P_level": sample["P_level"], "alt": sample["alt"]}
                for var in base_vars:
                    arr_list = [mhs[m][var] for m in months_used if var in mhs[m]]
                    if not arr_list:
                        continue
                    stack = np.stack(arr_list, axis=0).astype(float)  # (nmonths, lev, lat, lon)
                    annual_map[y][var] = np.nanmean(stack, axis=0)
                    
                    missing = [mm for mm in months_used if var not in mhs[mm]]
                    if missing:
                        print(f"[INFO] year {y}: '{var}' averaged over {stack.shape[0]}/{len(months_used)} months; missing: {missing}")
            
            # Years average
            years_used = sorted(annual_map.keys())

            sample_y = annual_map[years_used[-1]]
            output_data = {"lat": sample_y["lat"], "lon": sample_y["lon"]}
            output_data["P_level"] = sample_y["P_level"]
            output_data["alt"] = sample_y["alt"]

            axis_keys = {"lat", "lon", "P_level", "alt"}
            all_keys_years = set().union(*[d.keys() for d in annual_map.values()])
            base_vars = [k for k in all_keys_years if (k not in axis_keys) and (not k.endswith("_std"))]

            for var in base_vars:
                arr_list = [annual_map[y][var] for y in years_used if var in annual_map[y]]
                if not arr_list:
                    continue
                stack = np.stack(arr_list, axis=0).astype(float)  # (nyears, lev, lat, lon)
                output_data[var] = np.nanmean(stack, axis=0)
                output_data[f"{var}_std"] = np.nanstd(stack, axis=0)
                
                missing = [yy for yy in years_used if var not in annual_map[yy]]
                if missing:
                    print(f"[INFO] cyly: '{var}' averaged over {stack.shape[0]}/{len(years_used)} years; missing: {missing}")

            CAMS._saveNPZCAMS(self, data=output_data, dataType='cyly', y=years) # !!! [redudant]
        
        print("The data is combined.")
    
    def selectRegionCAMS(self, data, bbox): # !!! (redundant)
        """
        !!! BE AWARE OF BBOX DIFFERENCE THAN MERRA-2 !!!
        
        Select specific region of Earth of downloaded data.

        Parameters
        ----------
        data : DICT
            Data in dictionary form.
        bbox : LIST
            Requested coordinates, the bounding box.
            (upper_right_latitude, lower_left_longitude, lower_left_latitude, upper_right_longitude)

        Returns
        -------
        dataSel : DICT
            Data from requested region.
        """
        # Selected coordinates
        if len(bbox) == 2:
            bbox = bbox + bbox
            uniqueCoor = True
        elif len(bbox) == 4:
            bbox = bbox
            uniqueCoor = False
        else:
            print('\n!EcoSysEM.Error: boundaries() requires 2 `(lon, lat)` or 4 `(lower_left_lon, lower_left_lat, upper_right_lon, upper_right_lat)` positional arguments.')
            sys.exit()
        # BBox from data
        lonR = data['lon']
        latR = data['lat']
        bboxData = (lonR[0], latR[0], lonR[-1], latR[-1])
        # Get closest coordinates
        bbox =  MERRA2._closestCoord(self, lonR, latR, bbox) # !!! (redundant)
        # Check if requested region is in data
        cBBOX = [(bboxData[i] <= bbox[i] and bboxData[i+2] >= bbox[i+2]) for i in np.arange(2)]
        if all(cBBOX) == True:
            # Get indices
            # (lower_left_lon, lower_left_lat, upper_right_lon, upper_right_lat)
            idx = (int(np.argwhere(lonR == bbox[0])),
                   int(np.argwhere(latR == bbox[1])),
                   int(np.argwhere(lonR == bbox[2])) + 1,
                   int(np.argwhere(latR == bbox[3])) + 1)
            # Initialize dictionary of selected data
            dataSel = {}
            for var in data:
                if var == 'lon':
                    if uniqueCoor:
                        dataSel[var] = data[var][idx[0]]
                    else:
                        if idx[0] == idx[2]:
                            dataSel[var] = data[var][idx[0]]
                        else:
                            dataSel[var] = data[var][idx[0]:idx[2]]
                elif var == 'lat':
                    if uniqueCoor:
                        dataSel[var] = data[var][idx[1]]
                    else:
                        if idx[1] == idx[3]:
                            dataSel[var] = data[var][idx[1]]
                        else:
                            dataSel[var] = data[var][idx[1]:idx[3]]
                else:
                    if uniqueCoor:
                        dataSel[var] = data[var][idx[1],idx[0]]
                    else:
                        if idx[0] == idx[2]:
                            dataSel[var] = data[var][idx[1]:idx[3],idx[0]]
                        elif idx[1] == idx[3]:
                            dataSel[var] = data[var][idx[1],idx[0]:idx[2]]
                        else:
                            dataSel[var] = data[var][idx[1]:idx[3],idx[0]:idx[2]]
                if not isinstance(dataSel[var], (list, np.ndarray)): dataSel[var] = [dataSel[var]]
            return dataSel
        else:
            print('\n!EcoSysEM.Error: selected region is outside the data boundaries.')
            sys.exit()

    def _saveNPZCAMS(self, data, dataType, y, m = None, d = None): # !!! (redundant)
        """
        Create .npz file with downladed data.

        Parameters
        ----------
        date : DICT
            Data in dictionary form.
        dataType   : STR ('mly' or 'dly' or 'cmly')
            Type of data.
        y : INT or LIST of INT
            Year(s) of data.
        m : INT or LIST of INT
            Month(s) of data. 
        d : INT or LIST of INT
            Day of month of data.
        """
        path = f'data/CAMS/{dataType}/'
        if not os.path.isdir(path):
            os.makedirs(path)
        # File name (based on dataType)
        if dataType == 'dly':
            if not isinstance(y, int):
                print('\n!EcoSysEM.Error: argument \'y\' must be a integer.')
                sys.exit()
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer.')
                sys.exit()
            if not isinstance(d, int):
                print('\n!EcoSysEM.Error: argument \'d\' must be a integer.')
                sys.exit()
            file = f'{y}_{m}_{d}_day.npz'
        elif dataType == 'mly':
            if not isinstance(y, int):
                print('\n!EcoSysEM.Error: argument \'y\' must be a integer.')
                sys.exit()
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer.')
                sys.exit()
            file = f'{y}_{m}_month.npz'
        elif dataType == 'cmly':
            if not isinstance(y, list):
                print('\n!EcoSysEM.Error: argument \'y\' must be a list: [start_year, end_year]')
                sys.exit()
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
                sys.exit()
            file = f'{y[0]}_{y[-1]}_{m}.npz'
        elif dataType == 'yly':
            if not isinstance(y, int):
                print('\n!EcoSysEM.Error: argument \'y\' must be an integer')
                sys.exit()
            file = f'{y}_year.npz'
        elif dataType == 'cyly':
            if not isinstance(y, list):
                print('\n!EcoSysEM.Error: argument \'y\' must be a list: [start_year, end_year]')
                sys.exit()
            file = f'{y[0]}_{y[-1]}.npz'
        else:
            print('\n!EcoSysEM.Error: argument \'dataType\' must be \'mly\' (to generate monthly data), \'cmly\' (to generate combined monthly data),  \'yly\' (to generate annual data) or \'cyly\' (to generate combined annual data).')
            sys.exit()
        # Path generation
        pathfile = path + file
        # Save .npz file
        np.savez(pathfile, **data)   
        return pathfile
    
    def _openNPZCAMS(self, dataType, y, m = None, d = None): # !!! (redundant)
        """
        Open .npz file with downladed data.

        Parameters
        ----------
        date : DICT
            Data in dictionary form.
        dataType   : STR ('dly' or 'mly' or 'cmly')
            Type of data.
        y : INT or LIST of INT
            Year(s) of data.
        m : INT or LIST of INT
            Month(s) of data. 
        d : INT or LIST of INT
            Day of month of data.
        """
        path = f'data/CAMS/{dataType}/'
        if dataType == 'dly':
            if not isinstance(y, int):
                print('\n!EcoSysEM.Error: argument \'y\' must be a integer')
                sys.exit()
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
                sys.exit()
            if not isinstance(d, int):
                print('\n!EcoSysEM.Error: argument \'d\' must be a integer')
                sys.exit()
            file = f'{y}_{m}_{d}_day.npz'
        if dataType == 'mly':
            if not isinstance(y, int):
                print('\n!EcoSysEM.Error: argument \'y\' must be a integer')
                sys.exit()
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
                sys.exit()
            file = f'{y}_{m}_month.npz'
        elif dataType == 'cmly':
            if not isinstance(y, list):
                print('\n!EcoSysEM.Error: argument \'y\' must be a list: [start_year, end_year]')
                sys.exit()
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
                sys.exit()
            file = f'{y[0]}_{y[-1]}_{m}.npz'
        elif dataType == 'yly':
            if not isinstance(y, int):
                print('\n!EcoSysEM.Error: argument \'y\' must be an integer')
                sys.exit()
            file = f'{y}_year.npz'
        elif dataType == 'cyly':
            if not isinstance(y, list):
                print('\n!EcoSysEM.Error: argument \'y\' must be a list: [start_year, end_year]')
                sys.exit()
            file = f'{y[0]}_{y[-1]}.npz'
        else:
            print('\n!EcoSysEM.Error: argument \'dataType\' must be \'mly\' (to generate monthly data), \'cmly\' (to generate combined monthly data),  \'yly\' (to generate annual data) or \'cyly\' (to generate combined annual data).')
            sys.exit()
        return np.load(path + file)
    
    def dictCAMS(self, dataType, y, m = None, d = None, keys = 'All'): # !!! (redundant)
        """
        Get data in dictionary form.

        Parameters
        ----------
        dataType : STR ('mly' or 'dly' or 'cmly')
            Type of data.
        y : INT or LIST of INT
            Year(s) of data.
        m : INT or LIST of INT
            Month of data.
        d : INT
            Day of month of data.
        keys : LIST of STR
            List of requested variables. (Default: 'All')
        
        Returns
        -------
        dictVar : DICT
            Dictionary with requested variables.
        """
        npz = CAMS._openNPZCAMS(self, dataType, y, m, d) # !!! (redundant)
        if keys == 'All':
            keys = CAMS.keysCAMS(self, dataType, y, m, d) # !!! (redundant)
            dictVar = {key: npz[key] for key in keys}
        else:
            coor = []
            if not np.any(np.char.find('lat', keys) > -1):
                coor += ['lat']
            if not np.any(np.char.find('lon', keys) > -1):
                coor += ['lon']
            keys = coor + keys
            dictVar = {key: npz[key] for key in keys}
        npz.close()
        return dictVar
    
    def keysCAMS(self, dataType, y, m = None, d = None): # !!! (redundant)
        """
        Get variable list of data.

        Parameters
        ----------
        dataType : STR ('mly' or 'dly' or 'cmly')
            Type of data.
        y : INT or LIST of INT
            Year(s) of data.
        m : INT or LIST of INT
            Month of data. 
        d : INT
            Day of month of data.
        Returns
        -------
        keys  : LIST of STR
            Variable list of data.
        """
        npz = CAMS._openNPZCAMS(self, dataType, y, m, d) # !!! (redundant)
        keys = npz.files
        npz.close()
        return keys
    
    def deleteKeyCAMS(self, keys, dataType, y, m = None, d = None): # !!! (redundant)
        """
        Delete variable(s) from data.

        Parameters
        ----------
        keys : LIST of STR
            List of variables to be deleted.
        dataType : STR ('mly' or 'dly' or 'cmly')
            Type of data.
        y : INT or LIST of INT
            Year(s) of data.
        m : INT or LIST of INT
            Month of data.
        d : INT
            Day of month of data.
        """
        dictVar = CAMS.dictCAMS(self, dataType, y, m, d) # !!! (redundant)
        for key in keys:
            dictVar.pop(key, None)
        # Save new .npz file
        CAMS._saveNPZCAMS(self, dictVar, dataType, y, m, d) # !!! (redundant)
    
    def _deleteTempDataCAMS(self, filename):
        """
        Remove .nc file in "temporary" folder.
        
        Parameters
        ----------
        filename : STR
            Name of the file.
        """
        if os.path.exists(filename):
            os.remove(filename)
        else:
            print("The file does not exist")

class ISAMERRA2(ISA, MERRA2):
    """
    Combination of International Standard Atmosphere (ISA) model and Modern-Era 
    Retrospective analysis for Research and Applications Version 2 (MERRA-2).
    """
    def __init__(self, layers = 0, H2O = 0.0, pH = 8.0, resolution = 1000):
        ISA.__init__(self, layers = layers, H2O = H2O, pH = pH, resolution = resolution)
        MERRA2.__init__(self)
    
    def getConcISAMERRA2(self, phase, dataType, y, m, d = None, compound = None, bbox = (-180, -90, 180, 90), altArray = None, num = 50, surftrop = None):
        """
        Computation of vertical profiles of compounds (parcial pressure, Pi;
        gas concentration, Ci_G; liquid concentration in fresh water, Ci_L-FW;
        and liquid concentration in sea water, Ci_L-SW).
        Gas concentrations (Ci_G) are calculated using Dalton's law and the 
        ideal gas law, and liquid concentration (Ci_LFW and Ci_LSW) with 
        Henry's law.

        Parameters
        ----------
        phase : STR ('G', 'L-FW', 'L-SW', 'L' or 'All')
            Selection of phase of vertical profile.
                'G' - Gas.
                'L-FW' - Liquid fresh water.
                'L-SW' - Liquid sea water.
                'L' - Both liquid phases (L-FW, L-SW).
                'All' - All phases (G, L-FW, L-SW).
        dataType : STR ('mly', 'cmly', 'dly')
            Type of data
        y : INT or LIST of INT
            Year(s) of data
        m : INT or LIST of INT
            Month of data
        d : INT or LIST of INT
            Day(s) of data
        compound : STR or LIST, optional
            Interested compounds. Default: None -> All compounds.
        bbox : TUPLE, optional
            Earths region of data, the bounding box.
            (lower_left_longitude, lower_left_latitude, upper_right_longitude, upper_right_latitude)
        altArray : LIST or np.ndarray, optional
            List of altitudes
        num : INT, optional
            Number of altitude steps to generate.
        surftrop : STR ('surface', 'tropopause'), optional
            Get concentration from 2-meters air following topography (surftrop='surface') or tropopause height (surftrop='tropopause').
        Returns
        -------
        Dictionaries with pressures/concentrations.
        dict_Pi, dict_Ci_G : DICT (if phase='G')
        dict_Ci_LFW : DICT (if phase='L-FW')
        dict_Ci_LSW : DICT (if phase='L-SW')
        dict_Ci_LFW, dictCi_LSW : DICT (if phase='L')
        dict_Pi, dict_Ci_G, dict_Ci_LFW, dict_Ci_LSW : DICT (if phase='All')
             dict_Pi : Parcial pressure of desired compounds.
             dict_Ci_G : Concentration in gas of desired compounds.
             dict_Ci_LFW : Concentration in liquid (freshwater) of desired compounds.
             dict_Ci_LSW : Concentration in liquid (seawater) of desired compounds.

        """
        # Selection of compound and composition
        compounds = self.compounds
        compositions = np.array(list(self.compositions.values()))
        if compound:
            if type(compound) is str: compound = [compound]
            findC = compounds.reset_index().set_index('Compounds').loc[compound].reset_index().set_index('index').index
            compositions = compositions[findC]
            compounds = compounds[findC]
        # Temperature [K], pressure [Pa], altitude [m]
        if not surftrop:
            t, p, _ = MERRA2.getTPAlt(self, dataType, y, m, d, bbox, altArray, num)
        else:
            data = MERRA2.loadDataMERRA2(self, dataType, y, m, d)
            data = MERRA2.selectRegion(self, data, bbox)
            if surftrop == 'surface':
                t = np.array(data['T2M'])
                p = np.array(data['PS'])
            elif surftrop == 'tropopause':
                t = np.array(data['TROPT'])
                p = np.array(data['TROPPB'])
            else:
                print('\n!EcoSysEM.Error: argument \'surftrop\' must `surface` or `tropopause`.')
                sys.exit()
        # Constants
        R_g = 8314.46261815324  # Universal gas constant [(L·Pa)/(K·mol)]
        Hs_FW, notNaN_HsFW = eQ.solubilityHenry(compounds, 'FW', t)
        Hs_SW, notNaN_HsSW = eQ.solubilityHenry(compounds, 'SW', t)
        # Dictionaries initialization
        dict_Pi = {}
        dict_Ci_G = {}
        dict_Ci_LFW = {}
        dict_Ci_LSW = {}
        compounds = compounds.values
        for id_, composition in enumerate(compositions):
            # Gas phase - Partial pressure (Pi)
            Pi = p * composition # [Pa]
            # Gas phase - Gas concentration (Ci_G)
            Ci_G = (Pi / (R_g * t))
            # Liquid phase - Freshwater (Ci_LFW)
            if notNaN_HsFW[id_]:
                Ci_LFW = Pi * Hs_FW[..., id_] * (1/1000) # [mol/L]
            else:
                Ci_LFW = None
            # Liquid phase - Seawater (Ci_LSW)
            if notNaN_HsSW[id_]:
                Ci_LSW = Pi * Hs_SW[..., id_] * (1/1000) # [mol/L]
            else:
                Ci_LSW = None
            # Save data in dictionary
            dict_Pi[compounds[id_]] = Pi
            dict_Ci_G[compounds[id_]] = Ci_G
            dict_Ci_LFW[compounds[id_]] = Ci_LFW
            dict_Ci_LSW[compounds[id_]] = Ci_LSW
        if phase == 'G':
            return dict_Pi, dict_Ci_G
        elif phase == 'L-FW':
            return dict_Ci_LFW
        elif phase == 'L-SW':
            return dict_Ci_LSW
        elif phase == 'L':
            return dict_Ci_LFW, dict_Ci_LSW
        elif phase == 'All':
            return dict_Pi, dict_Ci_G, dict_Ci_LFW, dict_Ci_LSW
        else:
            print('!EcosysEM.Error: No phase selected. Use one of the following string:\n'+
                  '                             \'G\'       - Gas.\n'+
                  '                             \'L-FW\'    - Liquid fresh water.\n'+
                  '                             \'L-SW\'    - Liquid sea water.\n'+
                  '                             \'L\'       - Both liquid phases (L-FW, L-SW).\n'+
                  '                             \'All\'     - All phases (G, L-FW, L-SW).')
            sys.exit()

class CAMSMERRA2(CAMS, MERRA2):
    """
    Combination of Modern-Era Retrospective analysis for Research and 
    Applications Version 2 (MERRA-2) and Copernicus Atmosphere Monitoring
    Service (CAMS) database.
    """
    def __init__(self):
        CAMS.__init__(self)
        MERRA2.__init__(self)
    
    def interpolateCAMS(self, dataType, year, month, day=None, loc=None,
                            molecules = ('CO', 'CO2', 'CH4'), 
                            target_lats = np.arange(-90, 90.1, 0.5),
                            target_lons = np.arange(-180,  179.375 + 1e-3, 0.625),
                            method='linear'):
        """
        Interpolate CAMS .npz files onto target MERRA2 grid.
    
        Parameters
        ----------
        dataType : STR
            Name of the subfolder under `data/CAMS/` containing .npz files.
        year : INT
            Year of the desired dataset.
        month : INT
            Month of the desired dataset.
        day : INT, optional
            Day of the desired dataset.
        loc : STR
            Get concentration from 2-meters air following topography (loc='surface') or tropopause height (loc='tropopause').
        molecules : TUPLE of STR, optional
            Variable names to process (default: ('CO', 'CO2', 'CH4')).
        target_lats : 1D array, optional
            Desired latitudes for the CAMS grid.
            (default: np.arange(-90, 90.1, 0.5)).
        target_lons : 1D array, optional.
            Desired longitudes for the CAMS grid
            (default: np.arange(-180, 179.375+0.001, 0.625)).
        method : STR, optional
            Method of interpolation (default: 'linear').
    
        Returns
        -------
        result: DICT
            Interpolated data.
        """
        folder = f'data/CAMS/{dataType}/'
        if day is None:
            fname = f"{year}_{month}_month.npz"
        else:
            fname = f"{year}_{month}_{day}_day.npz"

        path = os.path.join(folder, fname)
        if not os.path.isfile(path):
            raise FileNotFoundError(f"CAMS file not found: {path}")

        npz = np.load(path)
        orig_lats = npz['lat']
        orig_lons = npz['lon']
        orig_alt  = npz['alt']
        orig_plev = npz['P_level']

        # Ensure ascending latitude order
        if orig_lats[0] > orig_lats[-1]:
            orig_lats = orig_lats[::-1]
            flip_lat = True
        else:
            flip_lat = False

        mesh_lat, mesh_lon = np.meshgrid(target_lats, target_lons, indexing='ij')
        points = np.stack([mesh_lat.ravel(), mesh_lon.ravel()], axis=-1)

        result = {}
        
        for var in list(molecules) + [f"{m}_std" for m in molecules]:
            if var not in npz:
                continue
            data3d = npz[var]
            if flip_lat:
                data3d = data3d[:, ::-1, :]
    
            n_alt = data3d.shape[0]
            out = np.empty((n_alt, target_lats.size, target_lons.size), dtype=data3d.dtype)
            for k in range(n_alt):
                interp = RegularGridInterpolator(
                    (orig_lats, orig_lons),
                    data3d[k],
                    method=method,
                    bounds_error=False,
                    fill_value=np.nan
                )
                out[k] = interp(points).reshape(target_lats.size, target_lons.size)
    
            result[var] = out
        
        result['lat'] = target_lats
        result['lon'] = target_lons
        result['alt'] = orig_alt
        result['P_level'] = orig_plev
    
        return result
    
    def getConcCAMS(self, phase, data, dataType, year, month, day=None, bbox = (-180, -90, 180, 90), altArray=None, loc=None, num=None):
        """
        Converts the mass ratio (kg/kg) to concentration (mol/L).

        Parameters
        ----------
        phase : STR ('G', 'L-FW', 'L-SW', 'L' or 'All')
            Selection of phase of vertical profile.
                'G' - Gas.
                'L-FW' - Liquid fresh water.
                'L-SW' - Liquid sea water.
                'L' - Both liquid phases (L-FW, L-SW).
                'All' - All phases (G, L-FW, L-SW).
        data : DICT
            Data in dictionary.
        dataType : STR
            Name of the subfolder under `data/CAMS/` containing .npz files.
        year : INT
            Year of the desired dataset.
        month : INT
            Month of the desired dataset.
        day : INT, optional
            Day of the desired dataset.
        altArray : LIST or np.ndarray, optional
            List of altitudes in m.
        loc : STR
            Get concentration from 2-meters air following topography (loc='surface') or tropopause height (loc='tropopause').
        num : INT, optional
            Number of altitude steps to generate.
        """
        from pyatmos import coesa76
        
        cams_molecules = ('CO', 'CO2', 'CH4')
        molecule_data = {}
        for mol in cams_molecules:
            if mol in data:
                molecule_data[mol] = data[mol]
        if not molecule_data:
            raise ValueError("No valid molecule keys found. Expected one of: CO, CO2, CH4")
    
        cams_plev = np.asarray(data['P_level'], dtype=float)  # (L_cams,)
        cams_alt  = np.asarray(data['alt'], dtype=float)      # (L_cams,)
    
        # MERRA2
        merra2 = CAMSMERRA2.loadDataMERRA2(self, dataType, year, month, day)
        PS = np.array(merra2['PS'])
        TS = np.array(merra2['T2M'])
        TROPPB = np.array(merra2['TROPPB'])
        TROPT = np.array(merra2['TROPT'])
    
        # Target T/P by loc
        if loc == 'surface':
            t_target, p_target = TS, PS
        elif loc == 'tropopause':
            t_target, p_target = TROPT, TROPPB
        else:
            t_target, p_target, z_m = CAMSMERRA2.getTPAlt(self, dataType, year, month, day, bbox, altArray, num)
    

        h_km = cams_alt * 1e-3 # km
        rho_kg_m3 = coesa76(h_km).rho # kg/m3
        rho_kg_L  = rho_kg_m3 * 1e-3 # kg/L
        rho = rho_kg_L[:, None, None]
    
        # Constants
        R_g = 8314.46261815324  # Universal gas constant [(L·Pa)/(K·mol)]
        
        # Dictionaries initialization
        dict_Pi = {}
        dict_Ci_G = {}
        dict_Ci_LFW = {}
        dict_Ci_LSW = {}
    
        for molecule, cams_array in molecule_data.items():
            cams_array = np.array(cams_array)
            M_kg_per_mol = Formula(molecule).mass * 1e-3    # g/mol to kg/mol

            conc_cams = (cams_array * rho) / M_kg_per_mol
            conc = self._reshapeAltCAMS(conc_cams, p_target, cams_plev)
            
            Pi = conc * (R_g * t_target)
            
            Hs_FW, notNaN_HsFW = eQ.solubilityHenry(molecule, 'FW', t_target)
            Hs_SW, notNaN_HsSW = eQ.solubilityHenry(molecule, 'SW', t_target)
            if Hs_FW.ndim >= 3 and Hs_FW.shape[-1] == 1: 
                Hs_FW = Hs_FW[..., 0]
            if Hs_SW.ndim >= 3 and Hs_SW.shape[-1] == 1: 
                Hs_SW = Hs_SW[..., 0]
            if notNaN_HsFW.ndim >= 3 and notNaN_HsFW.shape[-1] == 1: 
                notNaN_HsFW = notNaN_HsFW[..., 0]
            if notNaN_HsSW.ndim >= 3 and notNaN_HsSW.shape[-1] == 1: 
                notNaN_HsSW = notNaN_HsSW[..., 0]
    
            Ci_LFW = np.where(notNaN_HsFW, Pi * Hs_FW / 1000.0, np.nan)  # mol/L
            Ci_LSW = np.where(notNaN_HsSW, Pi * Hs_SW / 1000.0, np.nan)  # mol/L
    
            # Save data in dictionary
            dict_Pi[molecule]     = Pi
            dict_Ci_G[molecule]   = conc
            dict_Ci_LFW[molecule] = Ci_LFW
            dict_Ci_LSW[molecule] = Ci_LSW
    
        if phase == 'G':
            return dict_Pi, dict_Ci_G
        elif phase == 'L-FW':
            return dict_Ci_LFW
        elif phase == 'L-SW':
            return dict_Ci_LSW
        elif phase == 'L':
            return dict_Ci_LFW, dict_Ci_LSW
        elif phase == 'All':
            return dict_Pi, dict_Ci_G, dict_Ci_LFW, dict_Ci_LSW
        else:
            print('!EcosysEM.Error: No phase selected. Use one of the following string:\n'+
                  '                             \'G\'       - Gas.\n'+
                  '                             \'L-FW\'    - Liquid fresh water.\n'+
                  '                             \'L-SW\'    - Liquid sea water.\n'+
                  '                             \'L\'       - Both liquid phases (L-FW, L-SW).\n'+
                  '                             \'All\'     - All phases (G, L-FW, L-SW).')
            sys.exit()
    
    def _reshapeAltCAMS(self, orig_data, targ_data, plev):
        """
        Reshape altitude levels.

        Parameters
        ----------
        orig_data : 3D array
            CAMS concentration data.
        targ_data : 3D array
            MERRA2 pressure data.
        plev : 1D array
            CAMS pressure levels.

        Returns
        -------
        new_data : 3D array
            Concentration data in the shape of target data.

        """
        if orig_data.ndim == 2:
            orig_data = orig_data[np.newaxis, ...]
       
        if targ_data.ndim == 2:
            n_lat, n_lon = targ_data.shape
            new_data = np.full((n_lat, n_lon), np.nan, dtype=np.float64)
            tp = targ_data
    
            diffs = np.abs(plev[:, None, None] - tp[None, :, :])
            diffs[:, np.isnan(tp)] = np.inf
            idx_cell = diffs.argmin(axis=0)
            idx_cell = np.clip(idx_cell, 0, orig_data.shape[0] - 1)
            
            for i in range(n_lat):
                for j in range(n_lon):
                    if np.isnan(tp[i, j]):
                        continue
                    val = orig_data[idx_cell[i, j], i, j]
                    if not np.isnan(val):
                        new_data[i, j] = val
        elif targ_data.ndim == 3:
            n_target, n_lat, n_lon = targ_data.shape
            new_data = np.full((n_target, n_lat, n_lon), np.nan, dtype=np.float64)  # float64 for safe NaN handling
        
            for k in range(n_target):
                tp = targ_data[k]
        
                # Compute absolute pressure difference between each CAMS level and target pressure
                diffs = np.abs(plev[:, None, None] - tp[None, :, :])  # shape: (n_plev, n_lat, n_lon)
        
                # Prevent NaNs in targ_data from biasing argmin — set their diffs to infinity
                diffs[:, np.isnan(tp)] = np.inf
        
                # For each (i,j), take CAMS level index closest to target pressure
                idx_cell = diffs.argmin(axis=0)  # shape (n_lat, n_lon)
                idx_cell = np.clip(idx_cell, 0, orig_data.shape[0] - 1)
        
                # Assign concentration values only for valid target pressure points
                for i in range(n_lat):
                    for j in range(n_lon):
                        if np.isnan(tp[i, j]):
                            continue
                        val = orig_data[idx_cell[i, j], i, j]
                        if not np.isnan(val):
                            new_data[k, i, j] = val
    
        return new_data
