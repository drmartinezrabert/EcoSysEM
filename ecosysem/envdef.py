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
            return None
            
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
        fig, ax1 = plt.subplots(figsize = (3, 4))
        ax1.set_ylabel('Altitude (km)')
        ax1.set_xlabel('Temperature (K)', color = 'tab:red')
        ax1.set_xlim([150, 300])
        ax1.set_ylim([0, alt[-1]])
        ax1.plot(t, alt, color = 'tab:red')
        # Pressure
        ax2 = ax1.twiny()
        ax2.set_xlabel('Pressure (atm)', color = 'tab:blue')
        ax2.set_xlim([0, 1])
        ax2.plot(p, alt, color = 'tab:blue')
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
        bbox = (lonR[np.abs(lonR - Coor[0]).argmin()],
                latR[np.abs(latR - Coor[1]).argmin()],
                lonR[np.abs(lonR - Coor[2]).argmin()],
                latR[np.abs(latR - Coor[3]).argmin()])
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
        dataType   : STR ('mly' or 'cmly')
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
                return None
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
                return None
            if not isinstance(d, int):
                print('\n!EcoSysEM.Error: argument \'d\' must be a integer')
                return None
            file = f'{y}_{m}_{d}_day.npz'
        elif dataType == 'mly':
            if not isinstance(y, int):
                print('\n!EcoSysEM.Error: argument \'y\' must be a integer')
                return None
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
                return None
            file = f'{y}_{m}_month.npz'
        elif dataType == 'cmly':
            if not isinstance(y, list):
                print('\n!EcoSysEM.Error: argument \'y\' must be a list: [start_year, end_year]')
                return None
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
                return None
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
        dataType   : STR ('mly' or 'cmly')
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
                return None
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
                return None
            if not isinstance(d, int):
                print('\n!EcoSysEM.Error: argument \'d\' must be a integer')
                return None
            file = f'{y}_{m}_{d}_day.npz'
        if dataType == 'mly':
            if not isinstance(y, int):
                print('\n!EcoSysEM.Error: argument \'y\' must be a integer')
                return 0
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
                return 0
            file = f'{y}_{m}_month.npz'
        elif dataType == 'cmly':
            if not isinstance(y, list):
                print('\n!EcoSysEM.Error: argument \'y\' must be a list: [start_year, end_year]')
                return None
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
                return None
            file = f'{y[0]}_{y[-1]}_{m}.npz'
        return np.load(path + file)
    
    def getDataMERRA2_v1(self, years, months, 
                         days = 'All',
                         product = 'M2I1NXASM',
                         version = '5.12.4',
                         bbox = (-180, -90, 180, 90),
                         var = ['PS', 'T2M', 'TROPT', 'TROPPB'], 
                         update = False,
                         daily = False):
        """
        Download data from MERRA2 database.
        
        Parameters
        ----------
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
        # Coordinates (bbox)
        if len(bbox) == 2:
            bbox = bbox + bbox
        elif len(bbox) == 4:
            bbox = bbox
        else:
            print('\n!EcoSysEM.Error: boundaries() requires 2 `(lon, lat)` or 4 `(lower_left_lon, lower_left_lat, upper_right_lon, upper_right_lat)` positional arguments.')
            return None
        # Check arguments
        if not isinstance(years, list): years = [years]
        if not isinstance(months, (list, np.ndarray)): months = [months]
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
            start_date = datetime.date(y, m, 1)
            if days == 'All':
                last_day = calendar.monthrange(y, m)[1]
            else:
                if not isinstance(days, list): days = [days]
                last_day = max(days)
            end_date = datetime.date(y, m, last_day)
            delta = datetime.timedelta(days=1)
            date = start_date
            # Initialize dictionaries
            hourData = {}
            av_dayData = {}
            dayData = {}
            ## Day matrices
            while (date <= end_date):
                print(f"> {date.strftime('%Y-%m-%d')}")
                # Open granules to local phat
                results = earthaccess.search_data(
                    short_name = product, 
                    version = version,
                    temporal = (date, date),
                    bounding_box = bbox
                )
                # Open granules using Xarray
                fs = earthaccess.open(results)
                ds = xr.open_mfdataset(fs)
                lonR = ds.lon.to_numpy()
                latR = ds.lat.to_numpy()
                # Get closest coordinates
                bbox =  MERRA2._closestCoord(self, lonR, latR, bbox)
                lat_slice = slice(bbox[1], bbox[3])
                lon_slice = slice(bbox[0], bbox[2])
                ds = ds.sel(lat = lat_slice, lon = lon_slice)
                lon = ds.lon.to_numpy()
                lat = ds.lat.to_numpy()
                for key in var:
                    try:
                        ds[key]
                    except:
                        print(f'  {key} estimated.')
                    else:
                        hourData[key] = ds[key].values
                        print(f'  {key} done.')
                # Atmospheric altitudes
                varH = ['H', 'TROPH']
                varReqH = ['PS', 'TROPPB', 'TROPPT', 'TROPPV']
                if np.any([np.char.find(var, iVar) > -1 for iVar in varH]):
                    c = 0
                    for iV in var:
                        if np.any(np.char.find(iV, varReqH) > -1):
                            varNames = {'PS': 'H',
                                        'TROPPB': 'TROPH',
                                        'TROPPT': 'TROPH',
                                        'TROPPV': 'TROPH'}
                            hourData[varNames[iV]] = MERRA2._HfromP(self, hourData[iV])
                            c += 1
                    if c == 0:
                        print('\n!EcoSysEM.Error: missing required variable missing for atmospheric altitude(s) - \'PS\', \'TROPPB\', \'TROPPT\', or \'TROPPV\'.')
                        return None
                # Lapse rate
                varReqLR = ['T2M', 'TROPT', 'H', 'TROPH']
                if np.any(np.char.find(var, 'LR') > -1):
                    if np.all([np.any(np.char.find(var, iVarReq) > -1) for iVarReq in varReqLR]):
                        hourData['LR'] = MERRA2._LR(self, hourData['T2M'], hourData['TROPT'], hourData['H'], hourData['TROPH'])
                    else:
                        print('\n!EcoSysEM.Error: missing required variable missing for lapse rate - \'T2M\', \'TROPT\', \'H\', or \'TROPH\'.')
                        return None
                for iV in var:
                    # Daily averages and std
                    av_dayData[iV] = np.average(hourData[iV], axis = 0)
                    av_dayData[f'{iV}_std'] = np.std(hourData[iV], axis = 0)
                    # Day matrices
                    if date == start_date:
                        dayData[iV] = av_dayData[iV]
                    else:
                        dayData[iV] = np.dstack((dayData[iV], av_dayData[iV]))
                # Save daily data
                if daily == True:
                    MERRA2._saveNPZMERRA2(self, data = av_dayData, dataType = 'dly', y = y, m = m, d = int(date.day))
                # Next date
                date += delta
            # Month matrices
            monthData['lat'] = lat
            monthData['lon'] = lon
            for iV in var:
                if last_day != 1:
                    monthData[iV] = np.average(dayData[iV], axis = -1)
                    monthData[f'{iV}_std'] = np.std(dayData[iV], axis = -1)
                else:
                    monthData[iV] = dayData[iV]
                    monthData[f'{iV}_std'] = np.array([0.0])
            # Save numpy matrices in .npz format
            MERRA2._saveNPZMERRA2(self, data = monthData, dataType = 'mly', y = y, m = m)
        print("--- %s seconds ---" % (time.time() - start_time))
    
    def getDataMERRA2(self, dataType, years, months,
                      days = 'All',
                      product = 'M2I1NXASM',
                      version = '5.12.4',
                      bbox = (-180, -90, 180, 90),
                      var = ['PS', 'T2M', 'TROPT', 'TROPPB'], 
                      update = False,
                      daily = False):
        """
        Download data from MERRA2 database.
        
        Parameters
        ----------
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
        
        """
        # Coordinates (bbox)
        if len(bbox) == 2:
            bbox = bbox + bbox
        elif len(bbox) == 4:
            bbox = bbox
        else:
            print('\n!EcoSysEM.Error: boundaries() requires 2 `(lon, lat)` or 4 `(lower_left_lon, lower_left_lat, upper_right_lon, upper_right_lat)` positional arguments.')
            return None
        # Check arguments
        if not isinstance(years, list): years = [years]
        if not isinstance(months, (list, np.ndarray)): months = [months]
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
            start_date = datetime.date(y, m, 1)
            if days == 'All':
                last_day = calendar.monthrange(y, m)[1]
            else:
                if not isinstance(days, list): days = [days]
                last_day = max(days)
            end_date = datetime.date(y, m, last_day)
            delta = datetime.timedelta(days=1)
            date = start_date
            # Initialize dictionaries
            hourData = {}
            av_dayData = {}
            dayData = {}
            ## Day matrices
            while (date <= end_date):
                print(f"> {date.strftime('%Y-%m-%d')}")
                # Open granules to local phat
                results = earthaccess.search_data(
                    short_name = product, 
                    version = version,
                    temporal = (date, date),
                    bounding_box = bbox
                )
                # Open granules using Xarray
                fs = earthaccess.open(results)
                ds = xr.open_mfdataset(fs)
                lonR = ds.lon.to_numpy()
                latR = ds.lat.to_numpy()
                # Get closest coordinates
                bbox =  MERRA2._closestCoord(self, lonR, latR, bbox)
                lat_slice = slice(bbox[1], bbox[3])
                lon_slice = slice(bbox[0], bbox[2])
                ds = ds.sel(lat = lat_slice, lon = lon_slice)
                lon = ds.lon.to_numpy()
                lat = ds.lat.to_numpy()
                for key in var:
                    try:
                        ds[key]
                    except:
                        print(f'  {key} estimated.')
                    else:
                        hourData[key] = ds[key].values
                        print(f'  {key} done.')
                # Atmospheric altitudes
                varH = ['H', 'TROPH']
                varReqH = ['PS', 'TROPPB', 'TROPPT', 'TROPPV']
                if np.any([np.char.find(var, iVar) > -1 for iVar in varH]):
                    c = 0
                    for iV in var:
                        if np.any(np.char.find(iV, varReqH) > -1):
                            varNames = {'PS': 'H',
                                        'TROPPB': 'TROPH',
                                        'TROPPT': 'TROPH',
                                        'TROPPV': 'TROPH'}
                            hourData[varNames[iV]] = MERRA2._HfromP(self, hourData[iV])
                            c += 1
                    if c == 0:
                        print('\n!EcoSysEM.Error: missing required variable missing for atmospheric altitude(s) - \'PS\', \'TROPPB\', \'TROPPT\', or \'TROPPV\'.')
                        return None
                # Lapse rate
                varReqLR = ['T2M', 'TROPT', 'H', 'TROPH']
                if np.any(np.char.find(var, 'LR') > -1):
                    if np.all([np.any(np.char.find(var, iVarReq) > -1) for iVarReq in varReqLR]):
                        hourData['LR'] = MERRA2._LR(self, hourData['T2M'], hourData['TROPT'], hourData['H'], hourData['TROPH'])
                    else:
                        print('\n!EcoSysEM.Error: missing required variable missing for lapse rate - \'T2M\', \'TROPT\', \'H\', or \'TROPH\'.')
                        return None
                for iV in var:
                    # Daily averages and std
                    av_dayData[iV] = np.average(hourData[iV], axis = 0)
                    av_dayData[f'{iV}_std'] = np.std(hourData[iV], axis = 0)
                    # Day matrices
                    if date == start_date:
                        dayData[iV] = av_dayData[iV]
                    else:
                        dayData[iV] = np.dstack((dayData[iV], av_dayData[iV]))
                # Save daily data
                if daily == True:
                    MERRA2._saveNPZMERRA2(self, data = av_dayData, dataType = 'dly', y = y, m = m, d = int(date.day))
                # Next date
                date += delta
            # Month matrices
            monthData['lat'] = lat
            monthData['lon'] = lon
            for iV in var:
                if last_day != 1:
                    monthData[iV] = np.average(dayData[iV], axis = -1)
                    monthData[f'{iV}_std'] = np.std(dayData[iV], axis = -1)
                else:
                    monthData[iV] = dayData[iV]
                    monthData[f'{iV}_std'] = np.array([0.0])
            # Save numpy matrices in .npz format
            MERRA2._saveNPZMERRA2(self, data = monthData, dataType = 'mly', y = y, m = m)
        print("--- %s seconds ---" % (time.time() - start_time))
    
    def combDataMERRA2(self, years, month, dataType, keys = 'All'):
        """
        Get average and standard deviation from a group of data.
        
        Parameters
        ----------
        years : INT or LIST of INT
            Year(s) of data.
        month : INT or LIST of INT
            Month of data
        dataType   : STR ('mly')
            Type of data.
        keys  : LIST of STR
            List of requested variables. (Default: 'All')     
        
        """
        if not isinstance(month, int):
            print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
            return None
        if isinstance(years, int): years = [years]
        if isinstance(years, float): years = [years]
        if len(years) <= 1:
            print('Error: Introduce at least 2 years.')
            return None
        # Get all files from data\npz
        folder = f'data/MERRA2/{dataType}/'
        allFiles = np.array(os.listdir(folder))
        # Test elements
        tEl = np.empty((0))
        for y in years:
            if dataType == 'mly':
                el = f'{y}_{month}_month.npz'
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
                keys = MERRA2.keysMERRA2(self, dataType, years[0], month)
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
        # Monthly average and std
        for key in keys:
            if np.char.find(key, '_std') == -1:
                monthData[key] = np.average(combData[key], axis = -1)
            else:
                monthData[key] = np.std(combData[key], axis = -1)
            if key == 'lat' or key == 'lon':
                monthData[key] = np.squeeze(monthData[key])
        # Save numpy matrices in .npz format (v2)
        MERRA2._saveNPZMERRA2(self, data = monthData, dataType = 'cmly', y = [years[0], years[-1]], m = month) 
        
    def selectRegion(self, data, bbox):
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
            return None
        # Check user coordinates
        if bbox[0] > bbox[2]:
            print('\n!EcoSysEM.Error: `upper_right_longitude` (bbox[2]) must be higher than `lower_left_longitude` (bbox[0])')
            return None
        if bbox[1] > bbox[3]:
            print('\n!EcoSysEM.Error: `upper_right_latitude` (bbox[3]) must be higher than `lower_left_latitude` (bbox[1])')
            return None
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
                   int(np.argwhere(lonR == bbox[2])),
                   int(np.argwhere(latR == bbox[3])))
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
            return None
    
    def _saveNPZMERRA2(self, data, dataType, y, m, d = None):
        """
        Create .npz file with downladed data.

        Parameters
        ----------
        date : DICT
            Data in dictionary form.
        dataType   : STR ('mly' or 'cmly')
            Type of data.
        y : INT or LIST of INT
            Year(s) of data.
        m : INT or LIST of INT
            Month of data  
        
        """
        # Path check and folder creation if this does not exist
        path = f'data/MERRA2/{dataType}/'
        if not os.path.isdir(path):
            os.makedirs(path)
        # File name (based on dataType)
        if dataType == 'dly':
            if not isinstance(y, int):
                print('\n!EcoSysEM.Error: argument \'y\' must be a integer')
                return None
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
                return None
            if not isinstance(d, int):
                print('\n!EcoSysEM.Error: argument \'d\' must be a integer')
                return None
            file = f'{y}_{m}_{d}_day.npz'
        elif dataType == 'mly':
            if not isinstance(y, int):
                print('\n!EcoSysEM.Error: argument \'y\' must be a integer')
                return None
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
                return None
            file = f'{y}_{m}_month.npz'
        elif dataType == 'cmly':
            if not isinstance(y, list):
                print('\n!EcoSysEM.Error: argument \'y\' must be a list: [start_year, end_year]')
                return None
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
                return None
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
        dataType   : STR ('mly' or 'cmly')
            Type of data.
        y : INT or LIST of INT
            Year(s) of data.
        m : INT or LIST of INT
            Month of data  
        
        """
        path = f'data/MERRA2/{dataType}/'
        if dataType == 'dly':
            if not isinstance(y, int):
                print('\n!EcoSysEM.Error: argument \'y\' must be a integer')
                return None
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
                return None
            if not isinstance(d, int):
                print('\n!EcoSysEM.Error: argument \'d\' must be a integer')
                return None
            file = f'{y}_{m}_{d}_day.npz'
        if dataType == 'mly':
            if not isinstance(y, int):
                print('\n!EcoSysEM.Error: argument \'y\' must be a integer')
                return 0
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
                return 0
            file = f'{y}_{m}_month.npz'
        elif dataType == 'cmly':
            if not isinstance(y, list):
                print('\n!EcoSysEM.Error: argument \'y\' must be a list: [start_year, end_year]')
                return None
            if not isinstance(m, int):
                print('\n!EcoSysEM.Error: argument \'m\' must be a integer')
                return None
            file = f'{y[0]}_{y[-1]}_{m}.npz'
        return np.load(path + file)
    
    def dictMERRA2(self, dataType, y, m, d = None, keys = 'All'):
        """
        Get data in dictionary form.

        Parameters
        ----------
        dataType : STR ('mly' or 'cmly')
            Type of data.
        y : INT or LIST of INT
            Year(s) of data.
        m : INT or LIST of INT
            Month of data  
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
    
    def keysMERRA2(self, dataType, y, m, d = None):
        """
        Get variable list of data.

        Parameters
        ----------
        dataType : STR ('mly' or 'cmly')
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
        dataType : STR ('mly' or 'cmly')
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
    Description of class CAMS.
    """
    def _testCAMS(self):
        print('CAMS class works.')

class ISAMERRA2(ISA, MERRA2):
    """
    Combination of International Standard Atmosphere (ISA) model and Modern-Era 
    Retrospective analysis for Research and Applications Version 2 (MERRA-2).
    """
    def __init__(self, layers = 0, H2O = 0.0, pH = 8.0, resolution = 1000):
        ISA.__init__(self, layers = layers, H2O = H2O, pH = pH, resolution = resolution)
        MERRA2.__init__(self)
    
    def computeTandP(self, PS, TS, LR, HS, TROPH, num = 50):
        """
        Compute the change of temperature and pressure of the Earth's
        atmosphere over the range of altitudes. Based on ISA (ISO 2533:1975).
        
        Parameters
        ----------
        PS : FLOAT, LIST or ndarray
            Surface pressure. [Pa]
        TS : FLOAT, LIST or ndarray
            Surface temperature. [K]
        LR : FLOAT, LIST or ndarray
            Atmospheric lapse rate. [K/km]
        HS : FLOAT, LIST or ndarray
            Surface altitude. [m]
        TROPH : FLOAT, LIST or ndarray
            Tropopause altitude. [m]
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
        if not isinstance(PS, np.ndarray): PS = np.asarray(PS)
        if PS.ndim == 1: PS = np.array([PS])
        if not isinstance(TS, np.ndarray): TS = np.asarray(TS)
        if TS.ndim == 1: TS = np.array([TS])
        if not isinstance(LR, np.ndarray): LR = np.asarray(LR)
        if LR.ndim == 1: LR = np.array([LR])
        if not isinstance(HS, np.ndarray): HS = np.asarray(HS)
        if HS.ndim == 1: HS = np.array([HS])
        if not isinstance(TROPH, np.ndarray): TROPH = np.asarray(TROPH)
        if TROPH.ndim == 1: TROPH = np.array([TROPH])
        # Constants
        R = 8.3144598                   # Universal gas constant [J/mol/K]
        g0 = 9.80665                    # Gravitational acceleration [m/s^2]
        M0 = 0.0289644                  # Molar mass of Earth's air
        # Altitude. Shape: (alt, lat, lon)
        H = np.linspace(start = HS, stop = TROPH, num = num)                    # [m]
        # 3D matrix creation (TS, LR, HS)
        TS = np.repeat(TS[np.newaxis, :, :], H.shape[0], axis = 0)              # [K]
        LR = np.repeat(LR[np.newaxis, :, :], H.shape[0], axis = 0) / 1000       # [K/m]
        HS = np.repeat(HS[np.newaxis, :, :], H.shape[0], axis = 0)              # [m]
        # Temperature profile. Shape: (alt, lat, lon)
        T = TS + LR * (H - HS)
        # Pressure profile. Shape: (alt, lat, lon)
        P = PS * (1 + ((LR) / (TS)) * (H - HS)) ** (-(g0 * M0) / (R * LR))
        return np.squeeze(T), np.squeeze(P), np.squeeze(H)

class CAMSMERRA2(CAMS, MERRA2):
    """
    Combination of Modern-Era Retrospective analysis for Research and 
    Applications Version 2 (MERRA-2) and Copernicus Atmosphere Monitorinc
    Service (CAMS) database.
    """
    def __init__(self):
        CAMS.__init__(self)
        MERRA2.__init__(self)
