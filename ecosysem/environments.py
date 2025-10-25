# -*- coding: utf-8 -*-
"""
Created on Fri Aug 29 12:00:06 2025

@author: emartinez

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

class Environment:
    def _checkData(self, dataType, y, m = None, d = None):
        """
        Check data.

        Parameters
        ----------
        dataType   : STR ('dly', 'mly', 'cmly', 'yly', 'cyly')
            Type of data.
        y : INT or LIST of INT
            Year(s) of data.
        m : INT or LIST of INT, optional
            Month of data. The default is None.
        d : INT, optional
            Day of data. The default is None.

        Returns
        -------
        None.

        """
        path = f'data/{self.model}/{dataType}/'
        if dataType == 'dly':
            if not isinstance(y, int): raise ValueError('Argument \'y\' must be an integer.')
            if not isinstance(m, int): raise ValueError('Argument \'m\' must be an integer.')
            if not isinstance(d, int): raise ValueError('Argument \'d\' must be an integer.')
            fileName = f'{y}_{m}_{d}_day.npz'
        elif dataType == 'mly':
            if not isinstance(y, int): raise ValueError('Argument \'y\' must be an integer.')
            if not isinstance(m, int): raise ValueError('Argument \'m\' must be an integer.')
            fileName = f'{y}_{m}_month.npz'
        elif dataType == 'yly':
            if not isinstance(y, int): raise ValueError('Argument \'y\' must be an integer.')
            fileName = f'{y}_year.npz'
        elif dataType == 'cmly':
            if not isinstance(y, list): raise ValueError('Argument \'y\' must be a list: [start_year, end_year].')
            if len(y) != 2: raise ValueError('Argument \'y\' must be a list: [start_year, end_year].')
            if not isinstance(m, int):  raise ValueError('Argument \'m\' must be an integer.')
            fileName = f'{y[0]}_{y[1]}_{m}.npz'
        elif dataType == 'cyly':
            if not isinstance(y, list): raise ValueError('Argument \'y\' must be a list: [start_year, end_year].')
            if len(y) != 2: raise ValueError('Argument \'y\' must be a list: [start_year, end_year].')
            fileName = f'{y[0]}_{y[1]}.npz'
        else: raise ValueError(f"Unknown dataType ({dataType}). Existing dataType: 'dly', 'mly', 'yly', 'cmly', 'cyly'.")
        path += fileName
        if not os.path.isfile(path):
            raise OSError(f"File {path} not found. Please check dataType or download requested data (see README document).")

    def _openNPZ(self, dataType, y, m = None, d = None):
        """
        Open .npz file with downladed data.

        Parameters
        ----------
        dataType   : STR ('dly', 'mly', 'cmly', 'yly', 'cyly')
            Type of data.
        y : INT or LIST of INT
            Year(s) of data.
        m : INT or LIST of INT, optional
            Month of data. The default is None.
        d : INT, optional
            Day of data. The default is None.
        
        Returns
        -------
        None.
        
        """
        path = f'data/{self.model}/{dataType}/'
        if dataType == 'dly':
            file = f'{y}_{m}_{d}_day.npz'
        elif dataType == 'mly':
            file = f'{y}_{m}_month.npz'
        elif dataType == 'cmly':
            file = f'{y[0]}_{y[-1]}_{m}.npz'
        elif dataType == 'yly':
            file = f'{y}_year.npz'
        elif dataType == 'cyly':
            file = f'{y[0]}_{y[-1]}.npz'
        else:
            raise ValueError('Unknown dataType ({dataType}). Existing dataType: \'dly\', \'mly\', \'yly\', \'cmly\', \'cyly\'.')
        return np.load(path + file)
    
    def _saveNPZ(self, data, dataType, y, m = None, d = None):
        """
        Create .npz file with downladed data.

        Parameters
        ----------
        data : DICT
            Data in dictionary form.
        dataType : STR ('dly', 'mly', 'cmly', 'yly', 'cyly')
            Type of data.
        y : INT or LIST of INT
            Year(s) of data.
        m : INT or LIST of INT, optional
            Month of data.  The default is None.
        d : INT, optional
            Day of data. The default is None.
        
        Returns
        -------
        None.
        
        """
        # Path check and folder creation if this does not exist
        path = f'data/{self.environment}/{self.model}/{dataType}/'
        if not os.path.isdir(path):
            os.makedirs(path)
        # File name (based on dataType)
        if dataType == 'dly':
            if not isinstance(y, int): raise ValueError('Argument \'y\' must be an integer.')
            if not isinstance(m, int): raise ValueError('Argument \'m\' must be an integer.')
            if not isinstance(d, int): raise ValueError('Argument \'d\' must be an integer.')
            file = f'{y}_{m}_{d}_day.npz'
        elif dataType == 'mly':
            if not isinstance(y, int): raise ValueError('Argument \'y\' must be an integer.')
            if not isinstance(m, int): raise ValueError('Argument \'m\' must be an integer.')
            file = f'{y}_{m}_month.npz'
        elif dataType == 'cmly':
            if not isinstance(y, list): raise ValueError('Argument \'y\' must be a list: [start_year, end_year].')
            if len(y) != 2: raise ValueError('Argument \'y\' must be a list: [start_year, end_year].')
            if not isinstance(m, int):  raise ValueError('Argument \'m\' must be an integer.')
            file = f'{y[0]}_{y[-1]}_{m}.npz'
        elif dataType == 'yly':
            if not isinstance(y, int): raise ValueError('Argument \'y\' must be an integer.')
            file = f'{y}_year.npz'
        elif dataType == 'cyly':
            if not isinstance(y, list): raise ValueError('Argument \'y\' must be a list: [start_year, end_year].')
            if len(y) != 2: raise ValueError('Argument \'y\' must be a list: [start_year, end_year].')
            file = f'{y[0]}_{y[-1]}.npz'
        else:
            raise ValueError('Unknown dataType ({dataType}). Existing dataType: \'dly\', \'mly\', \'yly\', \'cmly\', \'cyly\'.')
        # Path generation
        pathfile = path + file
        # Save .npz file
        np.savez(pathfile, **data)
    
    def _reshapeData(self, data):
        """
        Reshape data with ascending latitude and longitude order.

        Parameters
        ----------
        data : DICT
            Data to be rehsaped.

        Returns
        -------
        data : DICT
            Reshaped data.

        """
        # Check latitude
        try: data['lat']
        except: raise KeyError('Latitude (lat) is missing in data.')
        # Check longitude
        try: data['lon']
        except: raise KeyError('Longitude (lon) is missing in data.')
        lats = data['lat']
        lons = data['lon']
        # Ensure ascending latitude order
        if lats[0] > lats[-1]:
            lats = lats[::-1]
            data['lat'] = lats
            for key in data:
                if isinstance(data[key], list): data[key] = np.array(data[key])
                if not isinstance(data[key], np.ndarray): continue
                if data[key].ndim == 3:
                    data[key] = data[key][:, ::-1, :]
                elif data[key].ndim == 2:
                    data[key] = data[key][::-1, :]
                else: continue
        # Ensure ascending lontidue order
        if lons[0] > lons[-1]:
            lons = lons[::-1]
            data['lon'] = lons
            for key in data:
                if isinstance(data[key], list): data[key] = np.array(data[key])
                if not isinstance(data[key], np.ndarray): continue
                if data[key].ndim == 3:
                    data[key] = data[key][:, :, ::-1]
                elif data[key].ndim == 2:
                    data[key] = data[key][:, ::-1]
                else: continue
        return data
    
    def getAttributeNames(self):
        """
        Return attribute names of an Environment object as a list.
        
        """
        return list(self.__dict__.keys())
    
    def loadData(self, dataType, y, m = None, d = None, keys = 'All'):
        """
        Get data in dictionary form.

        Parameters
        ----------
        dataType : STR ('dly', 'mly', 'cmly', 'yly', 'cyly')
            Type of data.
        y : INT or LIST of INT
            Year(s) of data.
        m : INT or LIST of INT, optional
            Month of data. The default is None
        d : INT or LIST of INT, optional
            Day(s) of data. The default is None.
        keys : LIST of STR
            List of requested variables. The default is 'All'.
        
        Returns
        -------
        dictVar : DICT
            Dictionary with requested variables.

        """
        npz = Environment._openNPZ(self, dataType, y, m, d)
        latlon_models = ['MERRA2', 'CAMS', 'ISAMERRA2', 'CAMSMERRA2']
        if (self.environment == 'Atmosphere') and (self.model in latlon_models): must_have_lat, must_have_lon = True, True
        if keys == 'All':
            keys = Environment.keys(self, dataType, y, m, d)
            dictVar = {key: npz[key] for key in keys}
        else:
            coor = []
            if not np.any(np.char.find('lat', keys) > -1) and must_have_lat is True:
                coor += ['lat']
            if not np.any(np.char.find('lon', keys) > -1) and must_have_lon is True:
                coor += ['lon']
            keys = coor + keys
            dictVar = {key: npz[key] for key in keys}
        npz.close()
        return dictVar    
    
    def combData(self, dataType, year, month, days = None, keys = 'All', dataDelete = False, method = 'linear', target_lats = None, target_lons = None):
        """
        Get average and standard deviation from a group of data.
        
        Parameters
        ----------
        dataType : STR ('cmly', 'mly', 'yly', 'cyly')
            Type of data.
        years : INT or LIST of INT
            Year(s) of data.
        month : INT or LIST of INT
            Month of data.
        days : INT, optional
            Last day of data.
        keys : LIST of STR, optional
            List of requested variables. The default is 'All'.
        dataDelete : BOOL
            Delete daily or monthly data after the average calculation. The default is False.
        
        Returns
        -------
        None
        
        """
        if self.model == 'CAMS':
            CAMS._combDataCAMS(self, dataType=dataType, 
                               years = year, 
                               months = month, 
                               method = method,
                               target_lats = target_lats, 
                               target_lons = target_lons)
        else:
            # Get all files from `data\`
            if dataType == 'cmly' or dataType == 'yly' or dataType == 'cyly':
                folder = 'data/{self.model}/mly/'
            elif dataType == 'mly':
                folder = 'data/{self.model}/dly/'
            allFiles = np.array(os.listdir(folder))
            # Test elements
            tEl = np.empty((0))
            if dataType == 'cmly':
                if not isinstance(year, list): raise ValueError('Argument \'y\' must be a list with length > 1')
                if len(year) <= 1: raise ValueError('Argument \'y\' must be a list with length > 1.')
                years = np.arange(year[0], year[-1]+1, 1)
                for y in years:    
                    el = f'{y}_{month}_month.npz'
                    tEl = np.append(tEl, el)
            elif dataType == 'mly':
                if not isinstance(year, int): raise ValueError('Argument \'year\' must be an integer')
                if not isinstance(month, int): raise ValueError('Argument \'month\' must be an integer')
                y = year
                m = month
                if not days: raise ValueError('Argument \'days\' must be defined: last day as integer or \'All\'')
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
            elif dataType == 'yly':
                if not isinstance(month, list): raise ValueError('Argument \'month\' must be a list: [start_month, end_month]')
                if len(month) != 2: raise ValueError('Argument \'month\' must be a list: [start_month, end_month]')
                if not isinstance(year, int): raise ValueError('Argument \'year\' must be an integer.')
                months = np.arange(month[0], month[-1]+1, 1)
                y = year
                for m in months:
                    el = f'{y}_{m}_month.npz'
                    tEl = np.append(tEl, el)
            elif dataType == 'cyly':
                if not isinstance(year, list): raise ValueError('Argument \'year\' must be a list: [start_year, end_year]')
                if not isinstance(month, list): raise ValueError('Argument \'month\' must be a list: [start_month, end_month]')
                years = np.arange(year[0], year[-1]+1, 1)
                months = np.arange(month[0], month[-1]+1, 1)            
                for y in years:
                    for m in months:
                        el = f'{y}_{m}_month.npz' 
                        tEl = np.append(tEl, el)
            else:
                raise ValueError('Unknown dataType ({dataType}). Existing dataType: \'mly\' (to generate monthly data), \'cmly\' (to generate combined monthly data),  \'yly\' (to generate annual data) or \'cyly\' (to generate combined annual data).')
            selFiles = allFiles[np.isin(allFiles, tEl)]
            # Stack matrices
            combData = {}
            resultData = {}
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
                if '_std' in key:
                    key_ = key.replace('_std', '')
                    resultData[key] = np.nanstd(combData[key_], axis = -1)
                else:
                    resultData[key] = np.nanmean(combData[key], axis = -1)
                if key == 'lat' or key == 'lon':
                    resultData[key] = np.squeeze(resultData[key])
            # Save numpy matrices in .npz format (v2)
            if dataType == 'cmly':
                Environment._saveNPZ(self, data = resultData, dataType = 'cmly', y = [year[0], year[-1]], m = month)
            elif dataType == 'mly':
                Environment._saveNPZ(self, data = resultData, dataType = 'mly', y = y, m = m)
            elif dataType == 'yly':
                Environment._saveNPZ(self, data = resultData, dataType = 'yly', y = y, m = None)
            elif dataType == 'cyly':
                Environment._saveNPZ(self, data = resultData, dataType = 'cyly', y = [year[0], year[-1]], m = None)
            # Delete the data used (if necessary)
            if dataDelete:
                for file in selFiles:
                    path = folder + file
                    os.remove(path)
    
    def keys(self, dataType, y, m = None, d = None):
        """
        Get variable list of data.

        Parameters
        ----------
        dataType : STR ('dly', 'mly', 'cmly', 'yly', 'cyly')
            Type of data.
        y : INT or LIST of INT
            Year(s) of data.
        m : INT or LIST of INT, optional
            Month of data. The default is None.
        d : INT, optional
            Day of data. The default is None.
        
        Returns
        -------
        keys  : LIST of STR
            Variable list of data.

        """
        npz = Environment._openNPZ(self, dataType, y, m, d)
        keys = npz.files
        npz.close()
        return keys

# Atmosphere ------------------------------------------------------------------
class Atmosphere(Environment):
    def _plotAtmosphere():
        pass
    
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
    
    def _selectRegion(self, data, bbox):
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
        else: raise ValueError('Argument bbox requires 2 `(lon, lat)` or 4 `(lower_left_lon, lower_left_lat, upper_right_lon, upper_right_lat)` positional values.')
        # Check user coordinates
        if bbox[0] > bbox[2]: raise ValueError(f'Argument bbox - `upper_right_longitude` ({bbox[2]}) must be higher than `lower_left_longitude` ({bbox[0]}).')
        if bbox[1] > bbox[3]: raise ValueError(f'Argument bbox - `upper_right_latitude` ({bbox[3]}) must be higher than `lower_left_latitude` ({bbox[1]}).')
        # BBox from data
        lonR = data['lon']
        latR = data['lat']
        bboxData = (lonR[0], latR[0], lonR[-1], latR[-1])
        # Get closest coordinates
        bbox =  Atmosphere._closestCoord(self, lonR, latR, bbox)
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
        else: raise ValueError('Argument bbox - selected region ({bbox}) is outside the data boundaries ({bboxData}).')
    
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
        return bbox

class ISA(Atmosphere):
    """
    International Standard Atmosphere (ISA). Static atmospheric model of how
    pressure, temperature, density and viscosity of the Earth's atmosphere
    change over a wide range of altitudes.
    --------------------------------------------------------------------------
    References: - ISO 2533:1975
                - National Oceanic and Atmospheric Administration (NOAA)
                - Goody & Yung (1989), doi: 10.1093/oso/9780195051346.001.0001
    
    """
    def __init__(self, layers = 'All', phase = 'All', H2O = 0.0, pH = 7.0, selCompounds = None, 
                 selAlt = None, resolution = 1000, showMessage = True):
        if showMessage:
            print('  > Creating ISA instance...')
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
        self.environment = 'Atmosphere'
        self.model = 'ISA'
        if selAlt:
            self._selectAltitude(selAlt)
        self._getConcISA(phase, selCompounds)
        if showMessage:
            print('  > Done.')
    
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
        if isinstance(layers, str):
            if layers == 'All':
                layers = range(8)
        elif isinstance(layers, int):
            layers = range(layers, layers+1)
        lapse_rate = (dISA.loc[layers]['Lapse rate']).to_numpy()
        base_T = (dISA.loc[layers]['Base temperature']).to_numpy()
        base_P = (dISA.loc[layers]['Base pressure']).to_numpy()
        start_alt = (dISA.loc[layers]['Start altitude']).to_numpy()
        end_alt = (dISA.loc[layers]['End altitude']).to_numpy()
        # Temperature calculation (ISA)
        alt = t = p = []
        for i in range(lapse_rate.size):
            # alt_aux = np.array(range(int(start_alt[i]), int(end_alt[i]), resolution))
            alt_aux = np.arange(start_alt[i], end_alt[i], resolution)
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
        self.altitude = alt # [m]
        self.temperature = t + 273.15 # [K]
        self.pressure = p # [Pa]
    
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
            self.compositions = pd.Series(newComp, index = dDC.Compounds).to_dict()
    
    def _getConcISA(self, phase, compound = None):
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
                        'All' - All phaes (G, L-FW, L-SW).
        compound : STR or LIST, optional
                Interested compounds. The default is None. (i.e., all compounds are considered).
        """
        # Data
        p = self.pressure
        t = self.temperature   # [K]
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
            self.Pi = dict_Pi
            self.Ci_G = dict_Ci_G
            self.Ci_LFW = None
            self.Ci_LSW = None
        elif phase == 'L-FW':
            self.Pi = None
            self.Ci_G = None
            self.Ci_LFW = dict_Ci_LFW
            self.Ci_LSW = None
        elif phase == 'L-SW':
            self.Pi = None
            self.Ci_G = None
            self.Ci_LFW = None
            self.Ci_LSW = dict_Ci_LSW
        elif phase == 'L':
            self.Pi = None
            self.Ci_G = None
            self.Ci_LFW = dict_Ci_LFW
            self.Ci_LSW = dict_Ci_LSW
        elif phase == 'All':
            self.Pi = dict_Pi
            self.Ci_G = dict_Ci_G
            self.Ci_LFW = dict_Ci_LFW
            self.Ci_LSW = dict_Ci_LSW
        else: raise ValueError(f'Unknown phase ({phase}). Existing phase: \'G\' (gas), \'L-FW\' (Liquid fresh water), \'L-SW\' (Liquid sea water), \'L\' (L-FW, L-SW), \'All\'- All phases.')
    
    def _selectAltitude(self, selAlt):
        """
        Select a specific region of atmosphere based on a minimum and maximum 
        altitud value.
        
        """
        if isinstance(selAlt, list):
            if len(selAlt) > 2: raise ValueError('Argument \'selAlt\' must be an integer, float or list [min_alt, max_alt].')
            elif len(selAlt) == 1:
                selAlt = [0, selAlt[0]]
            minAlt = selAlt[0]
            maxAlt = selAlt[1]
        elif isinstance(selAlt, (int, float)):
            minAlt = 0
            maxAlt = selAlt
        # Previous altitude, temperature and pressure
        prevAlt = self.altitude
        prevT = self.temperature
        prevP = self.pressure
        # Correct min and max altitude, if out of previous range
        minAlt = max(minAlt, min(prevAlt))
        maxAlt = min(maxAlt, max(prevAlt))
        # Indexes of min and max altitudes
        imAlt = int(np.argwhere(prevAlt <= minAlt)[-1])
        iMAlt = int(np.argwhere(prevAlt >= maxAlt)[0]) + 1
        # New altitude, temperature and pressure
        self.altitude = prevAlt[imAlt:iMAlt]
        self.temperature = prevT[imAlt:iMAlt]
        self.pressure = prevP[imAlt:iMAlt]
    
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
        
    ## Plotting functions 
    def plotTandP(self):
        """
        Plotting of pressure (in atm) and temperature (in Kelvin) along the
        atmosphere (altitude in km).
        
        """
        # Variables
        alt = self.altitude / 1000      # [km]
        t = self.temperature            # [K]
        p = self.pressure / 101325      # [atm]
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
        alt = self.altitude / 1000     # [km]
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

class MERRA2(Atmosphere):
    pass

class CAMS(Atmosphere):
    pass

class ISAMERRA2(ISA, MERRA2):
    pass

class CAMSMERRA2(CAMS, MERRA2):
    pass

# Hydrosphere -----------------------------------------------------------------

class Hydrosphere(Environment):
    pass

class Ocean(Hydrosphere):
    pass

class Sediments(Hydrosphere):
    pass

class River(Hydrosphere):
    pass

# Cryosphere ------------------------------------------------------------------

class Cryosphere(Environment):
    pass

class Glacier(Cryosphere):
    pass

# Lithosphere -----------------------------------------------------------------

class Lithosphere(Environment):
    pass

class Crust(Lithosphere):
    pass
