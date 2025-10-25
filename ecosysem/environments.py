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
    pass

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
