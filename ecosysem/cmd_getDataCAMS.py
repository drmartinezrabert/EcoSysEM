# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 06:55:03 2026

@author: emartinez
"""

from environments import CAMS
import argparse
import numpy as np

parser = argparse.ArgumentParser(prefix_chars='_',
    description="Execute getDataCAMS() via Command Line Interface (CLI).")
parser.add_argument('_dataType', nargs = '+', type = str,
                    help="[str] Type(s) of data ('mly', 'dly').")
parser.add_argument('_y', nargs = '+', type = int,
                    help="[int or list of int] Year(s) of requested data.")
parser.add_argument('_m', nargs = '+', type = int,
                    help="[int or list of int] Month(s) of requested data.")
parser.add_argument('_d', default='All', nargs = '+',
                    help="[int or list of int or str ('All'))] (Default: 'All') Day(s) of month of requested data. With 'All' get the whole month.")
parser.add_argument('_hr', nargs = '+', default = [0, 12],
                    help="[int or list of int] (Default: '0 12') Hour(s) of requested data.")
parser.add_argument('_dset', type = str,
                    help="[str] Name of dataset ('cams-global-greenhouse-gas-forecasts', 'cams-global-ghg-reanalysis-egg4', 'cams-global-atmospheric-composition-forecasts').")
parser.add_argument('_pressure', default=[50, 100, 200, 400, 600, 800, 900, 1000], nargs = '+', type = int,
                    help="[int] (Default: '50 100 200 400 600 800 900 1000') Pressure levels to download.")
parser.add_argument('_variables', default= None, nargs = "+",
                    help = '[str or list of str] (Default: None) A list of variables to download (Allowed: "co", "co2", "ch4").')
parser.add_argument('_bbox', default=[90, -180, -90, 180], nargs = '+', type = int,
                    help="[list] (Default: '90 -180 -90 180') Earth's region of data, the bounding box `-bbox upper_right_lat lower_left_lon lower_left_lat upper_right_lon`.")
parser.add_argument('_mode', default = None, type = str,
                    help="[str] Mode of download ('add').")
parser.add_argument('_method', default='linear', nargs = '+', type = str,
                    help="[str] (Default: 'linear') Method of interpolation in 'add' mode.")
parser.add_argument("_dropvariables", default = None, nargs = "+",
                    help="[str or list of str] (Default: None) A variable or list of variables to exclude from being parsed from the dataset.")
# Argument definition
args = parser.parse_args()
dataType = args.dataType
if not np.all(np.isin(dataType, ['dly', 'mly'])): raise ValueError(f"Unaccepted dataType ({dataType}). Accepted dataType: 'dly', 'mly'.")
years = args.y; years.sort()
months = args.m; months.sort()
days = args.d
if not np.any(np.isin(days, 'All')):
    days = sorted([int(d) for d in days])
else:
    days = 'All'
hours = args.hr
dataset = args.dset
pressure_levels = args.pressure; pressure_levels.sort()
variables = args.variables
bbox = list(args.bbox)
mode = args.mode
method = args.method
drop_variables = args.dropvariables
# Call function
CAMS = CAMS()
CAMS.getDataCAMS(dataType = dataType,
                 years = years,
                 months = months,
                 days = days,
                 hours = hours,
                 dataset = dataset,
                 pressure_levels = pressure_levels,
                 variables = variables,
                 bbox = bbox,
                 mode = mode,
                 method = method,
                 drop_variables = drop_variables)