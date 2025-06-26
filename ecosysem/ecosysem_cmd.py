# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 12:23:23 2024

@author: 2424069M
"""

from envdef import MERRA2, CAMS, ISAMERRA2, CAMSMERRA2
import argparse
import numpy as np
import sys

# Parse command line input

print('> Available functions: getDataMERRA2, getDataCAMS')
function = input('>> Enter the function: ')

if function == 'getDataMERRA2':
    parser = argparse.ArgumentParser(
        description="Execute getDataMERRA2() via Command Line Interface (CLI).")
    parser.add_argument('-dataType', nargs = '+', type = str,
                        help="[str or list] Type of data (dly: daily; mly: monthly; cmly: combined monthly; All: dly mly cmly).")
    parser.add_argument('-y', nargs = '+', type = int,
                        help="[int] Year of requested data.")
    parser.add_argument('-m', nargs = '+', type = int,
                        help="[int or list] Month(s) of requested data.")
    parser.add_argument('-d', default='All', nargs = '+',
                        help="[str or int] (Default: 'All') Last day of month of requested data. With 'All' get the whole month.")
    parser.add_argument('-product', default='M2I1NXASM',
                        help="[str] (Default: 'M2I1NXASM') Product of data (section of MERRA2 database).")
    parser.add_argument('-version', default='5.12.4',
                        help="[str] (Default: '5.12.4') Version of data.")
    parser.add_argument('-bbox', default=(-180, -90, 180, 90), nargs = '+', type = float,
                        help="[tuple] (Default: '-180 -90 180 90') Earths region of data, the bounding box `-bbox lower_left_lon lower_left_lat upper_right_lon upper_right_lat`.")
    parser.add_argument('-var', default="PS T2M TROPT TROPPB", nargs = '+', type = str,
                        help="[list of str] (Default: PS T2M TROPT TROPPB) List of requested variables.")
    # Argument definition
    args = parser.parse_args()
    dataType = args.dataType
    if not np.all(np.isin(dataType, ['dly', 'mly', 'cmly', 'All'])):
        print('\n!EcoSysEM.Error: dataType not found. Data type must be "dly", "mly", "cmly", list of data types or "All".')
        sys.exit()
    years = args.y; years = sorted(years)
    months = args.m; months = sorted(months)
    days = args.d
    if not np.any(np.isin(days, 'All')):
        days = sorted([int(d) for d in days])
    else:
        days = 'All'
    product = args.product
    version = args.version
    bbox = args.bbox
    var = args.var
    if type(var) == str:
        var = var.split()
    # Call function
    MERRA2 = MERRA2()
    MERRA2.getDataMERRA2(dataType = dataType,
                         years = years,
                         months = months,
                         days = days, 
                         product = product, 
                         version = version, 
                         bbox = bbox, 
                         var = var)

elif function == 'getDataCAMS':
    parser = argparse.ArgumentParser(
        description="Execute getDataCAMS() via Command Line Interface (CLI).")
    parser.add_argument('-type', nargs = '+', type = str,
                        help="[str] Type(s) of data ('mly', 'dly', 'cmly').")
    parser.add_argument('-y', nargs = '+', type = int,
                        help="[int or list of int] Year(s) of requested data.")
    parser.add_argument('-m', nargs = '+', type = int,
                        help="[int or list of int] Month(s) of requested data.")
    parser.add_argument('-d', default='All', nargs = '+',
                        help="[int or list of int or str ('All'))] (Default: 'All') Day(s) of month of requested data. With 'All' get the whole month.")
    parser.add_argument('-pressure', default="100 200 300 400 500 600 700 800 850 900 950 1000", nargs = '+', type = int,
                        help="[int] (Default: '100 200 300 400 500 600 700 800 850 900 950 1000') Pressure levels to download.")
    parser.add_argument('-bbox', default="90 -180 -90 180", nargs = '+', type = float,
                        help="[list] (Default: '90 -180 -90 180') Earth's region of data, the bounding box `-bbox upper_right_lat lower_left_lon lower_left_lat upper_right_lon`.")
    # Argument definition
    args = parser.parse_args()
    dataType = args.type
    years = args.y; years.sort()
    months = args.m; months.sort()
    days = args.d
    pressure_levels = args.pressure; pressure_levels.sort()
    bbox = list(args.bbox)
    # Call function
    CAMS = CAMS()
    CAMS.getDataCAMS(dataType = dataType,
                     years = years,
                     months = months,
                     days = days,
                     pressure_levels = pressure_levels,
                     bbox = bbox)
