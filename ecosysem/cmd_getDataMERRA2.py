# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 06:52:18 2026

@author: emartinez
"""

from environments import MERRA2
import argparse
import numpy as np

parser = argparse.ArgumentParser(prefix_chars='_',
    description="Execute getDataMERRA2() via Command Line Interface (CLI).")
parser.add_argument('_dataType', nargs = '+', type = str,
                    help="[str or list] Type of data (dly: daily; mly: monthly; cmly: combined monthly; All: dly mly cmly).")
parser.add_argument('_y', nargs = '+', type = int,
                    help="[int] Year of requested data.")
parser.add_argument('_m', nargs = '+', type = int,
                    help="[int or list] Month(s) of requested data.")
parser.add_argument('_d', default='All', nargs = '+',
                    help="[str, int or list] (Default: 'All') Last day of month of requested data. With 'All' get the whole month.")
parser.add_argument('_product', default='M2I1NXASM',
                    help="[str] (Default: 'M2I1NXASM') Product of data (section of MERRA2 database).")
parser.add_argument('_version', default='5.12.4',
                    help="[str] (Default: '5.12.4') Version of data.")
parser.add_argument('_bbox', default=(-180, -90, 180, 90), nargs = '+', type = float,
                    help="[tuple] (Default: '-180 -90 180 90') Earths region of data, the bounding box `-bbox lower_left_lon lower_left_lat upper_right_lon upper_right_lat`.")
parser.add_argument('_var', default="PS TROPPB T2M TROPT TROPH LR", nargs = '+', type = str,
                    help="[list of str] (Default: PS TROPPB T2M TROPT TROPH LR) List of requested variables.")
# Argument definition
args = parser.parse_args()
dataType = args.dataType
if not np.all(np.isin(dataType, ['dly', 'mly', 'cmly', 'All'])): raise ValueError(f"Unaccepted dataType ({dataType}). Accepted dataType: 'dly', 'mly', 'cmly', 'All'.")
years = args.y; years = sorted(years)
months = args.m; months = sorted(months)
days = args.d
if not np.any(np.isin(days, 'All')):
    days = sorted([int(d) for d in days])
else:
    days = 'All'
product = args.product
version = args.version
bbox = tuple(args.bbox)
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