# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 12:23:23 2024

@author: 2424069M
"""

from envdef import MERRA2, ISAMERRA2, CAMSMERRA2
import argparse

# Parse command line input

print('> Available functions: getDataMERRA2')
function = input('>> Enter the function: ')

if function == 'getDataMERRA2':
    parser = argparse.ArgumentParser(
        description="Execute getDataMERRA2() via Command Line Interface (CLI).")
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
    parser.add_argument('-bbox', default='(-180, -90, 180, 90)', nargs = '+', type = float,
                        help="[tuple] (Default: '(-180, -90, 180, 90)') Earths region of data, the bounding box `-bbox lower_left_lon lower_left_lat upper_right_lon upper_right_lat`.")
    parser.add_argument('-var', default="PS T2M TROPT TROPPB", nargs = '+', type = str,
                        help="[list of str] (Default: ['PS', 'T2M', 'TROPT', 'TROPPB']) List of requested variables.")
    parser.add_argument('--daily', default = False, action = 'store_true',
                        help="[bool] (Default: False) Daily data is saved.")
    # Argument definition
    args = parser.parse_args()
    years = args.y; years.sort()
    months = args.m; months.sort()
    days = args.d
    if days != 'All':
        days = [int(days[0])]
    product = args.product
    version = args.version
    bbox = tuple(args.bbox)
    var = args.var
    if type(var) == str:
        var = var.split()
    daily = args.daily
    # Call function
    MERRA2 = MERRA2()
    MERRA2.getDataMERRA2(years = years,
                         months = months,
                         days = days, 
                         product = product, 
                         version = version, 
                         bbox = bbox, 
                         var = var, 
                         daily = daily)