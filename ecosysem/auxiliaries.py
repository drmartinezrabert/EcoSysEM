# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 08:50:57 2026

@author: emartinez
"""

import re
import numpy as np
from geokernels.distance import geodist

def _grid_weights(lon, lat, alt = None, NaN_values = None):
    """
    Compute the contribution (weights) of each grid (2D or 3D) to the statistic (e.g., quantile, mean, median...).

    Parameters
    ----------
    lon : LIST or np.ndarray
        Array of longitudes.
    lat : LIST or np.ndarrray
        Array of latitudes.
    alt : LIST, optional
        Array of altitudes. The default is None.
    NaN_values : LIST or np.ndarray, optional
        Array or matrix to set NaN-value grids. The default is None.

    Returns
    -------
    weights : np.ndarray
        Array or matrix with grid weights.

    """
    #-Grid sizes
    d_lon = np.diff(lon); d_lon = np.append(d_lon, d_lon[-1])
    d_lat = np.diff(lat); d_lat = np.append(d_lat, d_lat[-1])
    if not isinstance(alt, (list, np.ndarray)):
        alt = np.array([0.0]); d_alt = np.array([1000.0, 1000.0])
    elif isinstance(alt, (list, np.ndarray)) and len(alt) == 1:
        alt = np.array(alt); d_alt = np.array([1000.0, 1000.0])
    else:
        d_alt = np.diff(alt); d_alt = np.append(d_alt, d_alt[-1])
    #-Matrices
    lons, lats = np.meshgrid(lon, lat)
    shape = lons.shape
    #-Grid coordinates
    i_west = lons - d_lon/2; i_west = np.where(i_west < -180, -180, i_west); i_west = i_west.flatten()
    i_east = lons + d_lon/2; i_east = np.where(i_east > 180, 180, i_east); i_east = i_east.flatten()
    j_south = lats.T - d_lat/2; j_south = np.transpose(np.where(j_south < -90, -90, j_south)); j_south = j_south.flatten()
    j_north = lats.T + d_lat/2; j_north = np.transpose(np.where(j_north > 90, 90, j_north)); j_north = j_north.flatten()
    if len(alt) > 1:
        k_down = alt - d_alt/2; k_down = np.where(k_down < 0, 0, k_down)
        k_up = alt + d_alt/2; k_up = np.where(k_up > max(alt), max(alt), k_up)
    else:
        k_down = alt - d_alt/2
        k_up = alt + d_alt/2
    #-Grid distances (Geodesic distances at Earth's surface, h=0)
    a = geodist(np.array(list(zip(j_south, i_east))), np.array(list(zip(j_north, i_east))), metric = 'meter') / 1000
    b = geodist(np.array(list(zip(j_north, i_east))), np.array(list(zip(j_north, i_west))), metric = 'meter') / 1000
    c = geodist(np.array(list(zip(j_north, i_west))), np.array(list(zip(j_south, i_west))), metric = 'meter') / 1000
    d = geodist(np.array(list(zip(j_south, i_west))), np.array(list(zip(j_south, i_east))), metric = 'meter') / 1000
    #-Matrix initialization
    V_grid = np.zeros((len(alt), len(lat), len(lon)))
    for id_alt, _ in enumerate(alt):
        i_kdown = k_down[id_alt] / 1000
        i_kup = k_up[id_alt] / 1000
        h = i_kup - i_kdown
        # Elevation factors
        EF_down = (6371 + i_kdown) / 6371
        EF_up = (6371 + i_kup) / 6371
        #-Grid distance for bases (A1 and A2)
        a1 = a * EF_down;  b1 = b * EF_down;  c1 = c * EF_down;  d1 = d * EF_down
        a2 = a * EF_up;    b2 = b * EF_up;    c2 = c * EF_up;    d2 = d * EF_up
        #-Areas of A1 and A2
        s1 = (a1+b1+c1+d1) / 2; A1 = np.sqrt((s1-a1)*(s1-b1)*(s1-c1)*(s1-d1))
        s2 = (a2+b2+c2+d2) / 2; A2 = np.sqrt((s2-a2)*(s2-b2)*(s2-c2)*(s2-d2))
        #-Grid volumes
        V_grid_ = (h/3) * (A1 + np.sqrt(A1*A2) + A2)
        V_grid_ = np.reshape(V_grid_, shape)
        V_grid[id_alt, ...] = V_grid_
    if isinstance(NaN_values, (list, np.ndarray)):
        V_grid = np.where(~np.isnan(NaN_values), V_grid, 0)
    return np.squeeze(V_grid / np.sum(V_grid))

def _elements(compound, minimum_elements = None):
    """
    Get elements of a compound.

    Parameters
    ----------
    compound : STR
        Molecular formula of compound.
    minimum_elements : LIST, optional
        Elements that must be check, in addition to elements of `compounds`. The default is None.

    Returns
    -------
    c_i : DICT
        Dictionary with elements of compound. Shape: {'Element': float}

    """
    def _is_decimal(string):
        dot = r'\.'
        number = r'[0-9]'
        check = bool(re.match(dot, string)) or bool(re.match(number, string))
        return check
    
    ind_v = [i for i, c in enumerate(compound) if c.isupper()]
    c_i = {}
    elements = []
    if '+' in compound:
        z_ind = compound.index('+')
        ind_v += [z_ind]
        sign = 1
    elif '-' in compound:
        z_ind = compound.index('-')
        ind_v += [z_ind]
        sign = -1
    else:
        ind_v += [len(compound)]
        sign = 0
    for ind in range(len(ind_v)-1):
        ind1, ind2 = int(ind_v[ind]), int(ind_v[ind+1])
        element = ''.join(list(filter(str.isalpha, compound[ind1:ind2])))
        elements += [element]
        ci_ = ''.join(list(filter(_is_decimal, compound[ind1:ind2])))
        if len(ci_) == 0: ci_ = 1.0
        c_i[element] = float(ci_)
    if isinstance(minimum_elements, (list, np.ndarray)):
        for other_element in minimum_elements:
            if not other_element in elements:
                c_i[other_element] = 0.0
    if sign != 0:
        if len(compound[z_ind+1:]) == 0:
            c_i['z'] = sign * 1.0
        else:
            c_i['z'] = sign * float(compound[z_ind+1:])
    else:
        c_i['z'] = 0.0
    return c_i

def _check_equation(equation, tol_decimals = 5, print_balance = False):
    """
    Check if an equation is well balanced.

    Parameters
    ----------
    equation : DICT
        Dictionary with stoichiometric parameters of equations. Negative values for substrates,
        and positive values for products. Shape: {'v_Compound': float}.
    tol_decimals : int, optional
        Maximum number of decials considered in stoichometries. The default is 5.
    print_balance : BOOL, optional
        Set whether resultance element balance is printed in console. The default is False.

    Returns
    -------
    check_balance : BOOL
        True if equation is well balanced and False when not.

    """
    if not isinstance(equation, dict): raise TypeError("Equation must be given as a dictionary: {'v_A': float, 'v_B': float}")
    check_equation = {}
    for v_i in equation:
        i = v_i.replace('v_', '')
        c_i = _elements(i)
        c_i.update((x, equation[v_i] * y) for x, y in c_i.items())
        check_equation = {i: round(check_equation.get(i, 0) + c_i.get(i, 0), tol_decimals) 
                          for i in set(check_equation).union(c_i)}
    if print_balance:
        print(f'· Balance: {check_equation}')
    check_balance = sum(list(check_equation.values())) == 0.0
    return check_balance