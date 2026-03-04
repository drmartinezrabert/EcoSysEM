# -*- coding: utf-8 -*-
"""
Created on Mon Nov 10 17:12:16 2025

@author: emartinez
"""

from reactions import Reactions as Rxn
from thermodynamics import ThEq
from environments import Atmosphere, ISAMERRA2, CAMSMERRA2
from auxiliaries import _grid_weights as wt

import datetime
from dateutil.relativedelta import relativedelta
from mpl_toolkits.basemap import Basemap
import pylab as pl
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap
import matplotlib.ticker as tkr
import os

def _savePlot(file, cPlots):
    if os.path.isfile(f'{file}.tiff'):
        while os.path.isfile(f'{file}_{cPlots}.tiff'):
            cPlots += 1
        plt.savefig(f'{file}_{cPlots}.tiff', bbox_inches='tight')
    else:
        plt.savefig(f'{file}.tiff', bbox_inches='tight')

def plot_seasonality(model, dataType, start_date, end_date, variable, delta_time = 1, pH = 7.0,
                     bbox = (-180, -90, 180, 90), compound = None, phase = 'All', altArray = None,
                     numAlt = 50, surftrop = None, figsize = None, lw = 2.0, color = 'deepskyblue', 
                     marker = None, ms = 5, logScale = False, alpha_fillbtw = 0.3, date_format = None,
                     yticks_format = None, xlabel = 'Time', ylabel = 'Variable (Units)', xlabel_rotation = 0.0, 
                     yLim = [None, None], show_right_labels = True, drl = 1.005, fsrl = 10.0, fontFamily = 'Arial', 
                     fontSize = 12.0, title = None, cf = 1.0, lines = [True, True, True, True, True], 
                     showMessage = False, fillMissing = False, weights = False, savePlot = False):
    """
    Plot seasonal variability of a variable, including all quartiles (0.0, 0.25, 0.5, 0.75, 1.0).

    Parameters
    ----------
    model : STR
        Environmental model. Available models: 'ISAMERRA2', 'CAMSMERRA2'.
    dataType : STR
        Type of data. Available types: 'dly', 'mly', 'yly', 'cmly'
    start_date : STR
        Start date of seasonality plot. Format: dly - 'yyyy-mm-dd'; mly - 'yyyy-mm'; cmly - 'yyyy-mm'; yly - 'yyyy'.
    end_date : STR
        Last date of seasonality plot. Format: dly - 'yyyy-mm-dd'; mly - 'yyyy-mm'; cmly - 'yyyy-mm'; yly - 'yyyy'.
    variable : STR
        Variable to be plotted. Some variables are dictionaries, and compound has to be given: attributeName-compound (e.g., Ci_LSW-CH4).
    delta_time : INT, optional
        Set time interval. The default is 1.
    pH : INT or FLOAT, optional
        pH value. The default is 7.0.
    bbox : Tuple
        Earths region of data, the bounding box. The default is (-180, -90, 180, 90).
        (lower_left_longitude, lower_left_latitude, upper_right_longitude, upper_right_latitude).
    compound : STR or LIST, optional
        Interested compounds. The default is None.
    phase : STR, optional
        Selection of phase. The default is 'All'.
            'G' - Gas phase.
            'L' - Liquid phase.
            'L-FW' - Liquid freshwater.
            'L-SW' - Liquid seawater.
            'All' - All phases.
    altArray : LIST or np.ndarray, optional
        Requested altitudes. The default is None.
    numAlt : INT, optional
        Number of altitude steps from 0.0 m to maximum tropopause altitude, if list of altitude is not given by the user with `altArray`. The default is 50.
    surftrop : STR, optional
        Get concentration from 2-meters air following topography ('surface') or tropopause height ('tropopause').. The default is None.
    figsize : (FLOAT, FLOAT), optional
        Figure size in inches. (Width, Height). The default is None.
    lw : FLOAT, optional
        Line width. The default is 2.0.
    color : STR, optional
        Color of lines and fill between. The default is 'deepskyblue'.
    marker : STR, optional
        Set the line marker. The default is None.
    ms : FLOAT or INT, optional
        Set the marker size in points. The default is 5.
    logScale : BOOL, optional
        Set whether coordinate-y is plotted in logarithmic scale. The default is False.
    alpha_fillbtw : FLOAT, optional
        Set alpha blending value, between 0.0 (transparent) and 1.0 (opaque). The default is 0.3.
    date_format : STR, optional
        Set date format shown in coordinate x. The default is None.
    yticks_format : STR, optional
        Set format of the y-coordinate ticks. The default is None.
    xlabel : STR, optional
        Set x-axis title. The default is 'Time'.
    ylabel : STR, optional
        Set y-axis title. The default is 'Variable (Units)'.
    xlabel_rotation : FLOAT, optional
        Set x-label rotation. The default is 0.0.
    yLim : [FLOAT, FLOAT], optional
        Set limits of y-coordinate [bottom, top]. The default is [None, None].
    show_right_labels : BOOL, optional
        Set whether quartile lables are shown. The default is True.
    drl : FLOAT, optional
        Set distance of right labels (Max, Q3, Q2, Q1, Min). The default is 1.005.
    fsrl : FLOAT, optional
        Font size of right labels. The default is 10.0.
    fontFamily : STR, optional
        Set font family. The default is 'Arial'.
    fontSize : FLOAT, optional
        Set font size. The default is 12.0.
    title : STR, optional
        Set plot title. The default is None.
    cf : FLOAT, optional
        Set conversion factor for variable. The default is 1.0.
    lines : [BOOL, BOOL, BOOL, BOOL, BOOL], optional
        Set what quartiles are plotted: [Max, Q3, Q2, Q1, Min]. The default is [True, True, True, True, True].
    showMessage : BOOL, optional
        Boolean to set whether informative messages are displayed in Console. The default is False.
    fillMissing : BOOL, optional
        Set whether missing data of air composition is filled with air composition from ISA (True) or leave data as np.nan (False). The default is False.
    weights : BOOL, optional
        Set whether grid weight contributions (grid defined by longitude, latitude and altitude arrays) are considered in statistics. The default is False.
    savePlot : BOOL, optional
        Save resultant plot in `results/` folder. The default is False.

    Returns
    -------
    Spyder plot or Plot in 'results/' folder.

    """
    dict_attributes = {'Pi', 'MolPct_G', 'Ci_G', 'Ci_LFW', 'Ci_LSW', 'DGr'}
    check_dict_attribute = '-' in variable
    variable_name = variable
    if check_dict_attribute:
        comp = variable[variable.find('-')+1:]
        variable = variable[0:variable.find('-')]
    # Check valid models and data types
    validModels = {'ISAMERRA2', 'CAMSMERRA2'}
    if not model in validModels:
        raise ValueError(f'Invalid model ({model}) to plot seasonality of {variable}. Valid models: {validModels}.')
    validDataTypes = {'dly', 'mly', 'yly', 'cmly'}
    if not dataType in validDataTypes:
        raise ValueError(f'Invalid data type ({dataType}) to plot seasonality of {variable}. Valid data types: {validDataTypes}.')
    print(f'> Plotting seasonality of {variable_name}.')
    # Define dates
    if dataType == 'dly':
        start_date = datetime.datetime.strptime(start_date, '%Y-%m-%d')
        end_date = datetime.datetime.strptime(end_date, '%Y-%m-%d')
        delta_date = datetime.timedelta(days = delta_time)
    elif dataType == 'mly':
        start_date = datetime.datetime.strptime(start_date, '%Y-%m')
        end_date = datetime.datetime.strptime(end_date, '%Y-%m')
        delta_date = relativedelta(months = delta_time)
    elif dataType == 'cmly':
        start_date = datetime.datetime.strptime(start_date, '%Y-%m')
        end_date_ = datetime.datetime.strptime(end_date, '%Y-%m')
        y_ed = int(end_date_.year)
        y_sd = int(start_date.year)
        m_ed = int(end_date_.month)
        end_date = datetime.datetime.strptime(f'{y_sd}-{m_ed}', '%Y-%m')
        delta_date = relativedelta(months = delta_time)
    elif dataType == 'yly':
        start_date = datetime.datetime.strptime(start_date, '%Y')
        end_date = datetime.datetime.strptime(end_date, '%Y')
        delta_date = relativedelta(years = delta_time)
    date = start_date
    maximum = []
    Q3 = []
    Q2 = []
    Q1 = []
    minimum = []
    x_ = 0
    x = []
    labels = []
    while date <= end_date:
        x_ += 1
        x += [x_]
        if not isinstance(date_format, str):
            if dataType == 'dly':   date_format = '%d-%m-%Y'
            elif dataType == 'mly': date_format = '%m-%Y'
            elif dataType == 'yly': date_format = '%Y'
            elif dataType == 'cmly': date_format = '%b'
        labels += [date.strftime(date_format)]
        if dataType == 'cmly':
            y = [int(date.year), y_ed]
        else:
            y = int(date.year)
        m = int(date.month)
        d = int(date.day)
        if dataType == 'cmly':
            print(f'    > Date: {date.strftime(date_format)} {y_sd}_{y_ed}')
        else:
            print(f'    > Date: {date.strftime(date_format)}')
        if model == 'ISAMERRA2':
            model_inst = ISAMERRA2(dataType, y, m, d, 
                                   pH = pH, 
                                   bbox = bbox,
                                   compound = compound,
                                   phase = phase,
                                   altArray = altArray,
                                   numAlt = numAlt, 
                                   surftrop = surftrop,
                                   keysAsAttributes = True,
                                   showMessage  = showMessage)
        elif model == 'CAMSMERRA2':
            model_inst = CAMSMERRA2(dataType, y, m, d,
                                    pH = pH,
                                    bbox = bbox,
                                    keys = 'All',
                                    phase = phase,
                                    altArray = altArray,
                                    numAlt = numAlt,
                                    surftrop = surftrop,
                                    showMessage  = showMessage,
                                    fillMissing = fillMissing)
        if weights:
            lon = model_inst.lon
            lat = model_inst.lat
            alt = model_inst.altitude
            if alt.ndim != 1: alt = None
            NaN_values = model_inst.temperature
            grid_weights = wt(lon, lat, alt = alt, NaN_values = NaN_values)
        else:
            grid_weights = None
        model_attributes = list(model_inst.__dict__.keys())
        try:
            var = getattr(model_inst, variable)
        except:
            raise ValueError(f'Variable {variable} not found in {model} model. Available attributes: {model_attributes}.')
        if variable in dict_attributes:
            dict_keys = list(var.keys())
            try:
                var = var[comp]
            except:
                raise ValueError(f'Compound/variable {comp} not found in .{variable} of {model}. Available compounds/variables: {dict_keys}.')
        var *= cf
        if lines[0]: maximum += [np.nanmax(var)]
        if lines[1]: Q3 += [np.nanquantile(var, 0.75, method = 'inverted_cdf', weights = grid_weights)]
        if lines[2]: Q2 += [np.nanquantile(var, 0.5, method = 'inverted_cdf', weights = grid_weights)]
        if lines[3]: Q1 += [np.nanquantile(var, 0.25, method = 'inverted_cdf', weights = grid_weights)]
        if lines[4]: minimum += [np.nanmin(var)]
        date += delta_date
    #-Plotting
    plt.rcParams["font.family"] = fontFamily
    plt.rcParams["font.size"] = fontSize
    plt.rcParams['figure.dpi'] = 400
    plt.rcParams['savefig.dpi'] = 400
    fig, ax = plt.subplots(figsize = figsize)
    if logScale:
        if lines[0]: plt.semilogy(x, maximum, ls = '-', lw = lw*0.5, color = color, marker = marker, ms = ms)
        if lines[1]: plt.semilogy(x, Q3, ls = '--', lw = lw*0.5, color = color, marker = marker, ms = ms)
        if lines[2]: plt.semilogy(x, Q2, ls = '-', lw = lw, color = color, marker = marker, ms = ms)
        if lines[3]: plt.semilogy(x, Q1, ls = '--', lw = lw*0.5, color = color, marker = marker, ms = ms)
        if lines[4]: plt.semilogy(x, minimum, ls = '-', lw = lw*0.5, color = color, marker = marker, ms = ms)
    else:
        if lines[0]: plt.plot(x, maximum, ls = '-', lw = lw*0.5, color = color, marker = marker, ms = ms)
        if lines[1]: plt.plot(x, Q3, ls = '--', lw = lw*0.5, color = color, marker = marker, ms = ms)
        if lines[2]: plt.plot(x, Q2, ls = '-', lw = lw, color = color, marker = marker, ms = ms)
        if lines[3]: plt.plot(x, Q1, ls = '--', lw = lw*0.5, color = color, marker = marker, ms = ms)
        if lines[4]: plt.plot(x, minimum, ls = '-', lw = lw*0.5, color = color, marker = marker, ms = ms)
    plt.fill_between(x, minimum, maximum, alpha = alpha_fillbtw, color = color)
    if show_right_labels:
        if lines[0]: plt.text(x[-1]*drl, maximum[-1], 'Max', va = 'center', size = fsrl)
        if lines[1]: plt.text(x[-1]*drl, Q3[-1], 'Q3', va = 'center', size = fsrl)
        if lines[2]: plt.text(x[-1]*drl, Q2[-1], 'Q2', va = 'center', size = fsrl)
        if lines[3]: plt.text(x[-1]*drl, Q1[-1], 'Q1', va = 'center', size = fsrl)
        if lines[4]: plt.text(x[-1]*drl, minimum[-1], 'Min', va = 'center', size = fsrl)
    plt.margins(x=0)
    plt.title(title)
    ax.set_xticks(x, labels, rotation = xlabel_rotation)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if yticks_format:
        ax.yaxis.set_major_formatter(yticks_format)
    plt.ylim(yLim[0], yLim[1])
    if savePlot:
        cPlots = 1
        file = f'results/{variable_name}_plotSeasonality'
        _savePlot(file, cPlots)
    plt.show()

def plotVarMap2D(data, varName, varUnits, cmap, bbox, vmin = None, vmax = None, numlevels = 100, fontFamily = 'Arial',  
                 figsize = (5.0, 3.6), formatColorbar = '{:0.1f}', fix_aspect = False, fontsize = 12, fwtl = 'normal',
                 drawCoastLines = True, clw = 0.5, drawParallels = True, plw = 0.5, cpl = 'k', title = None,
                 drawMeridians = True, mlw = 0.5, cml = 'k', colorbar = True, parallelsLabels = [1,0,0,0], 
                 parallels = [-90, -60, -30, 0, 30, 60, 90], meridians = [0., 60., 120., 180., 240., 300.], 
                 meridiansLabels = [0,0,0,1], continentColor = 'darkgrey', lakeColor = 'darkgrey', projection = 'cyl',
                 logColorbar = False, cb_minor_ticks = False, cb_ticks = None, num_cb_ticks = 8, cb_labels_rotation = 0.0,
                 colorbarSize = (10, 4), cbOrientation = 'horizontal', cbFontSize = 12, savePlot = False):
                 
    """
    Plot variable on world map.

    Parameters
    ----------
    data : np.ndarray
        Data to be plotted.
    varName : STR
        Name of variable.
    varUnits : STR
        Units of variable.
    cmap : STR or LIST
        Set color-mapping.
    bbox : TUPLE
        Earth's region of data, the bounding box.
        (lower_left_longitude, lower_left_latitude, upper_right_longitude, upper_right_latitude)
    title : STR
        Set title of plot. The default is None.
    vmin : FLOAT, optional
        Set minimum value that the plot covers. The default is None.
    vmax : FLOAT, optional
        Set maximum value that the plot covers. The default is None.
    numlevels : INT, optional
        Set the number and positions of the contour lines / regions. The default is 100.
    fontfamily : STR, optional
        Set which font family is picked up, either by specifying family names of fonts installed on user's system, or generic-families (e.g., 'serif', 'sans-serif', 'monospace', 'fantasy' or 'cursive'), or a combination of both.
    figsize : (FLOAT, FLOAT), optional
        Figure size. (Width, Height) in inches. The default is (5.0, 3.6).
    fwtl : STR, optional
        Font weight of title. The default is 'normal'.
    formatColorbar : STR, optional
        Set format colorbar. The default is '{:0.1f}'.
    fix_aspect : BOOL, optional
        Fix aspect ratio of plot to match aspect ratio of map projection region. The default is False.
    fontsize : FLOAT, optional
        Set font size. The default is 12.
    drawCoastLines : BOOL, optional
        Set whether drawing coast lines. The default is True.
    clw : FLOAT, optional
        Coast line width. The default is 0.5.
    drawParallels : BOOL, optional
        Set whether drawing parallels lines. The default is True.
    plw : FLOAT, optional
        Parallel line width. The default is 0.5.
    cpl : STR, optional
        Parellel line color. The default is 'k'.
    drawMeridians : BOOL, optional
        DESCRIPTION. The default is True.
    mlw : FLOAT, optional
        Meridian line width. The default is 0.5.
    cml : STR, optional
        Meridian line color. The default is 'k'.
    colorbar : BOOL, optional
        Set whether colorbar is displayed. The default is True.
    parallelsLabels : LIST, optional
        Set whether parallels are labelled where they intersect as a list [left, right, top, bottom]. The default is [1,0,0,0].
    parallels : LIST, optional
        Set where parallels are drawn. The default is [-90, -60, -30, 0, 30, 60, 90].
    meridians : LIST, optional
        Set where meridians are drawn. The default is [0., 60., 120., 180., 240., 300.].
    meridiansLabels : LIST, optional
        Set whether meridians are labelled where they intersect as a list [left, right, top, bottom]. The default is [0,0,0,1].
    continentColor : STR, optional
        Color of continents. The default is 'darkgrey'.
    lakeColor : STR, optional
        Color of water bodies (lakes and seas). The default is 'darkgrey'.
    projection : STR, optional
        Map projection. The default is 'cly' (Cylindrical Equidistant Projection).
    logColorbar : BOOL, optional
        Set logarithmic scale on contour colormap. The default is False.
    cb_minor_ticks : BOOL, optional
        Set whether minor ticks are shown or not. The default is False.
    cb_ticks : LIST, optional
        Set colorbar ticks. The default is None.
    num_cb_ticks : INT, optional
        Set number of ticks of colorbar (if these are not defined by 'cb_ticks'). The default is 8.
    cb_labels_rotation : FLOAT, optional
        Rotation of colorbar labels. The default is 0.0.
    colorbarSize : (FLOAT, FLOAT), optional
        Colorbar size. (Width, Height) in inches. The default is (10, 4).
    cbOrientation : STR, optional
        Colorbar orientation. The default is 'horizontal'.
    cbFontSize : FLOAT, optional
        Set font size of colorbar. The default is 12.
    savePlot : BOOL, optional
        Set whether the plot is saved in `/results` folder. The default is False.

    Returns
    -------
    Spyder plot.

    """
    savePath = 'results/'
    plt.rcParams["font.family"] = fontFamily
    # Plotting colour
    if isinstance(cmap, (list, np.ndarray)):
        cmap = ListedColormap(cmap, name = 'plot_cmap')
    if not vmin:
        vmin = np.nanmin(data)
    if not vmax:
        vmax = np.nanmax(data)
    # Plotting properties
    ny = data.shape[0]
    nx = data.shape[1]
    fig, ax = plt.subplots(figsize = figsize)
    m = Basemap(projection=projection, resolution='c',
                llcrnrlon=bbox[0], llcrnrlat=bbox[1],
                urcrnrlon=bbox[2], urcrnrlat=bbox[3], fix_aspect = fix_aspect)
    lons, lats = m.makegrid(nx, ny)
    x, y = m(lons, lats)
    if drawCoastLines:
        m.drawcoastlines(linewidth = clw)
    # Draw parallels and meridians.
    if drawParallels:
        m.drawparallels(parallels, color = cpl, labels=parallelsLabels, fontsize=fontsize, linewidth=plw)
    if drawMeridians:
        m.drawmeridians(meridians, color = cml, labels=meridiansLabels, fontsize=fontsize, linewidth=mlw)
    m.fillcontinents(color = continentColor, lake_color = lakeColor)
    levels = np.linspace(vmin, vmax, int(numlevels))
    if logColorbar:
        mcont = m.contourf(x, y, data, 
                           levels=levels,
                           norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax),
                           cmap=cmap)
    else:
        mcont = m.contourf(x, y, data, 
                           levels=levels,
                           norm=plt.Normalize(vmin=vmin, vmax=vmax),
                           cmap=cmap)
    levels_array = np.array(mcont.levels).reshape(1, -1)
    ax.set_facecolor(lakeColor)
    if title:
        plt.title(title, fontweight = fwtl, fontsize = fontsize+2)
    if savePlot:
        cPlots = 1
        file = f'{savePath}{varName}_varMap2D'
        _savePlot(file, cPlots)
    plt.show()
    # Standalone colorbar
    if colorbar:
        if logColorbar:
            norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
            if not isinstance(cb_ticks, (list, np.ndarray)):
                ticks_colorbar = np.logspace(np.log10(np.squeeze(levels_array)[0]), 
                                             np.log10(np.squeeze(levels_array)[-1]), num_cb_ticks)
            else:
                ticks_colorbar = cb_ticks
        else:
            norm = plt.Normalize(vmin=vmin, vmax=vmax)
            if not isinstance(cb_ticks, (list, np.ndarray)):
                ticks_colorbar = np.linspace(np.squeeze(levels_array)[0], np.squeeze(levels_array)[-1], num_cb_ticks)
            else:
                ticks_colorbar = cb_ticks
        pl.figure(figsize=colorbarSize)
        pl.imshow(levels_array, cmap = cmap, norm = norm)
        pl.gca().set_visible(False)
        clb = pl.colorbar(orientation = cbOrientation)
        clb.set_ticks(ticks_colorbar)
        labels_colorbar = [formatColorbar.format(x) for x in ticks_colorbar]
        clb.set_ticklabels(labels_colorbar, rotation = cb_labels_rotation)
        if logColorbar:
            clb.set_label(f'{varName} ({varUnits})\n[log scale]', fontsize = cbFontSize)
        else:
            clb.set_label(f'{varName} ({varUnits})', fontsize = cbFontSize)
        clb.ax.tick_params(labelsize = cbFontSize)
        if cb_minor_ticks:
            clb.minorticks_on()
        else:
            clb.minorticks_off()
        if savePlot:
            cPlots = 1
            file = f'{savePath}{varName}_varMap2D_cbar'
            _savePlot(file, cPlots)
        plt.show()
    
def plotZonalMean(altitude, data, color, varName, varUnits, zone, pH = None, T = None, compSpec = None, fontFamily = 'Arial',
                  fwtl = 'normal', fillBetween = True, semiLog = False, title = None, figsize = (5.0, 3.6), lw = 2.0, 
                  fontsize = 12, alpha = 0.4, nticks = 10, legend = True, ncol = 1, legendOrientation = None, latitude = None, 
                  longitude = None, vmin = None, vmax = None, colorbar = None, cbFontSize = 12, xTicks = None, 
                  colorbarSize = None, savePlot = False, cbOrientation = 'horizontal', formatColorbar = '{:0.1f}'):
    """
    Plot zonal mean of data.

    Parameters
    ----------
    altitude :  LIST or np.ndarray
        List of altitude values.
    data : DICT
        Data to be plotted.
    color : STR, LIST or np.ndarray
        Color map of plot.
    varName : STR
        Name of variable.
    varUnits : STR
        Units of variable.
    zone : STR
        Selection of zone where mean is computed.
            'lat': Zonal mean on latitude.
            'lon': Zonal mean on longitude.
            'latlon': Zonal mean on latitude and longitude.
    pH : FLOAT, optional
        Set pH to calculate pH speciation (if necessary). The default is None.
    T : FLOAT, LIST or np.ndarray, optional
        Set temperature to calculate pH speciation (if necessary). The default is None.
    compSpec : LIST, optional
        Selection of chemical species for pH speciation. The default is None.
    fillBetween : BOOL, optional
        Plot maximum and minimum range as shaded region. The default is True.
    semiLog : BOOL, optional
        Set whether plotting variable in logarithmic scale (only for zone = 'latlon'). The default is False.
    title : STR, optional
        Set plot title. The default is None.
    figsize : (FLOAT, FLOAT), optional
        Set figure size (width, height). The default is (5.0, 3.6).
    lw : FLOAT, optional
        Set line width of plot (only for zone = 'latlon'). The default is 2.0.
    fontsize : FLOAT, optional
        Set plot font size. The default is 12.
    alpha : FLOAT, optional
        Set alpha blending value, between 0 (transparent) and 1 (opaque). The default is 0.4.
    nticks : INT, optional
        Set number of ticks for x-coordinate (only for zone = 'latlon'). The default is 10.
    legend : bool, optional
        Set whether plot legend is displayed (only for zone = 'latlon'). The default is True.
    ncol : int, optional
        Set number of legend columns (only for zone = 'latlon'). The default is 1.
    legendOrientation : STR, optional
        Set legend orientation: 'horizontal' or 'vertical' (only for zone = 'latlon'). The default is None.
    latitude : LIST or np.ndarray, optional
        Set list of latitude values (only for zone = 'lon'). The default is None.
    longitude : LIST or np.ndarray, optional
        Set list of longitude values (only for zone = 'lat'). The default is None.
    vmin : FLOAT, optional
        Set minimum value that the plot covers. The default is None.
    vmax : FLOAT, optional
        Set maximum value that the plot covers. The default is None.
    colorbar : BOOL, optional
        Set whether the colorbar is displayed (only for zone = 'lat' and 'lon'). The default is None.
    cbFontSize : FLOAT, optional
        Set font size of colorbar. The default is 12.
    xTicks : LIST, optional
        Set latitude and longitude tick values (only for zone = 'lat' and 'lon'). The default is None.
    colorbarSize : (FLOAT, FLOAT), optional
        Set size of colorbar (width, height). Only for zone = 'lat' and 'lon'. The default is None.
    cbOrientation : STR, optional
        Set colorbar orientation: 'horizontal' or 'vertical' (only for zone = 'lat' and 'lon'). The default is 'horizontal'.
    formatColorbar : STR, optional
        Set format of tick values of colorbar. The default is '{:0.1f}'.
    savePlot : BOOL, optional
        Set whether the plot is saved in `/results` folder. The default is False.

    Returns
    -------
    Spyder plot.

    """
    savePath = 'results/' 
    font = {'size': fontsize}
    plt.rc('font', **font)
    plt.rcParams["font.family"] = fontFamily
    if zone == 'latlon':
        axis_ = (1, 2)
    elif zone == 'lon':
        axis_ = 2
        if not isinstance(latitude, (list, np.ndarray)): raise ValueError('Argument \'latitude\' must be given as a list or np.ndarray.')
        if xTicks:
            positions = xTicks
            labels = []
            for position in positions:
                if position > 0:
                    labels += [f'{abs(int(position))}°N']
                elif position == 0:
                    labels += [f'{abs(int(position))}°']
                else:
                    labels += [f'{abs(int(position))}°S']
        else:
            positions = [-90, -60, -30, 0, 30, 60, 90]
            labels = ['90°S', '60°S', '30°S', '0°', '30°N', '60°N', '90°N']
        split_points = [0, 6.0, np.max(altitude)]
        backgroundColor = ['black', 'lightgrey']
    elif zone == 'lat':
        axis_ = 1
        if not isinstance(longitude, (list, np.ndarray)): raise ValueError('Argument \'longitude\' must be given as a list or np.ndarray.')
        if xTicks:
            positions = xTicks
            labels = []
            for position in positions:
                if position > 0:
                    labels += [f'{abs(int(position))}°W']
                elif position == 0:
                    labels += [f'{abs(int(position))}°']
                else:
                    labels += [f'{abs(int(position))}°E']
        else:
            positions = [-180, -90, 0, 90, 179.375]
            labels = ['180°', '90°W', '0°', '90°E', '180°']
        split_points = [0, 6.0, np.max(altitude)]
        backgroundColor = ['black', 'lightgrey']
    else:
        raise ValueError(f'Invalid zone ({zone}). Available zones: \'latlon\' (latitude and longitude), \'lon\' (longitude), \'lat\' (latitude)')
    vars_ = tuple(data.keys())
    dataPlot = {}
    legend_elements = []
    fig, ax = plt.subplots(figsize = figsize)
    for idVar, var in enumerate(data):
        color_ = color[var]
        # Plotting colour
        if isinstance(color_, (list, np.ndarray)):
            levels = len(color_)
            color_ = ListedColormap(color_, name = 'plot_cmap')
        if zone == 'latlon':
            legend_elements += [Line2D([0], [0], color = color_, lw = lw, label = var)]
        # pH speciation
        if pH:
            rComp, _, _ = Rxn.getRxnpH(var)
            if rComp:
                comp = np.intersect1d(rComp, compSpec)
                if not comp: comp = [var]
                if len(comp) > 1: raise ValueError(f'Multiple chemical species have been selected for {var}: {comp}.')
                dataPlot = ThEq.pHSpeciation(comp[0], pH, T, data[var])
            else:
                dataPlot = data[var]
        else:
            dataPlot = data[var]
        # Data stats
        data_min = np.nanquantile(dataPlot, 0.0, axis = axis_)
        data_av = np.nanquantile(dataPlot, 0.5, axis = axis_)
        data_max = np.nanquantile(dataPlot, 1.0, axis = axis_)
        #-Plotting
        if zone == 'latlon':
            if semiLog:
                plt.semilogx(data_av, altitude, color_, lw = lw)
                maj_loc = matplotlib.ticker.LogLocator(subs='all', numticks=nticks)
                min_loc = matplotlib.ticker.LogLocator(subs='all', numticks=nticks)
                ax.xaxis.set_major_locator(maj_loc)
                ax.xaxis.set_minor_locator(min_loc)
            else:
                plt.plot(data_av, altitude, color_, lw = lw)
            if fillBetween:
                plt.fill_betweenx(altitude, data_min, data_max, color = color_, alpha = alpha)
            if vmin and vmax:
                plt.xlim(left = vmin, right = vmax)
            ax.set_ylabel('Altitude (km)')
            ax.set_xlabel(f'{varName} ({varUnits})')
            plt.minorticks_on()
            fig.tight_layout()
            plt.margins(y=0)
            if title:
                plt.title(title, fontweight = fwtl, fontsize = fontsize+2)
        elif zone == 'lon':
            ax.contourf(latitude, altitude, data_av, levels = levels, cmap = color_, vmin = vmin, vmax = vmax)
            ax.set_xlabel('Latitude (°)')
            ax.set_ylabel('Altitude (km)')
            ax.xaxis.set_major_locator(tkr.FixedLocator(positions))
            ax.xaxis.set_major_formatter(tkr.FixedFormatter(labels))
            plt.minorticks_on()
            for i in range(len(split_points) - 1):
                plt.axhspan(split_points[i], split_points[i+1], facecolor = backgroundColor[i], zorder = 0)
        elif zone == 'lat':
            ax.contourf(longitude, altitude, data_av, levels = levels, cmap = color_, vmin = vmin, vmax = vmax)
            ax.set_xlabel('Longitude (°)')
            ax.set_ylabel('Altitude (km)')
            ax.xaxis.set_major_locator(tkr.FixedLocator(positions))
            ax.xaxis.set_major_formatter(tkr.FixedFormatter(labels))
            plt.minorticks_on()
            for i in range(len(split_points) - 1):
                plt.axhspan(split_points[i], split_points[i+1], facecolor = backgroundColor[i], zorder = 0)
        if savePlot:
            cPlots = 1
            file = f'{savePath}{varName}_{zone}Mean'
            _savePlot(file, cPlots)
    plt.show()
    if zone == 'latlon':
        # Standalone legend
        if legend:
            fig, ax = plt.subplots()
            if legendOrientation == 'vertical':
                ax.legend(handles = legend_elements, loc='center')
            elif legendOrientation == 'horizontal':
                ax.legend(handles = legend_elements, loc='center', ncol = len(vars_))
            if not legendOrientation:
                ax.legend(handles = legend_elements, loc='center', ncol = ncol)
            ax.set_axis_off()
            fig.tight_layout()
            if savePlot:
                cPlots = 1
                file = f'{savePath}{varName}_{zone}Mean_legend'
                _savePlot(file, cPlots)
            plt.show()
    elif zone == 'lon' or zone == 'lat':
        # Standalone colorbar
        if colorbar:
            pl.figure(figsize=colorbarSize)
            limitVar = np.array([[vmin, vmax]])
            pl.imshow(limitVar, cmap = color_)
            pl.gca().set_visible(False)
            clb = pl.colorbar(orientation = cbOrientation)
            ticks_colorbar = np.linspace(np.squeeze(limitVar)[0], np.squeeze(limitVar)[1], 8)
            clb.set_ticks(ticks_colorbar)
            labels_colorbar = [formatColorbar.format(x) for x in ticks_colorbar]
            clb.set_ticklabels(labels_colorbar)
            clb.set_label(f'{varName} ({varUnits})', fontsize = cbFontSize)
            clb.ax.tick_params(labelsize = cbFontSize)
            if savePlot:
                cPlots = 1
                file = f'{savePath}{varName}_{zone}Mean_cbar'
                _savePlot(file, cPlots)
            plt.show()
    
def plotCrossSections(data2D, data3D, varName, varUnits, cmap, altitude, bbox = (-180, -90, 180, 90), fontFamily = 'Arial',
                      fwtl = 'normal', sections = None, depthArray = [0], fontsize = 8, vmin = None, vmax = None, title = None,
                      colorbar = True, xylabels = True, levels = 100, sectionFigSize = None, mapsize = (5.8, 4.5), 
                      fix_aspect = False, numTicks = None, clw = 0.5, continentColor = 'darkgrey', 
                      lakeColor = 'darkgrey', savePlot = False):
    """
    Plot three dimensional data on a world map (2D data) and different section plots (meridians and parallels; 3D)

    Parameters
    ----------
    data2D : np.ndarray
        Two-dimensional data to be plotted. (latitude and longitude).
    data3D : np.ndarray
        Three-dimensional data to be plotted (altitude, latitude and longitude).
    varName : STR
        Name of variable.
    varUnits : STR
        Units of variable.
    cmap : STR or LIST
        Set color-mapping.
    altitude : LIST or np.ndarray
        List of altitude values.
    bbox : TUPLE, optional
        Earths region of data, the bounding box. The default is (-180, -90, 180, 90).
        (lower_left_longitude, lower_left_latitude, upper_right_longitude, upper_right_latitude). 
    sections : DICT, optional
        Set of sections (meridians and parallels; {'A-B': (-180, 32, 180, 32), 'C-D': (90.0, -90, 90.0, 90)}). The default is None.
    depthArray : LIST or np.ndarray, optional
        List of depth values. The default is [0].
    fontsize : FLOAT, optional
        Set font size. The default is 8.
    vmin : FLOAT, optional
        Set minimum value that the plot covers. The default is None.
    vmax : FLOAT, optional
        Set maximum value that the plot covers. The default is None.
    title : STR
        Set title of world map plot. The default is None.
    colorbar : BOOL, optional
        Set whether colorbar is displayed.. The default is True.
    xylabels : BOOL, optional
        Set whether section labels are displayed in world map (2D data plotting). The default is True.
    levels : INT, optional
        Set the number and positions of the contour lines / regions. The default is 100.
    sectionFigSize : DICT, optional
        Section map sizes (3D data plotting; {'A-B': (7.2, 2.5)}). The default is None.
    mapsize : (FLOAT, FLOAT), optional
        World map size (2D data plotting). The default is (5.8, 4.5).
    fix_aspect : BOOL, optional
        Fix aspect ratio of plot to match aspect ratio of map projection region. The default is False.
    numTicks : DICT, optional
        Set latitude and longitude tick values of section plots. {'A-B': 9}. The default is None.
    clw : FLOAT, optional
        Coast line width. The default is 0.5.
    continentColor : STR, optional
        Color of continents. The default is 'darkgrey'.
    lakeColor : STR, optional
        Color of water bodies (lakes and seas). The default is 'darkgrey'.
    savePlot : BOOL, optional
        Set whether the plot is saved in `/results` folder. The default is False.

    Returns
    -------
    Spyder plot.

    """
    savePath = 'results/'
    font = {'size': fontsize}
    plt.rc('font', **font)
    plt.rcParams["font.family"] = fontFamily
    # Ocean depth
    npzDepth = np.load('data/GEBCO_LR.npz')
    keys = npzDepth.keys()
    dictDepth = {key: npzDepth[key] for key in keys}
    globalDepth = dictDepth['elevation']
    x = np.shape(globalDepth)[1]
    y = np.shape(globalDepth)[0]
    depth3D = np.repeat(globalDepth[np.newaxis, :, :], len(depthArray), axis = 0)
    # Aerosol concentration
    npzPHIS = np.load('data/MERRA2/PHIS.npz')
    lonR = npzPHIS['lon']
    latR = npzPHIS['lat']
    H = npzPHIS['PHIS']
    H3D = np.repeat(H[np.newaxis, :, :], len(altitude), axis = 0)
    # Closest coordinates
    for section in sections:
        section_coor = sections[section]
        ATM = Atmosphere()
        new_coor = ATM._closestCoord(lonR, latR, section_coor)
        sections[section] = new_coor
    labels = {}
    locusType = {}
    if sections:
        for section in sections:
            section_coor = sections[section]
            if section_coor[0] == -180 and section_coor[2] == 179.375 and section_coor[1] == section_coor[3]:
                # Parallel
                locusType[section] = 'parallel'
                if section_coor[1] > 0:
                    labels[section] = f'{"{:.1f}".format(section_coor[1])}°N'
                elif section_coor[1] == 0:
                    labels[section] = f'{"{:.1f}".format(section_coor[1])}°'
                else:
                    labels[section] = f'{"{:.1f}".format(section_coor[1])}°S'
            elif section_coor[1] == -90 and section_coor[2] == 90 and section_coor[0] == section_coor[2]:
                # Meridian
                locusType[section] = 'meridian'
                if section_coor[0] > 0:
                    labels[section] = f'{"{:.1f}".format(section_coor[0])}°E'
                elif section_coor[0] == 0:
                    labels[section] = f'{"{:.1f}".format(section_coor[0])}°'
                else:
                    labels[section] = f'{"{:.1f}".format(section_coor[0])}°W'
            else:
                locusType[section] = None
                labels[section] = None
    # Get vmin and vmax if not given
    if not vmin:
        vmin = min(np.nanmin(data2D), np.nanmin(data3D))
    if not vmax:
        vmax = max(np.nanmax(data2D), np.nanmax(data3D))
    backgroundColor = ['black', 'lightgrey']
    split_points = [np.min(depthArray)/1000, np.nanmax(H)/1000, np.max(altitude)/1000]
    # Plotting colour
    if isinstance(cmap, (list, np.ndarray)):
        plot_cmap = ListedColormap(cmap, name = 'plot_cmap')
        levels = len(cmap)
    else:
        plot_cmap = cmap
        levels = levels
    #-Section plots
    if sections:
        for section in sections:
            if numTicks:
                numTick = numTicks[section]
            else:
                if locusType[section] == 'parallel':
                    numTick = 9
                else:
                    numTick = 5
            coor = sections[section]
            locus = locusType[section]
            idx = (int(np.argwhere(lonR == coor[0])),
                   int(np.argwhere(latR == coor[1])),
                   int(np.argwhere(lonR == coor[2])),
                   int(np.argwhere(latR == coor[3])))
            if locus:
                if locus == 'parallel':
                    sectionData = data3D[:, idx[1], :]
                    matrixDepthArray = np.repeat(depthArray[:, np.newaxis], sectionData.shape[-1], axis = -1)
                    matrixAltArray = np.repeat(altitude[:, np.newaxis], sectionData.shape[-1], axis = -1)
                    depth = depth3D[:, idx[1], :]
                    depth = np.where(depth < matrixDepthArray, 9.99e99, -9.99e99)
                    topog = H3D[:, idx[1], :]
                    topog = np.where(topog >= matrixAltArray, 9.99e99, np.nan)
                    nanDepth = np.empty(depth.shape)
                    nanDepth.fill(np.nan)
                    xCoor = np.hstack((lonR[:], 180))
                    xLabel = 'Longitude (°)'
                    tickLoc = np.linspace(bbox[0], bbox[2], numTick)
                    tickLabels = []
                    for loc in tickLoc:
                        loc = int(loc)
                        if loc == 0 or loc == 180 or loc == -180:
                            tickLabels += [f'{abs(loc)}°']
                        elif loc > 0:
                            tickLabels += [f'{loc}°E']
                        else:
                            tickLabels += [f'{abs(loc)}°W']
                elif locus == 'meridian':
                    sectionData = data3D[:, :, idx[0]]
                    matrixDepthArray = np.repeat(depthArray[:, np.newaxis], sectionData.shape[-1], axis = -1)
                    matrixAltArray = np.repeat(altitude[:, np.newaxis], sectionData.shape[-1], axis = -1)
                    depth = depth3D[:, :, idx[0]]
                    depth = np.where(depth < matrixDepthArray, 9.99e99, -9.99e99)
                    topog = H3D[:, :, idx[0]]
                    topog = np.where(topog >= matrixAltArray, 9.99e99, np.nan)
                    nanDepth = np.empty(depth.shape)
                    nanDepth.fill(np.nan)
                    xCoor = latR[:]
                    xLabel = 'Latitude (°)'
                    tickLoc = np.linspace(bbox[1], bbox[3], numTick)
                    tickLabels = []
                    for loc in tickLoc:
                        loc = int(loc)
                        if loc > 0:
                            tickLabels += [f'{loc}°N']
                        elif loc == 0:
                            tickLabels += [f'{loc}°']
                        else:
                            tickLabels += [f'{abs(loc)}°S']
            else:
                sectionData = data3D[:, idx[1], idx[0]:idx[2]+1]
                matrixDepthArray = np.repeat(depthArray[:, np.newaxis], sectionData.shape[-1], axis = -1)
                matrixAltArray = np.repeat(altitude[:, np.newaxis], sectionData.shape[-1], axis = -1)
                depth = depth3D[:, idx[1], idx[0]:idx[2]+1]
                depth = np.where(depth < matrixDepthArray, 9.99e99, -9.99e99)
                topog = H3D[:, idx[1], idx[0]:idx[2]+1]
                topog = np.where(topog >= matrixAltArray, 9.99e99, np.nan)
                nanDepth = np.empty(depth.shape)
                nanDepth.fill(np.nan)
                xCoor = lonR[idx[0]:idx[2]+1]
                xLabel = 'Longitude (°)'
                if coor[1] == coor[3]:
                    # Parallel section
                    tickLoc = np.linspace(coor[0], coor[2], numTick)
                    tickLabels = []
                    for loc in tickLoc:
                        if loc > 0:
                            tickLabels += [f'{int(loc)}°E']
                        elif loc == 0 or loc == 180 or loc == -180:
                            tickLabels += [f'{int(loc)}°']
                        else:
                            tickLabels += [f'{int(abs(loc))}°W']
                elif coor[0] == coor[2]:
                    # Meridian section
                    tickLoc = np.linspace(coor[1], coor[3], numTick)
                    for loc in tickLoc:
                        if loc > 0:
                            tickLabels += [f'{loc}°N']
                        elif loc == 0:
                            tickLabels += [f'{loc}°']
                        else:
                            tickLabels += [f'{abs(loc)}°S']
            sectionData = np.vstack((sectionData, nanDepth))
            if locus == 'parallel':
                sectionData = np.column_stack((sectionData, sectionData[:,0]))
                depth = np.column_stack((depth, depth[:,0]))
            yCoor = np.hstack((altitude, depthArray))
            fig, ax = plt.subplots(figsize=sectionFigSize[section])
            ax.contourf(xCoor, yCoor/1000, sectionData, levels = 78, cmap = plot_cmap, vmin = vmin, vmax = vmax)
            ax.contourf(xCoor, depthArray/1000, depth, levels = 2, colors = ('k', 'cornflowerblue'))
            ax.set_ylabel('Altitude (km)', fontsize = 10)
            ax.set_xlabel(xLabel, fontsize = 10)
            ax.tick_params(labelsize = 10)
            ax.xaxis.set_major_locator(tkr.FixedLocator(tickLoc))
            ax.xaxis.set_major_formatter(tkr.FixedFormatter(tickLabels))
            ax.set_xlim(left = tickLoc[0], right = tickLoc[-1])
            plt.minorticks_on()
            for i in range(len(split_points) - 1):
                plt.axhspan(split_points[i], split_points[i+1], facecolor = backgroundColor[i], zorder = 0)
            if savePlot:
                cPlots = 1
                file = f'{savePath}{varName}_{section}Section'
                _savePlot(file, cPlots)
            plt.show()
    #-Main plot (global map)
    ny, nx = data2D.shape
    fig, ax = plt.subplots(figsize = mapsize)
    m = Basemap(projection='cyl', resolution='c',
                llcrnrlon=bbox[0], llcrnrlat=bbox[1],
                urcrnrlon=bbox[2], urcrnrlat=bbox[3], fix_aspect = fix_aspect)
    lons, lats = m.makegrid(nx, ny)
    x, y = m(lons, lats)
    m.drawcoastlines(linewidth=clw)
    m.fillcontinents(color = continentColor, lake_color = lakeColor)
    if xylabels:
        ax.set_xlabel('Longitude (°)', size = 10)
        ax.xaxis.set_label_coords(0.10, 1.05)
        ax.set_ylabel('Latitude (°)', size = 10)
        ax.yaxis.set_label_coords(-0.010, 0.12)
    m.drawmapboundary(fill_color='darkgrey')
    m.contourf(x, y, data2D, levels = levels, cmap = plot_cmap, vmin = vmin, vmax = vmax)
    if title:
        plt.title(title, fontweight = fwtl, fontsize = fontsize+2)
    if sections:
        for section in sections:
            locus = locusType[section]
            coor = sections[section]
            x1 = coor[0]
            x2 = coor[2]
            y1 = coor[1]
            y2 = coor[3]
            letterLabel1 = section[0]
            letterLabel2 = section[-1]
            if locus == 'parallel':
                parallelLabel = labels[section]
                parallel = coor[1]
                m.drawparallels([parallel], labels=[0,0,0,0], fontsize=10, linewidth=1.5, dashes=[4, 2])
                plt.plot(x1, y1, x2, y2, marker = '|', color = 'k', ms=12.0, mew=2.0, zorder = 3.5)
                ax.text(-232, parallel - 2.0, parallelLabel, size = 10)
                ax.text(-195, parallel - 3.0, letterLabel1, size = 12, weight = 'bold')
                ax.text(185, parallel - 3.0, letterLabel2, size = 12, weight = 'bold')
                        # bbox = dict(boxstyle=f"circle,pad={0.25}", fc = 'w'))
            elif locus == 'meridian':
                meridianLabel = labels[section]
                meridian = coor[0]
                m.drawmeridians([meridian], labels=[0,0,0,0], fontsize=10, linewidth=1.5, dashes=[4, 2])
                plt.plot(x1, y1, x2, y2, marker = '_', color = 'k', ms=12.0, mew=2.0, zorder = 3.5)
                ax.text(meridian - 15, 102, meridianLabel, size = 10)
                ax.text(meridian - 5, 92, letterLabel1, size = 12, weight = 'bold')
                ax.text(meridian - 5, -100, letterLabel2, size = 12, weight = 'bold')
            else:
                plt.plot(x1, y1, x2, y2, marker = '|', color = 'k', ms=12.0, mew=2.0, zorder = 3.5, linestyle = '-')
                ax.text(x1 - 5, y1 + 6, letterLabel1, size = 12, weight = 'bold')
                ax.text(x2 - 5, y1 + 6, letterLabel2, size = 12, weight = 'bold')
    if savePlot:
        cPlots = 1
        file = f'{savePath}{varName}_worldMap2D'
        _savePlot(file, cPlots)
    plt.show()
    if colorbar:
        # Standalone colorbar
        limitData = np.array([[vmin, vmax]])
        pl.imshow(limitData, cmap = plot_cmap)
        pl.gca().set_visible(False)
        ticks_colorbar = np.linspace(np.squeeze(limitData)[0], np.squeeze(limitData)[1], 8)
        labels_colorbar = ["{:0.1f}".format(x) for x in ticks_colorbar]
        clb = pl.colorbar(orientation = 'horizontal')
        clb.set_label(f'{varName} ({varUnits})')
        clb.set_ticks(ticks_colorbar)
        clb.set_ticklabels(labels_colorbar)
        if savePlot:
            cPlots = 1
            file = f'{savePath}{varName}_worldMap2D_cbar'
            _savePlot(file, cPlots)