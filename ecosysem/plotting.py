# -*- coding: utf-8 -*-
"""
Created on Mon Nov 10 17:12:16 2025

@author: emartinez
"""

from reactions import Reactions as Rxn
from thermodynamics import ThEq
from environments import Atmosphere

from mpl_toolkits.basemap import Basemap
import pylab as pl
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap
import matplotlib.ticker as tkr

def plotVarMap2D(data, varName, varUnits, cmap, bbox, vmin = None, vmax = None, numlevels = 100, 
                 figsize = (5.0, 3.6), formatColorbar = '{:0.1f}', fix_aspect = False, fontsize = 12, 
                 drawCoastLines = True, clw = 0.5, drawParallels = True, plw = 0.5, cpl = 'k', title = None,
                 drawMeridians = True, mlw = 0.5, cml = 'k', colorbar = True, parallelsLabels = [1,0,0,0], 
                 parallels = [-90, -60, -30, 0, 30, 60, 90], meridians = [0., 60., 120., 180., 240., 300.], 
                 meridiansLabels = [0,0,0,1], continentColor = 'darkgrey', lakeColor = 'darkgrey', 
                 colorbarSize = (10, 4), cbOrientation = 'horizontal', cbFontSize = 12):
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
        Earths region of data, the bounding box.
        (lower_left_longitude, lower_left_latitude, upper_right_longitude, upper_right_latitude)
    title : STR
        Set title of plot. The default is None.
    vmin : FLOAT, optional
        Set minimum value that the plot covers. The default is None.
    vmax : FLOAT, optional
        Set maximum value that the plot covers. The default is None.
    numlevels : INT, optional
        Set the number and positions of the contour lines / regions. The default is 100.
    figsize : (FLOAT, FLOAT), optional
        Figure size. (Width, Height) in inches. The default is (5.0, 3.6).
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
        Set whether meridians are labelled where they intersect as a list [left, right, top, bottom].. The default is [0,0,0,1].
    continentColor : STR, optional
        Color of continents. The default is 'darkgrey'.
    lakeColor : STR, optional
        Color of water bodies (lakes and seas). The default is 'darkgrey'.
    colorbarSize : (FLOAT, FLOAT), optional
        Colorbar size. (Width, Height) in inches. The default is (10, 4).
    cbOrientation : STR, optional
        Colorbar orientation. The default is 'horizontal'.
    cbFontSize : FLOAT, optional
        Set font size of colorbar. The default is 12.

    Returns
    -------
    Spyder plot.

    """
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
    m = Basemap(projection='cyl', resolution='c',
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
    m.contourf(x, y, data, 
               levels=levels,
               norm=plt.Normalize(vmin=vmin, vmax=vmax),
               cmap=cmap)
    ax.set_facecolor(lakeColor)
    if title:
        plt.title(title, fontweight = 'bold')
    plt.show()
    # Standalone colorbar
    if colorbar:
        pl.figure(figsize=colorbarSize)
        limitVar = np.array([[vmin, vmax]])
        pl.imshow(limitVar, cmap = cmap)
        pl.gca().set_visible(False)
        clb = pl.colorbar(orientation = cbOrientation)
        ticks_colorbar = np.linspace(np.squeeze(limitVar)[0], np.squeeze(limitVar)[1], 8)
        clb.set_ticks(ticks_colorbar)
        labels_colorbar = [formatColorbar.format(x) for x in ticks_colorbar]
        clb.set_ticklabels(labels_colorbar)
        clb.set_label(f'{varName} ({varUnits})', fontsize = cbFontSize)
        clb.ax.tick_params(labelsize = cbFontSize)
        plt.show()

def plotZonalMean(altitude, data, color, varName, varUnits, zone, pH = None, T = None, compSpec = None, 
                  fillBetween = True, semiLog = False, title = None, figsize = (5.0, 3.6), lw = 2.0, 
                  fontsize = 12, alpha = 0.4, nticks = 10, legend = True, ncol = 1, legendOrientation = None, 
                  latitude = None, longitude = None, vmin = None, vmax = None, colorbar = None, cbFontSize = 12,
                  xTicks = None, colorbarSize = None, cbOrientation = 'horizontal', formatColorbar = '{:0.1f}'):
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
        Units of variables.
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
        Set maximum value that the plot covers.. The default is None.
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

    Returns
    -------
    Spyder plot.

    """
    font = {'size': fontsize}
    plt.rc('font', **font)
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
            labels = ['180°W', '90°W', '0°', '90°E', '180°E']
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
                plt.title(title, fontweight = 'bold')
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
            plt.show()
    
