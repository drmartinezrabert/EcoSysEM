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

