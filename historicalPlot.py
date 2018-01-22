#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 16:07:38 2018

@author: taylorm
"""
import matplotlib as mpl
mpl.use('agg')
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as ft
import gc
from datetime import datetime, timedelta

def climoMEPlot(date, ensMean, ensStd, datefhour, dateArr):
    lons = np.arange(260,320.1,1)
    lats = np.arange(25,65.1,1)[::-1]
    datefhour = int(datefhour)
    colorsdiv = ((0.000, 0.655, 0.655, 1),
                 (0.314, 0.773, 0.773, 1),
                 (0.675, 0.875, 0.875, 1),
                 (0.890, 0.957, 0.957, 1),
                 (1.000, 1.000, 1.000, 0),
                 (0.984, 0.925, 0.961, 1),
                 (0.953, 0.792, 0.886, 1),
                 (0.910, 0.620, 0.796, 1),
                 (0.847, 0.396, 0.690, 1))

    colors = ((255/255, 255/255, 255/255, 0),
              (236./255., 192./255., 0/255., 1),
              (232./255., 133./255., 58./255., 1),
              (210./255., 78./255., 113./255., 1),
              (171./255., 20./255., 136./255., 1),
              (114./255., 0., 141./255., 1),
              (10./255., 45./255., 110./255., 1))
    states = ft.NaturalEarthFeature(category='cultural', scale='50m',
                                    facecolor='none',
                                    name='admin_1_states_provinces_lines')
    extent = [-100, -40, 25, 50]
    fig = plt.figure(figsize=(20, 20))
    ax3 = plt.axes(projection=ccrs.PlateCarree())
    ax3.set_extent(extent)
    ax3.coastlines(resolution=('10m'), zorder=4)
    ax3.add_feature(ft.BORDERS, alpha=0.7, zorder=3)
    ax3.add_feature(ft.OCEAN,facecolor='w', zorder=1)
    ax3.add_feature(ft.LAND,  facecolor='w', zorder=1)
    ax3.add_feature(ft.RIVERS, alpha=0.8)
    ax3.add_feature(ft.LAKES)
    x, y = np.meshgrid(lons, lats)
    sprdspace = [0, 1, 2, 4, 8, 12, 15, 20]
    mslspace = range(900, 1100, 2)
    ax3.add_feature(states, edgecolor='gray', zorder=3)

    cf = ax3.contourf(x, y, ensStd/100, sprdspace,
                      transform=ccrs.PlateCarree(), colors=colors, zorder=2)

    c1 = ax3.contour(x, y, ensMean/100, levels=mslspace, colors='k',
                     linewidths=0.5)
    cbar = plt.colorbar(cf, aspect=50, ticks=sprdspace,
                        orientation='horizontal',
                        pad=0.01, extend='max', extendrect=True)
    cbar.set_label('Spread (hPa)')
    plt.clabel(c1, fontsize=10, inline_spacing=-0.5, fmt='%3.0f')
    plt.title('GEFS Ensemble Mean MSLP, Spread ' +
              dateArr.strftime('%Y/%m/%d %Hz') + '\n' +
              str(datefhour) + 'h Forecast valid ' +
              date.strftime('%Y/%m/%d %Hz '), loc='left')
    plt.savefig('/home/taylorm/mcli/plots/'+date.strftime('%Y%m%d%H')
    +str(datefhour)+'me.png', bbox_inches='tight')
    plt.close(fig)
    
if __name__ == "__main__":
    mean = np.load('/home/taylorm/hourallgefs.npy')
    sprd = np.load('/home/taylorm/hourallgefssprd.npy')
    with open('/home/taylorm/allhoursdates.txt') as f:
        cont = f.readlines()
    cont = [x.strip() for x in cont]
    datearray = [datetime.fromordinal(int(q)) for q in cont]
    
    [climoMEPlot(datearray[n], mean[n,q-1], sprd[n,q-1], q*24, datearray[n]-timedelta(days=q)) for n in range(0,len(mean)) for q in np.arange(1,7.1,1)]