#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 14:33:27 2017

@author: taylorm
"""

import matplotlib as mpl
mpl.use('agg')
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as ft
import gc
#Globals
proj = ccrs.LambertConformal(central_longitude = -98,standard_parallels=(49,77))
mslspace = range(900, 1100, 4)
pwatspace = [0,1,2,5,10,20,30,40,50,60]
# Sea level pressure plotting function


def slpplotMaker(date, ensMean, ensStd, datefhour, dateArr, ssaAnom,
                 subsetPerc, totalPerc, pmm, lats, lons):
    extents=([-180, -50, 20, 65],'na')
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
    # Default Variables
    states = ft.NaturalEarthFeature(category='cultural', scale='50m',
                                    facecolor='none',
                                    name='admin_1_states_provinces_lines')

    x, y = np.meshgrid(lons, lats)
    sprdspace = [0, 1, 2, 4, 8, 12, 15, 20]
    extent = extents[0]
    loc = extents[1]
    fig = plt.figure(figsize=(20, 20))
    ax3 = plt.axes(projection=proj)
    ax3.set_extent(extent,crs=ccrs.PlateCarree())
    ax3.coastlines(resolution=('50m'), zorder=4)
    ax3.add_feature(ft.BORDERS,alpha=0.7,zorder=3)
    ax3.add_feature(ft.RIVERS,alpha=0.4,edgecolor='gray', zorder=3)
    ax3.add_feature(states, edgecolor='gray', zorder=3)
    ax3.add_feature(ft.LAKES,alpha=0.5,facecolor='gray', zorder=3)
    cf = ax3.contourf(x, y, ensStd/100, sprdspace,
                      transform=ccrs.PlateCarree(), colors=colors, zorder=2)

    c1 = ax3.contour(x, y, ensMean/100, levels=mslspace, colors='k',
                     linewidths=0.5,transform=ccrs.PlateCarree())
    cbar = plt.colorbar(cf, aspect=50, ticks=sprdspace,
                        orientation='horizontal',
                        pad=0.01, extend='max', extendrect=True)
    cbar.set_label('Spread (hPa)')
    clab=plt.clabel(c1, fontsize=10, inline_spacing=-0.5, fmt='%3.0f')
    plt.title('GEFS Ensemble Mean MSLP, Spread ' +
              date.strftime('%Y/%m/%d %Hz') + '\n' +
              str(datefhour) + 'h Forecast valid ' +
              dateArr.strftime('%Y/%m/%d %Hz '), loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/mslp/me/'+loc+'/init.png',
                    bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/mslp/me/'+loc+'/'+date.strftime('%Y%m%d%H')
    +str(datefhour)+'me'+loc+'.png', bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()

    cf = ax3.contourf(x,y,ensStd/100,sprdspace,
                  transform=ccrs.PlateCarree(),colors=colors,zorder=2)
    cf.cmap.set_over('#0a2d84')
    c1=ax3.contour(x,y,pmm/100,levels=mslspace,colors='k',linewidths=0.5,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,ticks=sprdspace, orientation='horizontal',pad=0.01,extend='max',extendrect=True)
    cbar.set_label('Spread (hPa)')
    clab = plt.clabel(c1,fontsize=10,inline_spacing=-0.5,fmt='%3.0f')
    plt.title('GEFS Probability Matched Mean MSLP, Spread '+'\n'+date.strftime('%Y/%m/%d %Hz ')+ str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz'),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/mslp/pmm/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/mslp/pmm/'+loc+'/'+date.strftime('%Y%m%d%H')+
                str(datefhour)+'pmm'+loc+'.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()


    cf = ax3.contourf(x,y,ssaAnom,[-8,-5,-3,-2,-1,1,2,3,5,8],
                      transform=ccrs.PlateCarree(),colors=colorsdiv,zorder=2)
    c1=ax3.contour(x,y,ensMean/100,levels=mslspace,colors='k',linewidths=0.5,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,ticks=[-8,-5,-3,-2,-1,1,2,3,5,8],
                      orientation='horizontal',pad=0.01,aspect=50)
    cbar.set_label('Standardized Anomaly (sigma)')
    clab = plt.clabel(c1,fontsize=10,inline_spacing=-0.5,fmt='%3.0f')
    plt.title('Standardized Spread Anomaly, MSLP '+'\n'+date.strftime('%Y/%m/%d %Hz ')+
              str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/mslp/ssa/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/mslp/ssa/'+loc+'/'+date.strftime('%Y%m%d%H')+
                str(datefhour)+'ssa'+loc+'.png',bbox_inches='tight')

    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()

    colors3 = ((0.000, 0.094, 0.537,0),
            (0.855, 1.000, 0.278,1),
            (0.925, 0.753, 0.000,1),
            (0.910, 0.522, 0.227,1),
            (0.824, 0.306, 0.443,1),
            (0.671, 0.078, 0.533,1))
    cf = ax3.contourf(x,y,subsetPerc,[0,90,95,99,99.5,99.9,100],transform=ccrs.PlateCarree(),colors=colors3,zorder=2)
    c1=ax3.contour(x,y,ensMean/100,levels=mslspace,colors='k',linewidths=0.5,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,ticks=[0,90,95,99,99.5,99.9,100],pad=0.01, orientation='horizontal')
    cbar.set_label('Restricted Anomaly Percentile')
    clab = plt.clabel(c1,fontsize=10,inline_spacing=-0.5,fmt='%3.0f')
    plt.title('Spread Percentile (Restricted M-Climate), '+date.strftime('%Y/%m/%d %Hz ')+'\n'+str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/mslp/subsetperc/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/mslp/subsetperc/'+loc+'/'+date.strftime('%Y%m%d%H')+
                str(datefhour)+'sP'+loc+'.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()

    cf = ax3.contourf(x,y,totalPerc,[0,90,95,99,99.5,99.9,100],transform=ccrs.PlateCarree(),colors=colors3,zorder=2)
    c1=ax3.contour(x,y,ensMean/100,levels=mslspace,colors='k',linewidths=0.5,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,pad=0.01,ticks=[0,90,95,99,99.5,99.9,100], orientation='horizontal')
    cbar.set_label('All Anomaly Percentile')
    clab = plt.clabel(c1,fontsize=10,inline_spacing=-0.5,fmt='%3.0f')
    plt.title('Spread Percentile (All M-Climate), '+date.strftime('%Y/%m/%d %Hz ')+'\n'+str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/mslp/totalperc/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/mslp/totalperc/'+loc+'/'+date.strftime('%Y%m%d%H')+
                str(datefhour)+'tP'+loc+'.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()
    plt.close(fig)

    gc.collect()

# 850mb temperature plot maker


def tmpplotMaker(date, ensMean, ensStd, datefhour, dateArr, ssaAnom,
                 subsetPerc, totalPerc, pmm, lats, lons):
    extents=([-180, -50, 20, 65],'na')
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
    colors = ((255/255, 255/255, 255/255, 0), (236./255.,192./255.,0/255.,1), (232./255.,133./255.,58./255.,1),(210./255.,78./255.,113./255.,1),
              (171./255.,20./255.,136./255.,1), (114./255.,0.,141./255.,1), (10./255.,45./255.,110./255.,1))
     ##Default Variables
    states = ft.NaturalEarthFeature(category='cultural',scale='50m',facecolor='none',
                                    name='admin_1_states_provinces_lines')

    x,y=np.meshgrid(lons,lats)
    sprdspace=[0,0.5,1,2,4,8,10,15]
    tmpspace = np.arange(-80,50,1)
    extent = extents[0]
    loc = extents[1]
    fig=plt.figure(figsize=(15,15))
    ax3 = plt.axes(projection=proj)
    ax3.set_extent(extent,crs=ccrs.PlateCarree())
    ax3.coastlines(resolution=('50m'), zorder=4)
    ax3.add_feature(ft.BORDERS,alpha=0.7,zorder=3)
    ax3.add_feature(ft.RIVERS,alpha=0.4,edgecolor='gray', zorder=3)
    ax3.add_feature(states, edgecolor='gray', zorder=3)
    ax3.add_feature(ft.LAKES,alpha=0.5,facecolor='gray', zorder=3)

    cf = ax3.contourf(x,y,ensStd,sprdspace,
                  transform=ccrs.PlateCarree(),colors=colors,zorder=2)
    c1=ax3.contour(x,y,ensMean,levels=tmpspace,colors='k',
                   linewidths=0.6,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,ticks=sprdspace, orientation='horizontal',
                      pad=0.01,extend='max',extendrect=True)
    c1.collections[np.where(c1.levels==0)[0][0]].set_linewidth(0.9)
    c1.collections[np.where(c1.levels==0)[0][0]].set_color('#d62919')
    c1.collections[np.where(c1.levels==0)[0][0]].set_linestyle('dotted')
    [c1.collections[np.where(c1.levels<0)[0][q]].set_color('#00798e') for q in range(0,len(np.where(c1.levels<0)[0][:]))]
    [c1.collections[np.where(c1.levels<0)[0][q]].set_linewidth(0.55) for q in range(0,len(np.where(c1.levels<0)[0][:]))]
    cbar.set_label('Spread (C)')
    clab = plt.clabel(c1,fontsize=12,inline_spacing=-0.5,fmt='%3.0f')
    plt.title('GEFS Ensemble Mean 850mb Temperature, Spread '+'\n'+date.strftime('%Y/%m/%d %Hz ')+
              str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/tmp/me/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/tmp/me/'+loc+'/'+date.strftime('%Y%m%d%H')+
                str(datefhour)+'me'+loc+'.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()

    cf = ax3.contourf(x,y,ensStd,sprdspace,
                  transform=ccrs.PlateCarree(),colors=colors,zorder=2)
    c1=ax3.contour(x,y,pmm,levels=tmpspace,colors='k',
                   linewidths=0.6,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,ticks=sprdspace, orientation='horizontal',
                      pad=0.01,extend='max',extendrect=True)
    c1.collections[np.where(c1.levels==0)[0][0]].set_linewidth(0.9)
    c1.collections[np.where(c1.levels==0)[0][0]].set_color('#d62919')
    c1.collections[np.where(c1.levels==0)[0][0]].set_linestyle('dotted')
    [c1.collections[np.where(c1.levels<0)[0][q]].set_color('#00798e') for q in range(0,len(np.where(c1.levels<0)[0][:]))]
    [c1.collections[np.where(c1.levels<0)[0][q]].set_linewidth(0.55) for q in range(0,len(np.where(c1.levels<0)[0][:]))]
    cbar.set_label('Spread (mm)')
    clab = plt.clabel(c1,fontsize=12,inline_spacing=-0.5,fmt='%3.0f')
    plt.title('GEFS Probability Matched Mean 850mb Temperature, Spread '+'\n'+date.strftime('%Y/%m/%d %Hz ')+
              str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/tmp/pmm/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/tmp/pmm/'+loc+'/'+date.strftime('%Y%m%d%H')+
                str(datefhour)+'pmm'+loc+'.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()


    cf = ax3.contourf(x,y,ssaAnom,[-8,-5,-3,-2,-1,1,2,3,5,8],
                      transform=ccrs.PlateCarree(),colors=colorsdiv,zorder=2)
    c1=ax3.contour(x,y,ensMean,levels=tmpspace,colors='k',linewidths=0.5,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,ticks=[-8,-5,-3,-2,-1,1,2,3,5,8],
                      orientation='horizontal',pad=0.01,aspect=50)
    cbar.set_label('Standardized Anomaly (sigma)')
    c1.collections[np.where(c1.levels==0)[0][0]].set_linewidth(0.9)
    c1.collections[np.where(c1.levels==0)[0][0]].set_color('#d62919')
    c1.collections[np.where(c1.levels==0)[0][0]].set_linestyle('dotted')
    [c1.collections[np.where(c1.levels<0)[0][q]].set_color('#00798e') for q in range(0,len(np.where(c1.levels<0)[0][:]))]
    [c1.collections[np.where(c1.levels<0)[0][q]].set_linewidth(0.55) for q in range(0,len(np.where(c1.levels<0)[0][:]))]
    clab = plt.clabel(c1,fontsize=10,inline_spacing=-0.5,fmt='%3.0f')
    plt.title('Standardized Spread Anomaly, 850 Temp '+'\n'+date.strftime('%Y/%m/%d %Hz ')+
              str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/tmp/ssa/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/tmp/ssa/'+loc+'/'+date.strftime('%Y%m%d%H')+
                str(datefhour)+'ssa'+loc+'.png',bbox_inches='tight')

    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()

    colors3 = ((0.000, 0.094, 0.537,0),
            (0.855, 1.000, 0.278,1),
            (0.925, 0.753, 0.000,1),
            (0.910, 0.522, 0.227,1),
            (0.824, 0.306, 0.443,1),
            (0.671, 0.078, 0.533,1))
    cf = ax3.contourf(x,y,subsetPerc,[0,90,95,99,99.5,99.9,100],transform=ccrs.PlateCarree(),colors=colors3,zorder=2)
    c1=ax3.contour(x,y,ensMean,levels=tmpspace,colors='k',linewidths=0.5,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,ticks=[0,90,95,99,99.5,99.9,100],pad=0.01, orientation='horizontal')
    c1.collections[np.where(c1.levels==0)[0][0]].set_linewidth(0.9)
    c1.collections[np.where(c1.levels==0)[0][0]].set_color('#d62919')
    c1.collections[np.where(c1.levels==0)[0][0]].set_linestyle('dotted')
    [c1.collections[np.where(c1.levels<0)[0][q]].set_color('#00798e') for q in range(0,len(np.where(c1.levels<0)[0][:]))]
    [c1.collections[np.where(c1.levels<0)[0][q]].set_linewidth(0.55) for q in range(0,len(np.where(c1.levels<0)[0][:]))]
    cbar.set_label('Restricted Anomaly Percentile')
    clab = plt.clabel(c1,fontsize=10,inline_spacing=-0.5,fmt='%3.0f')
    plt.title('Spread Percentile (Restricted M-Climate), '+date.strftime('%Y/%m/%d %Hz ')+'\n'+str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/tmp/subsetperc/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/tmp/subsetperc/'+loc+'/'+date.strftime('%Y%m%d%H')+
                str(datefhour)+'sP'+loc+'.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()

    cf = ax3.contourf(x,y,totalPerc,[0,90,95,99,99.5,99.9,100],transform=ccrs.PlateCarree(),colors=colors3,zorder=2)
    c1=ax3.contour(x,y,ensMean,levels=tmpspace,colors='k',linewidths=0.5,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,pad=0.01,ticks=[0,90,95,99,99.5,99.9,100], orientation='horizontal')
    c1.collections[np.where(c1.levels==0)[0][0]].set_linewidth(0.9)
    c1.collections[np.where(c1.levels==0)[0][0]].set_color('#d62919')
    c1.collections[np.where(c1.levels==0)[0][0]].set_linestyle('dotted')
    [c1.collections[np.where(c1.levels<0)[0][q]].set_color('#00798e') for q in range(0,len(np.where(c1.levels<0)[0][:]))]
    [c1.collections[np.where(c1.levels<0)[0][q]].set_linewidth(0.55) for q in range(0,len(np.where(c1.levels<0)[0][:]))]

    cbar.set_label('All Anomaly Percentile')
    clab = plt.clabel(c1,fontsize=10,inline_spacing=-0.5,fmt='%3.0f')
    plt.title('Spread Percentile (All M-Climate), '+date.strftime('%Y/%m/%d %Hz ')+'\n'+str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/tmp/totalperc/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/tmp/totalperc/'+loc+'/'+date.strftime('%Y%m%d%H')+
                str(datefhour)+'tP'+loc+'.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()
    plt.close(fig)

    gc.collect()

# 500mb heights plot maker


def hgtplotMaker(date, ensMean, ensStd, datefhour, dateArr, ssaAnom,
                 subsetPerc, totalPerc, pmm, lats, lons):
    extents=([-180, -50, 20, 65],'na')
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

    colors = ((255/255, 255/255, 255/255, 0), (236./255.,192./255.,0/255.,1), (232./255.,133./255.,58./255.,1),(210./255.,78./255.,113./255.,1),
              (171./255.,20./255.,136./255.,1), (114./255.,0.,141./255.,1), (10./255.,45./255.,110./255.,1))
    ##Default Variables
    states = ft.NaturalEarthFeature(category='cultural',scale='50m',facecolor='none',
                                    name='admin_1_states_provinces_lines')

    x,y=np.meshgrid(lons,lats)
    sprdspace=[0,1,2,4,8,12,15,20]
    hgtspace = range(460,600,6)
    extent = extents[0]
    loc = extents[1]
    fig=plt.figure(figsize=(15,15))
    ax3 = plt.axes(projection=proj)
    ax3.set_extent(extent,crs=ccrs.PlateCarree())
    ax3.coastlines(resolution=('50m'), zorder=4)
    ax3.add_feature(ft.BORDERS,alpha=0.7,zorder=3)
    ax3.add_feature(ft.RIVERS,alpha=0.4,edgecolor='gray', zorder=3)
    ax3.add_feature(states, edgecolor='gray', zorder=3)
    ax3.add_feature(ft.LAKES,alpha=0.5,facecolor='gray', zorder=3)

    cf = ax3.contourf(x,y,ensStd/10,sprdspace,
                  transform=ccrs.PlateCarree(),colors=colors,zorder=2)

    c1=ax3.contour(x,y,ensMean/10,levels=hgtspace,colors='k',
                   linewidths=0.5,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,ticks=sprdspace, orientation='horizontal',
                      pad=0.01,extend='max',extendrect=True)
    cbar.set_label('Spread (dam)')
    clab = plt.clabel(c1,fontsize=10,inline_spacing=-0.5,fmt='%3.0f')
    plt.title('GEFS Ensemble Mean 500mb Hgts, Spread '+date.strftime('%Y/%m/%d %Hz')+'\n'+
              str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/hgt/me/'+loc+'/init.png',
                    bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/hgt/me/'+loc+'/'+date.strftime('%Y%m%d%H')
    +str(datefhour)+'me'+loc+'.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()


    cf = ax3.contourf(x,y,ensStd/10,sprdspace,
                  transform=ccrs.PlateCarree(),colors=colors,zorder=2)
    cf.cmap.set_over('#0a2d84')
    c1=ax3.contour(x,y,pmm/10,levels=hgtspace,colors='k',linewidths=0.5,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,ticks=sprdspace, orientation='horizontal',pad=0.01,extend='max',extendrect=True)
    cbar.set_label('Spread (dam)')
    clab = plt.clabel(c1,fontsize=10,inline_spacing=-0.5,fmt='%3.0f')
    plt.title('GEFS Probability Matched Mean 500mb Hgts, Spread '+'\n'+date.strftime('%Y/%m/%d %Hz ')+ str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz'),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/hgt/pmm/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/hgt/pmm/'+loc+'/'+date.strftime('%Y%m%d%H')+
                str(datefhour)+'pmm'+loc+'.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()


    cf = ax3.contourf(x,y,ssaAnom,[-8,-5,-3,-2,-1,1,2,3,5,8],
                      transform=ccrs.PlateCarree(),colors=colorsdiv,zorder=2)
    c1=ax3.contour(x,y,ensMean/10,levels=hgtspace,colors='k',linewidths=0.5,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,ticks=[-8,-5,-3,-2,-1,1,2,3,5,8],
                      orientation='horizontal',pad=0.01,aspect=50)
    cbar.set_label('Standardized Anomaly (sigma)')
    clab = plt.clabel(c1,fontsize=10,inline_spacing=-0.5,fmt='%3.0f')
    plt.title('Standardized Spread Anomaly, 500mb Hgts '+'\n'+date.strftime('%Y/%m/%d %Hz ')+
              str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/hgt/ssa/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/hgt/ssa/'+loc+'/'+date.strftime('%Y%m%d%H')+
                str(datefhour)+'ssa'+loc+'.png',bbox_inches='tight')

    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()


    colors3 = ((0.000, 0.094, 0.537,0),
            (0.855, 1.000, 0.278,1),
            (0.925, 0.753, 0.000,1),
            (0.910, 0.522, 0.227,1),
            (0.824, 0.306, 0.443,1),
            (0.671, 0.078, 0.533,1))
    cf = ax3.contourf(x,y,subsetPerc,[0,90,95,99,99.5,99.9,100],transform=ccrs.PlateCarree(),colors=colors3,zorder=2)
    c1=ax3.contour(x,y,ensMean/10,levels=hgtspace,colors='k',linewidths=0.5,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,ticks=[0,90,95,99,99.5,99.9,100],pad=0.01, orientation='horizontal')
    cbar.set_label('Restricted Anomaly Percentile')
    clab = plt.clabel(c1,fontsize=10,inline_spacing=-0.5,fmt='%3.0f')
    plt.title('Spread Percentile (Restricted M-Climate), '+date.strftime('%Y/%m/%d %Hz ')+'\n'+str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/hgt/subsetperc/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/hgt/subsetperc/'+loc+'/'+date.strftime('%Y%m%d%H')+
                str(datefhour)+'sP'+loc+'.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()


    cf = ax3.contourf(x,y,totalPerc,[0,90,95,99,99.5,99.9,100],transform=ccrs.PlateCarree(),colors=colors3,zorder=2)
    c1=ax3.contour(x,y,ensMean/10,levels=hgtspace,colors='k',linewidths=0.5,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,pad=0.01,ticks=[0,90,95,99,99.5,99.9,100], orientation='horizontal')
    cbar.set_label('All Anomaly Percentile')
    clab = plt.clabel(c1,fontsize=10,inline_spacing=-0.5,fmt='%3.0f')
    plt.title('Spread Percentile (All M-Climate), '+date.strftime('%Y/%m/%d %Hz ')+'\n'+str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/hgt/totalperc/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/hgt/totalperc/'+loc+'/'+date.strftime('%Y%m%d%H')+
                str(datefhour)+'tP'+loc+'.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()
    plt.close(fig)

    gc.collect()

# Precipitable water plot maker


def pwatplotMaker(date, ensMean, ensStd, datefhour, dateArr, ssaAnom,
                 subsetPerc, totalPerc, pmm, lats, lons):
    extents=([-180, -50, 20, 65],'na')

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
    colors = ((255/255, 255/255, 255/255, 0), (236./255.,192./255.,0/255.,1), (232./255.,133./255.,58./255.,1),(210./255.,78./255.,113./255.,1),
              (171./255.,20./255.,136./255.,1), (114./255.,0.,141./255.,1), (10./255.,45./255.,110./255.,1))
     ##Default Variables
    states = ft.NaturalEarthFeature(category='cultural',scale='50m',facecolor='none',
                                    name='admin_1_states_provinces_lines')

    x,y=np.meshgrid(lons,lats)
    sprdspace=[0,0.5,1,2,4,8,10,15]
    pwatspace = [0,1,2,5,10,20,30,40,50,60]
    extent = extents[0]
    loc = extents[1]
    fig=plt.figure(figsize=(15,15))
    ax3 = plt.axes(projection=proj)
    ax3.set_extent(extent,crs=ccrs.PlateCarree())
    ax3.coastlines(resolution=('50m'), zorder=4)
    ax3.add_feature(ft.BORDERS,alpha=0.7,zorder=3)
    ax3.add_feature(ft.RIVERS,alpha=0.4,edgecolor='gray', zorder=3)
    ax3.add_feature(states, edgecolor='gray', zorder=3)
    ax3.add_feature(ft.LAKES,alpha=0.5,facecolor='gray', zorder=3)

    cf = ax3.contourf(x,y,ensStd,sprdspace,
                  transform=ccrs.PlateCarree(),colors=colors,zorder=2)
    c1=ax3.contour(x,y,ensMean,levels=pwatspace,colors='k',linewidths=1.2,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,ticks=sprdspace, orientation='horizontal',
                      pad=0.01,extend='max',extendrect=True)
    cbar.set_label('Spread (mm)')
    clab = plt.clabel(c1,fontsize=14,inline_spacing=-0.5,fmt='%3.0f',color='k')
    plt.title('GEFS Ensemble Precipitable Water Mean, Spread '+'\n'+date.strftime('%Y/%m/%d %Hz ')+
              str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/pwat/me/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/pwat/me/'+loc+'/'+date.strftime('%Y%m%d%H')+
                str(datefhour)+'me'+loc+'.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()

    cf = ax3.contourf(x,y,ensStd,sprdspace,
                  transform=ccrs.PlateCarree(),colors=colors,zorder=2)
    c1=ax3.contour(x,y,pmm,levels=pwatspace,colors='k',linewidths=1.2,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,ticks=sprdspace, orientation='horizontal',
                      pad=0.01,extend='max',extendrect=True)
    cbar.set_label('Spread (mm)')
    for l in c1:
        l.set_rotation(0) 
    clab = plt.clabel(c1,fontsize=14,inline_spacing=-0.5,fmt='%3.0f',color='k')
    plt.title('GEFS Probability Matched Mean Precipitable Water, Spread '+'\n'+date.strftime('%Y/%m/%d %Hz ')+
              str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/pwat/pmm/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/pwat/pmm/'+loc+'/'+date.strftime('%Y%m%d%H')+
                str(datefhour)+'pmm'+loc+'.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()




    cf = ax3.contourf(x,y,ssaAnom,[-8,-5,-3,-2,-1,1,2,3,5,8],
                      transform=ccrs.PlateCarree(),colors=colorsdiv,zorder=2)
    c1=ax3.contour(x,y,ensMean,levels=pwatspace,colors='k',linewidths=1.2,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,ticks=[-8,-5,-3,-2,-1,1,2,3,5,8],
                      orientation='horizontal',pad=0.01,aspect=50)
    cbar.set_label('Standardized Anomaly (sigma)')
    clab = plt.clabel(c1,fontsize=14,inline_spacing=-0.5,fmt='%3.0f',color='k')
    plt.title('Standardized Spread Anomaly, Precipitable Water '+'\n'+date.strftime('%Y/%m/%d %Hz ')+
              str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/pwat/ssa/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/pwat/ssa/'+loc+'/'+date.strftime('%Y%m%d%H')+
                str(datefhour)+'ssa'+loc+'.png',bbox_inches='tight')

    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()

    colors3 = ((0.000, 0.094, 0.537,0),
            (0.855, 1.000, 0.278,1),
            (0.925, 0.753, 0.000,1),
            (0.910, 0.522, 0.227,1),
            (0.824, 0.306, 0.443,1),
            (0.671, 0.078, 0.533,1))
    cf = ax3.contourf(x,y,subsetPerc,[0,90,95,99,99.5,99.9,100],transform=ccrs.PlateCarree(),colors=colors3,zorder=2)
    c1=ax3.contour(x,y,ensMean,levels=pwatspace,colors='k',linewidths=1.2,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,ticks=[0,90,95,99,99.5,99.9,100],pad=0.01, orientation='horizontal')
    cbar.set_label('Restricted Anomaly Percentile')
    clab = plt.clabel(c1,fontsize=14,inline_spacing=-0.5,fmt='%3.0f',color='k')
    plt.title('Spread Percentile (Restricted M-Climate), '+date.strftime('%Y/%m/%d %Hz ')+'\n'+str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/pwat/subsetperc/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/pwat/subsetperc/'+loc+'/'+date.strftime('%Y%m%d%H')+
                str(datefhour)+'sP'+loc+'.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()


    cf = ax3.contourf(x,y,totalPerc,[0,90,95,99,99.5,99.9,100],transform=ccrs.PlateCarree(),colors=colors3,zorder=2)
    c1=ax3.contour(x,y,ensMean,levels=pwatspace,colors='k',linewidths=1.2,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,pad=0.01,ticks=[0,90,95,99,99.5,99.9,100], orientation='horizontal')
    cbar.set_label('All Anomaly Percentile')
    clab = plt.clabel(c1,fontsize=14,inline_spacing=-0.5,fmt='%3.0f',color='k')
    plt.title('Spread Percentile (All M-Climate), '+date.strftime('%Y/%m/%d %Hz ')+'\n'+str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/pwat/totalperc/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/pwat/totalperc/'+loc+'/'+date.strftime('%Y%m%d%H')+
                str(datefhour)+'tP'+loc+'.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()
    plt.close(fig)

    gc.collect()


def slpplotMakerY(date, ensMean, ensStd, datefhour, dateArr, pmm, lats, lons):
    extents=([-180, -50, 20, 65],'na')
    datefhour = int(datefhour)

    colors = ((255/255, 255/255, 255/255, 0),
              (236./255., 192./255., 0/255., 1),
              (232./255., 133./255., 58./255., 1),
              (210./255., 78./255., 113./255., 1),
              (171./255., 20./255., 136./255., 1),
              (114./255., 0., 141./255., 1),
              (10./255., 45./255., 110./255., 1))
    # Default Variables
    states = ft.NaturalEarthFeature(category='cultural', scale='50m',
                                    facecolor='none',
                                    name='admin_1_states_provinces_lines')

    x, y = np.meshgrid(lons, lats)
    sprdspace = [0, 1, 2, 4, 8, 12, 15, 20]
    extent = extents[0]
    loc = extents[1]
    fig = plt.figure(figsize=(20, 20))
    ax3 = plt.axes(projection=proj)
    ax3.set_extent(extent,crs=ccrs.PlateCarree())
    ax3.coastlines(resolution=('50m'), zorder=4)
    ax3.add_feature(ft.BORDERS,alpha=0.7,zorder=3)
    ax3.add_feature(ft.RIVERS,alpha=0.4,edgecolor='gray', zorder=3)
    ax3.add_feature(states, edgecolor='gray', zorder=3)
    ax3.add_feature(ft.LAKES,alpha=0.5,facecolor='gray', zorder=3)
    cf = ax3.contourf(x, y, ensStd/100, sprdspace,
                      transform=ccrs.PlateCarree(), colors=colors, zorder=2)

    c1 = ax3.contour(x, y, ensMean/100, levels=mslspace, colors='k',
                     linewidths=0.5,transform=ccrs.PlateCarree())
    cbar = plt.colorbar(cf, aspect=50, ticks=sprdspace,
                        orientation='horizontal',
                        pad=0.01, extend='max', extendrect=True)
    cbar.set_label('Spread (hPa)')
    clab=plt.clabel(c1, fontsize=10, inline_spacing=-0.5, fmt='%3.0f')
    plt.title('GEFS Ensemble Mean MSLP, Spread ' +
              date.strftime('%Y/%m/%d %Hz') + '\n' +
              str(datefhour) + 'h Forecast valid ' +
              dateArr.strftime('%Y/%m/%d %Hz '), loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/mslp/me/'+loc+'/init.png',
                    bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/mslp/me/'+loc+'/'+str(datefhour)+'me'+loc+'.png', bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()

    cf = ax3.contourf(x,y,ensStd/100,sprdspace,
                  transform=ccrs.PlateCarree(),colors=colors,zorder=2)
    cf.cmap.set_over('#0a2d84')
    c1=ax3.contour(x,y,pmm/100,levels=mslspace,colors='k',linewidths=0.5,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,ticks=sprdspace, orientation='horizontal',pad=0.01,extend='max',extendrect=True)
    cbar.set_label('Spread (hPa)')
    clab = plt.clabel(c1,fontsize=10,inline_spacing=-0.5,fmt='%3.0f')
    plt.title('GEFS Probability Matched Mean MSLP, Spread '+'\n'+date.strftime('%Y/%m/%d %Hz ')+ str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz'),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/mslp/pmm/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/mslp/pmm/'+loc+'/'+str(datefhour)+'pmm'+loc+'.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()
    plt.close(fig)

    gc.collect()

# 850mb temperature plot maker


def tmpplotMakerY(date, ensMean, ensStd, datefhour, dateArr, pmm, lats, lons):
    extents=([-180, -50, 20, 65],'na')
    datefhour = int(datefhour)

    colors = ((255/255, 255/255, 255/255, 0), (236./255.,192./255.,0/255.,1), (232./255.,133./255.,58./255.,1),(210./255.,78./255.,113./255.,1),
              (171./255.,20./255.,136./255.,1), (114./255.,0.,141./255.,1), (10./255.,45./255.,110./255.,1))
     ##Default Variables
    states = ft.NaturalEarthFeature(category='cultural',scale='50m',facecolor='none',
                                    name='admin_1_states_provinces_lines')

    x,y=np.meshgrid(lons,lats)
    sprdspace=[0,0.5,1,2,4,8,10,15]
    tmpspace = np.arange(-80,50,1)
    extent = extents[0]
    loc = extents[1]
    fig=plt.figure(figsize=(15,15))
    ax3 = plt.axes(projection=proj)
    ax3.set_extent(extent,crs=ccrs.PlateCarree())
    ax3.coastlines(resolution=('50m'), zorder=4)
    ax3.add_feature(ft.BORDERS,alpha=0.7,zorder=3)
    ax3.add_feature(ft.RIVERS,alpha=0.4,edgecolor='gray', zorder=3)
    ax3.add_feature(states, edgecolor='gray', zorder=3)
    ax3.add_feature(ft.LAKES,alpha=0.5,facecolor='gray', zorder=3)

    cf = ax3.contourf(x,y,ensStd,sprdspace,
                  transform=ccrs.PlateCarree(),colors=colors,zorder=2)
    c1=ax3.contour(x,y,ensMean,levels=tmpspace,colors='k',
                   linewidths=0.6,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,ticks=sprdspace, orientation='horizontal',
                      pad=0.01,extend='max',extendrect=True)
    c1.collections[np.where(c1.levels==0)[0][0]].set_linewidth(0.9)
    c1.collections[np.where(c1.levels==0)[0][0]].set_color('#d62919')
    c1.collections[np.where(c1.levels==0)[0][0]].set_linestyle('dotted')
    [c1.collections[np.where(c1.levels<0)[0][q]].set_color('#00798e') for q in range(0,len(np.where(c1.levels<0)[0][:]))]
    [c1.collections[np.where(c1.levels<0)[0][q]].set_linewidth(0.55) for q in range(0,len(np.where(c1.levels<0)[0][:]))]
    cbar.set_label('Spread (C)')
    clab = plt.clabel(c1,fontsize=12,inline_spacing=-0.5,fmt='%3.0f')
    plt.title('GEFS Ensemble Mean 850mb Temperature, Spread '+'\n'+date.strftime('%Y/%m/%d %Hz ')+
              str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/tmp/me/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/tmp/me/'+loc+'/'+str(datefhour)+'me'+loc+'.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()

    cf = ax3.contourf(x,y,ensStd,sprdspace,
                  transform=ccrs.PlateCarree(),colors=colors,zorder=2)
    c1=ax3.contour(x,y,pmm,levels=tmpspace,colors='k',
                   linewidths=0.6,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,ticks=sprdspace, orientation='horizontal',
                      pad=0.01,extend='max',extendrect=True)
    c1.collections[np.where(c1.levels==0)[0][0]].set_linewidth(0.9)
    c1.collections[np.where(c1.levels==0)[0][0]].set_color('#d62919')
    c1.collections[np.where(c1.levels==0)[0][0]].set_linestyle('dotted')
    [c1.collections[np.where(c1.levels<0)[0][q]].set_color('#00798e') for q in range(0,len(np.where(c1.levels<0)[0][:]))]
    [c1.collections[np.where(c1.levels<0)[0][q]].set_linewidth(0.55) for q in range(0,len(np.where(c1.levels<0)[0][:]))]
    cbar.set_label('Spread (mm)')
    clab = plt.clabel(c1,fontsize=12,inline_spacing=-0.5,fmt='%3.0f')
    plt.title('GEFS Probability Matched Mean 850mb Temperature, Spread '+'\n'+date.strftime('%Y/%m/%d %Hz ')+
              str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/tmp/pmm/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/tmp/pmm/'+loc+'/'+str(datefhour)+'pmm'+loc+'.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()
    plt.close(fig)

    gc.collect()
# 500mb heights plot maker


def hgtplotMakerY(date, ensMean, ensStd, datefhour, dateArr, pmm, lats, lons):
    extents=([-180, -50, 20, 65],'na')
    datefhour = int(datefhour)


    colors = ((255/255, 255/255, 255/255, 0), (236./255.,192./255.,0/255.,1), (232./255.,133./255.,58./255.,1),(210./255.,78./255.,113./255.,1),
              (171./255.,20./255.,136./255.,1), (114./255.,0.,141./255.,1), (10./255.,45./255.,110./255.,1))
    ##Default Variables
    states = ft.NaturalEarthFeature(category='cultural',scale='50m',facecolor='none',
                                    name='admin_1_states_provinces_lines')

    x,y=np.meshgrid(lons,lats)
    sprdspace=[0,1,2,4,8,12,15,20]
    hgtspace = range(460,600,6)
    extent = extents[0]
    loc = extents[1]
    fig=plt.figure(figsize=(15,15))
    ax3 = plt.axes(projection=proj)
    ax3.set_extent(extent,crs=ccrs.PlateCarree())
    ax3.coastlines(resolution=('50m'), zorder=4)
    ax3.add_feature(ft.BORDERS,alpha=0.7,zorder=3)
    ax3.add_feature(ft.RIVERS,alpha=0.4,edgecolor='gray', zorder=3)
    ax3.add_feature(states, edgecolor='gray', zorder=3)
    ax3.add_feature(ft.LAKES,alpha=0.5,facecolor='gray', zorder=3)

    cf = ax3.contourf(x,y,ensStd/10,sprdspace,
                  transform=ccrs.PlateCarree(),colors=colors,zorder=2)

    c1=ax3.contour(x,y,ensMean/10,levels=hgtspace,colors='k',
                   linewidths=0.5,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,ticks=sprdspace, orientation='horizontal',
                      pad=0.01,extend='max',extendrect=True)
    cbar.set_label('Spread (dam)')
    clab = plt.clabel(c1,fontsize=10,inline_spacing=-0.5,fmt='%3.0f')
    plt.title('GEFS Ensemble Mean 500mb Hgts, Spread '+date.strftime('%Y/%m/%d %Hz')+'\n'+
              str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/hgt/me/'+loc+'/init.png',
                    bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/hgt/me/'+loc+'/'+str(datefhour)+'me'+loc+'.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()


    cf = ax3.contourf(x,y,ensStd/10,sprdspace,
                  transform=ccrs.PlateCarree(),colors=colors,zorder=2)
    cf.cmap.set_over('#0a2d84')
    c1=ax3.contour(x,y,pmm/10,levels=hgtspace,colors='k',linewidths=0.5,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,ticks=sprdspace, orientation='horizontal',pad=0.01,extend='max',extendrect=True)
    cbar.set_label('Spread (dam)')
    clab = plt.clabel(c1,fontsize=10,inline_spacing=-0.5,fmt='%3.0f')
    plt.title('GEFS Probability Matched Mean 500mb Hgts, Spread '+'\n'+date.strftime('%Y/%m/%d %Hz ')+ str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz'),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/hgt/pmm/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/hgt/pmm/'+loc+'/'+str(datefhour)+'pmm'+loc+'.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()
    plt.close(fig)

    gc.collect()
# Precipitable water plot maker


def pwatplotMakerY(date, ensMean, ensStd, datefhour, dateArr, pmm, lats, lons):
    extents=([-180, -50, 20, 65],'na')

    datefhour = int(datefhour)

    colors = ((255/255, 255/255, 255/255, 0), (236./255.,192./255.,0/255.,1), (232./255.,133./255.,58./255.,1),(210./255.,78./255.,113./255.,1),
              (171./255.,20./255.,136./255.,1), (114./255.,0.,141./255.,1), (10./255.,45./255.,110./255.,1))
     ##Default Variables
    states = ft.NaturalEarthFeature(category='cultural',scale='50m',facecolor='none',
                                    name='admin_1_states_provinces_lines')

    x,y=np.meshgrid(lons,lats)
    sprdspace=[0,0.5,1,2,4,8,10,15]
    pwatspace = [0,1,2,5,10,20,30,40,50,60]
    extent = extents[0]
    loc = extents[1]
    fig=plt.figure(figsize=(15,15))
    ax3 = plt.axes(projection=proj)
    ax3.set_extent(extent,crs=ccrs.PlateCarree())
    ax3.coastlines(resolution=('50m'), zorder=4)
    ax3.add_feature(ft.BORDERS,alpha=0.7,zorder=3)
    ax3.add_feature(ft.RIVERS,alpha=0.4,edgecolor='gray', zorder=3)
    ax3.add_feature(states, edgecolor='gray', zorder=3)
    ax3.add_feature(ft.LAKES,alpha=0.5,facecolor='gray', zorder=3)

    cf = ax3.contourf(x,y,ensStd,sprdspace,
                  transform=ccrs.PlateCarree(),colors=colors,zorder=2)
    c1=ax3.contour(x,y,ensMean,levels=pwatspace,colors='k',linewidths=1.2,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,ticks=sprdspace, orientation='horizontal',
                      pad=0.01,extend='max',extendrect=True)
    cbar.set_label('Spread (mm)')
    clab = plt.clabel(c1,fontsize=14,inline_spacing=-0.5,fmt='%3.0f',color='k')
    for l in clab:
        	l.set_rotation(0)
    plt.title('GEFS Ensemble Precipitable Water Mean, Spread '+'\n'+date.strftime('%Y/%m/%d %Hz ')+
              str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/pwat/me/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/pwat/me/'+loc+'/'+str(datefhour)+'me'+loc+'.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()

    cf = ax3.contourf(x,y,ensStd,sprdspace,
                  transform=ccrs.PlateCarree(),colors=colors,zorder=2)
    c1=ax3.contour(x,y,pmm,levels=pwatspace,colors='k',linewidths=1.2,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,ticks=sprdspace, orientation='horizontal',
                      pad=0.01,extend='max',extendrect=True)
    cbar.set_label('Spread')  

    clab = plt.clabel(c1,fontsize=14,inline_spacing=-0.5,fmt='%3.0f',color='k')
    for l in clab:
        l.set_rotation(0) 
    plt.title('GEFS Probability Matched Mean Precipitable Water, Spread '+'\n'+date.strftime('%Y/%m/%d %Hz ')+
              str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/pwat/pmm/'+loc+'/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/pwat/pmm/'+loc+'/'+str(datefhour)+'pmm'+loc+'.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()
    plt.close(fig)

    gc.collect()

def qpfplotMakerY(date, ensMean, ensStd, datefhour, dateArr, pmm, lats, lons):
    extents=([-180, -50, 20, 65],'na')

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
    colors = ((255/255, 255/255, 255/255, 0), (236./255.,192./255.,0/255.,1), (232./255.,133./255.,58./255.,1),(210./255.,78./255.,113./255.,1),
              (171./255.,20./255.,136./255.,1), (114./255.,0.,141./255.,1), (10./255.,45./255.,110./255.,1))
     ##Default Variables
    states = ft.NaturalEarthFeature(category='cultural',scale='50m',facecolor='none',
                                    name='admin_1_states_provinces_lines')

    x,y=np.meshgrid(lons,lats)
    sprdspace=[0,0.5,1,2,4,8,10,15]
    pwatspace = [0.1,.5,1,2,5,10,20]
    extent = extents[0]
    fig=plt.figure(figsize=(15,15))
    ax3 = plt.axes(projection=proj)
    ax3.set_extent(extent,crs=ccrs.PlateCarree())
    ax3.coastlines(resolution=('50m'), zorder=4)
    ax3.add_feature(ft.BORDERS,alpha=0.7,zorder=3)
    ax3.add_feature(ft.RIVERS,alpha=0.4,edgecolor='gray', zorder=3)
    ax3.add_feature(states, edgecolor='gray', zorder=3)
    ax3.add_feature(ft.LAKES,alpha=0.5,facecolor='gray', zorder=3)

    cf = ax3.contourf(x,y,ensStd,sprdspace,
                  transform=ccrs.PlateCarree(),colors=colors,zorder=2)
    c1=ax3.contour(x,y,ensMean,levels=pwatspace,colors='k',linewidths=1.2,transform=ccrs.PlateCarree())

    cbar=plt.colorbar(cf,aspect=50,ticks=sprdspace, orientation='horizontal',
                      pad=0.01,extend='max',extendrect=True)
    cbar.set_label('Spread (mm)')
    clab = plt.clabel(c1,fontsize=14,inline_spacing=-0.5,fmt='%3.0f',color='k')
    for l in clab:
        l.set_rotation(0) 
    plt.title('GEFS Ensemble QPF Mean, Spread '+'\n'+date.strftime('%Y/%m/%d %Hz ')+
              str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/qpf/me/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/qpf/me/'+str(datefhour)+'me.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()

    cf = ax3.contourf(x,y,ensStd,sprdspace,
                  transform=ccrs.PlateCarree(),colors=colors,zorder=2)
    c1=ax3.contour(x,y,pmm,levels=pwatspace,colors='k',linewidths=1.2,transform=ccrs.PlateCarree())
    cbar=plt.colorbar(cf,aspect=50,ticks=sprdspace, orientation='horizontal',
                      pad=0.01,extend='max',extendrect=True)
    cbar.set_label('Spread (mm)')
    clab = plt.clabel(c1,fontsize=14,inline_spacing=-0.5,fmt='%3.0f',color='k')
    for l in clab:
        l.set_rotation(0) 
    plt.title('GEFS Probability Matched Mean QPF, Spread '+'\n'+date.strftime('%Y/%m/%d %Hz ')+
              str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    if datefhour == 6 or datefhour == 0:
        plt.savefig('/home/taylorm/ssa/images/qpf/pmm/init.png',bbox_inches='tight')
    plt.savefig('/home/taylorm/ssa/images/qpf/pmm/'+str(datefhour)+'pmm.png',bbox_inches='tight')
    for a in c1.collections:
        a.remove()
    for a in cf.collections:
        a.remove()
    for a in clab:
        a.remove()
    cbar.remove()
    plt.close(fig)

    gc.collect()
