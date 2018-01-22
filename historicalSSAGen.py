#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 15:46:56 2017

@author: taylorm
"""

import numpy as np
import scipy.io as si
from netCDF4 import Dataset
from datetime import datetime, timedelta
import pdb
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as ft
from math import isnan
import seaborn as sns
import math
from matplotlib.colors import Normalize
import matplotlib as mpl
from scipy.stats import percentileofscore,binned_statistic,rankdata
from bisect import bisect
import pandas as pd
import mclifuncs as mc
import gc 

def loadfunc(dates,time):

    for x in range(1,21):
        try:
            ncepmem = Dataset('/EDATA1/TIGGE/NCEP/'+dates.strftime('%Y')+dates.strftime('%m')+dates.strftime('%d')+'/NCEP_'+dates.strftime('%Y')+dates.strftime('%m')+dates.strftime('%d')+'_mem_%02d.nc' % x)
            slpncep = ncepmem['SLP'][time]
            slpncep = np.expand_dims(slpncep,axis=0)
            ncepmem.close()
        except IOError:
            print dates.strftime('%Y')+dates.strftime('%m')+dates.strftime('%d')+'/NCEP_'+dates.strftime('%Y')+dates.strftime('%m')+dates.strftime('%d')+'_mem_%02d.nc' % x

            pdb.set_trace()
        if x == 1:
            expslpncep = slpncep
        else:
            try:
                expslpncep = np.append(slpncep,expslpncep,axis=0)

            except ValueError:
                pdb.set_trace()
    for x in range(1,21):
        try:
            cmcmem = Dataset('/EDATA1/TIGGE/CMC/'+dates.strftime('%Y')+dates.strftime('%m')+dates.strftime('%d')+'/CMC_'+dates.strftime('%Y')+dates.strftime('%m')+dates.strftime('%d')+'_mem_%02d.nc' % x)
            slpcmc = cmcmem['SLP'][time]
            slpcmc = np.expand_dims(slpcmc,axis=0)
            cmcmem.close()
        except IOError:
            print dates.strftime('%Y')+dates.strftime('%m')+dates.strftime('%d')+'/CMC_'+dates.strftime('%Y')+dates.strftime('%m')+dates.strftime('%d')+'_mem_%02d.nc' % x

            pdb.set_trace()
        if x == 1:
            expslpcmc = slpcmc
        else:
            try:
                expslpcmc = np.append(slpcmc,expslpcmc,axis=0)

            except ValueError:
                pdb.set_trace()
    #ssaAnom = zscorefunc(expslpmean,expslpstd,dates)
    return np.expand_dims(expslpncep,axis=0), np.expand_dims(expslpcmc,axis=0)

def slpplotMaker(date, ensMean, ensStd, datefhour, dateArr, ssaAnom,
                 subsetPerc, totalPerc, lats, lons):
    extents=([-100, -40, 25, 50])
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
    mslspace = range(900, 1100, 2)

    fig = plt.figure(figsize=(20, 20))
    ax3 = plt.axes(projection=ccrs.PlateCarree())
    ax3.set_extent(extents)
    ax3.coastlines(resolution=('10m'), zorder=4)
    ax3.add_feature(ft.BORDERS, alpha=0.7, zorder=3)
    ax3.add_feature(ft.OCEAN, facecolor='#d2eef2', zorder=1)
    ax3.add_feature(ft.LAND,  facecolor='#f5ffe0', zorder=1)
    ax3.add_feature(ft.RIVERS, alpha=0.8)
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
              date.strftime('%Y/%m/%d %Hz') + '\n' +
              str(datefhour) + 'h Forecast valid ' +
              dateArr.strftime('%Y/%m/%d %Hz '), loc='left')
    plt.savefig('/home/taylorm/ssa/'+dateArr.strftime('%d%b%Y')+'/dprogdt/me/'+date.strftime('%Y%m%d%H')
    +str(datefhour)+'me.png', bbox_inches='tight')
    plt.close(fig)

    fig=plt.figure(figsize=(20,20))
    ax3 = plt.axes(projection=ccrs.PlateCarree())
    ax3.set_extent(extents)
    ax3.coastlines(resolution=('10m'),zorder=4)
    ax3.add_feature(ft.BORDERS,alpha=0.7,zorder=3)
    ax3.add_feature(ft.OCEAN,facecolor='#d2eef2',zorder=1)
    ax3.add_feature(ft.LAND, facecolor='#f5ffe0',zorder=1)
    ax3.add_feature(ft.RIVERS,alpha=0.8)
    ax3.add_feature(states,edgecolor='gray',zorder=3)


    cf = ax3.contourf(x,y,ssaAnom,[-8,-5,-3,-2,-1,1,2,3,5,8],
                      transform=ccrs.PlateCarree(),colors=colorsdiv,zorder=2)
    c1=ax3.contour(x,y,ensMean/100,levels=mslspace,colors='k',linewidths=0.5)
    cbar=plt.colorbar(cf,ticks=[-8,-5,-3,-2,-1,1,2,3,5,8],
                      orientation='horizontal',pad=0.01,aspect=50)
    cbar.set_label('Standardized Anomaly (sigma)')
    plt.clabel(c1,fontsize=10,inline_spacing=-0.5,fmt='%3.0f')
    plt.title('Standardized Spread Anomaly, MSLP '+'\n'+date.strftime('%Y/%m/%d %Hz ')+
              str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    plt.savefig('/home/taylorm/ssa/'+dateArr.strftime('%d%b%Y')+'/dprogdt/ssa/'+date.strftime('%Y%m%d%H')+
                str(datefhour)+'ssa.png',bbox_inches='tight')

    plt.close(fig)

    fig=plt.figure(figsize=(20,20))
    ax3 = plt.axes(projection=ccrs.PlateCarree())
    ax3.set_extent(extents)
    ax3.coastlines(resolution=('10m'),zorder=4)
    ax3.add_feature(ft.BORDERS,alpha=0.7,zorder=3)
    ax3.add_feature(ft.OCEAN,facecolor='#d2eef2',zorder=1)
    ax3.add_feature(ft.LAND, facecolor='#f5ffe0',zorder=1)
    ax3.add_feature(ft.RIVERS,alpha=0.8)
    ax3.add_feature(states,edgecolor='gray',zorder=3)
    colors3 = ((0.000, 0.094, 0.537,0),
            (0.855, 1.000, 0.278,1),
            (0.925, 0.753, 0.000,1),
            (0.910, 0.522, 0.227,1),
            (0.824, 0.306, 0.443,1),
            (0.671, 0.078, 0.533,1))
    cf = ax3.contourf(x,y,subsetPerc,[0,90,95,99,99.5,99.9,100],transform=ccrs.PlateCarree(),colors=colors3,zorder=2)
    c1=ax3.contour(x,y,ensMean/100,levels=mslspace,colors='k',linewidths=0.5)
    cbar=plt.colorbar(cf,aspect=50,ticks=[0,90,95,99,99.5,99.9,100],pad=0.01, orientation='horizontal')
    cbar.set_label('Restricted Anomaly Percentile')
    plt.clabel(c1,fontsize=10,inline_spacing=-0.5,fmt='%3.0f')
    plt.title('Spread Percentile (Restricted M-Climate), '+date.strftime('%Y/%m/%d %Hz ')+'\n'+str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    plt.savefig('/home/taylorm/ssa/'+dateArr.strftime('%d%b%Y')+'/dprogdt/sP/'+date.strftime('%Y%m%d%H')+
                str(datefhour)+'sP.png',bbox_inches='tight')
    plt.close(fig)

    fig=plt.figure(figsize=(20,20))
    ax3 = plt.axes(projection=ccrs.PlateCarree())
    ax3.set_extent(extents)
    ax3.coastlines(resolution=('10m'),zorder=4)
    ax3.add_feature(ft.BORDERS,alpha=0.7,zorder=3)
    ax3.add_feature(ft.OCEAN,facecolor='#d2eef2',zorder=1)
    ax3.add_feature(ft.LAND, facecolor='#f5ffe0',zorder=1)
    ax3.add_feature(ft.RIVERS,alpha=0.8)
    ax3.add_feature(states,edgecolor='gray',zorder=3)

    cf = ax3.contourf(x,y,totalPerc,[0,90,95,99,99.5,99.9,100],transform=ccrs.PlateCarree(),colors=colors3,zorder=2)
    c1=ax3.contour(x,y,ensMean/100,levels=mslspace,colors='k',linewidths=0.5)
    cbar=plt.colorbar(cf,aspect=50,pad=0.01,ticks=[0,90,95,99,99.5,99.9,100], orientation='horizontal')
    cbar.set_label('All Anomaly Percentile')
    plt.clabel(c1,fontsize=10,inline_spacing=-0.5,fmt='%3.0f')
    plt.title('Spread Percentile (All M-Climate), '+date.strftime('%Y/%m/%d %Hz ')+'\n'+str(datefhour) +'h Forecast valid ' + dateArr.strftime('%Y/%m/%d %Hz '),loc='left')
    plt.savefig('/home/taylorm/ssa/'+dateArr.strftime('%d%b%Y')+'/dprogdt/tP/'+date.strftime('%Y%m%d%H')+
                str(datefhour)+'tP.png',bbox_inches='tight')
    plt.close(fig)

    gc.collect()

    
def mcliLoad(var=None, time=None, ind=None):
    print 'loading mcli array'
    filenames = ['/home/taylorm/mcli/mclidata/'+var+'NHm.nc',
             '/home/taylorm/mcli/mclidata/'+var+'NHs.nc']
    if ind is not None:

        if var == 'mslp':
            varnc = [Dataset(fn).variables['Pressure'] for fn in filenames]
            mArr = np.array(varnc[0][ind,4::4,::-1,:])
            sArr = np.array(varnc[1][ind,4::4,::-1,:])
        elif var == '850tmp':
            filenames = ['/home/taylorm/mcli/mclidata/'+var+'NHmn.nc',
                         '/home/taylorm/mcli/mclidata/'+var+'NHsn.nc']
            varnc = [Dataset(fn).variables['Temperature'] for fn in filenames]
            mArr = np.squeeze(np.array(varnc[0][ind,...,::-1,:]))
            sArr = np.squeeze(np.array(varnc[1][ind,...,::-1,:]))
        elif var == '500hgt':
            filenames = ['/home/taylorm/mcli/mclidata/'+var+'NHm.nc',
                         '/home/taylorm/mcli/mclidata/'+var+'NHs.nc']
            varnc = [Dataset(fn).variables['Geopotential_height'] for fn in filenames]
            mArr = np.squeeze(np.array(varnc[0][ind,...,::-1,:]))
            sArr = np.squeeze(np.array(varnc[1][ind,...,::-1,:]))
        elif var == 'pwat':
            filenames = ['/home/taylorm/mcli/mclidata/'+var+'NHm.nc',
                         '/home/taylorm/mcli/mclidata/'+var+'NHs.nc']
            varnc = [Dataset(fn).variables['Precipitable_water'] for fn in filenames]
            mArr = np.squeeze(np.array(varnc[0][ind,...,::-1,:]))
            sArr = np.squeeze(np.array(varnc[1][ind,...,::-1,:]))
        # Reforecast and GEFS have flipped lat axes,
        # so I'm flipping the Reforecast's

        print 'loaded mcli array'
        gc.collect()
        return mArr,sArr
    else:
        dataArr = 0
        return dataArr
    

def slpFuncHist(mslpMean, mslpStd, mslppmm, datefhour, dateArr,
            lats, lons, date, ind):
    subsetPerc = np.ones_like(mslpMean)
    totalPerc = np.ones_like(mslpMean)
    ssaAnom = np.ones_like(mslpMean)
    mArr,sArr = mc.mcliLoad(var='mslp', ind=ind)
    for i in range(0, len(mslpMean)):
        subsetPerc[i], totalPerc[i], ssaAnom[i] = mc.subsetMCli(mslpMean[i], mslpStd[i], mArr[:,i], sArr[:, i])
    print 'starting slp plots'
    [slpplotMaker(date, mslpMean[i], mslpStd[i], datefhour[i],
     dateArr[i], ssaAnom[i], subsetPerc[i], totalPerc[i], mslppmm[i],
     lats, lons) for i in range(0, len(mslpMean))]
    gc.collect()

def ssaGen(validTime):
#    validTime = datetime(2011,1,27,0,0)
    forecastTimes = [(validTime-timedelta(hours=t*24)) for t in np.arange(1,7.1,1)]
    
    expslpncep = np.ones((7,20,41,61))
    expslpcmc = np.ones((7,20,41,61))
    #expAgg = np.concatenate((expslpncep,expslpcmc),axis=1)
    mCliDates = mc.mcliTimeArray(var='mslp')

    subsetPerc = np.ones((7,41,61))
    totalPerc = np.ones((7,41,61))
    ssaAnom = np.ones((7,41,61))
    lats = np.arange(25,65.1,1)[::-1]
    lons = np.arange(260,310.1,1)
        
    for t in np.arange(1,7.1,1):
        t=int(t)
        hour = t*24

        mslpStd = np.load()   
        ind = mc.mcliSub(mCliDates, forecastTimes[t-1])
        mArr,sArr = mcliLoad(var='mslp', ind=ind)
        subsetPerc[t-1], totalPerc[t-1], ssaAnom[t-1] = mc.subsetMCli(mslpMean, mslpStd, mArr[:,t-1], sArr[:, t-1])

#        slpplotMaker(forecastTimes[t-1], mslpMean, mslpStd, hour,
#        validTime, ssaAnom[t-1], subsetPerc[t-1], totalPerc[t-1],
#        lats, lons)
        gc.collect()      
    np.save('ssaAnomSub.npy',ssaAnom)

if __name__ == "__main__":
    sns.set(font_scale=1.65, style="whitegrid", color_codes=True)
    with open('allobsdates1.txt') as f:
        cont = f.readlines()
    cont = [x.strip() for x in cont]
    cont[8]='2008-12-22 00:00:00'
    cont[17]='2009-01-24 00:00:00'
    cont[20]='2008-12-26 00:00:00'
        
        
        	
        
    cont = [datetime.strptime(date, '%Y-%m-%d %H:%M:%S') for date in cont]
    cont = [date.replace(hour=00) for date in cont]
    [ssaGen(validTime) for validTime in 

    
    
    