#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import matplotlib as mpl
mpl.use('agg')
import numpy as np
from datetime import timedelta, datetime
import seaborn as sns
from matplotlib.colors import Normalize
import sys
import mclifuncs as mc
import plotter as pt
import gc

from bs4 import BeautifulSoup
import requests
"""
Created on Wed Jun  7 16:29:35 2017

@author: taylorm
"""


class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


def mainfunc(token=0):
    mc.GEFScheck()

    ind, date, dateArr, mslpMean, mslpStd, mslppmm, lats, lons, tmpMean, tmpStd, hgtMean, hgtStd, tmppmm, hgtpmm,qpfMean, qpfStd, qpfpmm, pwatMean, pwatStd, pwatpmm = mc.liveLoad()

    gc.collect()
    datefhour = dateArr
    dateArr = [date+timedelta(hours=n) for n in dateArr]

    if token==1:
        slpFuncY(mslpMean, mslpStd, mslppmm, datefhour, dateArr,
                 lats, lons, date, ind)
        gc.collect()
        tmpFuncY(tmpMean, tmpStd, tmppmm, datefhour, dateArr,
                 lats, lons, date, ind)
        gc.collect()
        hgtFuncY(hgtMean, hgtStd, hgtpmm, datefhour, dateArr,
                 lats, lons, date, ind)
        gc.collect()
        pwatFuncY(pwatMean, pwatStd, pwatpmm, datefhour, dateArr,
                 lats, lons, date, ind)
        gc.collect()
        qpfFuncY(qpfMean, qpfStd, qpfpmm, datefhour, dateArr,
                 lats, lons, date, ind)
        gc.collect()
    else:
        slpFunc(mslpMean, mslpStd, mslppmm, datefhour, dateArr,
                 lats, lons, date, ind)
        gc.collect()
        tmpFunc(tmpMean, tmpStd, tmppmm, datefhour, dateArr,
                 lats, lons, date, ind)
        gc.collect()
        hgtFunc(hgtMean, hgtStd, hgtpmm, datefhour, dateArr,
                 lats, lons, date, ind)
        gc.collect()
        pwatFunc(pwatMean, pwatStd, pwatpmm, datefhour, dateArr,
                 lats, lons, date, ind)
        gc.collect()
    logfi = open('/home/taylorm/ssa/run.txt', 'w')
    # write the details into the log to have a record of
    # what runs were ran through the code
    logfi.write(date.strftime('%Y%m%d%H'))
    # closes file

    logfi.close()
    url1 = 'http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GEFS/Global_1p0deg_Ensemble/members/latest.html'

    # Play games with the url to get the TDSCatalog for siphon
    response = requests.get(url1)
    page = str(BeautifulSoup(response.content, 'lxml'))
    start_link = page.find("a href")

    start_quote = page.find('"', start_link)
    end_quote = page.find('"', start_quote + 1)
    url = page[start_quote + 1: end_quote]

    url = url.split('html')
#    endurl = url[1]

#    # load in year and month as string
#    year = endurl[83:87]
#    month = endurl[87:89]
#
#    # Load in run and day as string
#    run = endurl[92:94]
#    day = endurl[89:91]

    # Path to log
    pathCheck = '/home/taylorm/mcli/logs/'

    with open(pathCheck+'liveLog', "w") as f:
        f.seek(0)
        f.write("done")

    # path = '/home/taylorm/ssa'


def slpFuncY(mslpMean, mslpStd, mslppmm, datefhour, dateArr,
            lats, lons, date, ind):
    subsetPerc = np.ones_like(mslpMean)
    totalPerc = np.ones_like(mslpMean)
    ssaAnom = np.ones_like(mslpMean)
    saAnom = np.ones_like(mslpMean)
    if datetime.now().month >= 3 and datetime.now().month<=5:
        mArr,sArr = mc.mcliLoad(var='mslp',ind=ind,notDJF='MAM')
    elif datetime.now().month >= 6 and datetime.now().month<=8:
        mArr,sArr = mc.mcliLoad(var='mslp',ind=ind,notDJF='JJA')
    else:
        mArr,sArr = mc.mcliLoad(var='mslp', ind=ind)

    for i in range(0, len(mslpMean)):
        subsetPerc[i], totalPerc[i], ssaAnom[i], saAnom[i] = mc.subsetMCli(mslpMean[i], mslpStd[i], mArr[:,i], sArr[:, i])
    
    print 'completed SSA'
    print 'starting slp plots'
    if datetime.now().month >= 3 and datetime.now().month<=5:     
        [pt.slpplotMaker(date, mslpMean[i], mslpStd[i], datefhour[i],
     dateArr[i], ssaAnom[i], subsetPerc[i], totalPerc[i], mslppmm[i],
     lats, lons) for i in range(0, len(mslpMean))]
    elif datetime.now().month >= 6 and datetime.now().month<=8:
        [pt.slpplotMaker(date, mslpMean[i], mslpStd[i], datefhour[i],
     dateArr[i], ssaAnom[i], subsetPerc[i], totalPerc[i], mslppmm[i],
     lats, lons) for i in range(0, len(mslpMean))]
    else:
        [pt.slpplotMakerY(date, mslpMean[i], mslpStd[i], datefhour[i],
         dateArr[i],mslppmm[i], lats, lons) for i in range(0, len(mslpMean))]
    gc.collect()


def tmpFuncY(tmpMean, tmpStd, tmppmm, datefhour, dateArr,
            lats, lons, date, ind):
    print 'starting tmp plots'
    [pt.tmpplotMakerY(date, tmpMean[i], tmpStd[i], datefhour[i],
     dateArr[i], tmppmm[i], lats, lons) for i in range(0, len(tmpMean))]
    gc.collect()


def hgtFuncY(hgtMean, hgtStd, hgtpmm, datefhour, dateArr,
            lats, lons, date, ind):
    print 'starting hgt plots'
    [pt.hgtplotMakerY(date, hgtMean[i], hgtStd[i], datefhour[i],
     dateArr[i], hgtpmm[i], lats, lons) for i in range(0, len(hgtMean))]
    gc.collect()


def pwatFuncY(pwatMean, pwatStd, pwatpmm, datefhour, dateArr,
             lats, lons, date, ind):
    print 'starting pwat plots'
    [pt.pwatplotMakerY(date, pwatMean[i], pwatStd[i], datefhour[i],
     dateArr[i], pwatpmm[i], lats, lons) for i in range(0, len(pwatMean))]

def qpfFuncY(qpfMean, qpfStd, qpfpmm, datefhour, dateArr,
            lats, lons, date, ind):
    print 'starting qpf plots'
    [pt.qpfplotMakerY(date, qpfMean[i], qpfStd[i], datefhour[i],
     dateArr[i], qpfpmm[i], lats, lons) for i in range(0, len(qpfMean))]
    gc.collect()

def slpFunc(mslpMean, mslpStd, mslppmm, datefhour, dateArr,
            lats, lons, date, ind):
    subsetPerc = np.ones_like(mslpMean)
    totalPerc = np.ones_like(mslpMean)
    ssaAnom = np.ones_like(mslpMean)
    saAnom = np.ones_like(mslpMean)
    if datetime.now().month >= 3 or datetime.now().month<=5:
        mArr,sArr = mc.mcliLoad(var='mslp',ind=ind,notDJF='MAM')
    mArr,sArr = mc.mcliLoad(var='mslp', ind=ind)

    for i in range(0, len(mslpMean)):
        subsetPerc[i], totalPerc[i], ssaAnom[i], saAnom[i] = mc.subsetMCli(mslpMean[i], mslpStd[i], mArr[:,i], sArr[:, i])
    print 'starting slp plots'
    [pt.slpplotMaker(date, mslpMean[i], mslpStd[i], datefhour[i],
     dateArr[i], ssaAnom[i], subsetPerc[i], totalPerc[i], mslppmm[i],
     lats, lons) for i in range(0, len(mslpMean))]
    gc.collect()


def tmpFunc(tmpMean, tmpStd, tmppmm, datefhour, dateArr,
            lats, lons, date, ind):
    subsetPerc = np.ones_like(tmpMean)
    totalPerc = np.ones_like(tmpMean)
    ssaAnom = np.ones_like(tmpMean)
    saAnom = np.ones_like(tmpMean)
    mArr,sArr = mc.mcliLoad(var='850tmp', ind=ind)

    for i in range(0, len(tmpMean)):
        subsetPerc[i], totalPerc[i], ssaAnom[i], saAnom[i] = mc.subsetMCli(tmpMean[i], tmpStd[i], mArr[:, i], sArr[:,i])
    print 'starting tmp plots'
    [pt.tmpplotMaker(date, tmpMean[i], tmpStd[i], datefhour[i], dateArr[i], ssaAnom[i], subsetPerc[i], totalPerc[i], tmppmm[i], lats, lons)
     for i in range(0, len(tmpMean))]
    gc.collect()


def hgtFunc(hgtMean, hgtStd, hgtpmm, datefhour, dateArr,
            lats, lons, date, ind):
    subsetPerc = np.ones_like(hgtMean)
    totalPerc = np.ones_like(hgtMean)
    ssaAnom = np.ones_like(hgtMean)
    saAnom = np.ones_like(hgtMean)
    mArr,sArr = mc.mcliLoad(var='500hgt', ind=ind)
    for i in range(0, len(hgtMean)):
        subsetPerc[i], totalPerc[i], ssaAnom[i],saAnom[i] = mc.subsetMCli(hgtMean[i], hgtStd[i], mArr[:, i], sArr[:,i])
    print 'starting hgt plots'
    [pt.hgtplotMaker(date, hgtMean[i], hgtStd[i], datefhour[i], dateArr[i], ssaAnom[i], subsetPerc[i], totalPerc[i], hgtpmm[i], lats, lons)
     for i in range(0, len(hgtMean))]
    gc.collect()


def pwatFunc(pwatMean, pwatStd, pwatpmm, datefhour, dateArr,
             lats, lons, date, ind):
    subsetPerc = np.ones_like(pwatMean)
    totalPerc = np.ones_like(pwatMean)
    ssaAnom = np.ones_like(pwatMean)
    saAnom = np.ones_like(pwatMean)
    mArr,sArr = mc.mcliLoad(var='pwat', ind=ind)
    for i in range(0, len(pwatMean)):
        subsetPerc[i], totalPerc[i], ssaAnom[i],saAnom[i] = mc.subsetMCli(pwatMean[i], pwatStd[i], mArr[:,i], sArr[:, i])
    print 'starting pwat plots'
    [pt.pwatplotMaker(date, pwatMean[i], pwatStd[i], datefhour[i],
     dateArr[i], ssaAnom[i], subsetPerc[i], totalPerc[i], pwatpmm[i],
     lats, lons) for i in range(0, len(pwatMean))]
    gc.collect()

if __name__ == "__main__":
    
    sns.set(font_scale=1.65, style="whitegrid", color_codes=True)
    pathCheck = '/home/taylorm/mcli/logs/'
    if datetime.now().month >= 3 or datetime.now.month() < 12:
        try:
            mainfunc(token=1)
            sys.exit()
        except (KeyboardInterrupt, SystemExit):
            with open(pathCheck+'liveLog', "w") as f:
                f.seek(0)
                f.write("done")
                sys.exit()
    else:
        try:
            mainfunc()
            sys.exit()
        except (KeyboardInterrupt, SystemExit):
            with open(pathCheck+'liveLog', "w") as f:
                f.seek(0)
                f.write("done")
                sys.exit()
