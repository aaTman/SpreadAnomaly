#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 18:19:24 2017

@author: taylorm
"""
import matplotlib as mpl
mpl.use('agg')
import numpy as np
import sys
from bs4 import BeautifulSoup
from datetime import datetime, timedelta
from siphon import catalog, ncss
import requests
from netCDF4 import Dataset
import pdb
import gc
from bisect import bisect


# The crontab killer when the run has already completed - the code inserts the
# run into a log file and this checks if it already completed.
def percentileofscore(a, score, kind='rank'):

    a = np.array(a)
    n = len(a)

    if kind == 'rank':
        if not np.any(a == score):
            a = np.append(a, score)
            a_len = np.array(list(range(len(a))))
        else:
            a_len = np.array(list(range(len(a)))) + 1.0

        a = np.sort(a)
        idx = [a == score]
        pct = (np.mean(a_len[idx]) / n) * 100.0
        return pct

    elif kind == 'strict':
        return np.sum(a < score) / float(n) * 100
    elif kind == 'weak':
        return np.sum(a <= score) / float(n) * 100
    elif kind == 'mean':
        return (np.sum(a < score) + np.sum(a <= score)) * 50 / float(n)
    else:
        raise ValueError("kind can only be 'rank', 'strict', 'weak' or 'mean'")

def GEFScheck(init=0):
    if init == 1:
        url1 = 'http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GEFS/Global_1p0deg_Ensemble/members/latest.html'

        # Play games with the url to get the TDSCatalog for siphon
        response = requests.get(url1)
        page = str(BeautifulSoup(response.content, 'lxml'))
        start_link = page.find("a href")

        start_quote = page.find('"', start_link)
        end_quote = page.find('"', start_quote + 1)
        url = page[start_quote + 1: end_quote]

        url = url.split('html')
        endurl = url[1]

        # load in year and month as string
        year = endurl[83:87]
        month = endurl[87:89]

        # Load in run and day as string
        run = endurl[92:94]
        day = endurl[89:91]

        # Path to log
        pathCheck = '/home/taylorm/mcli/logs/'

        logC = open(pathCheck+'log', 'a')
        logC.write(run+' '+day+' '+month+' '+year+'\n')
        logC.close()
    else:
        # Path to log
        pathCheck = '/home/taylorm/mcli/logs/'
        with open(pathCheck+'liveLog', "r+") as f:
            old = f.read()
            
            f.seek(0)
            if old == 'running':
                sys.exit()
            f.write("running")

            
        # load the thredds latest GEFS run
        url1 = 'http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GEFS/Global_1p0deg_Ensemble/members/latest.html'

        # Play games with the url to get the TDSCatalog for siphon
        response = requests.get(url1)
        page = str(BeautifulSoup(response.content, 'lxml'))
        start_link = page.find("a href")

        start_quote = page.find('"', start_link)
        end_quote = page.find('"', start_quote + 1)
        url = page[start_quote + 1: end_quote]

        url = url.split('html')
        endurl = url[1]

        # load in year and month as string
        year = endurl[83:87]
        month = endurl[87:89]

        # Load in run and day as string
        run = endurl[92:94]
        day = endurl[89:91]



        # Checks if the log file has been made. If not, code will run.
        try:
            logCheck = open(pathCheck+'log')
        except IOError:
            print 'log file not created yet'
            return

        # write the details into the log to have a record
        # of what runs were ran through the code

        # lastCheck is the log file text, testStr is the run from THREDDS
        lastCheck = logCheck.readlines()
        testStr = run+' '+day+' '+month+' '+year+'\n'

        print lastCheck[-1]
        print testStr
        # Quits python/the code if it is the same
        if lastCheck[-1] == testStr:
            sys.exit()

# Load in the timeArray from the m-climate file to prepare to subset the 21-day
# index over the 30 years (erratic index values due to leap years and invalid
# data).


def mcliTimeArray(time=None, var=None):
    # time is used to backtest if needed, var must be declared or else code
    # will error out
    if time is None:

        # location of my m-climate files

        if datetime.now().month == 11:
            varFile = Dataset('/home/taylorm/mcli/mclidata/mslpNovS.nc')
        elif datetime.now().month >=3 and datetime.now().month <= 5:
            varFile = Dataset('/home/taylorm/mcli/mclidata/mslpMAMS.nc')
        else:
            varFile = Dataset('/home/taylorm/mcli/mclidata/mslpNHm.nc')
        # gets the time variable from the netcdf
        varFile = varFile['time'][:]

        # converts array from ordinal time to python datetime instances
        timedt = [datetime(1800, 1, 1, 0, 0)+timedelta(hours=n)
                  for n in varFile]

        return timedt


def liveLoad(wlon=180, elon=310, slat=20, nlat=80):

    print 'loading GEFS file'
    dateArr, date1, mslpMean, mslpStd, mslppmm, lats, lons, tmpMean, tmpStd, hgtMean, hgtStd, tmppmm, hgtpmm,qpfMean, qpfStd, qpfpmm, pwatMean, pwatStd, pwatpmm = GEFSload()


    print 'loaded GEFS file'
    liveDate = datetime.strptime(date1, "%Y%m%d%H")

    dates = mcliTimeArray(var='mslp')
    ind = mcliSub(dates, liveDate)

    print 'returning GEFS values to main m-climate function'
    return ind, liveDate, dateArr, mslpMean, mslpStd, mslppmm, lats, lons, tmpMean, tmpStd, hgtMean, hgtStd, tmppmm, hgtpmm,qpfMean, qpfStd, qpfpmm, pwatMean, pwatStd, pwatpmm


def subsetMCli(fcstm, fcsts, m630, s630):

    mstdgrid = np.ones_like(fcstm)
    mmngrid = np.ones_like(fcstm)
    pgrid = np.ones_like(fcstm)
    bsgrid = np.ones_like(fcstm)

    for y in range(0, np.shape(fcstm)[-2]):
        for z in range(0, np.shape(fcstm)[-1]):
            marg = m630[:, y, z].argsort()
            mbins = m630[marg, y, z]
            sbins = s630[marg, y, z]
            leng = len(mbins)
            centRange = bisect(mbins.tolist(), fcstm[y, z])

            if centRange <= leng/10:
                subsetSpread = sbins[0:leng/10]
                pgrid[y, z] = percentileofscore(subsetSpread, fcsts[y, z])
                mstdgrid[y, z] = np.std(subsetSpread)
                mmngrid[y, z] = np.mean(subsetSpread)
            elif centRange >= leng-leng/10:
                subsetSpread = sbins[leng-leng/10:]
                pgrid[y, z] = percentileofscore(subsetSpread, fcsts[y, z])
                mstdgrid[y, z] = np.std(subsetSpread)
                mmngrid[y, z] = np.mean(subsetSpread)
            else:
                subsetSpread = sbins[centRange-leng/20:centRange+leng/20]
                pgrid[y, z] = percentileofscore(subsetSpread, fcsts[y, z])
                mstdgrid[y, z] = np.std(subsetSpread)
                mmngrid[y, z] = np.mean(subsetSpread)

            bsgrid[y, z] = percentileofscore(sbins, fcsts[y, z])

    ssaAnom = (fcsts - mmngrid)/mstdgrid
    gc.collect()
    return pgrid, bsgrid, ssaAnom


def mcliLoad(var=None, time=None, ind=None,notDJF=None):
    filenames = ['/home/taylorm/mcli/mclidata/'+var+'NHm.nc',
             '/home/taylorm/mcli/mclidata/'+var+'NHs.nc']   

    if ind is not None:

        if var == 'mslp':
            if notDJF is not None:
                filenames = ['/home/taylorm/mcli/mclidata/'+var+notDJF+'M.nc',
                 '/home/taylorm/mcli/mclidata/'+var+notDJF+'S.nc']    
            varnc = [Dataset(fn).variables['Pressure'] for fn in filenames]
            mArr = np.array(varnc[0][ind,...,::-1,:])
            sArr = np.array(varnc[1][ind,...,::-1,:])
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

        gc.collect()
        return mArr,sArr
    else:
        dataArr = 0
        return dataArr


def mcliSub(dateList, fDate):

    years = [n.year for n in dateList]
    years = np.unique(years)

    val = [q for q in dateList for year in years if q-timedelta(days=10) <=
           fDate.replace(year=year) <= q+timedelta(days=10)]

    valset = set(val)
    ind = [i for i, item in enumerate(dateList) if item in valset]
    if len(ind) < 630:
        ind = None
    return ind


def pmm(ens):
    # Takes the mean of the ensemble
    ensMean = np.mean(ens, axis=1)

    # Sets the sorted index array for each forecast hour's ensemble mean values low to high
    ensMeanSort = np.array([np.argsort(ensMean[i].flatten())
                            for i in range(0, len(ensMean))])

    # Sets the index array of the sorted ensemble mean index array; in simple
    # terms this can be used to rearrange the data back to its original order
    ensMeanSortSort = np.array([np.argsort(ensMeanSort[i])
                                for i in range(0, len(ensMeanSort))])

    # Flattens the full ensemble by forecast hour
    ensFlat = [ens[i].flatten() for i in range(0, len(ens))]

    # Takes every 21st value to equal the size of the flattened mean array
    ensFlat = np.array([ensFlat[i][0::21]for i in range(0, len(ens))])

    # Sets the sorted index array from low to high for each forecast hour
    ensSort = np.array([np.argsort(ensFlat[i]) for i in range(0, len(ens))])

    # Replaces the ensemble mean values with the ensemble values, then returns
    # the values to the original index of the ensemble mean.
    enspmm = np.array([ensFlat[i][ensSort[i]][ensMeanSortSort[i]]
                       for i in range(0, len(ensSort))])

    # Catch for errors, haven't had an issue so this might be useless.
    try:
        enspmm = enspmm.reshape((len(ensSort), len(ensMean[0]),
                                 len(ensMean[0, 0])))
    except ValueError:

        pdb.set_trace()

    return enspmm


def GEFSload(wlon=180, elon=310, slat=20, nlat=80):
    # Initialize current time variable
    now = datetime.now()

    # load the thredds latest GEFS run
    url1 = 'http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GEFS/Global_1p0deg_Ensemble/members/latest.html'

    # Play games with the url to get the TDSCatalog for siphon
    response = requests.get(url1)
    page = str(BeautifulSoup(response.content, 'lxml'))
    start_link = page.find("a href")

    start_quote = page.find('"', start_link)
    end_quote = page.find('"', start_quote + 1)
    url = page[start_quote + 1: end_quote]

    url = url.split('html')
    endurl = url[1]

    # get the TDSCatalog file
    gefs = catalog.TDSCatalog(url1+endurl)
    ds = list(gefs.datasets.values())[0]

    # Pull the netcdf subset
    ncssl = ncss.NCSS(ds.access_urls['NetcdfSubset'])

    query = ncssl.query()
    # print the run that was loaded
    print 'obtained ' + endurl[83:94] + 'z run'

    # load in year and month as string
    year = endurl[83:87]
    month = endurl[87:89]

    print 'taking time range from ' + str(now) + ' to ' + str(now+timedelta(days=6.75))

    # Load in run and day as string
    run = endurl[92:94]
    day = endurl[89:91]

    # take only the latlon box, netcdf4 file, and the mentioned variables,
    # then assign them to variables
    query.lonlat_box(wlon, elon, slat,
                     nlat).time_range(now, now+timedelta(days=7))

    query.accept('netcdf4')
    query.variables('Pressure_reduced_to_MSL_msl_ens',
                    'Temperature_isobaric_ens',
                    'Geopotential_height_isobaric_ens',
                    'Precipitable_water_entire_atmosphere_single_layer_ens',
                    'Total_precipitation_surface_6_Hour_Accumulation_ens')
    data = ncssl.get_data(query)

    try:
        dateArr = data.variables['time2'][:]
    except KeyError:
        try:
            dateArr = data.variables['time'][:]
        except KeyError:
            try:
                dateArr = data.variables['time1'][:]
            except KeyError:
                try:
                    dateArr = data.variables['time3'][:]
                except KeyError:
                    pdb.set_trace()
    qpf = data.variables['Total_precipitation_surface_6_Hour_Accumulation_ens'][:]
    pwat = data.variables['Precipitable_water_entire_atmosphere_single_layer_ens'][:]
    mslp = data.variables['Pressure_reduced_to_MSL_msl_ens'][:]
    lats = data.variables['lat'][:]
    lons = data.variables['lon'][:]


    tmps = data.variables['Temperature_isobaric_ens'][:]
    tmps = np.squeeze(tmps[:, :, np.where(data.variables['isobaric2'][:]
                           == 85000), ...])
    hgt = data.variables['Geopotential_height_isobaric_ens'][:]
    hgt = np.squeeze(hgt[:, :, np.where(data.variables['isobaric2'][:]
                         == 50000), ...])


    # pwat=data.variables['Precipitable_water_entire_atmosphere_single_layer_ens']
    # rh700 = data.variables['Relative_humidity_isobaric_ens']

    # take means and stdevs
    qpfMean = np.mean(qpf,axis=1)
    qpfStd = np.std(qpf,axis=1)
    qpfPMM = pmm(qpf)

    pwatMean = np.mean(pwat,axis=1)
    pwatStd = np.std(pwat,axis=1)
    pwatPMM = pmm(pwat)

    mslpMean = np.mean(mslp, axis=1)
    mslpStd = np.std(mslp, axis=1)
    mslpPMM = pmm(mslp)

    tmps = tmps-273.15
    tmpMean = np.mean(tmps, axis=1)
    tmpStd = np.std(tmps, axis=1)
    tmpPMM = pmm(tmps)

    hgtMean = np.mean(hgt, axis=1)
    hgtStd = np.std(hgt, axis=1)
    hgtPMM = pmm(hgt)

    # string of date and run
    date1 = year+month+day+run

    return dateArr, date1, mslpMean, mslpStd, mslpPMM, lats, lons, tmpMean, tmpStd, hgtMean, hgtStd, tmpPMM, hgtPMM, qpfMean, qpfStd, qpfPMM, pwatMean, pwatStd, pwatPMM
