#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 14:02:15 2018

@author: taylorm
"""

import numpy as np
from netCDF4 import Dataset

class netCDFSSA:
    
    def __init__(self, ssa, time):
        self.ssa = ssa
        self.time = time
        
    def detNewApp(self):
        try:
            ds=Dataset('/home/taylorm/mcli/verifData.nc','a')
            netCDFSSA.appendNC(self,ds)
        except IOError:
            netCDFSSA.createNew(self)
            
    def createNew(self):
        ds = Dataset('/home/taylorm/mcli/verifData.nc','w')
        latitude = ds.createDimension("lat", size=61)
        longitude = ds.createDimension("lon", size=141) 
        fh = ds.createDimension("fhour", size=7)
        init = ds.createDimension("initTime",size=None)

        inits = ds.createVariable('initTime',np.float32,('initTime',))
        lats = ds.createVariable('latitude',np.int32, ('lat',))
        lons = ds.createVariable('longitude',np.int32, ('lon',))    
        fhour = ds.createVariable('fhour',np.int32, ('fhour',))
        sspanom = ds.createVariable('ssa', np.float64,('initTime','fhour','lat','lon'))
        

        inits.units = 'days from 0001-1-1 00:00 (ordinal time)'
        lats[:] = np.linspace(80,20,61)
        lons[:] = np.linspace(180,320,141)
        fhour[:]= np.arange(24,168.1,24)
        inits[:]=self.time
        sspanom[:]=self.ssa
        ds.close()
        
    def appendNC(self,ds):
        times = ds.variables['initTime']
        if self.time == times[len(times)-1]:
            pass
        else:
            ssa = ds.variables['ssa']
            ssa[len(ssa)] = self.ssa
            times[len(times)] = self.time
        ds.close()


        
        
        
        
