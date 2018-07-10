#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 14:02:15 2018

@author: taylorm
"""

import numpy as np
from netCDF4 import Dataset
import os

class netCDFSSA:
    
    def __init__(self, ssa, time):
        self.ssa = np.expand_dims(ssa,axis=0)
        self.time = time
        
    def detNewApp(self):
        try:
            ds=Dataset('/home/taylorm/mcli/verifData.nc','a')
            netCDFSSA.appendNC(self,ds)
        except IOError:
            netCDFSSA.createNew(self)
            
    def createNew(self):
        try:
            ds = Dataset('/home/taylorm/mcli/verifData.nc','w')
        except IOError:
            os.remove('/home/taylorm/mcli/verifData.nc')
            ds = Dataset('/home/taylorm/mcli/verifData.nc','w')
        ds.createDimension("lat", size=61)
        ds.createDimension("lon", size=141) 
        ds.createDimension("fhour", size=7)
        ds.createDimension("initTime",size=None)

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
        sspanom[0,:,:,:]=self.ssa
        ds.close()
        
    def appendNC(self,ds):
        times = ds.variables['initTime']
        if self.time == times[len(times)-1]:
            pass
        else:
            ssaNC = ds.variables['ssa']
            ssaNC[len(times)] = self.ssa
            times[len(times)] = self.time
        ds.close()


        
        
        
        
