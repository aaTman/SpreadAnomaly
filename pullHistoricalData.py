#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 18:36:38 2017

@author: taylorm
"""

#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
from datetime import datetime, timedelta

server = ECMWFDataServer()
date = [datetime(year,12,01) for year in range(1985,2018)]
for x in date:
    y = x + timedelta(days=89)
    server.retrieve({
        
        "class": "ei",
        "dataset": "interim",
        "date": x.strftime("%Y-%m-%d")+"/to/"+y.strftime("%Y-%m-%d"),
        "expver": "1",
        "grid": "1.0/1.0",
        "levtype": "sfc",
        "param": "151.128",
        "step": "0",
        "stream": "oper",
        "time": "00:00:00",
        "type": "an",
        "target": "eraBacktest"+str(x.year),
        "area":  "80/-180/20/-50",
        "format" : "netcdf",

    })