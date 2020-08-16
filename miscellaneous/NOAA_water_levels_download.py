# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 21:01:12 2020

@author: Arslaan Khalid
"""


import pandas as pd
import requests
import json
from datetime import datetime
from collections import OrderedDict
import matplotlib.pyplot as plt




def Get_NOAA_Water_Levels(gage, start, end, datum, time_zone, units):
    #--NOAA API 'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter'    
    fmt       = "json"                            #Format #@ json
    url       = 'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter'
    product   = 'hourly_height'                     #Product
    
    noaa_time_step = '60' + 'Min' # 6,10,30,60
  
    noaa = pd.DataFrame()
    gages = dict()
    
    t0     = start.strftime('%Y%m%d %H:%M')
    t1     = stop.strftime('%Y%m%d %H:%M')
    api_params = {'begin_date': t0, 'end_date': t1,
                'station': gage,'product':product,'datum':datum,'interval':'h',
                'units':units,'time_zone':time_zone,
                'application':'George_Mason_University','format':fmt }
    #--- getting the data
    r = requests.get(url, params = api_params)
    jdata =r.json()
    if 'error' in jdata.keys(): 
        print(jdata['error'],'\n!!!Fix the parameters based on above shown error message!!!')
    else:
        #---- adding data to list
        pred=[];t=[]  
        for j in jdata['data']:t.append(str(j['t']));pred.append(str(j['v']))  
        
        colname = str(gage)    
        noaa[colname]= pred
        noaa[colname] = noaa[colname].astype(float)
          
        idx = pd.date_range(start,periods = len(noaa.index), freq=noaa_time_step)   
        noaa = noaa.set_index(idx)  
        print('NOAA data successful reterived for stn: ' , gage) 
        return noaa

# =============================================================================
# 
# =============================================================================
    
# Start date (year, month, day)
y0, m0 ,d0, h0 = 2017, 2, 6  ,0        
y1, m1 ,d1, h1 = 2017, 2, 22 ,0  
  
# Create Datetime Objects
start     = datetime(y0, m0, d0, h0)    
stop      = datetime(y1, m1 ,d1, h1)  

noaa_vals = Get_NOAA_Water_Levels('8594900',start,stop,'NAVD','gmt','metric')   
