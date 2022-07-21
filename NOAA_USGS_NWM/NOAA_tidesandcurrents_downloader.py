# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 10:52:18 2022

@author: Arslaan.Khalid
"""


import matplotlib.pyplot as plt
import fileinput
from matplotlib.dates import DayLocator, HourLocator, DateFormatter
import datetime
import pandas as pd
import os
import requests
from scipy.interpolate import interp1d
import csv
import numpy as np



#===================== User Enteries =================
rps = ['5yr','10yr','20yr','50yr','100yr']
aep_erdc = {8747437: [1.8,2.59,3.47,4.41,4.94],8741533:[1.37,1.89,2.43,3.07,3.5]}
strt_end_years_gages = {8747437: [1978,2023],8741533:[2004,2023]}


for StationID in strt_end_years_gages.keys():
    empty = pd.DataFrame()   
    
    startYear=strt_end_years_gages[StationID][0]
    EndYear=strt_end_years_gages[StationID][1]
    # StationID=8747437#8741533#
    
    #AEP1p = 4.94#3.5 ## in meters and NAVD88
    
    
    
    #--NOAA API https://tidesandcurrents.noaa.gov/api/
    datum     = "NAVD"   #"NAVD"                  #Datum
    units     = "metric"                         #Units
    time_zone = "gmt"                         #Time Zone
    fmt       = "json"                            #Format
    url       = 'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter'
    product   = 'hourly_height'                     #Product
    
    
    #https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?product=hourly_height&application=NOS.COOPS.TAC.WL&begin_date=20210720&end_date=20220720&datum=NAVD&station=8747437&time_zone=GMT&units=metric&format=json
    
    
    
    #======================== Python Code ====================
    print(f'-------------------\nStarting NOAA {product} data downloading for {StationID}\nfrom {startYear} to {2023}\n------------------')
    Totalyears = EndYear-startYear
    
    for k in range(0,Totalyears):
        #print(k)
        
        startYear=startYear
        
        EndYear=startYear+1   
        # print(startYear)
        # print(EndYear)
        
        start, freq = "01-{}-01 00:00".format(startYear),"3600s" #---Date Format: %m-%Y-%d %H:%M
        adc_idx =  pd.date_range('01/01/{} 00:00:00'.format(startYear),'01/01/{} 00:00:00'.format(EndYear),freq='H')
        
        
        #----------------------definitation of stations---
        nodesx = {'1':{'{}'.format(StationID):[]}}
        
        noaa_time_step = '6T'
        
        #--------- Add duration of data set in hours
        period = len(adc_idx)-1
        
        ## creating data set
        #l =[]        
        #l.extend(range(1, period+1))
        #  
              
        #---------------------Ping NOAA API for Validation Data,Create NOAA Dataframe
        
        noaa = pd.DataFrame()
        gages = dict()
        no_data = []
        
        first = datetime.datetime.strptime(start,"%m-%Y-%d %H:%M" )
        last =  pd.date_range(first,periods = period, freq=freq)[-1]
         
        
        for n in nodesx:
            for key in nodesx[n]:
                g = int(key)
               
            t0     = first.strftime('%Y%m%d %H:%M')
            t1     = last.strftime('%Y%m%d %H:%M')
            api_params = {'begin_date': t0, 'end_date': t1,
                        'station': g,'product':product,'datum':datum,
                        'units':units,'time_zone':time_zone,'format':fmt,
                        'application':'NOS.COOPS.TAC.WL' }
                
            pred=[];obsv=[];t=[]
    
            try:
                r = requests.get(url, params = api_params)
                jdata =r.json()
                if 'error' in list(jdata.keys()):
                    print(f"ERROR in Year {startYear}: {jdata['error'].plit('.')[0]}")
                    continue
            
                for j in jdata['data']:
                    t.append(str(j['t']))
                    obsv.append(str(j['v']))
                    pred.append(str(j['s']))
                colname = str(g)    
                noaa[colname]= obsv
                noaa['time']= t
                noaa[colname] = noaa[colname].astype(float)
                gages[jdata['metadata']['id']]=jdata['metadata']['name']
            except:
                #print(g,'No Data')
                no_data.append(str(g)) 
             
        if 'error' in list(jdata.keys()):
            print('Skipping extracting data for year:',startYear)
            startYear=startYear+1        
            continue
        else:
            #assigning index
            noaa = noaa.set_index(noaa['time'])
            
            noaa.index = pd.to_datetime(noaa.index)
            
            # drop extra time axix
            noaa = noaa.drop('time', axis=1) 
            
            #
            #noaa = noaa[colname].astype(float)
            # appending Data to a single Data set
            #noaa[noaa[colname] == ""] = 0
            empty = empty.append(noaa, ignore_index=False)
            
            #itearting in the code for start and end year
            print('Processed data for year:',startYear)
            startYear=startYear+1
    
    empty[str(StationID)].replace('', np.nan, inplace=True)
    
    # =============================================================================
    # plotting
    # =============================================================================
    
    #plotting dataset
    #    
    empty1 = empty.astype(float) 
    
    plt.figure()
    plt.plot(empty1.index, empty1[str(StationID)],c='b',label = product)
    plt.xlabel(f'DateTime ({(time_zone).upper()})')
    plt.ylabel(f'Water Levels ({units}) [above {datum}]')
    plt.title(f"Water Levels for :{jdata['metadata']['name']} [{StationID}]\nLat:{jdata['metadata']['lat']}, Lon: {jdata['metadata']['lon']}")
    
    
    
    
    
    dict_aep = dict(zip(rps, aep_erdc[StationID]))   
    clrs = ['k','green','orange','r','m']
    dict_clrs = dict(zip(rps, clrs)) 
    for ss in dict_aep.keys():
        AEP = dict_aep[ss]
        plt.axhline(AEP,color = dict_clrs[ss],label = ss)
    plt.xlim(empty1.index[0],empty1.index[-1])
    
    plt.legend()
    
    
    # df_store = pd.DataFrame(empty1)   
    
    
    # with open('PeakStagesNOAA_{}.csv'.format(g),'w') as f:
    #    df_store.to_csv(f, index=True,header = True)
        
     
    #
    #


















#
#df = pd.DataFrame(l)
#
#
#df = df.set_index(adc_idx)
#
#
#df = df.merge(noaa,left_index=True, right_index=True)
#
#
#df.drop([0],axis=1,inplace=True)

#
# first data set




# merging dataframes 

#
#bla = [df_final , df]
#
#
#df_final1 = pd.concat(bla)
#
#df_final1 = df_final1.sort_index(axis=0, level=None, ascending=True)
#
#
#
#df_final1.plot()
#
#
#
#
#
#    
#df_store = df_store.merge(df,left_index=True, right_index=True)