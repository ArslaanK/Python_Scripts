# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 19:24:37 2020

@author: Arslaan Khalid
"""

import cdsapi
import pandas as pd
import datetime


file='/home/admin/Work/Programs/FirstStreet_datadownload/FirstStreet/ERA5_reanalysis/noreaster.txt'
Noreasters = pd.read_csv(file,sep='\t')

#Noreasters['start'] = Noreasters['byr'].astype(str)+ (Noreasters['bmt']).astype(str).apply(lambda x: x.zfill(2))+(Noreasters['bdy']).astype(str).apply(lambda x: x.zfill(2))
#Noreasters['end'] = Noreasters['eyr'].astype(str)+ (Noreasters['emt']).astype(str).apply(lambda x: x.zfill(2))+(Noreasters['edy']).astype(str).apply(lambda x: x.zfill(2))


        
for k in range(0,len(Noreasters)):#

    info1=Noreasters['nbmt'].loc[k]
    
    #---- extracting the year
    year1 = [str(Noreasters['byr'].loc[k])]    
    #---- extracting the month list
    month1=[]
    if Noreasters['nbmt'].loc[k] == Noreasters['emt'].loc[k]:
        tmp=Noreasters['nbmt'].loc[k]
        tmp = f"{tmp:02}"
        month1.append(tmp)
    else:
        tmp1,tmp2=Noreasters['nbmt'].loc[k],Noreasters['emt'].loc[k];
        tmp1,tmp2= f"{tmp1:02}" ,  f"{tmp2:02}"
        month1=[tmp1,tmp2]   
    
    #---- extracting the month list
    
    
    start = datetime.datetime(Noreasters['byr'].loc[k] , Noreasters['nbmt'].loc[k] , Noreasters['nbdy'].loc[k] )
    end = datetime.datetime(Noreasters['eyr'].loc[k] , Noreasters['emt'].loc[k] , Noreasters['edy'].loc[k])
    
    delta = end - start
    
    print(f'start:{start} to end:{end} is days:{delta.days}')
    
        
    day1=[]
    for i in range(delta.days + 1):
        #print (start + datetime.timedelta(days=i))
        dy = str(start + datetime.timedelta(days=i))
        dy1 = dy[8:10]
        day1.append(dy1)
    print(year1,month1,day1)
          
    #------------ use the extracted data to as server for data
    c = cdsapi.Client()
    
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type':'reanalysis',
            'format':'netcdf',
            'variable':[
                '10m_u_component_of_wind','10m_v_component_of_wind','mean_sea_level_pressure'
            ],
            'year':year1,
            'month':month1,
            'day':day1,
            'time':[
                '00:00', '06:00', '12:00',  '18:00' ]
        },
        f'/home/admin/Work/Programs/FirstStreet_datadownload/FirstStreet/ERA5_reanalysis/downloads/Noreaster{k}.nc')