# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 11:46:44 2017

@author: Arslaan khalid
"""
# for queries please contact @ akhalid6@gmu.edu


# following script downloads the latest point based forecast coming out of ESTOFS system
# ETSS forecast from the following code is un-biased

import matplotlib.pyplot as plt
import fileinput
from datetime import datetime, timedelta
import pandas as pd
import os
import requests
from scipy.interpolate import interp1d
import html2text
import re
import numpy as np

# file fort.221 is read before intializing the code to give the starting of the forecast data

fort22 = '/home/fhrl/Documents/Surge/input/Domain2/fort22/fort.221'



'''
with open(fort22,'r') as f:
    for i in range(0,1):
        line = f.readline().strip().split()
        timestamp = line[3]
        start_date = datetime.strptime(timestamp,'%Y%m%d%H')
'''

from datetime import datetime

timestamp = '2018041806'                          # change the following date to your required forecast


start_date = datetime.strptime(timestamp,'%Y%m%d%H')
        

        
date = timestamp.strip().split()[0][0:8]
hour = timestamp.strip().split()[0][8:]

#with open('df.txt') as f:
    #text = f.read()

#----------------=================== ESTOFS every 30 min _ Datum is MLLW___________________=================

url = 'http://nomads.ncep.noaa.gov/pub/data/nccf/com/estofs/prod/estofs_atl.{}/estofs.atl.t{}z.points.cwl.shef'.format(date,hour)

res = requests.get(url)
text = res.text


# data formatting


textm = text.replace('.E1','/')
textm = textm.replace('.E2','/')
textm = textm.replace('.E','/')


for i in range (1,10):
    textm = textm.replace('/{}'.format(i),'/')


textm = textm.replace('/0','/')
textm = textm.replace('/DIN30/','/DIN30/m/')
textm = textm.replace('SX','/m/SX')
textm = textm.replace('-9999',' ')


# renaming stations for indexing

textm = textm.replace('SLIM2','m/SLIM2/')
textm = textm.replace('WASD2','m/WASD2/')
textm = textm.replace('KPTV2','m/KPTV2/')
textm = textm.replace('CBBV2','m/CBBV2/')
textm = textm.replace('LWTV2','m/LWTV2/')
textm = textm.replace('BLTM2','m/BLTM2/')
textm = textm.replace('APAM2','m/APAM2/')
textm = textm.replace('BISM2','m/BISM2/')
textm = textm.replace('SWPV2','m/SWPV2/')

#---- identifing the index of stations

#df.loc[df[0] == 'SLIM2']

# Date splitting


n = textm.strip().split('/')[0:]



# removing blank spaces in the list

n = list(filter(None, n))

df = pd.DataFrame(n)



# stripping time
if df.iloc[5].str.contains('yuji').any():
    r =41
    #print('true')
else:
    r=39
    #print('false')
# stripping time

# default is 41 if fujo info not ion estofs output we add +1    

Date = text.strip().split()[r]

day = text.strip().split()[r][6:]
month = text.strip().split()[r][4:6]
year = text.strip().split()[r][0:4]

H = text.strip().split()[r+2][2:4]
Min = text.strip().split()[r+2][4:6]


datetime_start = month + '/'+ day + '/'+ year +' '+H+':'+Min


datetime = datetime.strptime(datetime_start,'%m/%d/%Y %H:%M')


# spliting datas for each station

#--------------==========----Extracting idices -----===========--------
indices = []
b = []

for x in range (0,len(n)):
    if n[x]=='m': 
        b = x       
        indices.append(b)
    
#------=====================----Station-Wise ETSS Datasets-----========


#-------------------------SLIM2-----------------    
ss = []

SLIM_idx = list(range(indices[41],indices[41]+362))

for s in SLIM_idx:
    
    f = df.iloc[s]
    ss.append(f)

SLIM = pd.DataFrame(ss,dtype=float)

SLIM = SLIM.drop(SLIM.index[0])

idx = []
idx = list(range(0,len(SLIM)))

SLIM.index = idx

SLIM = pd.DataFrame(SLIM,dtype=float)
SLIM = SLIM.rename(columns = {0:'SLIM2'})

#SLIM = pd.to_numeric(SLIM[0],errors='coerce')


#-------------------------WASD2-----------------   
wa = []

WASD_idx = list(range(indices[42],indices[42]+362))

for s in WASD_idx:
    
    w = df.iloc[s]
    wa.append(w)

WASD = pd.DataFrame(wa,dtype=float)

WASD = WASD.drop(WASD.index[0])

idx = []
idx = list(range(0,len(WASD)))

WASD.index = idx

WASD = pd.DataFrame(WASD,dtype=float)
WASD = WASD.rename(columns = {0:'WASD2'})


#WASD = pd.to_numeric(WASD[0],errors='coerce')




#-------------------------KPTV2-----------------   
kp = []

KPTV_idx = list(range(indices[43],indices[43]+362))

for s in KPTV_idx:
    
    w = df.iloc[s]
    kp.append(w)

KPTV = pd.DataFrame(kp,dtype=float)

KPTV = KPTV.drop(KPTV.index[0])

idx = []
idx = list(range(0,len(KPTV)))

KPTV.index = idx

KPTV = pd.DataFrame(KPTV,dtype=float)
KPTV = KPTV.rename(columns = {0:'KPTV2'})

#KPTV.plot()
#KPTV = pd.to_numeric(KPTV[0],errors='coerce')
df_final = SLIM.merge(WASD,left_index=True, right_index=True).merge(KPTV,left_index=True, right_index=True)

#-------------------------CBBV2-----------------   
cb = []

CBBV_idx = list(range(indices[49],indices[49]+362))

for s in CBBV_idx:
    
    w = df.iloc[s]
    cb.append(w)

CBBV = pd.DataFrame(cb,dtype=float)

CBBV = CBBV.drop(CBBV.index[0])

idx = []
idx = list(range(0,len(CBBV)))

CBBV.index = idx

CBBV = pd.DataFrame(CBBV,dtype=float)
CBBV = CBBV.rename(columns = {0:'CBBV2'})


#CBBV = pd.to_numeric(CBBV[0],errors='coerce')

#-------------------------LWTV-----------------   
lw = []

LWTV_idx = list(range(indices[45],indices[45]+362))

for s in LWTV_idx:
    
    w = df.iloc[s]
    lw.append(w)

LWTV = pd.DataFrame(lw,dtype=float)

LWTV = LWTV.drop(LWTV.index[0])

idx = []
idx = list(range(0,len(LWTV)))

LWTV.index = idx

LWTV = pd.DataFrame(LWTV,dtype=float)
LWTV = LWTV.rename(columns = {0:'LWTV2'})

#LWTV.tail()
#LWTV = pd.to_numeric(LWTV[0],errors='coerce')

#-------------------------BLTM-----------------   
bl = []

BLTM_idx = list(range(indices[39],indices[39]+362))

for s in BLTM_idx:
    
    w = df.iloc[s]
    bl.append(w)

BLTM = pd.DataFrame(bl,dtype=float)

BLTM = BLTM.drop(BLTM.index[0])

idx = []
idx = list(range(0,len(BLTM)))

BLTM.index = idx

BLTM = pd.DataFrame(BLTM,dtype=float)
BLTM = BLTM.rename(columns = {0:'BLTM2'})

#BLTM.tail()
df_final = df_final.merge(CBBV,left_index=True, right_index=True).merge(LWTV,left_index=True, right_index=True).merge(BLTM,left_index=True, right_index=True)

#BLTM = pd.to_numeric(BLTM[0],errors='coerce')

#-------------------------APAM-----------------   
ap = []

APAM_idx = list(range(indices[40],indices[40]+362))

for s in APAM_idx:
    
    w = df.iloc[s]
    ap.append(w)

APAM = pd.DataFrame(ap,dtype=float)

APAM = APAM.drop(APAM.index[0])

idx = []
idx = list(range(0,len(APAM)))

APAM.index = idx

APAM = pd.DataFrame(APAM,dtype=float)
APAM = APAM.rename(columns = {0:'APAM2'})

APAM.tail()
#APAM = pd.to_numeric(APAM[0],errors='coerce')

#-------------------------BISM2-----------------   
bm = []

BISM_idx = list(range(indices[35],indices[35]+362))

for s in BISM_idx:
    
    w = df.iloc[s]
    bm.append(w)

BISM = pd.DataFrame(bm,dtype=float)

BISM = BISM.drop(BISM.index[0])

idx = []
idx = list(range(0,len(BISM)))

BISM.index = idx

BISM = pd.DataFrame(BISM,dtype=float)
BISM = BISM.rename(columns = {0:'BISM2'})

df_final = df_final.merge(APAM,left_index=True, right_index=True).merge(BISM,left_index=True, right_index=True)

#BISM = pd.to_numeric(BISM[0],errors='coerce')



#-------------------------SWPV-----------------      
sw = []

SWPV_idx = list(range(indices[48],indices[48]+362))

for s in SWPV_idx:
    
    w = df.iloc[s]
    sw.append(w)

SWPV = pd.DataFrame(sw,dtype=float)

SWPV = SWPV.drop(SWPV.index[0])

idx = []
idx = list(range(0,len(SWPV)))

SWPV.index = idx

SWPV = pd.DataFrame(SWPV,dtype=float)
SWPV = SWPV.rename(columns = {0:'SWPV2'})



#===========---- Merging Dataframes ===-----------


df_final = df_final.merge(SWPV,left_index=True, right_index=True)



freq = "3600s" #---Date Format: %m-%Y-%d %H:%M


#BLTM = BLTM.fillna(to_replace=-1e+04,value=0,inplace=True)

df_idx = pd.date_range(datetime,periods = len(df_final), freq='30Min') 
df_final = df_final.set_index(df_idx)   

df_final.index.names = ['Datetime(UTC)']



with open('ESTOFS.csv','w') as f:
    df_final.to_csv(f, index=True,header = True)