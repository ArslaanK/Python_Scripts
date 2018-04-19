# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 11:46:44 2017

@author: Arslaan khalid

"""

# for queries please contact @ akhalid6@gmu.edu


# following script downloads the latest point based forecast coming out of ETSS system
# ETSS forecast from the following code is biased corrected



import matplotlib.pyplot as plt
import fileinput
from datetime import datetime, timedelta
import pandas as pd
import os
import requests
from scipy.interpolate import interp1d
import html2text
import re


url = 'http://tgftp.nws.noaa.gov/data/raw/sr/srus70.kwno.tid.twe.txt'
res = requests.get(url)
text = res.text

texto = re.sub('/ ','\n',text)

# data formatting


for i in range (2,10):
    texto = texto.replace('.E1','\n  m\n')    
    texto = texto.replace('.E{}'.format(i),'\n')
    
texto = texto.replace('/.E','  m\n')
texto = texto.replace('\n','')
texto = texto.replace('     ','  ')
texto = texto.replace('  ',' ')
texto = texto.replace('/.E',' m ')

n = texto.strip().split(' ')[0:]

# removing blank spaces in the list

n = list(filter(None, n))

df = pd.DataFrame(n)

# stripping time

Date = texto.strip().split()[21]

day = texto.strip().split()[21][6:]
month = texto.strip().split()[21][4:6]
year = texto.strip().split()[21][0:4]

hour = texto.strip().split('/')[23][2:4]
minute = texto.strip().split('/')[23][4:]


datetime_start = month + '/'+ day + '/'+ year +' '+hour+':'+minute


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

SLIM_idx = list(range(indices[0],indices[1]))

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

WASD_idx = list(range(indices[42],indices[43]))

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

#WASD.plot()
#WASD = pd.to_numeric(WASD[0],errors='coerce')



#-------------------------KPTV2-----------------   
kp = []

KPTV_idx = list(range(indices[88],indices[89]))

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

CBBV_idx = list(range(indices[92],indices[93]))

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

LWTV_idx = list(range(indices[96],indices[97]))

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

#LWTV = pd.to_numeric(LWTV[0],errors='coerce')

#-------------------------BLTM-----------------   
bl = []

BLTM_idx = list(range(indices[106],indices[107]))

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


df_final = df_final.merge(CBBV,left_index=True, right_index=True).merge(LWTV,left_index=True, right_index=True).merge(BLTM,left_index=True, right_index=True)

#BLTM = pd.to_numeric(BLTM[0],errors='coerce')

#-------------------------APAM-----------------   
ap = []

APAM_idx = list(range(indices[128],indices[129]))

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

#APAM = pd.to_numeric(APAM[0],errors='coerce')

#-------------------------BISM2-----------------   
bm = []

BISM_idx = list(range(indices[134],indices[135]))

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

#-------------------------AXTV2-----------------   
ax = []

AXTV_idx = list(range(indices[174],indices[175]))

for s in AXTV_idx:
    
    w = df.iloc[s]
    ax.append(w)

AXTV = pd.DataFrame(ax,dtype=float)

AXTV = AXTV.drop(AXTV.index[0])

idx = []
idx = list(range(0,len(AXTV)))

AXTV.index = idx

AXTV = pd.DataFrame(AXTV,dtype=float)
AXTV = AXTV.rename(columns = {0:'AXTV2'})


#AXTV = pd.to_numeric(AXTV[0],errors='coerce')

#-------------------------NCDV2-----------------    
nc = []

NCDV_idx = list(range(indices[184],indices[185]))

for s in NCDV_idx:
    
    w = df.iloc[s]
    nc.append(w)

NCDV = pd.DataFrame(nc,dtype=float)

NCDV = NCDV.drop(NCDV.index[0])

idx = []
idx = list(range(0,len(NCDV)))

NCDV.index = idx

NCDV = pd.DataFrame(NCDV,dtype=float)
NCDV = NCDV.rename(columns = {0:'NCDV2'})




#-------------------------SWPV-----------------      
sw = []

SWPV_idx = list(range(indices[56],indices[57]))

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


df_final = df_final.merge(AXTV,left_index=True, right_index=True).merge(NCDV,left_index=True, right_index=True).merge(SWPV,left_index=True, right_index=True)



freq = "3600s" #---Date Format: %m-%Y-%d %H:%M




df_idx = pd.date_range(datetime,periods = len(df_final), freq='H') 
df_final = df_final.set_index(df_idx)   

df_final.index.names = ['Datetime(UTC)']

#df_final.plot()

with open('ETSS.csv','w') as f:
    df_final.to_csv(f, index=True,header = True)




