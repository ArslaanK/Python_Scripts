# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 17:24:14 2018

@author: Arslaan Khalid
"""


# for queries please contact @ akhalid6@gmu.edu


# following script downloads the latest point based forecast coming out of CBOFS system
# ETSS forecast from the following code is un-biased


import sys,getopt
#from netCDF4 import Dataset  # use scipy instead
from scipy.io import netcdf #### <--- This is the library to import.
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from datetime import datetime, timedelta
import pandas as pd
import wget
from datetime import datetime



# datetime

fort22 = '/home/fhrl/Documents/Surge/input/Domain2/fort22/fort.221' #for Domain2
#fort22 = '/home/fhrl/Documents/AHPS_Forced_ADCIRC/INPUT/ModFEMAMeshVersion0/RUN/fort.221' #for Domain2



'''
with open(fort22,'r') as f:
    for i in range(0,1):
        line = f.readline().strip().split()
        timestamp = line[3]
        start_date = datetime.strptime(timestamp,'%Y%m%d%H')
'''

timestamp = '2018041806'                          # change the following date to your required forecast

start_date = datetime.strptime(timestamp,'%Y%m%d%H')

        
date = timestamp.strip().split()[0][0:8]
hour = timestamp.strip().split()[0][8:]

# file dowloading
print('Beginning file download with wget module')

url = 'https://opendap.co-ops.nos.noaa.gov/netcdf/cbofs/{}/nos.cbofs.stations.forecast.{}.t{}z.nc'.format(timestamp[:6],date,hour)  
output_file = '/home/fhrl/Documents/Surge/scripts/f1.nc'                 # location to save netcdf file
wget.download(url,output_file )


# Open file in a netCDF reader

file_name = output_file

nc = netcdf.netcdf_file(file_name,'r')

#Look at the variables available
nc.variables

#Look at the dimensions
nc.dimensions


#Look at a specific variable's dimensions
nc.variables['zeta'].dimensions   ## output is ('Time', 'south_north', 'west_east')

#Look at a specific variable's units
nc.variables['zeta'].units        ## output is ('K')


print(nc.variables['zeta'])

test1 = nc.variables['zeta'][:,:]

#stations_I = nc.variables['Ipos'][0:166]
#stations_J = nc.variables['Jpos'][0:166]
stations_lat = nc.variables['lat_rho'][0:166]
stations_lon = nc.variables['lon_rho'][0:166]

test = pd.DataFrame(test1)

# converting to meters  () default datum in NAVD 88

test = test*3.281

#test.columns

#for i in range(0,166):
 #   test[i].plot()
 
 #print(df.loc[df['A'] == 'foo'])


###++++++++++++++++++++++++++++++ STATIONWISE DATA EXTRACTION and Adjustment to timezone
# CB is 54
CB = test[54]+(0.566*3.281)

# KP is 47
KP = test[47]+(0.575*3.281)

# SW is 59
SW = test[59]+(0.450*3.281)

# BI is 29
BI = test[29]+(0.300*3.281)

# AP is 16
AP = test[16]+(0.231*3.281)

# BL is 8
BL = test[8]+(0.254*3.281)

# SL is 87
SL = test[87]+(0.262*3.281)

# LW is 32
LW = test[32]+(0.303*3.281)

# DA is 26
DA = test[26]+(0.356*3.281)

# AL is 160
AL = test[160]+(0.403*3.281)

#=========================================

df = pd.DataFrame()
df['CBBV2'] = CB
df['SWPV2'] = SW
df['BLTM2'] = BL
df['BISM2'] = BI
df['SLIM2'] = SL
df['APAM2'] = AP
df['KPTV2'] = KP
df['LWTV2'] = LW



#df.plot()


with open('/home/fhrl/Documents/Surge/scripts/CBOFS.csv','w') as f:
    df.to_csv(f, sep='\t',index=False,header = True)

