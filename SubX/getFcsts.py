#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import numpy as np
import pandas as pd

import os
from datetime import datetime, timedelta, date
import time
import argparse

import matplotlib.pyplot as plt
import proplot as pplt

import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
import cartopy.feature as cfeature
import requests
from subx_utils import *

# Eliminate Warnings
import warnings
warnings.filterwarnings("ignore")

# Set xarray to keep attributes
xr.set_options(keep_attrs=True)

# Get Start Time
start_time = time.time()

# Parse commend line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--date",nargs='?',default=None,help="make subx forecasts based on this date")
args = parser.parse_args()

# ### Models and Forecast Settings

#_,subxclimos_list,_,_=initSubxModels()
#a,b,c,d,e=initSubxModels()

subxclimos_list,_,_,_,_=initSubxModels()
model_labels=[item['group']+'-'+item['model'] for item in subxclimos_list]
nweeks=4
interactive_vars=['uas','vas','psl']


#print(initSubxModels())
#print(initSubxModels())

#model_labels=[item['group']+'-'+item['model'] for item in subxclimos_list]
#nweeks=4


if (args.date):
    fcstdate,fcst_week_dates=getFcstDates(date=args.date)
else:
    fcstdate,fcst_week_dates=getFcstDates()


# ### File paths
timeout_seconds = 60
url='http://iridl.ldeo.columbia.edu/SOURCES/.Models/.SubX/'
datatype='forecast'
#hcstPath='/shared/subx/hindcast/'
hcstPath='/mcs/scratch/kpegion/subx/hindcast/'
hcstPath='/home/admin/Work/Operational/SubX/Pre-Process/Operational2023_subx/'
#outPath='/mcs/scratch/kpegion/subx/figs_test/'
outPath='/home/admin/Work/Operational/SubX/Pre-Process/Operational2023_subx/'

#print(subxclimos_list)
print('PROCESSING Forecasts: ')

# Loop over all the SubX Models
ds_models_list=[]

not_downloaded = []

for imodel,subx_model in enumerate(subxclimos_list):
    #print(str(date.toordinal(fcstdate)),fcstdate)

    # Get the model, group, variables, and levels from the dictionary 
    varnames=subx_model['varnames']
    plevstrs=subx_model['plevstrs']
    model=subx_model['model']
    group=subx_model['group']


    #print(model)
    if model in []: # skipping some of the models 'GEOS_V2p1'
        continue
    else:

        #print(plevstrs)
        
        print('===> '+group+'-'+model)
    
        
        # Loop over variables for this model
        ds_anoms_list=[]
        for varname,plevstr in zip(varnames,plevstrs):
            
            # Read Data
            baseURL=url+'.'+group+'/.'+model+'/.'+datatype+'/.'+varname
            inFname=baseURL+'/dods'  
            print(f'DataLink: {inFname}')
            try:
                
                # Use requests to open the dataset with a timeout
                response = requests.get(inFname, timeout=timeout_seconds)
                response.raise_for_status()  # Raise an HTTPError for bad responses
                
                # Open the dataset from the response conten
                ds_tmp=xr.open_dataset(inFname)
        
            	# getting model initializing dates
                startdates=ds_tmp['S'].values
        
        
                # Set files & paths for writing
                fcst_path=hcstPath+varname+plevstr+'/fcst/'+group+'-'+model
                print()
                #print(fcst_path)
                if not os.path.exists(fcst_path):
                   print("MAKING model forecast directory: ",fcst_path)
                   os.makedirs(fcst_path)
        
            	# restructing and slicing the file before downloading
                ds_tmp['S']=ds_tmp['S'].dt.floor('d')
        
                ds_tmp=ds_tmp.rename({'X':'lon','Y':'lat','L':'time'})  
                
                ds_tmp['time']=np.arange(len(ds_tmp['time']))
                ds_tmp['lon'].attrs['units']='degrees_east'
                ds_tmp['lat'].attrs['units']='degrees_north'
                ds_tmp=ds_tmp.squeeze(drop=True)
                
                # slicing the extent
                ds_forecast_recent=ds_tmp.sel(S=str(startdates[-1])[:10],lon=slice(245,330),lat=slice(-5,60)) 
                # reversing lats
                ds_forecast_recent=ds_forecast_recent.reindex(lat=list(reversed(ds_forecast_recent['lat'])))
        
                try:
                    if 'M' in ds_tmp.dims:
                            
                        for mms in range(1,len(ds_tmp[varname].M)+1):
                            ds_each_ens=ds_forecast_recent.sel(M=mms)
                            min_va = ds_each_ens[varname].min().data.tolist()
                            #print(min_va)
                            if np.isnan(min_va) == False:
                             f2 = varname+'_'+group+'-'+model+f'_{str(startdates[-1])[:10]}_e{mms}.nc'
                             fcst_fname=fcst_path+'/'+varname+'_'+group+'-'+model+f'_{str(startdates[-1])[:10]}_e{mms}.nc'
                             print("Downloading file:",f2)
                             ds_each_ens.to_netcdf(fcst_fname)
                    else:
                        ds_each_ens=ds_forecast_recent#.sel(M=mms)
                        min_va = ds_each_ens[varname].min().data.tolist()
                        #print(min_va)
                        if np.isnan(min_va) == False:
                         f2 = varname+'_'+group+'-'+model+f'_{str(startdates[-1])[:10]}_e{mms}.nc'
                         fcst_fname=fcst_path+'/'+varname+'_'+group+'-'+model+f'_{str(startdates[-1])[:10]}_e{mms}.nc'
                         print("Downloading file:",f2)
                         ds_each_ens.to_netcdf(fcst_fname)
    
                except:
                    f2 = varname+'_'+group+'-'+model+f'_{str(startdates[-1])[:10]}_e{mms+1}.nc'
                    print(f'[Warning]: Data not downloaded:',f2)
                    not_downloaded.append(f2)
            
            except requests.exceptions.Timeout:
                # Handle timeout exception
                print("Timeout occurred while fetching the dataset.")
            except requests.exceptions.RequestException as e:
                # Handle other request exceptions (e.g., network errors)
                print(f"Error: {e}")
            except Exception as e:
                # Handle other exceptions
                print(f"Error: {e}")        
                print(f'[Warning]: Online file not downloadable')  
 

print()
print("TIME INFO: --- %s seconds ---" % (time.time() - start_time))

