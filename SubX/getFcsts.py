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

# adding functions


def initSubxModels():
    
    #all_varnames=['ua','ua','rlut','tas','ts','zg','va','va','pr','zg','uas','vas','psl']
    #all_plevstrs=['850','200','toa','2m','sfc','500','200','850','sfc','200','10m','10m','msl']
    #all_units=['ms-1','ms-1','Wm-2','degC','degC','m','ms-1','ms-1','mmday-1','m','ms-1','ms-1','hPa']
    
    all_varnames=['tas','pr','zg']
    all_plevstrs=['2m','sfc','500']
    all_units=['degC','mmday-1','m']
   
    sub1_varnames=['tas','pr','zg']
    sub1_plevstrs=['2m','sfc','500']
    sub1_units=['degC','mmday-1','m']

    all_varnames=['uas','vas','psl']
    all_plevstrs=['10m','10m','msl']
    all_units=['ms-1','ms-1','hPa']
   
    sub1_varnames=['uas','vas','psl']
    sub1_plevstrs=['10m','10m','msl']
    sub1_units=['ms-1','ms-1','hPa']




    
    ccsm4_dict={'model':'CCSM4','group':'RSMAS','varnames': all_varnames, 'plevstrs': all_plevstrs, 'plot_loc':2}
    geos_dict={'model':'GEOS_V2p1','group':'GMAO','varnames': all_varnames, 'plevstrs': all_plevstrs,'plot_loc':4}
    fim_dict={'model':'FIMr1p1','group':'ESRL','varnames': all_varnames, 'plevstrs': all_plevstrs,'plot_loc':1}
    #geps_dict={'model':'GEPS6','group':'ECCC','varnames': all_varnames, 'plevstrs': all_plevstrs,'plot_loc':6}
    geps_dict={'model':'GEPS7','group':'ECCC','varnames': all_varnames, 'plevstrs': all_plevstrs,'plot_loc':6}
    nrl_dict={'model':'NESM','group':'NRL','varnames': all_varnames, 'plevstrs': all_plevstrs,'plot_loc':5}
    gefs_dict={'model':'GEFSv12_CPC','group':'EMC','varnames': sub1_varnames, 'plevstrs': sub1_plevstrs,'plot_loc':3}
    #gefs_dict={'model':'GEFSv12','group':'EMC','varnames': sub1_varnames, 'plevstrs': sub1_plevstrs,'plot_loc':3}
    cfsv2_dict={'model':'CFSv2','group':'NCEP','varnames': sub1_varnames, 'plevstrs': sub1_plevstrs,'plot_loc':7}

    subxclimos_list=[fim_dict,ccsm4_dict,geos_dict,nrl_dict,geps_dict,gefs_dict,cfsv2_dict]
    subxmodels_list=[fim_dict,ccsm4_dict,geos_dict,nrl_dict,geps_dict,gefs_dict,cfsv2_dict]
   
    return subxmodels_list, subxclimos_list,all_varnames, all_plevstrs, all_units


def getFcstDates(date=None):

    if (date):
        print("Using Input Date: ",date)
        currentdate=datetime.strptime(date,'%Y%m%d')
        print("Using Input Date: ",currentdate)
    else:
        currentdate=datetime.today().replace(microsecond=0,second=0,minute=0,hour=0)
        print("Using Current Date: ",currentdate)

    # How far are we from the most recent thurs?
    diffdate=(currentdate.weekday()-3) % 7

    # fcstdate is the most recent Thurs
    fcstdate=currentdate-timedelta(days=diffdate)

    # Identify the previous Friday (6 days earlier) as start of fcst week
    weekdate=fcstdate - timedelta(days=6)

    print("USING FCSTS ICS FROM: ",weekdate.strftime(('%Y%m%d')), " to ",fcstdate.strftime('%Y%m%d'))
    print("Output files and figures will be labeled as: ",fcstdate.strftime('%Y%m%d'))

    # Get the range of dates over the last week from the start of previous fcst week to fcstdate
    fcst_week_dates=pd.date_range(start=weekdate,end=fcstdate,freq='D') 
    
    return fcstdate,fcst_week_dates


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
hcstPath='/home/admin/Work/Operational/SubX/Pre-Process/Operational2023_subx/'
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

