# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 12:05:12 2020

@author: Arslaan Khalid
"""


#==============================================================================
# GFS grib data downloader and reader 
#==============================================================================
 
import sys,getopt
import pandas as pd
import numpy.core.multiarray 
import urllib
from urllib.request import urlretrieve # this is how it works for python3
#import geojsoncontour
import numpy as np
from datetime import datetime, timezone, timedelta
import subprocess
import os
import shutil
import netCDF4
import matplotlib.pyplot as plt


print('\nWelcome to meterological pre-processor for iFLOOD (version2)\n.\n.\n')
print('\nGFS 0.5 degree historical download\n')


#==================================================
#configin = '/home/admin/Work/Operational/iFLOOD/config.json'
forecasted_exports_dir='/home/fhrl/Documents/iFLOOD_Brazil/Scripts/historical'
#forecasted_exports_dir=r'C:\Users\Arslaan Khalid\Desktop\Papers_by_Arslaan\Paper10-iFLOODBrazil\iFLOOD_Brazil\Scripts\historical'
#
#
#with open(configin) as json_file:  
#    confin = json.load(json_file)
timestamp1 = input("Storm start date (YYYYMMDDHH)? ")
#timestamp1 = "2016102000"
timestamp2 = input("Storm end date (YYYYMMDDHH)? ")
#timestamp2 = "2016102200"

eventname = input("Name of storm?")
#eventname="storm1"


#==================================================

rootdir='{}'.format(forecasted_exports_dir)


start_date = datetime.strptime(timestamp1,'%Y%m%d%H')
end_date = datetime.strptime(timestamp2,'%Y%m%d%H')

len1 = ((end_date-start_date).days)*4

mydir = os.path.join(rootdir, eventname)


if os.path.isdir(mydir):
    shutil.rmtree(mydir) 
os.makedirs(mydir)

#--- getting the directory of downloaded files

os.chdir(mydir)
os.popen('cp /home/fhrl/Documents/iFLOOD_Brazil/Scripts/wgrib2 .')


print('\n\n______________________________\nStarting data download for:\n'.format(eventname))

for i in range(0,(len1+1),1): 
    #i=0
    file_number = '%03d'%i
    
    temp_date = start_date + timedelta(hours=6*i)
    temp_date0 = datetime.strftime(temp_date,'%Y%m%d%H')
    print('forecasted date {} | {}/{}'.format(temp_date,i,len1))

    url = 'https://nomads.ncdc.noaa.gov/data/gfsanl/{}/{}/gfsanl_4_{}_0000_000.grb2'.format(temp_date0[:6],temp_date0[:8],temp_date0[:8]) 
    output_file = '{}/{}.grb2'.format(mydir,i)
    
    urllib.request.urlretrieve(url,output_file )
    

print('\n\n___________\n')


#----------------- using the Wgrib2 to extract only the required datasets----
#cycletype = input("Model Cycle Type(Nowcast/Forecast/Hindcast)? ")

#if usrin =="Y" or usrin=="y":
#    cycletype = input("Model Cycle Type(Nowcast/Forecast/Hindcast)? ")
#else:
#    cycletype ='Forecast'

cycletype ='Forecast'


#if cycletype=='Nowcast':
# command="./wgrib2 0.grb2 -s | egrep '(UGRD:10 m above ground|VGRD:10 m above ground|PRMSL:mean sea level)' | ./wgrib2 -i 0.grb2 -netcdf uvpNowcast.nc"
# #os.popen('/usr/bin/cdo -s remapbil,LCC2GCS_latlongrid 0.grb2 nowcast.nc')
# p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)    
# (output, err) = p.communicate()
# #This makes the wait possible
# p_status = p.wait()    
# print('\nExtracting u10m,v10m and psl for nowcast file')
 
 
if cycletype=='Forecast': 
 for i in range(0,(len1+1),1): 
  file_number = '%0d'%i
  command="./wgrib2 {}.grb2 -s | egrep '(UGRD:10 m above ground|VGRD:10 m above ground|PRMSL:mean sea level)' | ./wgrib2 -i {}.grb2 -netcdf uvpForecast{}.nc".format(file_number,file_number,file_number)
  #os.popen('/usr/bin/cdo -s remapbil,LCC2GCS_latlongrid 0.grb2 nowcast.nc')
  p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)    
  (output, err) = p.communicate()
  #This makes the wait possible
  p_status = p.wait()
  print('Extracting u10m,v10m and psl for file: {}'.format(file_number))

#if cycletype=='Hindcast': 
# for i in range(0,37,3): 
#  file_number = '%0d'%i
#  command="./wgrib2 {}.grb2 -s | egrep '(UGRD:10 m above ground|VGRD:10 m above ground|PRMSL:mean sea level)' | ./wgrib2 -i {}.grb2 -netcdf uvpForecast{}.nc".format(file_number,file_number,file_number)
#  #os.popen('/usr/bin/cdo -s remapbil,LCC2GCS_latlongrid 0.grb2 nowcast.nc')
#  p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)    
#  (output, err) = p.communicate()
#  #This makes the wait possible
#  p_status = p.wait()
#  print('Extracting u10m,v10m and psl for file: {}'.format(file_number))




# =============================================================================
#  convert the data to adcirc NWS 6
# =============================================================================
    
for e in range (0,(len1+1),1):
 #e=0
 tmp_url = '{}/uvpForecast{}.nc'.format(mydir,e)
 #tmp_url =r'C:\Users\Arslaan Khalid\Desktop\Papers_by_Arslaan\Paper10-iFLOODBrazil\iFLOOD_Brazil\Scripts\historical\storm1\uvp1.nc'
 file1 = netCDF4.Dataset(tmp_url)
 #===== subseting the data ===========================
 lat  = file1.variables['latitude'][:]
    
 lon  = file1.variables['longitude'][:]

 u_data  = file1.variables['UGRD_10maboveground'][:,:,:]
    
 v_data  = file1.variables['VGRD_10maboveground'][:,:,:]
    
 P_data  = file1.variables['PRMSL_meansealevel'][:,:,:]
    
 #Time_data = netCDF4.num2date(time_var_re[res:],time_var_re.units)

 wind_mag_re = np.sqrt(u_data**2 + v_data**2)  

 #lat1 = lat[300:650] # for US
 #lon1 = lon[1100:] # for US

 lon1 = lon[900//2:1500//2] #brazil
 lat1 = lat[60//2:370//2] #brazil
    
    
 grid_space_x = round((lon[1]-lon[0]),4)
 grid_space_y = round((lat[1]-lat[0]),4)

 if e==0:print("Grid resolution is {} degree".format(grid_space_x))
    
 file_number = '%02d'%e
 # =============================================================================
 #     flipping the data
 # =============================================================================
 #========================================== u wind 10m ====================================
 x1,x2,y1,y2=60//2,370//2,900//2,1500//2 # the 2 is based on resolution, this is divisible of 1 degree resolution
 
 domain1_u = u_data[:,x1:x2,y1:y2]
 #domain1_u = u_data[:,60:370,900:1500]
 domain1_u_flip = np.flipud(domain1_u[0])

 domain1_v = v_data[:,x1:x2,y1:y2]
 domain1_v_flip = np.flipud(domain1_v[0]) 
 #    plt.imshow(domain1_u, cmap='coolwarm')
 
    
 domain1_p = P_data[:,x1:x2,y1:y2]
 domain1_p_flip = np.flipud(domain1_p[0])  

    
 domain1_m = wind_mag_re[:,x1:x2,y1:y2]
 domain1_m_flip = np.flipud(domain1_m[0])

 #plt.imshow(domain1_m_flip, cmap='jet')    
 #==============================================================
 container_u = []
    
 domain11_u = domain1_u_flip.data

 #code this number to read from len(lat1)
 lat10 = len(lat1)
 lon10 = len(lon1)


 for i in range(0,lat10):
  for j in range(0,lon10):
   #print(round(domain11_u[i][j],1))
   key1 = round(domain11_u[i][j],1)
   container_u.append(key1)
 #============================================================    
 container_v = []
    
 domain11_v = domain1_v_flip.data
    
 for i in range(0,lat10):
  for j in range(0,lon10):
   #print(round(domain11_v[i][j],1))
   key2 = round(domain11_v[i][j],1)
   container_v.append(key2)
 #=============================================================
 container_p = []

 domain11_p = domain1_p_flip.data
    
 for i in range(0,lat10):
  for j in range(0,lon10):
   #print(round(domain11_p[i][j],0))
   key3 = round(domain11_p[i][j],0)
   container_p.append(key3)
 #=============================================================
 winds_tmp = pd.DataFrame()
    
 winds_tmp['u']=container_u
 winds_tmp['v']=container_v
 winds_tmp['p']=container_p
    
 winds_tmp1 =  winds_tmp.round(1)
 winds_tmp1 = winds_tmp1.replace(9.999000260554009e+20,np.nan)
 winds_tmp1[winds_tmp1 > 201325] = np.nan
 winds_tmp1[winds_tmp1 < -99] = np.nan
 #defaults
 #winds_tmp1 = winds_tmp1.replace(9.999000260554009e+20,0)
 # have conditions here to buffer down the winds and pressure in the no data areas
    
 winds_tmp1['u'] = winds_tmp1['u'].replace(np.nan,(np.nanmean(winds_tmp1['u'])))
 winds_tmp1['v'] = winds_tmp1['v'].replace(np.nan,(np.nanmean(winds_tmp1['v'])))
 winds_tmp1['p'] = winds_tmp1['p'].replace(np.nan,101325)        
 #------------------
 print('processing file:',e)
 #print('processed ({}/84)'.format(e)

 # writing the data to file
 if e==0:                        
  with open('{}/{}_fort.22'.format(forecasted_exports_dir,eventname),'w') as f:
   winds_tmp1.to_csv(f, sep='\t',index=False,header = False)
 else:
  with open('{}/{}_fort.22'.format(forecasted_exports_dir,eventname),'a') as f:
   winds_tmp1.to_csv(f, sep='\t',index=False,header = False)


print('Finished writing NWS 6 file')

ctrl_file = "following is the fort22 domain extents\n{} {} {} {} {} {} 10800    ! NWLAT, NWLON, WLATMAX, WLONMIN, WLATINC, WLONINC, WTIMINC".format(lat10,lon10,round(lat1.max(),4),(round(lon1.min(),4)-360),round(grid_space_x,4),round(grid_space_y,4))

with open('{}/{}_ctrl_file_domain_desc.txt'.format(forecasted_exports_dir,eventname),'w') as f:
           f.write(ctrl_file) 
