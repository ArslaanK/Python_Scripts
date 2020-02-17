# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 23:43:27 2019

@author: Arslaan Khalid
"""



#==============================================================================
# ERA5 file read and Write for ADCIRC
#==============================================================================
 

#
import matplotlib.pyplot as plt # Import the Matplotlib package
#from osgeo import gdal # Import the GDAL library
#import sys,getopt
import pandas as pd
import numpy.core.multiarray 
#import wget
#import geojsoncontour
import numpy as np
from datetime import datetime, timezone, timedelta
import netCDF4
import numpy.ma as ma
import time

#=================

print('\n======================\nStarting Preprocessor for NWS6 winds from raw ERA5\n======================\nAuthor: ArslaanKhalid\nakhalid6@gmu.edu\n\n\n')
# storm name and date
stormname="joaquin"
date="2015091500"
durationrun = 30 #days
interval = 6 #hours

url1= '/home/admin/Work/Programs/FirstStreet_datadownload/FirstStreet/Batch1/ERA5/joaquin2015.nc'

 

file1 = netCDF4.Dataset(url1)

      

time_var_re = file1.variables['time']
Time_data_re = netCDF4.num2date(time_var_re[:],time_var_re.units)

# finding the startdate of the reanalysis date for 
# going one week back 

chk_date1 = datetime.strptime(date,'%Y%m%d%H')

if chk_date1 in Time_data_re:
    res = (Time_data_re == chk_date1).argmax()

    print('extracting from file the data starting:',Time_data_re[res])


#===== subseting the data ===========================
lat  = file1.variables['latitude'][:]

lon  = file1.variables['longitude'][:]

u_data  = file1.variables['u10'][res:,:,:]

v_data  = file1.variables['v10'][res:,:,:]

#time.sleep(2)
P_data  = file1.variables['msl'][res:,:,:]

Time_data = netCDF4.num2date(time_var_re[res:],time_var_re.units)

#wind_mag_re = np.sqrt(u_data**2 + v_data**2)  


#lat1 = lat[300:650] # for brazil
#lon1 = lon[1100:] # for brazil

lon1 = lon[900:1300] #x
lat1 = lat[120:350] #y


#lon1 = 360-lon[120:350]
#lat1 = lat


P_data[0,120:350,900:1300]

#lon1 = 360-lon[120:350]
#lat1 = lat


#P_data[0,120:350,900:1300]
#
#for g in range(0,len(u_data),1):
#    #m.drawcoastlines(linewidth=1.0,ax=ax1)
#    plt.imshow(P_data[g,120:350,900:1300], cmap='jet')
#    plt.title(str(Time_data[g]))
#    plt.pause(0.1)
#
#
#
#
#m = Basemap(projection='cyl',llcrnrlat=2.75,urcrnrlat=60,llcrnrlon=225,urcrnrlon=324.75,resolution='l')
#
#f, (ax1) = plt.subplots(1,figsize=(7,5)) #, sharex='col', sharey='row')
#
##m.bluemarble() 
#ax1.set_title('NCEP-NAM Gridded Winds')
#m.drawcoastlines(linewidth=1.0,ax=ax1)
#
#lat1 = np.arange(-2, 61, 0.25)
#lon1 = np.arange(225, 315, 0.25)
#
#x, y = m(*np.meshgrid(lon,lat))
#         
#data1=wind_mag_re[96,:,:]   # flipping the data for Contour plot
#rangecontour = np.arange(0, 25, 0.5)
#
#cs = m.contourf(x,y,data1,rangecontour, cmap='jet',ax=ax1)
#
##u = u_data[96,:,:] 
##v = v_data[96,:,:]    # streamline provide the direction understanding of the winds
##m.streamplot(x, y, u, v, color='w',density=2, linewidth=0.5,arrowsize=1.,ax=ax1) #cmap=plt.cm.coolwarm
#m.drawparallels(np.arange(-3,60,10),labels=[1,0,0,0],linewidth=0.25,ax=ax1) # draw parallels
#m.drawmeridians(np.arange(225,315,15),labels=[0,0,0,1],linewidth=0.25,ax=ax1) # draw meridians
#
#
#
#from mpl_toolkits.axes_grid1 import make_axes_locatable
##f.colorbar(cs, orientation="horizontal", ax=axlist)
#
#
#divider = make_axes_locatable(ax1)
#cax = divider.append_axes("right", size="3%", pad=0.05)
#plt.colorbar(cs, cax=cax,label='meters per second (m/s)',format = "%d", cmap='jet')
#
#plt.savefig(r'C:\Users\Arslaan Khalid\Documents\FirstStreet\Deliverables\Griddedwinds.png',dpi=100, bbox_inches = 'tight', pad_inches = 0.1)




for e in range (0,(durationrun*(24)),6):
    
    print(f'Processing data file {e}/{(durationrun*(24))} & timestamp:{Time_data[e]}')
    file_number = '%02d'%e

    domain1_u = u_data[e,120:350,900:1300]
    #plt.imshow(domain1_u, cmap='coolwarm')
    
    #========================================== v wind 10m ====================================
    domain1_v = v_data[e,120:350,900:1300]
    #plt.imshow(domain1_v, cmap='coolwarm')
    
    #========================================== MSLP ====================================
    domain1_p = P_data[e,120:350,900:1300]
    #plt.imshow(domain1_p, cmap='coolwarm')
    #==============================================================
    container_u = []      
    
    for i in range(0,domain1_u.shape[0]):
        for j in range(0,domain1_u.shape[1]):
            #print(round(domain1_u[i][j],1))
            key1 = round(domain1_u[i][j],1)
            container_u.append(key1)
    #============================================================    
    container_v = []
        
    for i in range(0,domain1_u.shape[0]):
        for j in range(0,domain1_u.shape[1]):
            #print(round(domain1_v[i][j],1))
            key2 = round(domain1_v[i][j],1)
            container_v.append(key2)
    #=============================================================
    container_p = []
        
    for i in range(0,domain1_u.shape[0]):
        for j in range(0,domain1_u.shape[1]):
            #print(round(domain1_p[i][j],0))
            key3 = round(domain1_p[i][j],0)
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
    #print(f'processed ({e}/{len(u_data)})')
    # writing the data to file                        
    with open('/home/admin/Work/Programs/FirstStreet_datadownload/FirstStreet/Batch1/ERA5/{}.22'.format(stormname),'a') as f:
       winds_tmp1.to_csv(f, sep='\t',index=False,header = False)   

 
##    # controlfile specs also calculate these varaibles as lon[0],lon[-1],lat[0],lat[-1]
ctrl_file = "following is the fort22 domain extents\n{} {} {} {} {} {} {}    ! NWLAT, NWLON, WLATMAX, WLONMIN, WLATINC, WLONINC, WTIMINC(doesnot depedn on timestep)".format(domain1_u.shape[0],domain1_u.shape[1],lat1.max(),(lon1.min()-360),lat1[0]-lat1[1],lon1[1]-lon1[0],(3600*interval))       
with open('/home/admin/Work/Programs/FirstStreet_datadownload/FirstStreet/Batch1/ERA5/ctrl_file_domain_desc.txt','w') as f:
           f.write(ctrl_file)        
