# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 12:24:17 2018

@author: Arslaan Khalid
"""


#==============================================================================
# S2S file read and Write for ADCIRC
#==============================================================================
 

#
import matplotlib.pyplot as plt # Import the Matplotlib package
#from osgeo import gdal # Import the GDAL library
#import sys,getopt
import pandas as pd
import numpy.core.multiarray 
#import wget
#import geojsoncontour
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np

from datetime import datetime, timezone, timedelta


import netCDF4

import numpy.ma as ma



# file specificatons

#Xinterval=1 
#Yinterval=1
#DX=90
#DY=63
#max60
#min-135
#=================


start = 239

#for t in range(1,12):   

    

    # set up the URL to access the data server.

    # See the NWW3 directory on NOMADS

    # for the list of available model run dates.

   

    #mydate='20180503'  # this needs to be a recent data since only recent models are available for download

url1= r'Z:\Project_S2S\3_subx_data\isabel\fort22\verified_winds\uas_NCEPR1-REAN_2003daily.nc'

url2= r'Z:\Project_S2S\3_subx_data\isabel\fort22\verified_winds\vas_NCEPR1-REAN_2003daily.nc'

url3= r'Z:\Project_S2S\3_subx_data\isabel\fort22\verified_winds\pslmsl_NCEPR1-REAN_2003daily.nc'

 

 

file1 = netCDF4.Dataset(url1)

lat  = file1.variables['lat'][:]

lon  = file1.variables['lon'][:]

u_data  = file1.variables['uas'][start:start+35,:,:]

 

file2 = netCDF4.Dataset(url2)

v_data  = file2.variables['vas'][start:start+35,:,:]

 

file3 = netCDF4.Dataset(url3)

P_data  = file3.variables['pslmsl'][start:start+35,:,:]

#lat3  = file3.variables['lat'][:]
#
#lon3  = file3.variables['lon'][:]

#########lon[lon < 180] = lon[lon < 180] + 360


u_flip = np.flipud(u_data[0])

v_flip = np.flipud(v_data[0])

p_flip = np.flipud(P_data[0])


 

#wind_mag = np.sqrt(u_data**2 + v_data**2)

'''
# plotting scheme =================
for i in range(0,35):
    plt.imshow(u_data[i], cmap='jet',origin='upper')
    plt.savefig('S2S{}.png'.format(i))
'''



### for checking the variabl
#gridvars = netCDF4.Dataset(gridfile).variables

#### GRID INFO=============================

# this data is every 1 degree 60 to -2 north and 225 to 315 degree east
# find the cut numbers

m = Basemap(projection='cyl',llcrnrlat=-2,urcrnrlat=60,llcrnrlon=225,urcrnrlon=314,resolution='l')
 


for i in range (0,35):
    print(i)
    #i=0
    file_number = '%02d'%i
    #========================================== u wind 10m ====================================
    u_data_flip = np.flipud(u_data[i])           # flipping the 
    domain1_u = u_data_flip[30:93]    # 150,333
    plt.imshow(domain1_u, cmap='jet')
    
    
    domain11_u = numpy.delete(domain1_u, numpy.s_[0:225], axis=1)# right of US cut
    domain11_u = numpy.delete(domain11_u, numpy.s_[90:], axis=1)   # left of US cut
    plt.imshow(domain11_u, cmap='jet')
    
    #========================================== v wind 10m ====================================
    v_data_flip = np.flipud(v_data[i])
    domain1_v = v_data_flip[30:93]    # 150,333
    #plt.imshow(domain1_v, cmap='jet')
    
    
    domain11_v = numpy.delete(domain1_v, numpy.s_[0:225], axis=1)# right of US cut
    domain11_v = numpy.delete(domain11_v, numpy.s_[90:], axis=1)   # left of US cut
    plt.imshow(domain11_v, cmap='jet')
    
    #========================================== P ====================================
    p_data_flip = np.flipud(P_data[i])
    domain1_p = p_data_flip[30:93]    # 150,333
    plt.imshow(domain1_p, cmap='jet')
    
    
    domain11_p = numpy.delete(domain1_p, numpy.s_[0:225], axis=1)# right of US cut
    domain11_p = numpy.delete(domain11_p, numpy.s_[90:], axis=1)   # left of US cut
    plt.imshow(domain11_p, cmap='jet')
    
    
    #========================================== Wind Magnitude ====================================
#    domain1_m = wind_mag[i][30:93]    # 150,333
#    #plt.imshow(domain1_m, cmap='jet')
#    
#    
#    domain11_m = numpy.delete(domain1_m, numpy.s_[0:225], axis=1)# right of US cut
#    domain11_m = numpy.delete(domain11_m, numpy.s_[90:], axis=1)   # left of US cut
#    plt.imshow(domain11_m, cmap='jet')
#    
#    
#    
#    
#    
#    # check using panoply
#    lat = np.arange(-2,61,1)   # domain of S2S cut
#    lon = np.arange(225,315,1)
#    rangecontour1 = np.arange(wind_mag.min(), wind_mag.max(), 0.5)
#    
        
    
    
    
    
    
    
    
    #======================= plotting scheme ==================================
    
    #======================= Wind ==================================
    
#    fig = plt.figure(num=None, figsize=(8, 6) )
#   
#    #m.bluemarble()    
#    m.drawcoastlines()
#    
#    m.imshow(domain11_m, cmap='jet',origin='upper',vmin=wind_mag.min(), vmax=wind_mag.max(), aspect='auto')
#    x, y = m(*np.meshgrid(lon,lat))
#    
#    cbar = plt.colorbar(cmap='jet',label='meters per second (m/s)',format = "%d")
#    plt.title('Wind Magnitude\n S2S (1 Degree Resolution)\n')
#     
#    data1=np.flip(domain11_m, axis=0)    # flipping the data for Contour plot
#    
#    cs = m.contourf(x,y,data1,rangecontour1)
#    u = np.flip(domain11_u, axis=0) 
#    v = np.flip(domain11_v, axis=0)   # streamline provide the direction understanding of the winds
#    m.streamplot(x, y, u, v, color='w',density=1, linewidth=0.5,arrowsize=3.) #cmap=plt.cm.coolwarm
#    #m.quiver(x, y, u, v, units='width')
#
#    x, y = m(283, 38.9)
#    plt.plot(x, y, '*', markersize=8,color='w')
#    x, y = m(275.0, 37.5)
#    plt.text(x, y, ' Washington,\n        D.C.\n', fontsize=7,color='w')
#
#    m.drawparallels(np.arange(-3,60,10),labels=[1,0,0,0],linewidth=0.25) # draw parallels
#    m.drawmeridians(np.arange(225,315,15),labels=[0,0,0,1],linewidth=0.25) # draw meridians

#    
#    plt.savefig(r'C:\Users\Arslaan Khalid\Desktop\W{}.png'.format(file_number),dpi=300, bbox_inches = 'tight', pad_inches = 0.1)
#    plt.close()
    
    #plt.clabel(cs, fmt='%d', fontsize=12, colors='k')
    
    # for directions (u and v)
    

    
#    #=============================== Pressure ==============
#    
#         
#        flipped_winds=np.flip(domain11_p, axis=0)    # flipping the data for Contour plot
#        
#        fig = plt.figure(num=None, figsize=(8, 6) )
#        m = Basemap(projection='cyl',llcrnrlat=-2,urcrnrlat=60,llcrnrlon=225,urcrnrlon=314,resolution='l')
#        
#        #m.bluemarble()    
#        m.drawcoastlines()
#        
#        m.imshow(flipped_winds, cmap='jet',origin='upper',vmin=domain11_p.min(), vmax=domain11_p.max(), aspect='auto')
#        x, y = m(*np.meshgrid(lon,lat))
#        
#        cbar = plt.colorbar(cmap='jet',label='Pressure (pa)',format = "%d")
#        plt.title('Mean Sea Level Pressure \n S2S (1 Degree Resolution)\n')
#         
#        data2=np.flip(domain11_p, axis=0)    # flipping the data for Contour plot
#        
#        rangecontour2 = np.arange(domain11_p.min(), domain11_p.max(), 50)
#        cs = m.contourf(x,y,data2,rangecontour2)
#    
#        x, y = m(283, 38.9)
#        plt.plot(x, y, '*', markersize=8,color='w')
#        x, y = m(275.0, 37.5)
#        plt.text(x, y, ' Washington,\n        D.C.\n', fontsize=7,color='k')
#        m.drawparallels(np.arange(-3,60,10),labels=[1,0,0,0],linewidth=0.25) # draw parallels
#        m.drawmeridians(np.arange(225,315,15),labels=[0,0,0,1],linewidth=0.25) # draw meridians
#    
##    
#        plt.savefig('S2S_P{}.png'.format(file_number),dpi=300, bbox_inches = 'tight', pad_inches = 0.1)
#        plt.close()
    
    #plt.clabel(cs, fmt='%d', fontsize=12, colors='k')
    
    
    
    ##### subsection of a area
    import pandas as pd
    
    #==============================================================
    container_u = []
    
    for i in range(0,63):
        for j in range(0,90):
            print(round(domain11_u[i][j],1))
            key1 = round(domain11_u[i][j],1)
            container_u.append(key1)
    #============================================================    
    container_v = []
    
    for i in range(0,63):
        for j in range(0,90):
            print(round(domain11_v[i][j],1))
            key2 = round(domain11_v[i][j],1)
            container_v.append(key2)
    
    #=============================================================
    container_p = []
    
    for i in range(0,63):
        for j in range(0,90):
            print(round(domain11_p[i][j],0))
            key3 = round(domain11_p[i][j],0)
            container_p.append(key3)
     
     
    #=============================================================
    
     
    #=============================================================
    S2S_data = pd.DataFrame()
    
    S2S_data['u']=container_u
    S2S_data['v']=container_v
    S2S_data['p']=container_p
    
    S2S_data1 =  S2S_data.round(1) 
    
    # writing the data to file
    
    with open(r'C:\Users\Arslaan Khalid\Desktop\fort.222','a') as f:
       S2S_data1.to_csv(f, sep='\t',index=False,header = False)   
    


