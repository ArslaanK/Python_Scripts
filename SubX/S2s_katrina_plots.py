# -*- coding: utf-8 -*-
"""
Created on Wed May 15 13:11:44 2019

@author: Arslaan Khalid
"""



#==============================================================================
# S2S file read and Write for ADCIRC
#==============================================================================
 

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
import os
import shutil
from string import Template
import pathlib as pl
path_dir = pl.Path(r'Z:\Project_S2S\3_subx_data\Katrina2005')
model_sets = ['RSMAS-CCSM4','EMC-GEFS','ESRL-FIMr1p1','GMAO-GEOS_V2p1','NRL-NESM','ECCC-GEM']
variable = ['pslmsl','uas10m','vas10m'] #for hindcast data psl is pslmsl
model_available_ensembles_hincast = {'ECCC-GEM':4,'EMC-GEFS':11,'ESRL-FIMr1p1':4,'GMAO-GEOS_V2p1':4,'NRL-NESM':1,'RSMAS-CCSM4':3,'NCEPR1-REAN':1}


StormNamenYear='Katrina2005'

# add YearMonthDate of the weeks under analysis
AnalysisDates = {'EMC-GEFS':[20050803,20050810,20050817,20050824],
                 'ECCC-GEM':[20050726,20050802,20050809,20050816],
                 'ESRL-FIMr1p1':[20050803,20050810,20050817,20050824],
                 'NRL-NESM':[20050731,20050807,20050814,20050821],
                 'GMAO-GEOS_V2p1':[20050730,20050809,20050819,20050824],
                 'RSMAS-CCSM4':[20050730,20050806,20050813,20050820],
                 'NCEPR1-REAN':[20050720,20050727,20050803,20050810,20050817,20050824]
                 } # reanalysis data should be subset according to the models initialization


# plotting reanalysis data ----

# creating datafolder for enesemble
#rootdir=str(path_dir/'{}').format(models)      
#tmp_re = str(path_dir/'Updated_2019_Analysis'/'{}'/'S2S_process_data'/ 'NCEPR1-REAN' / 'tmp' / 'weekbefore.22').format(StormNamenYear)



#============= generating the reanalysis data===============
# ------reading NCEP-Reanalysis files
url_psl_re = str(path_dir/ 'NCEPR1-REAN' / 'pslmsl_NCEPR1-REAN_{}daily.nc').format(StormNamenYear[-4:])
url_uwind_re = str(path_dir/'NCEPR1-REAN' / 'uas_NCEPR1-REAN_{}daily.nc').format(StormNamenYear[-4:])
url_vwind_re = str(path_dir/'NCEPR1-REAN' / 'vas_NCEPR1-REAN_{}daily.nc').format(StormNamenYear[-4:])
  
file4 = netCDF4.Dataset(url_uwind_re)
lat_re  = file4.variables['lat'][:]
lon_re  = file4.variables['lon'][:]
u_data_re  = file4.variables['uas'][:]

file5 = netCDF4.Dataset(url_vwind_re)
v_data_re  = file5.variables['vas'][:]

file6 = netCDF4.Dataset(url_psl_re)
P_data_re  = file6.variables['pslmsl'][:]            

time_var_re = file4.variables['time']
Time_data_re = netCDF4.num2date(time_var_re[:],time_var_re.units)

# finding the startdate of the reanalysis date for 
# going one week back 

date = '20050820'
chk_date = str(date) + '00'
chk_date1 = datetime.strptime(chk_date,'%Y%m%d%H')

if chk_date1 in Time_data_re:
    res = (Time_data_re == chk_date1).argmax()

    print('writing file for data:',Time_data_re[res])

#===== subseting the data ===========================
u_data_re  = file4.variables['uas'][res:][:] # cross check in panoply
v_data_re  = file5.variables['vas'][res:][:]
P_data_re  = file6.variables['pslmsl'][res:][:] 
wind_mag_re = np.sqrt(u_data_re**2 + v_data_re**2)       

#for i in range (0,15):
i = 8 # storm near landfall 28 august was 
print(i)
file_number = '%02d'%i
#========================================== u wind 10m ====================================

domain1_u_flip = np.flipud(u_data_re[i])
domain1_u = domain1_u_flip[30:93]
#domain1_u = u_data_re[i][30:93]    # 150,333
#plt.imshow(domain1_u, cmap='jet')        

domain11_u = numpy.delete(domain1_u, numpy.s_[0:225], axis=1)# right of US cut
domain11_u = numpy.delete(domain11_u, numpy.s_[90:], axis=1)   # left of US cut
#plt.imshow(domain11_u, cmap='jet')      
#========================================== v wind 10m ====================================
domain1_v_flip = np.flipud(v_data_re[i])
domain1_v = domain1_v_flip[30:93]
#domain1_v = v_data_re[i][30:93]    # 150,333
#plt.imshow(domain1_v, cmap='jet')
        
domain11_v = numpy.delete(domain1_v, numpy.s_[0:225], axis=1)# right of US cut
domain11_v = numpy.delete(domain11_v, numpy.s_[90:], axis=1)   # left of US cut
#plt.imshow(domain11_v, cmap='jet')       
#========================================== P ====================================
domain1_p_flip = np.flipud(P_data_re[i])
domain1_p = domain1_p_flip[30:93]
#plt.imshow(domain1_p, cmap='jet')

domain11_p = numpy.delete(domain1_p, numpy.s_[0:225], axis=1)# right of US cut
domain11_p = numpy.delete(domain11_p, numpy.s_[90:], axis=1)   # left of US cut
#plt.imshow(domain11_p, cmap='jet')
#print(i)
#========================================== Wind Magnitude ====================================
domain1_m_flip = np.flipud(wind_mag_re[i])
domain1_m = domain1_m_flip[30:93]

#plt.imshow(domain1_m, cmap='jet')       

domain11_m = numpy.delete(domain1_m, numpy.s_[0:225], axis=1)# right of US cut
domain11_m = numpy.delete(domain11_m, numpy.s_[90:], axis=1)   # left of US cut
#plt.imshow(domain11_m, cmap='jet')


#---------------- plotting scheme--------------------

f, (ax1) = plt.subplots(1,figsize=(7,5)) #, sharex='col', sharey='row')

m = Basemap(projection='cyl',llcrnrlat=-2,urcrnrlat=60,llcrnrlon=225,urcrnrlon=314,resolution='l')
#m.bluemarble() 
ax1.set_title('True (NCEP)')
m.drawcoastlines(linewidth=1.0,ax=ax1)

lat1 = np.arange(-2, 61, 1)
lon1 = np.arange(225, 315, 1)

x, y = m(*np.meshgrid(lon1,lat1))
#            
#            cbar = plt.colorbar(cmap='jet',label='meters per second (m/s)',format = "%d")
#            plt.title('Wind Magnitude\n S2S (1 Degree Resolution)\n')
#             
data1=np.flip(domain11_m, axis=0)    # flipping the data for Contour plot
rangecontour = np.arange(0, 25, 0.5)

cs = m.contourf(x,y,data1,rangecontour, cmap='jet',ax=ax1)

u = np.flip(domain11_u, axis=0) 
v = np.flip(domain11_v, axis=0)   # streamline provide the direction understanding of the winds
m.streamplot(x, y, u, v, color='w',density=3, linewidth=0.5,arrowsize=1.,ax=ax1) #cmap=plt.cm.coolwarm
m.drawparallels(np.arange(-3,60,10),labels=[1,0,0,0],linewidth=0.25,ax=ax1) # draw parallels
m.drawmeridians(np.arange(225,315,15),labels=[0,0,0,1],linewidth=0.25,ax=ax1) # draw meridians



from mpl_toolkits.axes_grid1 import make_axes_locatable
#f.colorbar(cs, orientation="horizontal", ax=axlist)


divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="3%", pad=0.05)
plt.colorbar(cs, cax=cax,label='meters per second (m/s)',format = "%d", cmap='jet')

plt.savefig(r'C:\Users\Arslaan Khalid\Desktop\S2S\Updated_2019_Analysis\Katrina2005\Figures\true.png',dpi=300, bbox_inches = 'tight', pad_inches = 0.1)









#now doing this for all the s2s models
    ##-------------------- Generating S2S datasets -------------
model_sets = ['RSMAS-CCSM4','EMC-GEFS','ESRL-FIMr1p1','GMAO-GEOS_V2p1','NRL-NESM','ECCC-GEM']
model_available_ensembles_hincast = {'ECCC-GEM':4,'EMC-GEFS':11,'ESRL-FIMr1p1':4,'GMAO-GEOS_V2p1':4,'NRL-NESM':1,'RSMAS-CCSM4':3,'NCEPR1-REAN':1}
model_plot = {'ECCC-GEM':4,'EMC-GEFS':10,'ESRL-FIMr1p1':3,'GMAO-GEOS_V2p1':4,'NRL-NESM':1,'RSMAS-CCSM4':3}



f, ((ax1, ax2,ax3),(ax4,ax5, ax6)) = plt.subplots(2, 3,figsize=(15,15)) #, sharex='col', sharey='row')

m = Basemap(projection='cyl',llcrnrlat=-2,urcrnrlat=60,llcrnrlon=225,urcrnrlon=314,resolution='l')
#m.bluemarble() 
#ax1.set_title('True (NCEP)')
#m.drawcoastlines(linewidth=1.0,ax=ax1)

lat1 = np.arange(-2, 61, 1)
lon1 = np.arange(225, 315, 1)

x, y = m(*np.meshgrid(lon1,lat1))
#            
#            cbar = plt.colorbar(cmap='jet',label='meters per second (m/s)',format = "%d")
#            plt.title('Wind Magnitude\n S2S (1 Degree Resolution)\n')



k = 1

for models in model_sets:    
    
    #models = 'NRL-NESM'
    date = str(AnalysisDates[models][3])   
    landfall = '20050828'
    start_date = datetime.strptime(date,'%Y%m%d')
    landfall_date = datetime.strptime(landfall,'%Y%m%d')

    diff = int((landfall_date-start_date).days)

    
#for k in range(1,model_available_ensembles_hincast[models]+1): 
#   print(k)
    #for ensemble in range(1,model_available_ensembles_hincast[models]+1):
    ensemble=model_plot[models] # uncomment this when doing for ncep renalysis
    ensembleno = "e{}".format(ensemble)
    stdout = "{}_{}".format(date,ensembleno)
    
    mydir = os.path.join(rootdir, stdout)
                
    # reading S2S model files
    url_psl = str(path_dir  / '{}' / 'psl_msl_{}_{}.e{}.daily.nc').format(models,models,date,ensemble)
    url_uwind = str(path_dir / '{}' / 'uas_10m_{}_{}.e{}.daily.nc').format(models,models,date,ensemble)
    url_vwind = str(path_dir  / '{}' / 'vas_10m_{}_{}.e{}.daily.nc').format(models,models,date,ensemble)

           
    #-------------------reading data from S2S files----
  
    file1 = netCDF4.Dataset(url_uwind)
    lat  = file1.variables['lat'][:]
    lon  = file1.variables['lon'][:]
    u_data  = file1.variables['uas'][:]
    
    file2 = netCDF4.Dataset(url_vwind)
    v_data  = file2.variables['vas'][:]
    
    file3 = netCDF4.Dataset(url_psl)
    P_data  = file3.variables['psl'][:]

    wind_mag = np.sqrt(u_data**2 + v_data**2)

    
    i = diff # use this yourself to match the 28th august
    #print(i)
    file_number = '%02d'%i
    #========================================== u wind 10m ====================================
    if models == "RSMAS-CCSM4" or models == "ESRL-FIMr1p1" or models == "GMAO-GEOS_V2p1" or models == "NCEPR1-REAN" or models == "NRL-NESM" or models == "ECCC-GEM":
        domain1_u_flip = np.flipud(u_data[i])
        domain1_u = domain1_u_flip[30:93]
    else:
        domain1_u = u_data[i][30:93]    # 150,333
    #plt.imshow(domain1_u, cmap='jet')        
    
    domain11_u = numpy.delete(domain1_u, numpy.s_[0:225], axis=1)# right of US cut
    domain11_u = numpy.delete(domain11_u, numpy.s_[90:], axis=1)   # left of US cut
    #plt.imshow(domain11_u, cmap='jet')      
    #========================================== v wind 10m ====================================
    if models == "RSMAS-CCSM4" or models == "ESRL-FIMr1p1" or models == "GMAO-GEOS_V2p1" or models == "NCEPR1-REAN" or models == "NRL-NESM" or models == "ECCC-GEM":
        domain1_v_flip = np.flipud(v_data[i])
        domain1_v = domain1_v_flip[30:93]
    else:
        domain1_v = v_data[i][30:93]    # 150,333
    #plt.imshow(domain1_v, cmap='jet')
            
    domain11_v = numpy.delete(domain1_v, numpy.s_[0:225], axis=1)# right of US cut
    domain11_v = numpy.delete(domain11_v, numpy.s_[90:], axis=1)   # left of US cut
    #plt.imshow(domain11_v, cmap='jet')       
    #========================================== P ====================================
    if models == "RSMAS-CCSM4" or models == "ESRL-FIMr1p1" or models == "GMAO-GEOS_V2p1" or models == "NCEPR1-REAN" or models == "NRL-NESM" or models == "ECCC-GEM":
        domain1_p_flip = np.flipud(P_data[i])
        domain1_p = domain1_p_flip[30:93]
    else:
        domain1_p = P_data[i][30:93]    # 150,333
    #plt.imshow(domain1_p, cmap='jet')
    
    domain11_p = numpy.delete(domain1_p, numpy.s_[0:225], axis=1)# right of US cut
    domain11_p = numpy.delete(domain11_p, numpy.s_[90:], axis=1)   # left of US cut
    #plt.imshow(domain11_p, cmap='jet')
    #========================================== Wind Magnitude ====================================
    if models == "RSMAS-CCSM4" or models == "ESRL-FIMr1p1" or models == "GMAO-GEOS_V2p1" or models == "NCEPR1-REAN" or models == "NRL-NESM" or models == "ECCC-GEM":
        domain1_m_flip = np.flipud(wind_mag[i])
        domain1_m = domain1_m_flip[30:93]
    else:                        
        domain1_m = wind_mag[i][30:93]    # 150,333
    #plt.imshow(domain1_m, cmap='jet')       
    
    domain11_m = numpy.delete(domain1_m, numpy.s_[0:225], axis=1)# right of US cut
    domain11_m = numpy.delete(domain11_m, numpy.s_[90:], axis=1)   # left of US cut
#    plt.imshow(domain11_m, cmap='jet')
#    plt.title(k)
#    plt.pause(0.1)
#    
    #--------------------------------------------------------
    
    exec('g=ax{}'.format(k))
    data1=np.flip(domain11_m, axis=0)    # flipping the data for Contour plot
    
    g.set_title(models)
    m.drawcoastlines(linewidth=1.0,ax=g)

    
    cs = m.contourf(x,y,data1,rangecontour, cmap='jet',ax=g)
#    if k ==3 or k ==6:
#        cbar = plt.colorbar(cs,cmap='jet',label='meters per second (m/s)',format = "%d",ax=g)
#    #plt.clabel(cs, cs.levels, fmt='%.1f', fontsize=8, colors='k')  
    
    u = np.flip(domain11_u, axis=0) 
    v = np.flip(domain11_v, axis=0)   # streamline provide the direction understanding of the winds
    m.streamplot(x, y, u, v, color='w',density=1.5, linewidth=0.5,arrowsize=1.,ax=g) #cmap=plt.cm.coolwarm
    m.drawparallels(np.arange(-3,60,10),labels=[1,0,0,0],linewidth=0.25,ax=g) # draw parallels
    m.drawmeridians(np.arange(225,315,15),labels=[0,0,0,1],linewidth=0.25,ax=g) # draw meridians

    k = k+1    
    

#cbar = m.colorbar(cmap='jet',label='meters per second (m/s)',format = "%d")

from mpl_toolkits.axes_grid1 import make_axes_locatable
#f.colorbar(cs, orientation="horizontal", ax=axlist)


divider = make_axes_locatable(ax6)
cax = divider.append_axes("right", size="3%", pad=0.05)
plt.colorbar(cs, cax=cax,label='meters per second (m/s)',format = "%d", cmap='jet')

divider = make_axes_locatable(ax3)
cax = divider.append_axes("right", size="3%", pad=0.05)
plt.colorbar(cs, cax=cax,label='meters per second (m/s)',format = "%d", cmap='jet')


plt.savefig(r'C:\Users\Arslaan Khalid\Desktop\S2S\Updated_2019_Analysis\Katrina2005\Figures\models2.png',dpi=300, bbox_inches = 'tight', pad_inches = 0.1)







#
#
#    
#        
#        
#        
#        #=============================== Pressure ==============
#        
#        
#        fig = plt.figure(num=None, figsize=(8, 6) )
#        m = Basemap(projection='cyl',llcrnrlat=-2,urcrnrlat=60,llcrnrlon=225,urcrnrlon=314,resolution='l')
#        
#        #m.bluemarble()    
#        m.drawcoastlines()
#        
#        #m.imshow(domain11_p, cmap='jet',origin='upper',vmin=domain11_p.min(), vmax=domain11_p.max(), aspect='auto')
#       
#        
#        #cbar = plt.colorbar(cmap='jet',label='Pressure (pa)',format = "%d")
#        #plt.title('Mean Sea Level Pressure \n S2S (1 Degree Resolution)\n')
#         
#        data2=np.flip(domain11_p, axis=0)    # flipping the data for Contour plot
#        
#        rangecontour2 = np.arange(P_data.min(), P_data.max(), 250)
#        cs = m.contourf(x,y,data2,rangecontour2)
#    
#        
#        m.drawparallels(np.arange(-3,60,10),labels=[1,0,0,0],linewidth=0.25) # draw parallels
#        m.drawmeridians(np.arange(225,315,15),labels=[0,0,0,1],linewidth=0.25) # draw meridians
#    
#    
#        plt.savefig('S2S_P{}.png'.format(file_number),dpi=300, bbox_inches = 'tight', pad_inches = 0.1)
#        plt.close()
#        
#        plt.clabel(cs, fmt='%d', fontsize=12, colors='k')
#        
#        
#
#
#






