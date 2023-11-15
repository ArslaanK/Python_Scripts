
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 12:10:16 2022

@author: Arslaan.Khalid
"""



import matplotlib.pyplot as plt # Import the Matplotlib package
import numpy as np
import pandas as pd
import netCDF4
import requests
from datetime import datetime, timedelta
import matplotlib.lines as mlines
import pyproj
import json
from shapely.geometry import Point, mapping
from functools import partial
from shapely.ops import transform
import time
from math import pi, atan
import matplotlib.dates as mdates
import h5py
import numpy as np
import tqdm
import time
import pyproj
import json
from shapely.geometry import Point, Polygon, mapping
from functools import partial
from shapely.ops import transform
import matplotlib.dates as mdates
myxlabeldatesFormat = mdates.DateFormatter('%m/%d')

# =============================================================================
# Loading Functions
# =============================================================================
def descend_obj(obj,sep='\t'):
    """
    Iterate through groups in a HDF5 file and prints the groups and datasets names and datasets attributes
    """
    if type(obj) in [h5py._hl.group.Group,h5py._hl.files.File]:
        for key in obj.keys():
            print (sep,'-',key,':',obj[key])
            descend_obj(obj[key],sep=sep+'\t')
    # elif type(obj)==h5py._hl.dataset.Dataset:
    #     for key in obj.attrs.keys():
    #         print (sep+'\t','-',key,'d:',obj.attrs[key])

def h5dump(path,group='/'):
    """
    print HDF5 file metadata

    group: you can give a specific group, defaults to the root group
    """
    with h5py.File(path,'r') as f:
         descend_obj(f[group])

def adcirc_nodalIdentifier(lat,lon,x2,y2):
    chosen_lat = x2  
    chosen_lon = y2
    min_distance = None
    best_index = 0
    
    for i in range(len(lat)):
        current_distance = (lat[i] - chosen_lat)**2 + (lon[i] - chosen_lon)**2
        if min_distance is None or current_distance < min_distance:
            best_index = i
            min_distance = current_distance    
    return best_index         



# =============================================================================
# Read Points file
# =============================================================================
stations_file = r'E:/Projects/LWI R6/RAS_outputs/Isaac_Plans/ReducedBoundary_v1.xlsx'

stns_ = pd.read_excel(stations_file)
stns_.index = stns_.Point

# =============================================================================
# model information
# =============================================================================
plan_file = r'E:\Projects\LWI R6\RAS_outputs\Isaac_Plans/SLaMM_Ray.p27.hdf'


data = h5py.File(plan_file, 'r')


available = data.keys()



model_info_name = data['/Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/2D Flow Areas/'].keys()
for mdl_inf_nm in model_info_name:mdl_inf_nm
    #print(ss)


  
p22=f'/Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/2D Flow Areas/{mdl_inf_nm}/Water Surface'
WSE_lnk = data[p22]
   

WSE_attributes = WSE_lnk.attrs.keys()
wse_units = WSE_lnk.attrs['Units'].decode("utf-8")

Results_summary = data.get('Results').get('Summary').keys()

#-----extracting Geometry datasets
geo = data.get('Geometry')
geo_vars = np.array(geo)


geom_xy = data.get('Geometry').get('2D Flow Areas').get(mdl_inf_nm).get('Cells Center Coordinate')[:,:]#.keys()

x,y = geom_xy[:,0],geom_xy[:,1]

print('some Ras units conversion, therefore have to divide x,y of cells by 3.281\n to convert to meters for display')
x,y=x/3.281,y/3.281


#-----extracting Results datasets
results = data.get('Results')
results.keys()

unsteady = results.get('Unsteady')
unsteady.keys()

ts_ = unsteady.get('Output').get('Output Blocks').get('Base Output').get('Unsteady Time Series').get('2D Flow Areas')

possible_attributes_outputs = ts_.keys()

output_variables = ts_.get(mdl_inf_nm).keys()

wse_grid = ts_.get(mdl_inf_nm).get('Water Surface')
wdepth_grid = ts_.get(mdl_inf_nm).get('Cell Cumulative Excess Depth')
wvel_grid = ts_.get(mdl_inf_nm).get('Face Velocity')



tsteps,cells = wse_grid.shape


wse_grid_out = wse_grid[:,:]
wvel_grid_out = wvel_grid[:,:]
# wdpth_grid_out = wdepth_grid[:,:]

df_wse_out = pd.DataFrame(wse_grid_out)
df_wvel_out = pd.DataFrame(wvel_grid_out)
# df_wdpth_out = pd.DataFrame(wdpth_grid_out)


WSE_tistmps = data.get('/Results/Unsteady/Output/Output Blocks/DSS Profile Output')['Unsteady Time Series']
# WSE_tistmps.keys()
wse_timesteps = WSE_tistmps['Time Date Stamp'][:]

ts_lst = []
for t in wse_timesteps:
    ss = t.decode('utf-8')
    ts_lst.append(ss)
      
df_wse_out.index = pd.to_datetime(ts_lst[:len(df_wse_out)])
df_wvel_out.index = pd.to_datetime(ts_lst[:len(df_wse_out)])
# df_wdpth_out.index = pd.to_datetime(ts_lst)


#%%

import pyproj
transformer = pyproj.Transformer.from_crs("esri:102039", "epsg:4326")
# transformer.transform(x,y)



x2,y2 = stns_['Latitude'].values,stns_['Longitude'].values

xlist,ylist = transformer.transform(x,y)

closet_nodes_RAS= {}



for ppt in tqdm.tqdm(range(0,len(x2))): #     
#for ppt in range(0,len(x1)):
    pp1,pp2 = x2[ppt],y2[ppt]
    
    #pp1,pp2 = transformer.transform(pp1,pp2)
    #print(pp1,pp2)

    idd = adcirc_nodalIdentifier(xlist,ylist,pp1,pp2)
    closet_nodes_RAS[stns_.iloc[ppt]['Point']]=idd
    #closet_nodes.append(idd)
    # else:






#%%

myxlabeldatesFormat = mdates.DateFormatter('%m/%d')
# =============================================================================
# 
# =============================================================================
# write code for user to choose from drop down      
timestep_to_display = 0
#stnx,stny = x[5000],y[5000]
if timestep_to_display!=0:
    z= df_wse_out.iloc[timestep_to_display].values
    plt_lb = f'{df_wse_out.index[50]}'
else:
    z= df_wse_out.max(axis=0).values
    plt_lb = 'Maximum of Simulation'


# =============================================================================
#     
# =============================================================================
# import os
# os.environ['PROJ_LIB'] = r'C:\Users\arslaan.khalid\Ana3\pkgs\proj-7.2.0-h3e70539_0\Library\share'
from mpl_toolkits.basemap import Basemap


llcrnrlat0,urcrnrlat0,llcrnrlon0,urcrnrlon0= 28.5,30.6,268,271.37
xlat,ylat = transformer.transform(x,y)

m = Basemap(projection='cyl',llcrnrlat=llcrnrlat0,urcrnrlat=urcrnrlat0,llcrnrlon=llcrnrlon0,\
            urcrnrlon=urcrnrlon0,resolution='h', epsg = 4269)
    

# =============================================================================
#     #create spatial plot as well  
# =============================================================================
    
fig = plt.figure(figsize = (10,6))
fig.subplots_adjust(wspace=0.1)
ax = plt.subplot(1,1,1) 
#ss1 = ax.scatter(x,y,z,c=z, cmap='jet',vmin=0,vmax=10)

#m.bluemarble(scale=0.2)   # full scale will be overkill
m.drawcoastlines(color='white', linewidth=0.2)  # add coastlines

x4, y4 = m(xlist, ylist)  # transform coordinates
ss1 =plt.scatter( 360+y4,x4, z,c=z, cmap='jet',vmin=0,vmax=10) 
#m.arcgisimage(service='World_Imagery', xpixels = 1500, verbose= False) #1500

plt.colorbar(ss1,label='Water Surface Elevation  ft_NAVD88',format = "%.0f",extend = 'max')


for u3 in stns_.index:   
   
    x2,y2 = stns_.loc[u3]['Latitude'],stns_.loc[u3]['Longitude']

    y2 = y2+360
    
    #ax.scatter(stnx,stny,s=8,c='k',marker='s')
    m.scatter(y2,x2,s=2,c='k',marker='s')
    #plt.text(y2,x2,f'{u3}',fontsize=11,c='m',fontweight='bold')#,f'LocID:{y1}'#f'Stn{y2}'


#ax.set_xlabel('X-coordinates');ax.set_ylabel('Y-coordinates')

ax.set_title(f'Water Surface Elevations [{wse_units}_NAVD88]\n[Timestep: {plt_lb}]')
ax.set_aspect('auto')

    
    
# =============================================================================
#     
# =============================================================================
import random

df_wse_out_clipped = df_wse_out.copy()
  
df_bound = pd.DataFrame(closet_nodes_RAS.items())  

df_ras_boundary = pd.DataFrame()

for colss in df_bound[1].values:
    
    if colss in df_ras_boundary.columns:
        num_gen = round(random.random(),2)
        df_ras_boundary[f'{colss}_{num_gen}']=df_wse_out_clipped[colss]
    else:
        df_ras_boundary[colss]=df_wse_out_clipped[colss]
    
    
bounda_cols = df_bound[1].values 


plt.figure()

k1=0
for k in bounda_cols[:100]:
    if k1==0:
        plt.plot(df_ras_boundary.index,df_ras_boundary[k],c='k',lw=0.75,label='West')    
    else:
        plt.plot(df_ras_boundary.index,df_ras_boundary[k],c='k',lw=0.75) 
    k1=k1+1
for k in bounda_cols[100:200]:
    if k1==100:
        plt.plot(df_ras_boundary.index,df_ras_boundary[k],c='m',lw=0.75,label='Middle')    
    else:
        plt.plot(df_ras_boundary.index,df_ras_boundary[k],c='m',lw=0.75)     
    k1=k1+1
for k in bounda_cols[200:]:
    if k1==200:
        plt.plot(df_ras_boundary.index,df_ras_boundary[k],c='b',lw=0.75,label='East')    
    else:  
        plt.plot(df_ras_boundary.index,df_ras_boundary[k],c='b',lw=0.75)    
    k1=k1+1
    
    
plt.gca().xaxis.set_major_formatter(myxlabeldatesFormat)
 
plt.legend()   





df_ras_boundary.columns = stns_.index.values



df_ras_boundary.to_csv(r'E:\Projects\LWI R6\RAS_outputs\Isaac_Plans\timeseries_along_requested_boundary_Isaac_v2.csv')




#%%


