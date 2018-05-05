
# -*- coding: utf-8 -*-
"""
Created on Thu May 04 17:24:14 2018
@author: Arslaan Khalid, akhalid6@gmu.edu
"""

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt # Import the Matplotlib package
from osgeo import gdal # Import the GDAL library
import sys,getopt
import pandas as pd
import numpy.core.multiarray 
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
from datetime import datetime, timedelta


from matplotlib import pylab
from matplotlib.dates import DayLocator, HourLocator, DateFormatter
from matplotlib.font_manager import FontProperties
import pandas as pd


from matplotlib.pyplot import text


import netCDF4
import numpy.ma as ma
import matplotlib.tri as tri
from django.utils import timezone


#---Get Start Date From fort.221
from datetime import datetime

fort22 = '/home/fhrl/Documents/Surge/input/Domain2/fort22/fort.221' #for Domain2


# this needs to be a recent data since only recent models are available for download

with open(fort22,'r') as f:
    for i in range(0,1):
        line = f.readline().strip().split()
        timestamp = line[3]
        start_date = datetime.strptime(timestamp,'%Y%m%d%H')



url='/home/fhrl/Documents/KalpanaTesting_DONE/ChesapeakeBay/fort.63.nc'
#Date_UTC = '' Incase cant't access the model start date(manual entry)
#start_date = datetime.strptime(Date_UTC ,'%Y%m%d%H')





# ==================================== Reading the NETCDF formatted fort.63 ====================
file1 = netCDF4.Dataset(url)
lat  = file1.variables['y'][:]
lon  = file1.variables['x'][:]

#Longitudinal Correction

lon[lon < 180] = lon[lon < 180] + 360

gridfile=url
gridvars = netCDF4.Dataset(gridfile).variables

var_element = 'element'

elems = gridvars[var_element][:,:]-1  # Move to 0-indexing by subtracting 1, elements indexing starts with '1' in netcdf file


# basemap 

m = Basemap(projection='cyl',llcrnrlat=36,urcrnrlat=40,llcrnrlon=282,urcrnrlon=286,resolution='h', epsg = 4269)
#US is 4269, and you can google the region you want


for i in range(0,83):
    
    i=i+1
   
    data = file1.variables['zeta'][i,:]
    z = data.data
    file_number = '%02d'%i
    
    import matplotlib.tri as tri
    print ('\nTriangulating ...\n')
    triang = tri.Triangulation(lon,lat, triangles=elems)

    if data.mask.any():
      # -99999 entries in 'data' array are usually masked, mask all corresponding triangles
      point_mask_indices = numpy.where(data.mask)
      tri_mask = numpy.any(numpy.in1d(elems, point_mask_indices).reshape(-1, 3), axis=1)
      triang.set_mask(tri_mask)

    # First create the x and y coordinates of the points.
  
    levels = np.arange(-1.5, 5, 0.1)
  
    #plt.scatter(lon, lat, s=0.0125, c=z,cmap='jet')
        
    print ('Making contours ...\n')
    plt.tricontourf(triang, data, levels=levels,alpha=0.9,vmin=-1.5, vmax=5, aspect='auto',cmap='jet')
   
    
    m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 800, verbose= False)
    #World_Street_Map,ESRI_Imagery_World_2D,NatGeo_World_Map ,World_Imagery (type of Data to use from)  

        
    m.drawcoastlines(color='k')

    
    plt.title('Water Levels\n Chesapeake Bay Flood Forecast System\n')

    # Scatte point for capital and specified cities

    plt.xlim([282, 286])
    plt.ylim([36, 40])
    
    plt.colorbar(cmap='jet',label='Water Levels (meters realtive to MSL)',format = "%.1f")
    
    x, y = m(283, 38.9)
    plt.plot(x, y, '*', markersize=12,color='white')
    plt.text(x, y, ' Washington,D.C.\n', fontsize=14,color='white')
    x, y = m(284, 37.2)
    plt.plot(x, y, 'ok', markersize=8,color='white')
    plt.text(x, y, ' Kiptopeke\n', fontsize=14,color='white')
    plt.xlabel('\n\nDate (UTC):{} \nForecast Hour: {} \n'.format(start_date+ timedelta(hours=i),file_number))
    
    plt.savefig('/home/fhrl/Documents/Surge/results/figs_WL/WL{}.png'.format(file_number),dpi=300, bbox_inches = 'tight', pad_inches = 0.1)
    plt.close() 
