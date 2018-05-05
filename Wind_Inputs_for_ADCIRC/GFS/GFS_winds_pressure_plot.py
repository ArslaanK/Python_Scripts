

# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 17:09:02 2018

@author: Arslaan Khalid
"""

#==============================================================================
# GFS grib data downloader and reader 
#==============================================================================
 
import matplotlib.pyplot as plt # Import the Matplotlib package
from osgeo import gdal # Import the GDAL library
import sys,getopt
import pandas as pd
import numpy.core.multiarray 
import wget
#import geojsoncontour
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
from datetime import datetime, timezone, timedelta


#Basemap?   for details
# file dowloading (for first 10 hours)

for i in range(10,85): 
    file_number = '%03d'%i
    print('Beginning file download with wget module: forecated hour {}'.format(file_number))
    

    url = 'http://ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/gfs.2018050112/gfs.t12z.pgrb2.0p25.f{}'.format(file_number) 
    output_file = r'C:\Users\Arslaan Khalid\Desktop\GFS\{}.grb2'.format(i)
    
    wget.download(url,output_file )


# creating a new file 
myfile = open(r'C:\Users\Arslaan Khalid\Desktop\GFS\fort.22', 'w').close()

#from datetime import datetime
#now = str(datetime.datetime.now('UTC'))

Date_UTC = url[52:62]

start_date = datetime.strptime(Date_UTC ,'%Y%m%d%H')

for i in range(0,85):
    # Read the GRIB file
    print(i)
    file_number = '%02d'%i
    print('Beginning file download with wget module: forecated hour {}'.format(file_number))
    

    #url = 'http://ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/gfs.2018050112/gfs.t12z.pgrb2.0p25.f{}'.format(file_number) 
   
    input_file = r'C:\Users\Arslaan Khalid\Desktop\GFS\{}.grb2'.format(i)

    print(input_file)
    grib = gdal.Open(input_file)         # one alternative is to download the file or else read through online, takes about the same time
    
    
    ds = grib 
    # Grid Information
    
    #lat-lon grid:(1440 x 721) units 1e-06 input WE:NS output WE:SN res 48
    #lat 90.000000 to -90.000000 by 0.250000
    #lon 0.000000 to 359.750000 by 0.250000 #points=1038240
    
    #////////////////
    
    if i == 0:
        u_value=266
        v_value=267
        p_value=352
    else:
        u_value=288
        v_value=289
        p_value=415
        
    #u_data of wind
    u_wind = grib.GetRasterBand(u_value)      # 288 if after 0 hour 266
    u_data = u_wind.ReadAsArray() 
    
    #v_data of wind
    v_wind = grib.GetRasterBand(v_value)     #289 if after 0 hour else 267
    v_data = v_wind.ReadAsArray() 
    
    
    #P_data of wind
    P = grib.GetRasterBand(p_value)         #415 if after 0 hour else 352
    P_data = P.ReadAsArray() 
    
    
    import numpy as np 
    
    wind_mag = np.sqrt(u_data**2 + v_data**2)
    
    '''
    # Read the number of bands inside the GRIB file
    number_of_bands = grib.RasterCount
    print(number_of_bands)
    
    
    # Get the band name and description
    for i in range(1,418):
        print(i)
        band = grib.GetRasterBand(i)
        metadata = band.GetMetadata()
        band_name = metadata['GRIB_COMMENT']
        band_description = band.GetDescription()
        print(band_name)
        print(band_description)
        print("-----------------------------")
    
    
    
    # Show the image
    plt.imshow(u_data, cmap='coolwarm')
    
    plt.imshow(v_data, cmap='jet')
    
    plt.imshow(P_data)
    
    plt.imshow(wind_mag, cmap='coolwarm')
    '''
    
    ##### subsection of a area
    
    
    
    #========================================== u wind 10m ====================================
    domain1_u = u_data[120:375]    # 150,333
    #plt.imshow(domain1_u, cmap='coolwarm')
    
    
    domain11_u = numpy.delete(domain1_u, numpy.s_[0:1000], axis=1)# right of US cut
    domain11_u = numpy.delete(domain11_u, numpy.s_[220:], axis=1)   # left of US cut
    #plt.imshow(domain11_u, cmap='coolwarm')
    #========================================== v wind 10m ====================================
    domain1_v = v_data[120:375]    # 150,333
    #plt.imshow(domain1_v, cmap='coolwarm')
    
    
    domain11_v = numpy.delete(domain1_v, numpy.s_[0:1000], axis=1)# right of US cut
    domain11_v = numpy.delete(domain11_v, numpy.s_[220:], axis=1)   # left of US cut
    #plt.imshow(domain11_v, cmap='coolwarm')
    #========================================== MSLP ====================================
    domain1_p = P_data[120:375]    # 150,333
    #plt.imshow(domain1_p, cmap='coolwarm')
    
    
    domain11_p = numpy.delete(domain1_p, numpy.s_[0:1000], axis=1)# right of US cut
    domain11_p = numpy.delete(domain11_p, numpy.s_[220:], axis=1)   # left of US cut
    #plt.imshow(domain11_p, cmap='coolwarm')
    #domain11[163] # represents US
    
    #domain11[2][0:]
    
    #==================================== wind_mag ====================
    
    domain1_m = wind_mag[120:375]    # 150,333
    #plt.imshow(domain1_p, cmap='coolwarm')
    
    
    domain11_m = numpy.delete(domain1_m, numpy.s_[0:1000], axis=1)# right of US cut
    domain11_m = numpy.delete(domain11_m, numpy.s_[220:], axis=1)   # left of US cut
    

    lat = np.arange(-3.5,60.25,0.25)   # domain of GFS cut
    lon = np.arange(250,305.0,0.25)
    rangecontour1 = np.arange(0, 30, 5.0)
    
    
    #
    #============================== Plotting schemes ================================
    #
    fig = plt.figure(num=None, figsize=(8, 6) )
    m = Basemap(projection='cyl',llcrnrlat=-3.5,urcrnrlat=60,llcrnrlon=250,urcrnrlon=305,resolution='l')
    
    m.bluemarble()    
    m.drawcoastlines()
    m.imshow(domain11_m, cmap='jet',origin='upper',alpha=0.7,vmin=0, vmax=30, aspect='auto')
    x, y = m(*np.meshgrid(lon,lat))
    cbar = plt.colorbar(cmap='jet',label='meters per second (m/s)',format = "%d")
    plt.title('Wind Magnitude\n Global Forecast System (0.25 Degree Resolution)\n')
 
    data1=np.flip(domain11_m, axis=0)    # flipping the data for Contour plot

    cs = m.contour(x,y,data1,rangecontour1)
    plt.clabel(cs, fmt='%d', fontsize=12, colors='k')
    
    # for directions (u and v)
    
    u = np.flip(domain11_u, axis=0) 
    v = np.flip(domain11_v, axis=0)     
    
    #========================
    #m.barbs(x,y,u,v)

    
    m.drawparallels(np.arange(-3.5,60.,15.),labels=[1,0,0,0]) # draw parallels
    m.drawmeridians(np.arange(250.,305.,15.),labels=[0,0,0,1]) # draw meridians

    plt.xlabel('\n\nDate (UTC):{} \nForecast Hour: {} \n'.format(start_date+ timedelta(hours=i),file_number))

    plt.savefig(r'C:\Users\Arslaan Khalid\Desktop\GFS\figs\W{}.png'.format(file_number),dpi=300, bbox_inches = 'tight', pad_inches = 0.25)
    plt.close()  


    #=======================================================================================================
    # for MSLP
    fig = plt.figure(num=None, figsize=(8, 6) )
    m = Basemap(projection='cyl',llcrnrlat=-3.5,urcrnrlat=60,llcrnrlon=250,urcrnrlon=305,resolution='l')
    m.bluemarble() 
    
    m.drawcoastlines()
    m.imshow(domain11_p,  cmap='jet',origin='upper',alpha=0.7,vmin=90900, vmax=105000, aspect='auto')
    x, y = m(*np.meshgrid(lon,lat))
 
    cbar = plt.colorbar(cmap='jet',label='Pressure reduced to MSL (Pa)',format = "%d")
    plt.title('Mean Sea Level Pressure\n Global Forecast System (0.25 Degree Resolution)\n')

    data=np.flip(domain11_p, axis=0)    # flipping the data for Contour plot
    
    rangecontour = np.arange(90900, 105000, 500.0)
    
    cs = m.contour(x,y,data,rangecontour)
    plt.clabel(cs, fmt='%d', fontsize=10, colors='k')
        
    m.drawparallels(np.arange(-3.5,60.,15.),labels=[1,0,0,0]) # draw parallels
    m.drawmeridians(np.arange(250.,305.,15.),labels=[0,0,0,1]) # draw meridians

    plt.xlabel('\n\nDate (UTC):{} \nForecast Hour: {} \n'.format(start_date+ timedelta(hours=i),file_number))

    plt.savefig(r'C:\Users\Arslaan Khalid\Desktop\GFS\figs\P{}.png'.format(file_number),dpi=300, bbox_inches = 'tight', pad_inches = 0.25)
    plt.close() 

    
    '''
    container_u = []
    
    for i in range(0,255):
        for j in range(0,220):
            print(round(domain11_u[i][j],1))
            key1 = round(domain11_u[i][j],1)
            container_u.append(key1)
    #============================================================    
    container_v = []
    
    for i in range(0,255):
        for j in range(0,220):
            print(round(domain11_v[i][j],1))
            key2 = round(domain11_v[i][j],1)
            container_v.append(key2)
    
    #=============================================================
    container_p = []
    
    for i in range(0,255):
        for j in range(0,220):
            print(round(domain11_p[i][j],0))
            key3 = round(domain11_p[i][j],0)
            container_p.append(key3)
     
     
    #=============================================================
    
     
    #=============================================================
    GFS_data = pd.DataFrame()
    
    GFS_data['u']=container_u
    GFS_data['v']=container_v
    GFS_data['p']=container_p
    
    '''
    
    # writing the data to file
    
    #with open(r'C:\Users\Arslaan Khalid\Desktop\GFS\fort.22','a') as f2:
       #GFS_data.to_csv(f2, sep='\t',index=False,header = False)   
       
       
    #open(r'C:\Users\Arslaan Khalid\Desktop\GFS\fort.22', 'a').close()

#open(r'C:\Users\Arslaan Khalid\Desktop\GFS\fort.22', 'a').close()





'''

#
#   wave watch 3 -significant wave height
#


import netCDF4

# set up the figure
plt.figure()

# set up the URL to access the data server.
# See the NWW3 directory on NOMADS 
# for the list of available model run dates.

mydate='20180503'  # this needs to be a recent data since only recent models are available for download
url='http://nomads.ncep.noaa.gov:9090/dods/wave/nww3/nww3'+ mydate+'/nww3'+mydate+'_00z'

file = netCDF4.Dataset(url)
lat  = file.variables['lat'][:]
lon  = file.variables['lon'][:]
data = file.variables['htsgwsfc'][1,:,:]
file.close()

m=Basemap(projection='mill',lat_ts=10,llcrnrlon=lon.min(), \
  urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(), \
  resolution='c')



# convert the lat/lon values to x/y projections.

x, y = m(*np.meshgrid(lon,lat))

# plot the field using the fast pcolormesh routine 
# set the colormap to jet.

m.pcolormesh(x,y,data,shading='flat',cmap=plt.cm.jet)
m.colorbar(location='right')
m.contourf(x,y,data)
# Add a coastline and axis values.

m.drawcoastlines()
m.fillcontinents()
m.drawmapboundary()
m.drawparallels(np.arange(-90.,90.,30.),labels=[1,0,0,0])
m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])















# Convert matplotlib contour to geojson
geojsoncontour.contour_to_geojson(
    contour=contour,
    geojson_filepath='gfs.geojson',
    min_angle_deg=10.0,
    ndigits=3,
    unit='m'
)
'''