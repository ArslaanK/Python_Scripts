# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 21:37:42 2020

@author: Arslaan Khalid
"""


import os
os.environ['PROJ_LIB'] = r'C:\Users\Arslaan Khalid\AppData\Local\conda\conda\pkgs\proj4-5.2.0-ha925a31_1\Library\share'

## loading libraries
import time
from netCDF4 import Dataset
import pandas as pd
import matplotlib.pyplot as plt
import netCDF4 as nc4
import fiona
#import gdal
from shapely.geometry import shape,mapping, Point, Polygon, MultiPolygon
import matplotlib.path as mpltPath
import numpy as np

#-- optional for basemap in plots
from mpl_toolkits.basemap import Basemap


#----------------- USER INPUTS ---------------------------

Index_extraction_type = 'AOI' # 'Statewide' or 'AOI' also called Area of Interest

State="Potomac1" # if statewide, please choose the Full Name of state, i.e. Virginia, Maryland, etc. 
                 # if AOI, please choose the name of specific polygon

# if AOI, please create a shapefile surrounding your AOI, and add two attributes, FID and NAME, That name will be used below
root_path=r'C:\Users\Arslaan Khalid\Desktop\Papers_by_Arslaan\Paper2\Potomac_version2\NWM\Summer_2019-master\Summer_2019-master\NWM\utils'


if Index_extraction_type =='Statewide' :
    Shapefile_url=r"{}\Political_Boundaries_Area\Political_Boundaries_Area.shp".format(root_path)
    Shapefile_name="Political_Boundaries_Area"


elif Index_extraction_type =='AOI' :
    Shapefile_url=r"{}\awein\AOI_potomac.shp".format(root_path)
    Shapefile_name="AOI_potomac" # own created file

else:
    print('Choose the right type')
    
NWM_geo_file = r'{}\nwm-v1.2-channel_spatial_index.nc'.format(root_path)



# =============================================================================
#  Extracting NWM points for User specified State
# =============================================================================

# reading the NWM provided raw geo file
nc_fid = Dataset(NWM_geo_file, 'r') 

# reading the important variables to georeference
feature_id0 = nc_fid.variables['feature_id'][:]
feature_lat0 = nc_fid.variables['latitude'][:]
feature_lon0 = nc_fid.variables['longitude'][:]


# Pandas dataframe for easy access to data (renaming data columns)
lats = pd.DataFrame(feature_lat0)
lats=lats.rename(columns={0:'lats'})

lons = pd.DataFrame(feature_lon0)
lons = lons.rename(columns={0:'lons'})

ids =  pd.DataFrame(feature_id0)
ids = ids.rename(columns={0:'ID'})

# merging all the 3 data columns with same index column
merged0 = pd.concat([lats,lons,ids], axis=1)

# merging lat and lon data columns with same index column
merged1 = pd.concat([lons,lats], axis=1)

# convert the column data into list
USall = merged1.values
list_coordinates = USall.tolist()

# for plotting all the NWM stream IDs (US)
#plt.plot(merged0['lons'],merged0['lats'],marker='o',linestyle='None',c='k',markersize=1)

# =============================================================================
#  Using shapefile for finding the NWM points in the state polygon
# =============================================================================


#---Reading the records in the shapefile
multipolUS = fiona.open(Shapefile_url)

#--- locating the polygon ID of the record
for i in range(len(list(multipolUS))):
    if multipolUS[i]['properties']['NAME'] == State:
        print('found State at polygon ID:',i)
        break
    

multi= multipolUS[i] # only one feature in the shapefile

print(f"Please verify if the state printed here is the one requested: {multi['properties']['NAME']}")

#---extracting the bounding coordiantes of the polygon
bound = multi['geometry']['coordinates'][0]

# using the bounding coordinates to find with NWM points exist in the state polygon 
path = mpltPath.Path(bound)
inside = path.contains_points(list_coordinates)
#print(inside)

#  converting the boolean array data to get the indexes of the Stream IDs in the polygon
StreamIDs = pd.DataFrame(np.where(inside)[0],index=None)
StreamIDs = StreamIDs.set_index(StreamIDs[0])
StreamIDs_indx = StreamIDs[0].values.tolist()

# writing the indexes to a file
if Index_extraction_type =='Statewide' :
    StreamIDs.to_csv(r'state_indexes.csv')
elif Index_extraction_type =='AOI':
    StreamIDs.to_csv(r'state_indexes_AOI.csv')


# concating the 3 initial dataframes again but now based on the state indexes found
# this yields the points for the states only
mergedfinal1 = pd.concat([lats,lons,ids], axis=1,join_axes=[StreamIDs.index])


# plotting the data for the state only
plt.figure()
plt.plot(mergedfinal1['lons'],mergedfinal1['lats'],marker='o',linestyle='None',c='k',markersize=1)

IDnames = mergedfinal1['ID'].values.tolist()

 # if you want to plot the Streams IDs
#for k in range(len(mergedfinal1)-1):
#    print(k)
#    plt.text(mergedfinal1['lons'].iloc[k],mergedfinal1['lats'].iloc[k],mergedfinal1['ID'].iloc[k],color='m')

plt.title(f'NWM stream ID locations for State:{State}\nTotal NWM points for state:{len(StreamIDs)}')

plt.show()
# =============================================================================
#  Initializing the basemap for the better visuals of the georeferenced plots
# =============================================================================

# uncomment if you want to improve the plotting of NWM points

#map = Basemap(projection='cyl', 
#              lat_0=0, lon_0=0,llcrnrlat=mergedfinal1['lats'].min()-1,urcrnrlat=mergedfinal1['lats'].max()+1,llcrnrlon=mergedfinal1['lons'].min()-1,urcrnrlon=mergedfinal1['lons'].max()+1,resolution='h', epsg = 4326)
#
#map.plot(mergedfinal1['lons'],mergedfinal1['lats'],marker='o',linestyle='None',c='k',markersize=0.5)
#
#map.bluemarble()
##map.drawcoastlines(linewidth=0.5)
#map.readshapefile(f'{Shapefile_url[:-4]}', f'{Shapefile_name}', drawbounds=True, zorder=None, linewidth=1.5, color='r', antialiased=1, ax=None, default_encoding='utf-8')

#
## =============================================================================
##  Geocoding the NWM raw files
## =============================================================================
#
## the list of indexes found for the user defined state
#
## reading the indexes file created earlier
#
#StreamIDs = pd.read_csv(r'state_indexes.csv', index_col=0)
#get_col_name = StreamIDs.columns[0]
#StreamIDs_indx = StreamIDs[get_col_name].values.tolist()
#
#
#
## reading the raw downloaded file
#import urllib.request
#import urllib.error
#
#
#
#
#fullfilename = 'C:\\Users\\arslaan.khalid\\Desktop\\Temp2\\NOMADs_Downloads\\Analysis_and_Assimilation\\nwm.20190709.t00z.analysis_assim.channel_rt.tm00.conus.nc'
#url = 'http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.20190709/analysis_assim/nwm.t00z.analysis_assim.channel_rt.tm00.conus.nc'
#urllib.request.urlretrieve(url, fullfilename)
#
#
#
#print(NWM_raw_file)
#
#start=time.time()
#
#NWM_raw_file = fullfilename
#nc_fid1 = Dataset(NWM_raw_file, 'r') 
#
## initializing a new netcdf file to write state only data
#
#f = nc4.Dataset('sample4.nc','w', format='NETCDF4') #'w' stands for write
#
## setting the file descriptive from the raw NWM file
#f.setncatts(nc_fid1.__dict__)
#
## creating the dimensions in the new netcdf
## note that the new size od the 1 d array is the same as the total number of NWM points in a State
#f.createDimension('feature_id', len(StreamIDs_indx))
##f.createDimension('longitude', len(mergedfinal1['lons']))
##f.createDimension('latitude', len(mergedfinal1['lats']))
#
## reading the variables from the raw NWM file
#ncvar = nc_fid1.variables['streamflow']
##ncvar1 = nc_fid1.variables['velocity']
#feature_crs1 = nc_fid1.variables['crs']
#
## reading the variables from the geo NWM file
#feature_id0 = nc_fid.variables['feature_id']
#feature_lat0 = nc_fid.variables['latitude']
#feature_lon0 = nc_fid.variables['longitude']
#
#
## using the State wide index extracted earlier to extract the state only data from the full NWM dataset
#feature_id1 = nc_fid1.variables['feature_id'][StreamIDs_indx]
#stream_flow1 = nc_fid1.variables['streamflow'][StreamIDs_indx]
##velocity1 = nc_fid1.variables['velocity'][StreamIDs_indx]
#
##way1
#lats1 = mergedfinal1['lats'].values.tolist()
#lons1 = mergedfinal1['lons'].values.tolist()
### alternate
##lats1 = nc_fid.variables['latitude'][StreamIDs_indx]
##lons1 = nc_fid.variables['longitude'][StreamIDs_indx]
#
## creating the new variables in the newly created georef file 
#var0 = f.createVariable('crs',feature_crs1.dtype,feature_crs1.dimensions)
#var1 = f.createVariable('streamflow',ncvar.dtype,ncvar.dimensions)
##var2 = f.createVariable('velocity',ncvar1.dtype,ncvar1.dimensions)
#var3 = f.createVariable('longitude',feature_lon0.dtype,feature_lon0.dimensions)
#var4 = f.createVariable('latitude',feature_lat0.dtype,feature_lat0.dimensions)
#var5 = f.createVariable('feature_id',feature_id0.dtype,feature_id0.dimensions)
#
## setting the variable description to the new netcdf using the raw NWM file
#var0.setncatts(feature_crs1.__dict__)
#var1.setncatts(ncvar.__dict__)
##var2.setncatts(ncvar1.__dict__)
#var3.setncatts(feature_lon0.__dict__)
#var4.setncatts(feature_lat0.__dict__)
#var5.setncatts(feature_id0.__dict__)
#
## assigining variable data to the newly created variables
#
##len(stream_flow1)
#var1[:]=stream_flow1[:]
##var2[:]=velocity1[:]
#var3[:]=lons1[:]
#var4[:]=lats1[:]
#var5[:]=feature_id1[:]
#
## closing the newly created file
#f.close()
#
#end=time.time()
#totaltime=round((end-start),2)
#
#
#print("\nProcessing time:",totaltime,"secs")
#
#
#
#
#current = 0.42
#previous = 1.86
#
#
#def get_change(current, previous):
#    if current == previous:
#        return 100.0
#    try:
#        return (abs(current - previous) / previous) * 100.0
#    except ZeroDivisionError:
#        return 0
#
#get_change(current,previous)
#
#
## =============================================================================
##  Opening the newly created netcdf file
## =============================================================================
#
#
#nc_out1 = Dataset(r'sample4.nc', 'r') 
#
#
#nc_out1.variables
#
#id_2 = nc_out1.variables['feature_id'][:]
#lat_2 = nc_out1.variables['latitude'][:]
#lon_2 = nc_out1.variables['longitude'][:]
#
#
#plt.figure()
#plt.plot(lon_2,lat_2.data,marker='o',linestyle='None',c='k',markersize=1)
#
#
#
#
## =============================================================================
##  Method 2 for georeferencing is sloppy
##  Note:  I am trying to just make the changes in the Raw NWM file rather creating a new file 
## =============================================================================
##
###--- cropping data from global to state --------------
## 1. load the NWM raw file
## 2. use the State wide index to keep the state wide data only
## 3. Somehow the shape of the 2.7 million needs to be changed to match the total number of points in
##    the state
## 4. create the lat and lon variable in the same file
##
##nc_fid1.variables['feature_id'][:] = nc_fid1.variables['feature_id'][StreamIDs_indx]
#
#
#
#
## =============================================================================
##  convert to function
## =============================================================================
#
#raw1=r'C:\Users\arslaan.khalid\Desktop\test1\NWMArslaanK\temp\NOMADs_Downloads\Analysis_and_Assimilation\nwm.20190701.t00z.analysis_assim.channel_rt.tm00.conus.nc'
#geo1 = r'sample4.nc'
#georef_streamflows(raw1,geo1)
#
#
#import time
#import netCDF4 as nc4
#
#
#
#
#import time
#import netCDF4 as nc4
#import pandas as pd
#
##-loading the index files
#
#StreamIDs = pd.read_csv(r'state_indexes.csv', index_col=0)
#get_col_name = StreamIDs.columns[0]
#StreamIDs_indx = StreamIDs[get_col_name].values.tolist()
#
#
#def georef_streamflows(raw,out,geo):
#
#    start=time.time()    
#
#    # # reading the NWM provided geofile
#   
#    nc_fid = nc4.Dataset(geo, 'r') 
#
#    # reading the raw file   
#    nc_fid1 = nc4.Dataset(raw, 'r')
#    # initializing a new netcdf file to write state only data
#    f = nc4.Dataset(out,'w', format='NETCDF4') #'w' stands for write
#    
#    # setting the file descriptive from the raw NWM file
#    f.setncatts(nc_fid1.__dict__)
#    
#    # creating the dimensions in the new netcdf
#    # note that the new size od the 1 d array is the same as the total number of NWM points in a State
#    f.createDimension('feature_id', len(StreamIDs_indx))
#    
#    # reading the variables from the raw NWM file
#    ncvar = nc_fid1.variables['streamflow']
#    #ncvar1 = nc_fid1.variables['velocity']
#    feature_crs1 = nc_fid1.variables['crs']
#    
#    # reading the variables from the geo NWM file
#    feature_id0 = nc_fid.variables['feature_id']
#    feature_lat0 = nc_fid.variables['latitude']
#    feature_lon0 = nc_fid.variables['longitude']
#    
#    # using the State wide index extracted earlier to extract the state only data from the full NWM dataset
#    feature_id1 = nc_fid1.variables['feature_id'][StreamIDs_indx]
#    stream_flow1 = nc_fid1.variables['streamflow'][StreamIDs_indx]
#    #velocity1 = nc_fid1.variables['velocity'][StreamIDs_indx]
#    lats1 = mergedfinal1['lats'].values.tolist()
#    lons1 = mergedfinal1['lons'].values.tolist()
#    
#    # creating the new variables in the newly created georef file 
#    var0 = f.createVariable('crs',feature_crs1.dtype,feature_crs1.dimensions)
#    var1 = f.createVariable('streamflow',ncvar.dtype,ncvar.dimensions)
#    #var2 = f.createVariable('velocity',ncvar1.dtype,ncvar1.dimensions)
#    var3 = f.createVariable('longitude',feature_lon0.dtype,feature_lon0.dimensions)
#    var4 = f.createVariable('latitude',feature_lat0.dtype,feature_lat0.dimensions)
#    var5 = f.createVariable('feature_id',feature_id0.dtype,feature_id0.dimensions)
#    
#    # setting the variable description to the new netcdf using the raw NWM file
#    var0.setncatts(feature_crs1.__dict__)
#    var1.setncatts(ncvar.__dict__)
#    #var2.setncatts(ncvar1.__dict__)
#    var3.setncatts(feature_lon0.__dict__)
#    var4.setncatts(feature_lat0.__dict__)
#    var5.setncatts(feature_id0.__dict__)
#    
#    # assigining variable data to the newly created variables
#    
#    #len(stream_flow1)
#    var1[:]=stream_flow1[:]
#    #var2[:]=velocity1[:]
#    var3[:]=lons1[:]
#    var4[:]=lats1[:]
#    var5[:]=feature_id1[:]
#    
#    # closing the newly created file
#    f.close()
#    
#    end=time.time()
#    totaltime=round((end-start),2)
#    
#    
#    print("\nProcessing time:",totaltime,"secs")
#    
#    return