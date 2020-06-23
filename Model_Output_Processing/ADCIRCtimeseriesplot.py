

# -*- coding: utf-8 -*-
"""
Created on Sat Jun 23 16:13:19 2020

@author: Arslaan Khalid
# for queries please contact @ akhalid6@gmu.edu
# following script prepare the timeseries plots from adcirc netcdf files
"""

# =============================================================================
# Loading Libraries
# =============================================================================


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

# =============================================================================
# Loading Functions
# =============================================================================

def adcirc_nodalIdentifier(keys,lat,lon):
    chosen_lat = database_lats[keys][0]   
    chosen_lon =database_lons[keys][0]
       
    min_distance = None;best_index = 0
    
    for i in range(len(lat)):
        current_distance = (lat[i] - chosen_lat)**2 + (lon[i] - chosen_lon)**2
        if min_distance is None or current_distance < min_distance:
            best_index = i
            min_distance = current_distance    
    return best_index

def adcirc_windDirection(u,v):
  #u = list
  #v = list
  windDir = []  
  if (v > 0):
      windDir = ((180 / pi) * np.arctan(u/v) + 180)
  if (u < 0 and v < 0):
      windDir =((180 / pi) * np.arctan(u/v) + 0)
  if (u > 0 and v < 0):
      windDir = ((180 / pi) * np.arctan(u/v) + 360)    

  return windDir


def convert_vdatum(latitude: float, longitude: float, elevation: float, ivert: str, ihorz: str = 'NAD83_2011', iunit: str = 'm', geoid: str = 'geoid12b', overt: str = 'NAVD88', ounit: str = 'm'):
    url_str = f'https://vdatum.noaa.gov/vdatumweb/api/tidal?lon={longitude}&lat={latitude}' \
              + f'&height={elevation}&s_h_frame={ihorz}&s_v_frame={ivert}&s_v_unit={iunit}' \
              + f'&t_v_frame={overt}&t_v_unit={ounit}'
    vd_value = requests.get(url = url_str).content.decode()
    new = pd.read_json(vd_value, lines = True)
    if 'errorCode' in new.columns:
    	new.insert(0,'tar_height','NaN')
    return new


# =============================================================================
#  USER INPUTS
# =============================================================================

start_timestamp="2017091619"
end_timestamp="2017092914"

Input_dir = 'path to input files'
Out_dir = 'path to processed plots'
your_stations_File = 'path to stations.txt file'
                     # this can be a txt file with 4 columns ['FullName','Abbreviaton','lat','lon','StationID'] , see below
                     #Fullname	lon	 lat		Abbreviaton	StationID
                     #Kiptopeke, VA	-75.99539591	37.16343353		KPTV	8632200	
                     #S1_t2	-75.940369	37.149697		S1T2	8111111                   
                     #Note: if your station has a NOAA ID use that, i.e. Kiptopeke, VA is 8632200 
                     # else give it some random number i.e. S1_t2 is 8111111    



plotting_stations_list=[8571421,8531680,8111111,8111112,8111113,8111114] # add the stations to view on a 3x3 panel [i.e. 3x3 will need 9 stations]
paneltype = [2,3] # [columns, rows] [i.e. 3x3 will need 9 stations]

plot_type= 'Winds' # choose from these [Winds,Pressure,WaterLevels,Currents]
convert_yourmodelreults_to_datum='LMSL' # if you want to change your model results datum , valid for water levels only


#plotting selections [for external data]
use_NOAA_observed=True # True or False
date_type="wind"      # list of products you might be able to get
                                # "water_level" for observed water levels, maximum download data limit is 30 days
                                # "wind" for observed wind speeds and direction
                                #	air_pressure" for observed air pressure
                                #	"hourly_height" for observed water levels, maximum download data limit is 365 days
                                #	"predictions" for NOAA predicted tides, maximum download data limit is x days
                                #	 "currents" for observed water currents
	      
col_name = 'Speed'       # this must be known in order to use NOAA data, use from below
                                #'Water Level' for observed water levels 
                                #'Prediction' for NOAA predicted tides 
                                #'Speed' for observed wind speeds and direction
                                #'Pressure' for observed air pressure

plt_ylims= [0,20]       # add the y limits of plots (examples below)
                        # for water may be [-2,2]
                        # for winds may be [0,20]
                        # for currents may be [0,2]
                        # for pressure may be [900,1050]
#default_NOAA_waterlevel_datum= 'MSL' # hard coded in the main code, can be changed there

myxlabeldatesFormat = mdates.DateFormatter('%m/%d')
Figur_title='Hurricane Maria & Jose [18Sep-29Sep 2017]'

# =============================================================================
#   MAIN CODE
# =============================================================================
# =============================================================================
#  Assigning input files for Timeseries plots 
# [fort.63.nc,fort.64.nc,fort.74.nc,fort.73.nc]
# =============================================================================

if plot_type =='Winds':
    stdin1 = Input_dir + '//fort.74.nc' # wind file
    variable_plot1,variable_plot2 = 'windx','windy'
elif plot_type =='Pressure':
    stdin1 = Input_dir + '//fort.73.nc' #  pressure file
    variable_plot1 = 'pressure' 
    variable_plot2 = 0
elif plot_type =='WaterLevels':
    stdin1 = Input_dir + '//fort.63.nc'
    variable_plot1 = 'zeta'
    variable_plot2 = 0
elif plot_type =='Current':
    stdin1 = Input_dir + '//fort.64.nc'
    variable_plot1,variable_plot2 = 'u-vel','v-vel'

# =============================================================================
# reading the stations file
# =============================================================================
try:
    stations_input = pd.read_csv(your_stations_File,sep='\t')
    
    database = stations_input.groupby('StationID')['Abbreviaton'].apply(list).to_dict()
    databaseFullnames = stations_input.groupby('StationID')['Fullname'].apply(list).to_dict()
    database_lats = stations_input.groupby('StationID')['lat'].apply(list).to_dict()
    database_lons = stations_input.groupby('StationID')['lon'].apply(list).to_dict()
    print(f'file read correctly!, first3 rows of file\n',stations_input.head(3))
except:
    print('File is not formatted correctly or doesnot have the required columns')    

# =============================================================================
# # read the netcdf parameters
# =============================================================================
    
file1 = netCDF4.Dataset(stdin1) # this read the file
print('File contains following variables:\n',file1.variables)
lat  = file1.variables['y'][:]
lon  = file1.variables['x'][:]    
depth  = file1.variables['depth'][:] 
elems = file1.variables['element'][:,:]-1  # Move to 0-indexing by subtracting 1, elements indexing starts with '1' in netcdf file         

time_var_re = file1.variables['time'] # read the datetime
# sometimes the format isn't good so did this
mn_num = datetime.strptime(f"{time_var_re.units.split('-')[1]}",'%b').month
yr_num = f"{time_var_re.units.split(' ')[2][-4:]}"
dt_num = f"{time_var_re.units.split(' ')[2][:2]}"

Time_data_re = netCDF4.num2date(time_var_re[:],f"seconds since {yr_num}-{mn_num}-{dt_num} 00:00:00 +00:00")


# =============================================================================
# looping through all the stationsID to extract node number
# =============================================================================
station_nodes = {}
start1 = time.time()
# this finds index of your stations from the numerical mesh file
for keys in database:
   best_index = adcirc_nodalIdentifier(keys,lat,lon) 
   print("best_index is :{} for station {} ".format(best_index,databaseFullnames[keys][0]))
   station_nodes['{}'.format(keys)] = best_index 

end1 = time.time()
print('Finished finding node numbers after {} seconds'.format(end1 - start1))
  
#--- loading only the data for the stations in my AOI
Timeseries = pd.DataFrame()
  
for keys in database:   
    Node_ID = station_nodes[str(keys)]
    print('loading station data : ',database[keys])               
    data1 = file1.variables[variable_plot1][:,Node_ID] 
    if variable_plot2 != 0:
        data2 = file1.variables[variable_plot2][:,Node_ID]
        magnitude = np.sqrt(np.square(data1)+np.square(data2))
        Timeseries['{}'.format(Node_ID)] = magnitude.data 
    else:        
        Timeseries['{}'.format(Node_ID)] = data1.data                                            
Timeseries.replace(to_replace=-99999.000000,value=np.nan,inplace=True)

# this sets the datetime index 
start_date = datetime.strptime(start_timestamp,'%Y%m%d%H')
end_date = datetime.strptime(end_timestamp,'%Y%m%d%H')
# this sets the date teime index (you may not may not need timedelta(hours=1))
time_index = pd.date_range(start_date+timedelta(hours=1),periods = len(Timeseries), freq='h') # use 'h' if hourly or '30Min' if every 30 mins
TimeseriesT = Timeseries.set_index(time_index)


# =============================================================================
# 
# Download the NOAA stations data (water level) [Note: can be tuned to download wind, pressure and other datasets on NOAA tides and currents website]
# 
# =============================================================================


if use_NOAA_observed == True:
    NOAA_data=pd.DataFrame()
    NOAA_data_otherdatum=pd.DataFrame()
    
    end_date1 = datetime.strftime(end_date,'%Y%m%d%H')
    start_date1 = datetime.strftime(start_date,'%Y%m%d%H')
                                    
    plotting_list=[]
    
    for keys in database:
        try:
            keys=int(keys)
            if plot_type =='WaterLevels':
                temp_observed_NOAA_url = 'https://tidesandcurrents.noaa.gov/api/datagetter?product={}&application=NOS.COOPS.TAC.WL&begin_date={}&end_date={}&datum=MSL&station={}&time_zone=GMT&units=metric&format=csv'.format(date_type,start_date1[:8],end_date1[:8],keys)
                datumcorrection0 = convert_vdatum(database_lats[keys][0], database_lons[keys][0], 0,'LMSL') # this converts from MSL to NAVD88(navd88 is hard coded in the function above, change there if need conversion to other formats)
                datumcorrection = datumcorrection0['tar_height'][0] 
            else:
                temp_observed_NOAA_url = 'https://tidesandcurrents.noaa.gov/api/datagetter?product={}&application=NOS.COOPS.TAC.WL&begin_date={}&end_date={}&station={}&time_zone=GMT&units=metric&format=csv'.format(date_type,start_date1[:8],end_date1[:8],keys)
    
            temp_observed_NOAA_reterival = pd.read_csv(temp_observed_NOAA_url)
            temp_observed_NOAA = pd.DataFrame(temp_observed_NOAA_reterival[f' {col_name}'])
            temp_observed_NOAA = temp_observed_NOAA.set_index(temp_observed_NOAA_reterival['Date Time'])
            temp_observed_NOAA.index = pd.to_datetime(temp_observed_NOAA.index)
            #temp_observed_NOAA_hourly = temp_observed_NOAA.resample('H').mean() 
    
       
            #exec("Node_ID = station_nodes1['{}']".format(keys))
            Node_ID=keys  # replacing Node_ID by keys
            
            NOAA_data['{}'.format(Node_ID)] = temp_observed_NOAA[f' {col_name}'] 
            if plot_type =='WaterLevels':NOAA_data_otherdatum['{}'.format(Node_ID)] = temp_observed_NOAA[f' {col_name}'] + datumcorrection
            print('NOAA data loaded for station:',keys, '[:)]')
            plotting_list.append(keys)
        except:
            
            Node_ID=keys  # replacing Node_ID by keys
            NOAA_data['{}'.format(Node_ID)] = [np.nan]*len(NOAA_data)
            if plot_type =='WaterLevels':NOAA_data_otherdatum['{}'.format(Node_ID)] = [np.nan]*len(NOAA_data)
            print('(X)NOAA data not available for station:',keys)        


# =============================================================================
# 
# =============================================================================
# initialize a figure

fig = plt.figure(figsize = (9,12))
fig.subplots_adjust(wspace=0.1)

# do the loop for each subplot
for y in range(0,(paneltype[0]*paneltype[1])):  
    
    keys=plotting_stations_list[y]


    plt.subplot(paneltype[0],paneltype[1],1+y) 
    
    # =============================================================================
    #  plotting scheme
    # =============================================================================

    Node_ID=str(station_nodes[f'{keys}'])

    if Node_ID in NOAA_data.columns:
        plt.plot(NOAA_data.index, NOAA_data[Node_ID],label='NOAA',c='b') # Note: use this if NOAA data is downloaded correctly
 
    plt.plot(TimeseriesT[Node_ID].index,TimeseriesT[Node_ID],label='modeled',c='g')

    plt.grid(True, which='both', linestyle='--',color='gray',linewidth=0.2)
    plt.xlabel('Date Time')
    plt.ylabel('Water level in meters[NAVD88]',fontsize=8)
    
    plt.title(f'{keys},{database[keys][0]}')

    plt.xlim(TimeseriesT.index[0],TimeseriesT.index[-1])
    
    plt.ylim(plt_ylims[0],plt_ylims[1])
    plt.gca().xaxis.set_major_formatter(myxlabeldatesFormat)
    #pltaxis.set_major_formatter(myFmt)".format(i))
    if y==0:
        plt.legend(ncol=2,fontsize=7)


# fit subplots and save fig
#fig.tight_layout()

#fig.subplots_adjust(
#top=0.88,
#bottom=0.11,
#left=0.11,
#right=0.9,
#hspace=0.28,
#wspace=0.0)


fig.suptitle(Figur_title)

fig.set_size_inches(w=12,h=12)

fig.savefig(r'{}\sample_timeseries_{}.png'.format(Out_dir,plot_type),  bbox_inches = 'tight',dpi = 100, pad_inches = 0.2)
print(f'Successfully created the Timeseries Figure for Type: {plot_type}')

