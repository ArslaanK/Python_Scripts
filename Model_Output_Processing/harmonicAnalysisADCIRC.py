# -*- coding: utf-8 -*-
"""
Created on Sun Aug 11 23:21:38 2024

@author: Arslaan.Khalid
"""

   
# =============================================================================
# loading libraries
# =============================================================================
    
import netCDF4 as nc
import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup
import matplotlib.pyplot as plt
import utide

print(utide.__version__)


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


def solve_tide(sample_tide,  col_name= ' Prediction', lat=38, method="ols", conf_int="MC", verbose=False):
    
    # Usage
    # result = solve_tide(sample_tide)
    # print(result)
    
    # Solve for tidal constituents
    coef = utide.solve(
        sample_tide.index,
        sample_tide[col_name],
        lat=lat,
        method=method,
        conf_int=conf_int,
        verbose=verbose,
    )

    # Extract names, amplitude, and phase
    names = coef['name']
    amplitude = coef['A']
    phase = coef['g']

    # Create DataFrame
    model_tide = pd.DataFrame([names, amplitude, phase]).T
    model_tide = model_tide.set_index(model_tide[0])
    del model_tide[0]

    # Rename columns and index
    model_tide.columns = ['Amplitude', 'Phase']
    model_tide.index.names = ['Constituent Name']

    # Convert to float and round values
    model_tide = model_tide.astype(float)
    model_tide.Amplitude = model_tide.Amplitude.round(3)
    model_tide.Phase = model_tide.Phase.round(1)

    return model_tide


def fetch_noaa_tide_data(GageID, units='metric'):
    # Usage
    # GageID = '8594900'
    # units = 'metric'
    # df = fetch_noaa_tide_data(GageID, units)
    # print(df)
    
    # Determine the unit code
    if units == 'english':
        code_ = 1
    elif units == 'metric':
        code_ = 0
    else:
        raise ValueError("Units must be either 'english' or 'metric'")

    # Construct the URL
    url_harmonic = f'https://tidesandcurrents.noaa.gov/harcon.html?unit={code_}&timezone=0&id={GageID}'

    # Fetch the HTML content
    response = requests.get(url_harmonic)
    html_content = response.text

    # Parse the HTML content using BeautifulSoup
    soup = BeautifulSoup(html_content, 'html.parser')

    # Extract the table containing harmonic constituents
    table = soup.find('table')

    # Extract table headers
    headers = [header.text.strip() for header in table.find_all('th')]

    # Extract table rows
    rows = []
    for row in table.find_all('tr')[1:]:  # Skip the header row
        cells = row.find_all('td')
        row_data = [cell.text.strip() for cell in cells]
        rows.append(row_data)

    # Create DataFrame
    df_noaa_reported = pd.DataFrame(rows)
    df_noaa_reported.columns = headers

    # Convert Amplitude and Phase to float
    df_noaa_reported['Amplitude'] = df_noaa_reported['Amplitude'].astype(float)
    df_noaa_reported['Phase'] = df_noaa_reported['Phase'].astype(float)

    # Set index to 'Name'
    df_noaa_reported = df_noaa_reported.set_index(df_noaa_reported['Name'])
    del df_noaa_reported['Name']

    return df_noaa_reported



#%%


your_stations_File = 'locations.txt'

# this can be a txt file with 4 columns ['FullName','Abbreviaton','lat','lon','StationID'] , see below
#Fullname	lon	 lat		Abbreviaton	StationID
#Kiptopeke, VA	-75.99539591	37.16343353		KPTV	8632200	
#S1_t2	-75.940369	37.149697		S1T2	8111111                   
#Note: if your station has a NOAA ID use that, i.e. Kiptopeke, VA is 8632200 
# else give it some random number i.e. S1_t2 is 8111111    


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

# ====================

#%%

# =============================================================================
# # Example usage for netcdf
# =============================================================================


file_path = r'path_to_tidal_fort.63.nc'


index_of_station = 55500 # locate these based on some code or SMS for all NOAA stations
dataset = nc.Dataset(file_path)

time = dataset.variables['time'][:]
time_axis = nc.num2date(time, units=time.units) 
# time = nc.num2date(time, units='seconds since 1970-01-01 00:00:00') # change this based on model data initialization


water_level = dataset.variables['zeta'][:,index_of_station].data


lat  = dataset.variables['y'][:]
lon  = dataset.variables['x'][:]    


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
    data1 = dataset.variables['zeta'][:,Node_ID] 
    
    Timeseries['{}'.format(Node_ID)] = data1.data    
                                        
Timeseries.replace(to_replace=-99999.000000,value=np.nan,inplace=True)
TimeseriesT = Timeseries.set_index(time_axis)


# now loop over each stationD in station_nodes and get the tidal constituents

for stationID in station_nodes.keys():
    water_ = TimeseriesT[stationID]
    
    # using utide for tide decomposition    
    model_tides = solve_tide (water_,  col_name=stationID, lat=38, method="ols", conf_int="MC", verbose=False)
        
        



#%%

# =============================================================================
# Since no adcirc model output yet, I used the NOAA tide itself to get the amplitude and tidal consitutents
# =============================================================================

units = 'metric' # english or metric
start_date = '20240101'
end_date  = '20241231'
GageID = '8594900'

url_noaa = f'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?product=predictions&application=NOS.COOPS.TAC.WL&begin_date={start_date}&end_date={end_date}&datum=MSL&station={GageID}&time_zone=GMT&units={units}&interval=h&format=csv'

sample_tide = pd.read_csv(url_noaa)

sample_tide = sample_tide.set_index(sample_tide['Date Time']);
sample_tide.index = pd.to_datetime(sample_tide.index)
del sample_tide['Date Time']

model_tide = solve_tide(sample_tide,  col_name= ' Prediction', lat=38, method="ols", conf_int="MC", verbose=False)
    
print(model_tide)




#%%

# =============================================================================
# extracting NOAA tides
# =============================================================================
GageID = '8594900'
units = 'metric' # english or metric

df_noaa_reported = fetch_noaa_tide_data(GageID, units)
print(df_noaa_reported)


#%%

# =============================================================================
# # comparing model and NOAA reported tidal constituents
# =============================================================================

# amplitude
df_amplitude = pd.DataFrame(df_noaa_reported['Amplitude'])
df_amplitude.columns = ['NOAA']
df_amplitude['Predicted'] = model_tide['Amplitude']


# phase 
df_phase = pd.DataFrame(df_noaa_reported['Phase'])
df_phase.columns = ['NOAA']
df_phase['Predicted'] = model_tide['Phase']

# dropping Nans
df_amplitude = df_amplitude.dropna()
df_phase = df_phase.dropna()




#%%


# =============================================================================
# create a scatter plot
# =============================================================================

predicted_amplitudes = df_amplitude.NOAA.values
model_amplitudes = df_amplitude.Predicted.values

predicted_phases  = df_phase.NOAA.values
model_phases  = df_phase.Predicted.values


# Calculate 10% error bars
amplitude_errors = 0.1 * np.array(predicted_amplitudes)
phase_errors = 0.1 * np.array(predicted_phases)

# scatter plot for amplitudes with error bars
plt.figure(figsize=(6, 6))
plt.errorbar(model_amplitudes, predicted_amplitudes, yerr=amplitude_errors, fmt='o', color='blue', label='Amplitudes')
plt.plot([0, max(model_amplitudes)], [0, max(predicted_amplitudes)], color='k', linestyle='-',lw=0.5)
plt.xlabel('Model Amplitudes')
plt.ylabel('Predicted Amplitudes')
plt.title('Comparison of Model and Predicted Amplitudes with\n 10% Error Bars')
plt.legend()
plt.grid(True)
plt.show()

# scatter plot for phases with error bars
plt.figure(figsize=(6, 6))
plt.errorbar(model_phases, predicted_phases, yerr=phase_errors, fmt='o', color='green', label='Phases')
plt.plot([0, 360], [0, 360], color='k', linestyle='-',lw=0.5)
plt.xlabel('Model Phases')
plt.ylabel('Predicted Phases')
plt.title('Comparison of Model and Predicted Phases with\n 10% Error Bars')
plt.legend()
plt.grid(True)
plt.show()







