
Documentation = '''

Map Visualization:
    The main dashboard includes a map (using Mapbox) where incidents reported by VDOT are displayed.
    Users can click on points on the map to get more information about the incidents.

Time Filtering:
    Users can select a specific time period using a range slider, allowing them to focus on incidents reported within a particular timeframe. 

Region Selection:
    The dashboard supports the selection of specific regions using buttons.
    Each button represents a region, and clicking on a button filters the displayed incidents based on the chosen region.

County Filtering:
    Users can select specific counties from a dropdown menu, narrowing down the incidents displayed on the map.

Data Label Selection:
    Users can choose different data labels (e.g., 'County', 'Category', 'Road Type') to categorize and visualize incidents.

Event Category and Type Filtering:
    Users can select specific event categories and types using dropdown menus, allowing them to filter incidents based on these criteria.

Road Network Visualization:
    The code includes the visualization of the road network near the selected incidents.
    Roads are displayed on the map, and users can click on them to get more information.

Graphical Plots:
    The dashboard includes various plots such as bar plots and pie charts displaying the distribution of incidents based on selected criteria.

USGS Stream and Meteorological Data:
    Users can synchronize USGS stream and meteorological data with the incidents on the map.
    Clicking on the checkbox allows users to display stream gages and meteorological gages on the map.

Download Data:
    Users can download the displayed data in CSV format using a download button.

Responsive Design:
    The dashboard is designed using Dash and Bootstrap, ensuring a responsive and user-friendly layout.

Tabs for Additional Information:
    The dashboard includes tabs for 'About,' 'Alerts,' and 'Contact Us,' providing additional information and functionality.

Alerts Signup Page:
    There is a tab dedicated to an "Alerts Signup Page," allowing users to select alerts for a given location or road.

Contact Us Page:
    A "Contact Us" tab provides information for users to get in touch with the dashboard administrators.

'''





from dash import Dash, html, dcc, Input, Output,callback_context
import plotly.express as px
import dash_bootstrap_components as dbc
import geopandas as gpd
import dash_ag_grid as dag
import pandas as pd
import json
import dataretrieval.nwis as nwis
from shapely.geometry import Point, LineString
import plotly.graph_objects as go
from dash.dependencies import Input, Output, State
import dash_mantine_components as dmc
from dash import html, Output, Input, callback
import dash
from dash.dependencies import ALL
import tqdm
import numpy as np
import shapely.geometry
import os
import datetime

import warnings
warnings.filterwarnings("ignore")


current_working_directory = os.getcwd()

# data_root = f'{current_working_directory}//Data_Clean'

data_root = 'https://data.iflood.vse.gmu.edu/VDOT_dataset'

output_dir = f'{current_working_directory}//Outputs'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)


# data_root = r'C:/Users/Arslaan Khalid/Desktop/Papers_by_Arslaan/VDOT_Dashboard/Data'
# data_root = current_working_directory


file_names = { 'VDOT closures': 'road_closures.csv',
               'VDOT Region': 'VDOT_regions.csv',
               'Gages': 'USGS_gagesVA.csv',
               'Roads': 'road_lines_simplified.csv', 
               'Obs Stage': 'usgs_stage_2019_2024_h.csv',
               'Obs Discharge': 'usgs_discharges_2019_2024_h.csv' ,
               'Obs Precipitation': 'usgs_precip_2019_2024_h.csv'  
    }

# loading observations


disc_df = pd.read_csv(f"{data_root}/{file_names['Obs Stage']}")       
stage_df = pd.read_csv(f"{data_root}/{file_names['Obs Discharge']}")          
met_df = pd.read_csv(f"{data_root}/{file_names['Obs Precipitation']}") 

# converting to datetime index
disc_df.index = pd.to_datetime(disc_df['datetime']);del disc_df['datetime']
stage_df.index = pd.to_datetime(stage_df['datetime']);del stage_df['datetime']
met_df.index = pd.to_datetime(met_df['datetime']);del met_df['datetime']


# loading vdot regions


vdot_region1 = pd.read_csv(f"{data_root}/{file_names['VDOT Region']}")
del vdot_region1['Unnamed: 0']

vdot_region1 = gpd.GeoDataFrame(vdot_region1, geometry=gpd.GeoSeries.from_wkt(vdot_region1['geometry']), crs='EPSG:26918')
vdot_region = vdot_region1.to_crs(4326)



# loading USGS Stream data
VA = pd.read_csv(f"{data_root}/{file_names['Gages']}", dtype ='str')
usgs_gages_all = gpd.GeoDataFrame(VA, geometry=gpd.points_from_xy(VA['dec_long_va'], VA['dec_lat_va']), crs='EPSG:4326')

usgs_gages_all['lat'] = usgs_gages_all['geometry'].y
usgs_gages_all['lon'] = usgs_gages_all['geometry'].x

usgs_gages_all['STAID'] = usgs_gages_all.site_no 

usgs_gages_all = usgs_gages_all.to_crs(26918)

usgs_gages = usgs_gages_all[usgs_gages_all['parm_cd'].isin(['00060', '00065'])]
usgs_gages.drop_duplicates(subset='STAID', keep='first', inplace=True)



# loading USGS Met 
usgs_gages_wind = usgs_gages_all[usgs_gages_all.parm_cd=='00035']

met_gages = usgs_gages_all[usgs_gages_all.parm_cd=='00045']

# Read the VDOT data
rd_file = pd.read_csv(f"{data_root}/{file_names['VDOT closures']}")

# rd_file = rd_file[:1000]
# convert to dataframe
gdf = gpd.GeoDataFrame(rd_file, geometry=gpd.points_from_xy(rd_file['x'], rd_file['y']), crs='EPSG:4326')


# loading Roads network
roads_in_2 = pd.read_csv(f"{data_root}/{file_names['Roads']}")
del roads_in_2['Unnamed: 0']
# roads_in_2 = roads_in_2[:1000]

# convert to dataframe
roads_in_2['geometry'] = roads_in_2['geometry'].apply(lambda x: x if pd.notnull(x) else None) 
roads = gpd.GeoDataFrame(roads_in_2, geometry=gpd.GeoSeries.from_wkt(roads_in_2['geometry']), crs='EPSG:4326')
# remove null geometry, if any
roads = roads[roads.geometry.notnull()]

# join all roads with the vdot data, to clean some road loading
roads_near_vdot = gpd.sjoin_nearest(roads, gdf,max_distance=0.0001)


# keep, lat lon columns
gdf['lat'] = gdf['geometry'].y
gdf['lon'] = gdf['geometry'].x

# Convert the geometry column to GeoJSON
gdf['geometry'] = gdf['geometry'].apply(lambda geom: geom.__geo_interface__)

# keeping interested columns
intersted_columns = ['fid', 'containing','containi_1', 'event_cate', 'event_subc','reason_for',  'road_syste',
       'route_name', 'type_event', 'update_tim', 'geometry',
       'lat', 'lon']
gdf = gdf[intersted_columns]

gdf.columns = ['FID', 'County', 'Region','Category','SubCategory', 'Reason Closure', 'Road Type',
       'Route Name', 'Event Type', 'Reported Time', 'geometry',
       'lat', 'lon']

# get all unique regions
available_regions = gdf['Region'].unique()

# convert to datetime
gdf['Reported Time2'] = pd.to_datetime(gdf['Reported Time'], errors='coerce')
gdf = gdf.dropna(subset=['Reported Time2'])
gdf = gdf.sort_values(['Reported Time2'])

# define a global variable for usage later
gdf_global = gdf.copy()

# get unique month and years
unique_dates = pd.to_datetime(gdf['Reported Time2']).dt.to_period('M').unique()

# create the labels for the slider
timewindows = {}

for date in unique_dates:
    timestamp = int(date.to_timestamp().timestamp())
    
    if date == unique_dates[0]:
        label = f'{date.month}' 
    elif date.month == 6:
        label = f'{date.month}' 
    elif date.month == 12:
        label = f'{date.year}' 
    else:
        label = ''
    timewindows[timestamp] = {'label': label}


# slider 1st and last
slider_start = list(timewindows.keys())[0]
slider_end = list(timewindows.keys())[-1]


# create region buttons

region_buttons = [
    dbc.Button(
        region,
        id={'type': 'region-button', 'index': i},
        color='secondary',  # You can change 'primary' to any Bootstrap color class
        style={'color': 'white'},  # Set the text color
        className='mr-1',
    )  
    for i, region in enumerate(available_regions)
]

# initialize a column
initial_column = 'Category'
initial_region = 'Northern'
# create DASH app
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = dbc.Container([
    html.H1("Interactive VDOT Dashboard", className='mb-2', style={'textAlign': 'center'}),
    html.Br(),

    dcc.Tabs([
            dcc.Tab(label='Dashboard', children=[

    dbc.Row([
        dbc.Col([
            dcc.Graph(id='basemap', clickData={'points': [{'pointNumber': 0}]}),  # Initial clickData
        ], width=8),

        dbc.Col([

            #------------ Define Region ---- 
            html.H3("Select Region", className='mb-2', style={'textAlign': 'center'}),

            # Add buttons to select regions
            dbc.Row([
                dbc.Col(button, width=2)  # Set the width to 2 for each button
                for button in region_buttons
            ], className='mt-2'),

            # ------------ time slider
            html.H3("Select Time Period", className='mb-2', style={'textAlign': 'center'}),

            dbc.Row([
                dbc.Col([
                    dcc.RangeSlider(
                        id='time-slider',
                        min=slider_start,
                        max=slider_end,
                        step=None,
                        marks=timewindows
                    )
                ], width=12)
            ]),
            #------------ Define County ---- 
            html.H3("Select County", className='mb-2', style={'textAlign': 'center'}),
            html.Br(),
            # Dropdown for counties
            dbc.Row([
                dbc.Col([
                    dcc.Dropdown(
                        id='county-dropdown',
                        multi=True,
                        placeholder='Showing all Counties as Default. Select Counties from Dropdown if interested',
                        clearable=False,
                    )
                ], width=12)
            ]),
            html.Br(),

            #------------ Define Category ---- 
            html.H4("Select Data Label", className='mb-2', style={'textAlign': 'center'}),

            dbc.Row([
                dbc.Col([
                    dcc.Dropdown(
                        id='category',
                        value=initial_column,
                        clearable=False,
                        options=[{'label': col, 'value': col} for col in ['County', 'Category', 'Road Type']]
                    )
                ], width=12)
            ]),
            #------------ Define Event Category ---- 
            html.H4("Select Event Category", className='mb-2', style={'textAlign': 'center'}),
            dbc.Row([
                dbc.Col([
                    dcc.Dropdown(
                        id='category1',
                        multi=False,
                        placeholder='Select Event Categories',
                        clearable=False,
                    )
                ], width=12)
            ]),

            #------------ Define Event Category ----
            html.H4("Select Event Type", className='mb-2', style={'textAlign': 'center'}),
            dbc.Row([
                dbc.Col([
                    dcc.Dropdown(
                        id='subcategory',
                        multi=True,
                        placeholder='Select Event Type',
                        clearable=False,
                    )
                ], width=12)
            ]),

            #------------ Define USGS Sync Button ----
            # html.Br(),
            # dbc.Row([
            #     dbc.Col([
            #         dmc.Switch(
            #             id='sync-usgs-data',
            #             checked=False,  # Set the initial state of the switch
            #             label='Sync USGS Data'
            #         )
            #     ], width=8)



            # ]),


        ], width=4)
    ]),

    #------------ Define Download Shapefile ----
    dbc.Row([
        dbc.Col([
        html.Button("Download Shapefile", id="shp-button", n_clicks=0),
    ], width=4),  # Corrected indentation
    ]),

    #------------ Define Stats Plots ----

    dbc.Row([
        dbc.Col([
            dcc.Graph(id='bar-plot')
        ], width=12, md=4),
        dbc.Col([
            dcc.Graph(id='pie-plot')
        ], width=12, md=4),
        dbc.Col([
            dcc.Graph(id='road-plot')
        ], width=12, md=4),

    ]),

    #------------ Define Timeseries Plots ----

    dbc.Row([
        dbc.Col([
            dcc.Graph(id='discharge-plot'),
            html.Div(style={'text-align': 'center'})  # Optional title for the plot
        ], width=4, md=4),  # Adjust the width as needed

        dbc.Col([
            dcc.Graph(id='stage-plot'),
            html.Div(style={'text-align': 'center'})  # Optional title for the plot
        ], width=4, md=4),  # Adjust the width as needed

        dbc.Col([
            dcc.Graph(id='precip-plot'),
            html.Div(style={'text-align': 'center'})  # Optional title for the plot
        ], width=4, md=4),  # Adjust the width as needed
    ]),

    # add weather report 

    #------------ Download CSV ----

    dbc.Row([
        dbc.Col([
            html.Button("Download CSV", id="csv-button", n_clicks=0),
            dag.AgGrid(
                id='grid',
                rowData=gdf.to_dict("records"),
                columnDefs=[{"field": i} for i in ['FID', 'County', 'Category', 'SubCategory', 'Reason Closure','Event Type',  'Road Type', 'Reported Time','lat','lon']],
                columnSize="auto",  # Adjusted column width
                defaultColDef={"filter": True, "sortable": True, "floatingFilter": True},
                suppressDragLeaveHidesColumns=True,
                selectedRows=[],
                csvExportParams={
                "fileName": "filtered_data.csv",
                },


            )
        ], width=24, md=12),
    ], className='mt-4'),
    ]),

    #------------ Here we define tabs for further info ----

        dcc.Tab(label='About', children=[
            dbc.Row([
                dbc.Col([
                    html.H1("About this Dashboard", className='mb-2', style={'textAlign': 'center'}),
                    html.P(Documentation),
                    # Add more information about the dashboard here...
                ], width=8),
            ]),
        ]),
        
        dcc.Tab(label='Alerts', children=[
            dbc.Container([  # Add this container
                dbc.Row([
                    dbc.Col([
                        html.H1("Alerts Signup Page", className='mb-2', style={'textAlign': 'center'}),
                        html.P("Here users should be able to select alerts on a given location or road."),
                        # Add more information about the dashboard here...
                    ], width=3),
                ]),
                # ... (Rest of your layout for Alerts)
            ])
        ]),

        dcc.Tab(label='Contact Us', children=[
            dbc.Container([  # Add this container
                dbc.Row([
                    dbc.Col([
                        html.H1("Contact Us", className='mb-2', style={'textAlign': 'center'}),
                        html.P("Contact Us"),
                        # Add more information about the dashboard here...
                    ], width=3),
                ]),
                # ... (Rest of your layout for Contact Us)
            ])
        ]),

    ]),

    #------------ Here tabs end ----

], fluid=True)




#------------ Function for Download csv ----

@callback(
    Output("grid", "exportDataAsCsv"),
    Input("csv-button", "n_clicks"),
    Input('grid', 'rowData'),
    )


def export_data_as_csv(n_clicks,rodata1):
    if n_clicks:
        # Define the file path for the CSV file
        csv_file_path = f"{output_dir}//filtered_data.csv"
        
        dd = pd.DataFrame.from_records(rodata1)
        #print(rodata1)
        dd.to_csv(csv_file_path)
        # Return the file path to trigger the export
        return False
    return False


#------------ Function for defining the dropdown of categories ----

@app.callback(
    Output('category1', 'options'),
    Input('category', 'value'),
)
def update_category_options(selected_category):
    if selected_category:
        # Assuming your data has a DataFrame named 'gdf'
        if selected_category == 'Category':
            
            category_options = [{'label': subcat, 'value': subcat} for subcat in gdf['Category'].unique()]
            category_options.insert(0, {'label': 'All', 'value': 'All'})

        else:
            category_options = []  # Add handling for other categories if needed
        return category_options
    else:
        return []


#------------ Function for defining the dropdown of sub-categories ----

@app.callback(
    Output('subcategory', 'options'),
    Input('category1', 'value'),
)
def update_subcategory_options(selected_category2):
    #print(selected_category2)

    if selected_category2:
        if selected_category2 == 'All':
            subcategory_options = [{'label': subcat, 'value': subcat} for subcat in
                                   gdf['SubCategory'].unique()]
        elif selected_category2 == 'Incident':
            subcategory_options = [{'label': subcat, 'value': subcat} for subcat in
                                   gdf[gdf['Category'] == 'Incident']['SubCategory'].unique()]
        elif selected_category2 == 'Weather':
            subcategory_options = [{'label': subcat, 'value': subcat} for subcat in
                                   gdf[gdf['Category'] == 'Weather']['SubCategory'].unique()]
        elif selected_category2 == 'Planned Event':
            subcategory_options = [{'label': subcat, 'value': subcat} for subcat in
                                   gdf[gdf['Category'] == 'Planned Event']['SubCategory'].unique()]  # Fixed this line

        if len(gdf['SubCategory'].unique()) > 1:
            return subcategory_options
        else:
            return []
    else:
        return []


#------------ Function interactivity between dropdown component, scatter plot on the basemap, and the bar plot ----

@app.callback(
    Output('basemap', 'figure'),
    Output('bar-plot', 'figure'),
    Output('pie-plot', 'figure'),
    Output('road-plot', 'figure'),
    Output('discharge-plot', 'figure'),   
    Output('stage-plot', 'figure'),
    Output('precip-plot', 'figure'),
    Output('county-dropdown', 'options'),
    Output('county-dropdown', 'value'),
    Output('grid', 'rowData'),
    Input('category', 'value'),
    Input('time-slider', 'value'),
    Input('basemap', 'clickData'),  # Add clickData as input
    # Input('sync-usgs-data', 'checked'),  # Add checkbox state as input
    Input('shp-button', 'n_clicks'),  # Add shp checkbox state as input
    Input('county-dropdown', 'value'),
    Input('category1', 'value'),
    Input('subcategory', 'value'),
    *[Input({'type': 'region-button', 'index': i}, 'n_clicks_timestamp') for i in range(len(available_regions))],
    *[State({'type': 'region-button', 'index': i}, 'n_clicks') for i in range(len(available_regions))]
)



def update_map(selected_yaxis, selected_time_range,click_data,shape_checked,selected_counties,selected_category, selected_category2,*region_click_timestamps_and_states): #sync_usgs_data

    global gdf  # Declare gdf as a global variable
   
    # Check if the callback is being triggered for the initial load
    ctx = callback_context
    #print(ctx.triggered)
    if not ctx.triggered:
        #print(callback_context)
        # Set the default region to 'Northern'
        clicked_region = 'Northern'

    else:
    # Find the index of the last clicked button with a valid timestamp
        valid_timestamps = [ts for ts in region_click_timestamps_and_states[:len(available_regions)] if ts is not None]
        
        # this may be causing to refresh for the county
        if not valid_timestamps:

            return dash.no_update

        last_clicked_timestamp = max(valid_timestamps)
        clicked_button_index = region_click_timestamps_and_states.index(last_clicked_timestamp)

        clicked_region = available_regions[clicked_button_index]

    # Filter the gdf DataFrame based on the selected region
    gdf2 = gdf[gdf['Region'] == clicked_region]

    # get vdot selected region
    vdot_region1 = vdot_region[vdot_region.Region == clicked_region]
    roads_in_vdot_reg = gpd.clip(roads_near_vdot,vdot_region1)

    # get values added to a dict
    def extract_coordinates(row):
        lats, lons, names = np.array([]), np.array([]), np.array([])
        feature, name = row['geometry'], row['ST_FULL']
        if isinstance(feature, (shapely.geometry.linestring.LineString, shapely.geometry.multilinestring.MultiLineString)):
            linestrings = [feature] if isinstance(feature, shapely.geometry.linestring.LineString) else feature.geoms
            for linestring in linestrings:
                x, y = linestring.xy
                lats = np.append(lats, y)
                lons = np.append(lons, x)
                names = np.append(names, [name] * len(y))
                lats = np.append(lats, None)
                lons = np.append(lons, None)
                names = np.append(names, None)
        return lats, lons, names

    # 
    result = roads_in_vdot_reg.apply(extract_coordinates, axis=1)

    lats, lons, names = np.concatenate(result.apply(lambda x: x[0]).values), np.concatenate(result.apply(lambda x: x[1]).values), np.concatenate(result.apply(lambda x: x[2]).values)

    road_lats = lats
    road_lons = lons
    road_names = names



    # finding the counties in the given region
    counties_in_region = gdf2[(gdf2['Region'] == clicked_region)]['County'].unique()
    
    plot_pi_bar = True
    if selected_counties:
        if selected_counties is None or not selected_counties:
            gdf2 = gdf[gdf['Region'] == clicked_region].copy()
        else:
            gdf2 = gdf[(gdf['Region'] == clicked_region) & (gdf['County'].isin(selected_counties))].copy()
            #print(selected_counties)
            if len(gdf2) < 2:
                plot_pi_bar = False
                gdf2 = gdf[gdf['Region'] == clicked_region].copy()
        #print(gdf2)

    # finding the counties in the given region
    #county_options = [{'label': county, 'value': county} for county in counties_in_region]
    county_options = [county for county in counties_in_region if len(gdf[(gdf['Region'] == clicked_region) & (gdf['County'] == county)]) > 1]

    # here we filter the dataframe

    #gdf3 = gdf2.copy()
    #print(selected_category,selected_category2)

    column_use =  'Category'

    if selected_yaxis == 'County':
        column_use =  'County'
    if selected_yaxis == 'Road Type':
        column_use =  'Road Type'
    if selected_yaxis == 'Category':
        column_use =  'Category'



    if selected_category is None or (selected_category == 'All' and (selected_category2 is None or not selected_category2)):
        gdf3 = gdf2.copy()


    elif selected_category != 'All' and (selected_category2 is None or not selected_category2): 
        gdf3 = gdf2[gdf2['Category'] == selected_category].copy()

        column_use =  'SubCategory'



    elif selected_category != 'All' and selected_category2 is not None: 
        gdf3 = gdf2[(gdf2['Category'] == selected_category) & (gdf2['SubCategory'].isin(selected_category2))].copy()
        if selected_category in ['Incident','Planned Event'] :  
            column_use =  'Event Type'
        elif selected_category in ['Weather']:
            column_use =  'Reason Closure'
        else:
            column_use =  'SubCategory'

    elif selected_category == 'All':
        gdf3 = gdf2.copy()

        column_use =  'Category'

        if selected_category2 is not None:
            gdf3 = gdf2[gdf2['SubCategory'].isin(selected_category2)].copy() 

            column_use =  'Reason Closure'
            
    else:
        gdf3 = gdf2.copy()


    # Initialize discharge_plot outside the conditional block
    discharge_plot = px.line(title='Discharge Observations')
    stage_plot = px.line(title='Stage Observations')
    precip_plot = px.line(title='Precipitation Observations')

    # Convert timestamp back to datetime
    if selected_time_range == None:
        selected_time_range = [slider_start,slider_end]

    start_date = pd.to_datetime(selected_time_range[0], unit='s')
    end_date = pd.to_datetime(selected_time_range[1], unit='s')

    filtered_gdf = gdf3[
        (gdf3['Reported Time2'] >= start_date) & (gdf3['Reported Time2'] <= end_date)
    ]
    selected_countieslabel = selected_counties
    if selected_counties == None or len(selected_counties)==0:
        selected_countieslabel = ''
    # Build the Plotly scatter plot on the basemap with hover information
    #print(filtered_gdf)
    if len(filtered_gdf)>2: 

        if selected_category == 'All': 
            # Concatenate the two Series
            combined_series = pd.concat([gdf3['Reason Closure'].dropna(), gdf3['Event Type'].dropna()])

            # Sample unique values for color to match the length of other parameters
            sampled_color_data = combined_series.sample(n=len(filtered_gdf))

            # Scatter Mapbox Plot
            fig = px.scatter_mapbox(filtered_gdf, 
                                    lat=filtered_gdf.lat, 
                                    lon=filtered_gdf.lon, 
                                    color=sampled_color_data,
                                    mapbox_style="carto-positron", 
                                    title=f'Incidents reported during {str(start_date)[:7]} & {str(end_date)[:7]}: {column_use} {selected_countieslabel}\n Total Incidents: {len(filtered_gdf)}',
                                    hover_name=filtered_gdf.index, 
                                    hover_data={column_use: True})
        else:


            fig = px.scatter_mapbox(filtered_gdf, lat=filtered_gdf.lat, lon=filtered_gdf.lon, color=column_use,
                                    mapbox_style="carto-positron", title=f'Incidents reported during {str(start_date)[:7]} & {str(end_date)[:7]}: {column_use} {selected_countieslabel}\n Total Incidents: {len(filtered_gdf)}',
                                    hover_name=filtered_gdf.index, hover_data={column_use: True})
        fig.add_trace(go.Scattermapbox(
                    lat=usgs_gages.lat,
                    lon=usgs_gages.lon,
                    mode='markers',
                    text=usgs_gages.STAID,
                    hoverinfo='text',
                    marker=dict(size=10, color='black'),  # Adjust the marker size and color as needed
                    opacity=0.3,
                    name='USGS stream gages', # Set the legend name
                    visible ='legendonly'
                ))

        fig.add_trace(go.Scattermapbox(
                    lat=met_gages.lat,
                    lon=met_gages.lon,
                    mode='markers',
                    text=met_gages.STAID,
                    hoverinfo='text',
                    marker=dict(size=10, color='blue'),  # Adjust the marker size and color as needed
                    opacity=0.3,
                    name='USGS met gages', # Set the legend name
                    visible ='legendonly'
                ))

        fig.add_trace(go.Scattermapbox(
            lat=road_lats,
            lon=road_lons,
            mode='lines',
            text=road_names,  # Road names in the 'ST_NAME' column
            hoverinfo='text',  # Display text on hover
            line=dict(width=1, color='black'),
            opacity=0.2,
            name='roads',
            showlegend=True,
            visible ='legendonly',
            below='Roads'  # Set to minimum display order

        ))





        fig.update_layout(height=800)  # Adjust the height value as needed
    else:
            # Create an empty scatter map with basemap and USGS gages
        fig = px.scatter_mapbox(filtered_gdf,lat=[], lon=[], color=[], title=f'Incidents reported during {str(start_date)[:7]} & {str(end_date)[:7]}: {column_use} {selected_countieslabel}',
                                mapbox_style="carto-positron", hover_data={column_use: True})

        # Add a text annotation to indicate no data
        fig.add_annotation(
            go.layout.Annotation(
                text="No data for selected filters",
                x=0.5,
                y=0.5,
                xref="paper",
                yref="paper",
                showarrow=False,
                font=dict(size=16),
            )
        )


        fig.update_layout(height=800) 

    selected_road = []

    if click_data:
        #print(click_data,type(click_data))

        if 'lon' in click_data['points'][0]:
         
            # Extract clicked point number
            #point_number = click_data['points'][0]['pointNumber']
            AOI = Point(click_data['points'][0]['lon'],click_data['points'][0]['lat'])
            print(click_data['points'][0]['lon'],click_data['points'][0]['lat'])

            aoi_gdf = gpd.GeoDataFrame(geometry=[AOI])

            aoi_gdf.crs = '4326'
            aoi_gdf = aoi_gdf.to_crs(26918)

            # Add a scattermapbox trace with a marker for the clicked point
            fig.add_trace(go.Scattermapbox(
                lat=[click_data['points'][0]['lat']],
                lon=[click_data['points'][0]['lon']],
                mode='markers+text',
                marker=dict(size=10, color='black'),  # Adjust the marker size and color as needed
                text='Selected Point',
                textposition='top center',
                name='Selected Point'
            ))


            # Identify the road or feature clicked using the geometry
            clicked_point1 = Point(click_data['points'][0]['lon'],click_data['points'][0]['lat'])
            clicked_point = gpd.GeoDataFrame(geometry=[clicked_point1])

            clicked_point.crs = '4326'
            selected_road = identify_selected_road(clicked_point, roads)

            selected_road = selected_road[0:1] # always select the first road


            # Update the selected road output
            

            # get the info on the selected road
            if len(selected_road)>0:
                selected_road_output = f"Selected Road: {selected_road.ST_FULL.tolist()}"
                #print(selected_road_output)
                lat1_se,lon1_se = [],[]

                for index, row in selected_road.iterrows():
                    for coords in list(row['geometry'].coords):
                           lat1,lon1 = coords[1],coords[0]                
                           lat1_se.append(lat1);lon1_se.append(lon1)
                           
                name_se=selected_road.ST_FULL.tolist()*len(lon1_se) 
                    
                lat1_se = np.append(lat1_se, None)
                lon1_se = np.append(lon1_se, None)
                name_se = np.append(name_se, None)   
                    
                 

                fig.add_trace(go.Scattermapbox(
                    lat=lat1_se,
                    lon=lon1_se,
                    mode='lines',
                    text=name_se,  # Road names in the 'ST_NAME' column
                    hoverinfo='text',  # Display text on hover
                    line=dict(width=1, color='red'),
                    opacity=0.8,
                    name='selected road',
                    showlegend=True,
                    #visible ='legendonly',
                    

                ))

                # Update map extents
                fig.update_layout(
                    mapbox=dict(
                        center=dict(lat=np.mean(lat1_se[:-1]), lon=np.mean(lon1_se[:-1])),
                        zoom=12
                    )
                )

                roads_selected_buff = selected_road.buffer(0.0001)

                filtered_gdf_cleaned = gpd.GeoDataFrame(filtered_gdf, geometry=gpd.points_from_xy(filtered_gdf.lon, filtered_gdf.lat), crs='EPSG:4326')
                # filtered_gdf_cleaned = filtered_gdf_cleaned.crs(4326)
                roads_incidents = gpd.clip(filtered_gdf_cleaned,roads_selected_buff)


                roads_incidents['Reported Time2'] = pd.to_datetime(roads_incidents['Reported Time'])

                roads_incidents.index = roads_incidents['Reported Time2'] 



        # Check if points are clicked on the map

   
    road_plot = px.line(title='Road Stats Over Time')
    if len(filtered_gdf)>2:
        # Build the Plotly bar plot for unique entries count
        #print(filtered_gdf,column_use)




        if selected_category == 'All':
            combined_series = pd.concat([filtered_gdf['Reason Closure'].dropna(), filtered_gdf['Event Type'].dropna()])
            bar_plot = px.bar(combined_series.value_counts(),
                              x=combined_series.value_counts().index,
                              y=combined_series.value_counts().values,
                              color=combined_series.dropna().unique(),
                              title=f'Unique Entries Count: {column_use}',
                              labels={'y': 'Number of occurrences'})

            pie_plot = px.pie(values=combined_series.value_counts(), names=combined_series.dropna().unique())
            #road_plot = px.line(title='Road Stats')
            if len(selected_road)>0:
                road_plot =px.scatter(
                roads_incidents,
                x=roads_incidents.index,
                y='Category',
                color='Category',  # Color by incident category
                title="Road Stats Over Time",
                labels={'Category': 'Event Category'}  # Custom label for y-axis
                )




        else:
            bar_plot = px.bar(filtered_gdf[column_use].value_counts(),
                              x=filtered_gdf[column_use].value_counts().index,
                              y=filtered_gdf[column_use].value_counts().values,
                              color=filtered_gdf[column_use].dropna().unique(),
                              title=f'Unique Entries Count: {column_use}',
                              labels={'y': 'Number of occurrences'})


            pie_plot = px.pie(values=filtered_gdf[column_use].value_counts(), names=filtered_gdf[column_use].dropna().unique())
            #road_plot = px.line(title='Road Stats')
            #road_plot = px.line(title='Road Stats')

            if len(selected_road)>0:
                road_plot =px.scatter(
                roads_incidents,
                x=roads_incidents.index,
                y='Category',
                color='Category',  # Color by incident category
                title="Road Stats Over Time",
                labels={'Category': 'Event Category'}  # Custom label for y-axis
                )


    else:
        # Initialize discharge_plot outside the conditional block
        bar_plot = px.line(title='Bar Plot')
        pie_plot = px.line(title='Pie Plot')
        road_plot = px.line(title='Road Stats Over Time')


    if shape_checked:
            filtered_shp = gpd.GeoDataFrame(filtered_gdf, geometry=gpd.points_from_xy(filtered_gdf['lon'], filtered_gdf['lat']), crs='EPSG:4326')
            del filtered_shp['Reported Time2']
            tmz = str(datetime.datetime.now())[:16]
            tmz =''
            filtered_shp.to_file(f'{output_dir}//filtered_shapefile_{tmz}.shp')
 
    # if sync_usgs_data:

    if click_data:
        #print(click_data,type(click_data))

        if 'lon' in click_data['points'][0]:
         
            # Extract clicked point number
            #point_number = click_data['points'][0]['pointNumber']
            AOI = Point(click_data['points'][0]['lon'],click_data['points'][0]['lat'])


            aoi_gdf = gpd.GeoDataFrame(geometry=[AOI])

            aoi_gdf.crs = '4326'
            aoi_gdf = aoi_gdf.to_crs(26918)

            # Add a circle around the clicked point
            circle_radius = 100  # Adjust the radius as needed
            circle = aoi_gdf.buffer(circle_radius)

            


            nearest_gage = gpd.sjoin_nearest(usgs_gages,aoi_gdf,max_distance=5000) #
            if len(nearest_gage)==0:
                nearest_gage = gpd.sjoin_nearest(usgs_gages,aoi_gdf,max_distance=20000) #
            print('Stream Gages')
            print(nearest_gage)

            nearest_met_gages = gpd.sjoin_nearest(met_gages,aoi_gdf,max_distance=5000)
            if len(nearest_met_gages)==0:
                nearest_met_gages = gpd.sjoin_nearest(met_gages,aoi_gdf,max_distance=20000) #
            print('Precipitation Gages')
            print(nearest_met_gages)

            


            if selected_time_range == None:
                selected_time_range = [slider_start,slider_end]

            start_date = pd.to_datetime(selected_time_range[0], unit='s')
            end_date = pd.to_datetime(selected_time_range[1], unit='s')

            disc_tmp = pd.DataFrame()
            stage_tmp = pd.DataFrame()
            print(start_date,end_date)
            
            # selecting preloaded data for nearest gages
            usgs_ = nearest_gage.STAID.tolist()
            usgs_disc,usgs_stg = [],[]

            for ii in usgs_:
                if ii in disc_df.columns:
                    usgs_disc.append(ii)
                if ii in disc_df.columns:
                    usgs_stg.append(ii)

            disc_tmp = disc_df[usgs_disc]
            stage_tmp = stage_df[usgs_stg]
            # applying time indexing
            disc_tmp = disc_tmp[str(start_date):str(end_date)]
            stage_tmp = stage_tmp[str(start_date):str(end_date)]

            # for stn in tqdm.tqdm(nearest_gage.STAID):

            #     #stn='01652500'
            #     site = str(stn)
            #     # get instantaneous values (iv), daily values (dv)
            #     df = nwis.get_record(sites=site, service='iv', start=f'{str(start_date)[:7]}-01', end=f'{str(end_date)[:7]}-01') # change this with slider

                
            #     if '00060' in df.columns: # discharge
            #         disc_tmp[site]=df['00060']
            #     if '00065' in df.columns: # stage
            #         stage_tmp[site]=df['00065']



            if len(disc_tmp)>1:
                disc_tmp[disc_tmp < -999] = np.nan
                discharge_plot = px.line(disc_tmp,
                                    x=disc_tmp.index,
                                    y=disc_tmp.columns,
                                    title='Discharge Observations',
                                    labels={'y': 'Discharge (CFS)','x': 'Date Time'})

                # Update the layout of the discharge plot to add ylabel
                discharge_plot.update_layout(
                    yaxis_title='Discharge (CFS)'
                )
            else:
                    # If no observations, create an empty plot with a message
                discharge_plot = go.Figure()
                discharge_plot.add_annotation(
                    text='No discharge observations available for the selected gage and time range.',
                    showarrow=False,
                    arrowhead=1,
                    arrowcolor='black',
                    arrowwidth=2,
                    ax=0,
                    ay=-40
                )


            if len(stage_tmp)>1:
                stage_tmp[stage_tmp > 100] = np.nan
                stage_tmp[stage_tmp < -99] = np.nan
                stage_plot = px.line(stage_tmp,
                                    x=stage_tmp.index,
                                    y=stage_tmp.columns,
                                    title='Stage Observations',
                                    labels={'y': 'Stage (ft)','x': 'Date Time'})

                # Update the layout of the discharge plot to add ylabel
                stage_plot.update_layout(
                    yaxis_title='Stage (ft)'
                )
            else:
                # If no observations, create an empty plot with a message
                stage_plot = go.Figure()
                stage_plot.add_annotation(
                    text='No stage observations available for the selected gage and time range.',
                    showarrow=False,
                    arrowhead=1,
                    arrowcolor='black',
                    arrowwidth=2,
                    ax=0,
                    ay=-40
                )


            # clipping the precip data
            met_tmp = pd.DataFrame()

            # selecting preloaded data for nearest gages
            usgs_p = nearest_met_gages.STAID.tolist()
            usgs_met = []

            for ii in usgs_p:
                if ii in met_df.columns:
                    usgs_met.append(ii)
                

            met_tmp = met_df[usgs_met]
            # applying time indexing
            met_tmp = met_tmp[str(start_date):str(end_date)]
            
            # for stn in tqdm.tqdm(nearest_met_gages.STAID):
                
            #     #stn='01652500'
            #     site = str(stn)
            #     # get instantaneous values (iv)
            #     df = nwis.get_record(sites=site, service='iv', start=f'{str(start_date)[:7]}-01', end=f'{str(end_date)[:7]}-01') # change this with slider

                
            #     if '00045' in df.columns: # discharge
            #         met_tmp[site]=df['00045']

            if len(met_tmp)>1 and len(usgs_met)>0:
                met_tmp[met_tmp > 100] = np.nan
                met_tmp[met_tmp < -99] = np.nan

                precip_plot = px.line(met_tmp,
                                    x=met_tmp.index,
                                    y=met_tmp.columns,
                                    title='Precipitation Observations',
                                    labels={'y': 'inches (in)','x': 'Date Time'})

                # Update the layout of the discharge plot to add ylabel
                precip_plot.update_layout(
                    yaxis_title='inches (in)'
                )
            else:

                # If no observations, create an empty plot with a message
                precip_plot = go.Figure()
                precip_plot.add_annotation(
                    text='No precipitation observations available for the selected gage and time range.',
                    showarrow=False,
                    arrowhead=1,
                    arrowcolor='black',
                    arrowwidth=2,
                    ax=0,
                    ay=-40
                )















    my_cellStyle = {
        "styleConditions": [
            {
                "condition": f"params.colDef.field == '{selected_yaxis}'",
                "style": {"backgroundColor": "white"},
            },
            {
                "condition": f"params.colDef.field != '{selected_yaxis}'",
                "style": {"color": "black"}
            },
        ]
    }


    

    return fig, bar_plot, pie_plot, road_plot, discharge_plot, stage_plot, precip_plot, county_options, selected_counties, filtered_gdf.to_dict("records")



def identify_selected_road(clicked_point, roads_gdf):
    selected_road = None

    # Check if the clicked point intersects with any road geometry
    
    selected_road = gpd.sjoin_nearest(roads_gdf,clicked_point,max_distance=0.00001)
    
    return selected_road

#------------ Function to update button color ---
previous_button_state = [0] * len(available_regions)

# Callback to update button color
@app.callback(
    [Output({'type': 'region-button', 'index': i}, 'color') for i in range(len(available_regions))],
    [Input({'type': 'region-button', 'index': i}, 'n_clicks') for i in range(len(available_regions))],
    [State({'type': 'region-button', 'index': i}, 'id') for i in range(len(available_regions))],
    prevent_initial_call=True,
)
def update_button_color(*n_clicks_values_and_ids):
    global previous_button_state
    
    # Identify the clicked button index
    clicked_button_index = None
    for i, (n_clicks, button_id) in enumerate(zip(n_clicks_values_and_ids, range(len(available_regions)))):
        if n_clicks and previous_button_state[i] != n_clicks:
            clicked_button_index = i
            previous_button_state[i] = n_clicks
    
    if clicked_button_index is None:
        # No button has been clicked, set all colors to 'secondary'
        colors = ['secondary'] * len(available_regions)
    else:
        # Create a list of colors for each button
        colors = ['secondary'] * len(available_regions)
    
        # Set the color of the clicked button to yellow ('warning' class)
        colors[clicked_button_index] = 'warning'
    
    return colors


#========================================================================================================

if __name__ == '__main__':
    app.run_server(debug=False, port=8002)
