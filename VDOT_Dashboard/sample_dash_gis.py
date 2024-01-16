from dash import Dash, html, dcc, Input, Output
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

# roads = gpd.read_file(r"C:/Users/Arslaan Khalid/Desktop/Papers_by_Arslaan/VDOT_Dashboard/Data/roads_testing_v2.shp")
 
# roads = roads.to_crs(4326)

# lats = []
# lons = []
# names = []
# # instead of this, just get this into np lists and then load them only, or dicts, so no need to redo on the fly
# for feature, name in zip(roads.geometry, roads.ST_FULL):
#     if isinstance(feature, shapely.geometry.linestring.LineString):
#         linestrings = [feature]
#     elif isinstance(feature, shapely.geometry.multilinestring.MultiLineString):
#         linestrings = feature.geoms
#     else:
#         continue
#     for linestring in linestrings:
#         x, y = linestring.xy
#         lats = np.append(lats, y)
#         lons = np.append(lons, x)
#         names = np.append(names, [name]*len(y))
#         lats = np.append(lats, None)
#         lons = np.append(lons, None)
#         names = np.append(names, None)

# road_lats = lats
# road_lons = lons
# road_names = names



# loadin vdot regions

vdot_region = gpd.read_file(r'C:/Users/Arslaan Khalid/Desktop/Papers_by_Arslaan/VDOT_Dashboard/Data/VDOT_regions.shp')
vdot_region = vdot_region.to_crs(4326)

# loading USGS Stream data

usgs_gages = gpd.read_file(r'C:/Users/Arslaan Khalid/Desktop/Papers_by_Arslaan/VDOT_Dashboard/Data/realstx_shp/realstx.shp')


usgs_gages['lat'] = usgs_gages['geometry'].y
usgs_gages['lon'] = usgs_gages['geometry'].x

usgs_gages = usgs_gages.to_crs(26918)

# loading USGS Met data
met_gages = gpd.read_file(r'C:\Users\Arslaan Khalid\Desktop\Papers_by_Arslaan\VDOT_Dashboard\Data\USGS_met_gages.shp')

met_gages['lat'] = met_gages['geometry'].y
met_gages['lon'] = met_gages['geometry'].x

met_gages = met_gages.to_crs(26918)



# Read the shapefile
#gdf = gpd.read_file(r"C:\Users\Arslaan Khalid\Desktop\Papers_by_Arslaan\VDOT_Dashboard\Data\sample_data.shp")
#gdf = gpd.read_file(r"C:\Users\Arslaan Khalid\Desktop\Papers_by_Arslaan\VDOT_Dashboard\Data\road_closures_sample.shp")
rd_file = pd.read_csv(r"C:\Users\Arslaan Khalid\Desktop\Papers_by_Arslaan\VDOT_Dashboard\Data\road_closures.csv")

# rd_file = rd_file[:1000]

gdf = gpd.GeoDataFrame(rd_file, geometry=gpd.points_from_xy(rd_file['x'], rd_file['y']), crs='EPSG:4326')



# loading Roads network

roads_in_2 = pd.read_csv(r'C:/Users/Arslaan Khalid/Desktop/Papers_by_Arslaan/VDOT_Dashboard/Data/road_lines_undissolved.csv')
del roads_in_2['Unnamed: 0']

# roads_in_2 = roads_in_2[:1000]


roads_in_2['geometry'] = roads_in_2['geometry'].apply(lambda x: x if pd.notnull(x) else None) 

roads = gpd.GeoDataFrame(roads_in_2, geometry=gpd.GeoSeries.from_wkt(roads_in_2['geometry']), crs='EPSG:4326')

roads = roads[roads.geometry.notnull()]


roads_near_vdot = gpd.sjoin_nearest(roads, gdf,max_distance=0.0001)


# working with the vdot incidents

gdf['lat'] = gdf['geometry'].y
gdf['lon'] = gdf['geometry'].x

# Convert the geometry column to GeoJSON
gdf['geometry'] = gdf['geometry'].apply(lambda geom: geom.__geo_interface__)

intersted_columns = ['fid', 'containing','containi_1', 'event_cate', 'event_subc','reason_for',  'road_syste',
       'route_name', 'type_event', 'update_tim', 'geometry',
       'lat', 'lon']

gdf = gdf[intersted_columns]

gdf.columns = ['FID', 'County', 'Region','Category','SubCategory', 'Reason Closure', 'Road Type',
       'Route Name', 'Event Type', 'Reported Time', 'geometry',
       'lat', 'lon']

available_regions = gdf['Region'].unique()


gdf['Reported Time2'] = pd.to_datetime(gdf['Reported Time'], errors='coerce')

gdf = gdf.dropna(subset=['Reported Time2'])

gdf = gdf.sort_values(['Reported Time2'])


gdf_global = gdf.copy()
unique_dates = pd.to_datetime(gdf['Reported Time2']).dt.to_period('M').unique()


# timewindows = {int(date.to_timestamp().timestamp()): {'label': f'{date.month}\n{date.year}' if date.month == 1 else str(date.month)} for date in unique_dates}

# timewindows = {
#     int(date.to_timestamp().timestamp()): {
#         'label': f'{date.month}\n{date.year}' if date.month in [6, 12] else ''
#     } for date in unique_dates
# }


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
    
    #label = f'{date.month}\n{date.year}' 
    timewindows[timestamp] = {'label': label}



slider_start = list(timewindows.keys())[0]
slider_end = list(timewindows.keys())[-1]



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


initial_column = 'Category'

app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.layout = dbc.Container([
    html.H1("Interactive VDOT Dashboard", className='mb-2', style={'textAlign': 'center'}),
    html.Br(), html.Br(), html.Br(),


    dcc.Tabs([
            dcc.Tab(label='Main Dashboard', children=[

    dbc.Row([
        dbc.Col([
            dcc.Graph(id='basemap', clickData={'points': [{'pointNumber': 0}]}),  # Initial clickData
        ], width=8),

        dbc.Col([
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

            #------------ Define Region ---- 
            html.H3("Select Region", className='mb-2', style={'textAlign': 'center'}),

            # Add buttons to select regions
            dbc.Row([
                dbc.Col(button, width=2)  # Set the width to 2 for each button
                for button in region_buttons
            ], className='mt-2'),

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
            html.Br(),
            dbc.Row([
                dbc.Col([
                    dmc.Switch(
                        id='sync-usgs-data',
                        checked=False,  # Set the initial state of the switch
                        label='Sync USGS Data'
                    )
                ], width=8)



            ]),
        ], width=4)
    ]),
    dbc.Row([
        dbc.Col([
        html.Button("Download Shapefile", id="shp-button", n_clicks=0),
    ], width=4),  # Corrected indentation
    ]),

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

    dbc.Row([
        dbc.Col([
            dcc.Graph(id='discharge-plot'),
            html.Div("Discharge Observations", style={'text-align': 'center'})  # Optional title for the plot
        ], width=4, md=4),  # Adjust the width as needed

        dbc.Col([
            dcc.Graph(id='stage-plot'),
            html.Div("Stage Observations", style={'text-align': 'center'})  # Optional title for the plot
        ], width=4, md=4),  # Adjust the width as needed

        dbc.Col([
            dcc.Graph(id='precip-plot'),
            html.Div("Precipitation Observations", style={'text-align': 'center'})  # Optional title for the plot
        ], width=4, md=4),  # Adjust the width as needed
    ]),

    # add a Precip and USGS gage data
    # add weather report

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

        dcc.Tab(label='About', children=[
            dbc.Row([
                dbc.Col([
                    html.H1("About this Dashboard", className='mb-2', style={'textAlign': 'center'}),
                    html.P("This dashboard provides interactive visualizations for VDOT data."),
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

], fluid=True)

@callback(
    Output("grid", "exportDataAsCsv"),
    Input("csv-button", "n_clicks"),
)
def export_data_as_csv(n_clicks):
    if n_clicks:
        return True
    return False


# Combine the two callback functions for 'subcategory' options into one
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



# Create interactivity between dropdown component, scatter plot on the basemap, and the bar plot
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
    Input('sync-usgs-data', 'checked'),  # Add checkbox state as input
    Input('county-dropdown', 'value'),
    Input('category1', 'value'),
    Input('subcategory', 'value'),
    *[Input({'type': 'region-button', 'index': i}, 'n_clicks_timestamp') for i in range(len(available_regions))],
    *[State({'type': 'region-button', 'index': i}, 'n_clicks') for i in range(len(available_regions))]
)




def update_map(selected_yaxis, selected_time_range,click_data,sync_usgs_data,selected_counties,selected_category, selected_category2,*region_click_timestamps_and_states):

    global gdf  # Declare gdf as a global variable
    # Determine which button was clicked

    

    # Find the index of the last clicked button with a valid timestamp
    valid_timestamps = [ts for ts in region_click_timestamps_and_states[:len(available_regions)] if ts is not None]
    
    # this may be causing to refresh for the county
    if not valid_timestamps:
        return dash.no_update

    last_clicked_timestamp = max(valid_timestamps)
    clicked_button_index = region_click_timestamps_and_states.index(last_clicked_timestamp)

    clicked_region = available_regions[clicked_button_index]
    #print(clicked_button_index,clicked_region)


    # Filter the gdf DataFrame based on the selected region
    gdf2 = gdf[gdf['Region'] == clicked_region]


    roads_near_vdot

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
            print(selected_counties)
            if len(gdf2) < 2:
                plot_pi_bar = False
                gdf2 = gdf[gdf['Region'] == clicked_region].copy()
        #print(gdf2)

    # finding the counties in the given region
    #county_options = [{'label': county, 'value': county} for county in counties_in_region]
    county_options = [county for county in counties_in_region if len(gdf[(gdf['Region'] == clicked_region) & (gdf['County'] == county)]) > 1]

    # here we filter the dataframe

    #gdf3 = gdf2.copy()
    print(selected_category,selected_category2)

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
    print(filtered_gdf)
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
                    text=met_gages.SiteID,
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

            # Update the selected road output
            

            # get the info on the selected road
            if len(selected_road)>0:
                selected_road_output = f"Selected Road: {selected_road.ST_FULL.tolist()}"
                print(selected_road_output)
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
                        center=dict(lat=np.mean(lat1_se), lon=np.mean(lon1_se)),
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



    if sync_usgs_data:
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

                


                nearest_gage = gpd.sjoin_nearest(usgs_gages,aoi_gdf,max_distance=5000)
                print('Stream Gages')
                print(nearest_gage)

                nearest_met_gages = gpd.sjoin_nearest(met_gages,aoi_gdf,max_distance=10000)
                print('Precipitation Gages')
                print(nearest_met_gages)

                disc_tmp = pd.DataFrame()


                if selected_time_range == None:
                    selected_time_range = [slider_start,slider_end]

                start_date = pd.to_datetime(selected_time_range[0], unit='s')
                end_date = pd.to_datetime(selected_time_range[1], unit='s')


                for stn in tqdm.tqdm(nearest_gage.STAID):

                    #stn='01652500'
                    site = str(stn)
                    # get instantaneous values (iv), daily values (dv)
                    df = nwis.get_record(sites=site, service='iv', start=f'{str(start_date)[:7]}-01', end=f'{str(end_date)[:7]}-01') # change this with slider

                    
                    if '00060' in df.columns: # discharge
                        disc_tmp[site]=df['00060']




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



                stage_tmp = pd.DataFrame()
                
                for stn in tqdm.tqdm(nearest_gage.STAID):
                    
                    #stn='01652500'
                    site = str(stn)
                    # get instantaneous values (iv)
                    df = nwis.get_record(sites=site, service='iv', start=f'{str(start_date)[:7]}-01', end=f'{str(end_date)[:7]}-01') # change this with slider

                    
                    if '00065' in df.columns: # discharge
                        stage_tmp[site]=df['00065']

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

                met_tmp = pd.DataFrame()
                
                for stn in tqdm.tqdm(nearest_met_gages.SiteID):
                    
                    #stn='01652500'
                    site = str(stn)
                    # get instantaneous values (iv)
                    df = nwis.get_record(sites=site, service='iv', start=f'{str(start_date)[:7]}-01', end=f'{str(end_date)[:7]}-01') # change this with slider

                    
                    if '00045' in df.columns: # discharge
                        met_tmp[site]=df['00045']

                if len(met_tmp)>1:
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




if __name__ == '__main__':
    app.run_server(debug=False, port=8002)
