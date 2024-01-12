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

usgs_gages = gpd.read_file(r'C:/Users/Arslaan Khalid/Desktop/Papers_by_Arslaan/VDOT_Dashboard/Data/realstx_shp/realstx.shp')


usgs_gages['lat'] = usgs_gages['geometry'].y
usgs_gages['lon'] = usgs_gages['geometry'].x

usgs_gages = usgs_gages.to_crs(26918)

# Read the shapefile
#gdf = gpd.read_file(r"C:\Users\Arslaan Khalid\Desktop\Papers_by_Arslaan\VDOT_Dashboard\Data\sample_data.shp")
gdf = gpd.read_file(r"C:\Users\Arslaan Khalid\Desktop\Papers_by_Arslaan\VDOT_Dashboard\Data\road_closures_sample.shp")

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

unique_dates = pd.to_datetime(gdf['Reported Time2']).dt.to_period('M').unique()


timewindows = {int(date.to_timestamp().timestamp()): {'label': f'{date.month}\n{date.year}' if date.month == 1 else str(date.month)} for date in unique_dates}

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
    
    html.Br(),html.Br(),html.Br(),
    
    html.H3("Select Time Period", className='mb-2', style={'textAlign': 'center'}),

    dbc.Row([
        dbc.Col([
            dcc.RangeSlider(
                id='time-slider',
                min=slider_start,
                max=slider_end,
                step=None,
                marks = timewindows


            )
        ], width=12)
    ]),

    #------------ Define Region ---- 
    html.H3("Select Region", className='mb-2', style={'textAlign': 'center'}),

    # Add buttons to select regions
    dbc.Row(region_buttons, className='mt-2'),

    #------------ Define County ---- 
    html.H3("Select County", className='mb-2', style={'textAlign': 'center'}),
    html.Br(),
    # Dropdown for counties
    dbc.Row([
        dbc.Col([
            dcc.Dropdown(
                id='county-dropdown',
                multi=True,
                placeholder='Showing all Counties as Default. Select one County from Dropdown if interested',
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
                options=[{'label': col, 'value': col} for col in ['County','Category','Road Type']]
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
    html.H4("Select Event Sub-Category", className='mb-2', style={'textAlign': 'center'}),
    dbc.Row([
        dbc.Col([
            dcc.Dropdown(
                id='subcategory',
                multi=True,
                placeholder='Select Event Subcategories',
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
        ], width=12)
    ]),


    # load the basemap

    dbc.Row([
        dbc.Col([
            dcc.Graph(id='basemap', clickData={'points': [{'pointNumber': 0}]}),  # Initial clickData
        ], width=14)
    ]),


    dbc.Row([
        dbc.Col([
            dcc.Graph(id='bar-plot')
        ], width=12, md=6),
        dbc.Col([
                    dcc.Graph(id='pie-plot')
                ], width=12, md=6),
    ]),


    dbc.Row([
        dbc.Col([
            dcc.Graph(id='discharge-plot')
        ], width=12, md=6),
        dbc.Col([
            dcc.Graph(id='stage-plot')
        ], width=12, md=6),
        
    ]),


    # add a Precip and USGS gage data
    # add weather report

    dbc.Row([

        dbc.Col([
            dag.AgGrid(
                id='grid',
                rowData=gdf.to_dict("records"),
                columnDefs=[{"field": i} for i in ['FID','County','Category', 'Road Status', 'Road Type','Event Type','Reported Time']],
                columnSize="auto",  # Adjusted column width
            )
        ], width=24, md=12),
    ], className='mt-4'),

])


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
        if selected_category2 == 'Incident':
            subcategory_options = [{'label': subcat, 'value': subcat} for subcat in
                                   gdf[gdf['Category'] == 'Incident']['SubCategory'].unique()]
        elif selected_category2 == 'Weather':
            subcategory_options = [{'label': subcat, 'value': subcat} for subcat in
                                   gdf[gdf['Category'] == 'Weather']['SubCategory'].unique()]
        elif selected_category2 == 'Planned Event':
            subcategory_options = [{'label': subcat, 'value': subcat} for subcat in
                                   gdf[gdf['Category'] == 'Planned Event']['SubCategory'].unique()]  # Fixed this line
        else:
            subcategory_options = []  # Add handling for other categories if needed
        return subcategory_options
    else:
        return []


# Create interactivity between dropdown component, scatter plot on the basemap, and the bar plot
@app.callback(
    Output('basemap', 'figure'),
    Output('bar-plot', 'figure'),
    Output('pie-plot', 'figure'),
    Output('discharge-plot', 'figure'),   
    Output('stage-plot', 'figure'),  
    Output('grid', 'defaultColDef'),
    Output('county-dropdown', 'options'),
    Output('county-dropdown', 'value'),
    Input('category', 'value'),
    Input('time-slider', 'value'),
    Input('basemap', 'clickData'),  # Add clickData as input
    Input('sync-usgs-data', 'checked'),  # Add checkbox state as input
    Input('county-dropdown', 'value'),
    *[Input({'type': 'region-button', 'index': i}, 'n_clicks_timestamp') for i in range(len(available_regions))],
    *[State({'type': 'region-button', 'index': i}, 'n_clicks') for i in range(len(available_regions))]
)




def update_map(selected_yaxis, selected_time_range,click_data,sync_usgs_data,selected_counties,*region_click_timestamps_and_states):

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


    

    if selected_counties:
        #print(selected_counties)
        gdf2 = gdf[(gdf['Region'] == clicked_region) & (gdf['County'].isin(selected_counties))]

    # finding the counties in the given region
    counties_in_region = gdf[(gdf['Region'] == clicked_region)]['County'].unique()
    county_options = [{'label': county, 'value': county} for county in counties_in_region]


    # Initialize discharge_plot outside the conditional block
    discharge_plot = px.line(title='Discharge Observations')
    stage_plot = px.line(title='Stage Observations')


    # Convert timestamp back to datetime
    if selected_time_range == None:
        selected_time_range = [slider_start,slider_end]

    start_date = pd.to_datetime(selected_time_range[0], unit='s')
    end_date = pd.to_datetime(selected_time_range[1], unit='s')

    filtered_gdf = gdf2[
        (gdf2['Reported Time2'] >= start_date) & (gdf2['Reported Time2'] <= end_date)
    ]

    # Build the Plotly scatter plot on the basemap with hover information
    fig = px.scatter_mapbox(filtered_gdf, lat=filtered_gdf.lat, lon=filtered_gdf.lon, color=selected_yaxis,
                            mapbox_style="carto-positron", title=f'Incidents reported during {str(start_date)[:7]} & {str(end_date)[:7]}: {selected_yaxis}',
                            hover_name=filtered_gdf.index, hover_data={selected_yaxis: True})
    fig.add_trace(go.Scattermapbox(
                lat=usgs_gages.lat,
                lon=usgs_gages.lon,
                mode='markers',
                text=usgs_gages.STAID,
                hoverinfo='text',
                marker=dict(size=10, color='black'),  # Adjust the marker size and color as needed
                opacity=0.3,
                name='USGS gages', # Set the legend name
            ))

    
    if click_data:
        #print(click_data,type(click_data))

        if 'lon' in click_data['points'][0]:
         
            # Extract clicked point number
            #point_number = click_data['points'][0]['pointNumber']
            AOI = Point(click_data['points'][0]['lon'],click_data['points'][0]['lat'])


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
        # Check if points are clicked on the map

   
   

    # Build the Plotly bar plot for unique entries count
    bar_plot = px.bar(filtered_gdf[selected_yaxis].value_counts(),
                      x=filtered_gdf[selected_yaxis].value_counts().index,
                      y=filtered_gdf[selected_yaxis].value_counts().values,
                      color=filtered_gdf[selected_yaxis].dropna().unique(),
                      title=f'Unique Entries Count: {selected_yaxis}',
                      labels={'y': 'Number of occurrences'})


    pie_plot = px.pie(values=filtered_gdf[selected_yaxis].value_counts(), names=filtered_gdf[selected_yaxis].dropna().unique())



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

                


                nearest_gage = gpd.sjoin_nearest(usgs_gages,aoi_gdf,max_distance=100)
                print(nearest_gage)


                disc_tmp = pd.DataFrame()

                if selected_time_range == None:
                    selected_time_range = [slider_start,slider_end]

                start_date = pd.to_datetime(selected_time_range[0], unit='s')
                end_date = pd.to_datetime(selected_time_range[1], unit='s')


                for stn in nearest_gage.STAID:
                    
                    #stn='01652500'
                    site = str(stn)
                    # get instantaneous values (iv), daily values (dv)
                    df = nwis.get_record(sites=site, service='iv', start=f'{str(start_date)[:7]}-01', end=f'{str(end_date)[:7]}-01') # change this with slider

                    
                    if '00060' in df.columns: # discharge
                        disc_tmp[site]=df['00060']

                if len(disc_tmp)>1:
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
                
                for stn in nearest_gage.STAID:
                    
                    #stn='01652500'
                    site = str(stn)
                    # get instantaneous values (iv)
                    df = nwis.get_record(sites=site, service='iv', start=f'{str(start_date)[:7]}-01', end=f'{str(end_date)[:7]}-01') # change this with slider

                    
                    if '00065' in df.columns: # discharge
                        stage_tmp[site]=df['00065']

                if len(stage_tmp)>1:
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



    my_cellStyle = {
        "styleConditions": [
            {
                "condition": f"params.colDef.field == '{selected_yaxis}'",
                "style": {"backgroundColor": "#d3d3d3"},
            },
            {
                "condition": f"params.colDef.field != '{selected_yaxis}'",
                "style": {"color": "black"}
            },
        ]
    }


    

    return fig, bar_plot, pie_plot, discharge_plot, stage_plot, {'cellStyle': my_cellStyle}, county_options, None




if __name__ == '__main__':
    app.run_server(debug=False, port=8002)
