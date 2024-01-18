# VDOT Dashboard

This repository contains the code for an interactive VDOT (Virginia Department of Transportation) Dashboard built using Dash and Plotly. The dashboard visualizes road closures and incidents reported in the state of Virginia, providing an intuitive interface for exploring the data.

## Features

- **Time Period Selection:** Users can select a specific time period using a range slider to filter incidents reported during a particular timeframe.

- **Region Selection:** The dashboard allows users to explore incidents by region, providing buttons for quick selection. Users can click on the region buttons to filter incidents for a specific region.

- **County Selection:** A dropdown menu enables users to choose specific counties for a more focused view of incidents.

- **Data Label Selection:** Users can choose the data label (Category, Road Type, etc.) to visualize on the map and other plots.

- **Event Category and Sub-Category Selection:** Users can further filter incidents based on event category and sub-category using dropdown menus.

- **Sync USGS Data:** A switch allows users to synchronize USGS (United States Geological Survey) data, providing additional information related to selected incidents.

- **Interactive Basemap:** The main map visualizes incidents on a mapbox with an option to click on specific points for detailed information.

- **Plots:** The dashboard includes various plots, such as bar plots, pie plots, discharge observations, and stage observations, providing insights into the selected incidents.

- **Data Grid:** A data grid displays detailed information about the selected incidents, allowing for easy exploration.

## Setup

1. Clone the repository:

```bash
git clone <repository-url>

python sample_dash_gis.py

The dashboard will be accessible at http://127.0.0.1:8002/ in your web browser.
```

## Requirements
Python 3.x
Dash
Plotly
Dash Bootstrap Components
Geopandas
Dash AG Grid
Pandas
Dataretrieval
Shapely
Dash Mantine Components

## Usage
Launch the application by running app.py.
Access the dashboard in your web browser.
Use the provided controls to interact with the data and explore incidents in Virginia.

## License
This project is licensed under the MIT License. Feel free to use and modify the code as needed.
