# -*- coding: utf-8 -*-
"""
Created on Thu Jun 19 16:39:55 2025

@author: akhalid
"""


import pandas as pd
import h5py
import numpy as np



def process_adcirc_spatial_timeseries(file_path, start_date=None, end_date=None):
    """
    Reads an ADCIRC Excel file, optionally clips the dataset to the specified time window,
    converts the TIME column to a fraction of a day, and computes a cumulative time that continues across days.

    Parameters:
    - file_path: Path to the Excel file.
    - start_date: (Optional) Start date in ISO format (e.g., '2005-10-04T00:00:00').
    - end_date: (Optional) End date in ISO format (e.g., '2005-10-06T23:59:59').

    Returns:
    - DataFrame with Cumulative_Seconds column added at first column, all other data stays same or clipped to timewindow.
    """
    
    df = pd.read_excel(file_path, engine='openpyxl')
    df['TIME'] = pd.to_datetime(df['TIME'])

    if start_date:
        df = df[df['TIME'] >= pd.to_datetime(start_date)]
    if end_date:
        df = df[df['TIME'] <= pd.to_datetime(end_date)]
        
        
    # finds fraction of the day in seconds
    df['Fraction_of_Day'] = (
        df['TIME'].dt.hour / 24 +
        df['TIME'].dt.minute / (24 * 60) +
        df['TIME'].dt.second / (24 * 3600)
    )

    # create a cumulative sum column that adds 1 after every next day
    df['Cumulative_Sum'] = 0.0
    previous_cumulative_sum = 0.0
    previous_day = df['TIME'].dt.date.iloc[0]

    for i in range(len(df)):
        current_day = df['TIME'].dt.date.iloc[i]
        if current_day != previous_day:
            previous_cumulative_sum += 1.0
            previous_day = current_day
        df.loc[i, 'Cumulative_Sum'] = previous_cumulative_sum + df.loc[i, 'Fraction_of_Day']

    df_out = df.drop(columns=['TIME', 'Fraction_of_Day'])
    
    # bring the cummulative sum column as the first column
    cols = ['Cumulative_Sum'] + [col for col in df_out.columns if col != 'Cumulative_Sum']
    
    
    df_out = df_out[cols]
    df_out = df_out.rename(columns={'Cumulative_Sum': 'Cumulative_Seconds'})
    
    print('Data in csv:')
    print(df_out.head())
    return df_out


def replace_stage_hydrograph_data(hdf5_path, new_data, bc_index=0, fill_value=-99999.0):
    """
    Replaces a Stage Hydrograph dataset in an HDF5 file with processed ADCIRC data.

    Parameters:
    - hdf5_path: Path to the HDF5 file.
    - df_processed: DataFrame containing the new data to insert.
    - bc_index: Index of the boundary condition to replace (default is 0).
    - fill_value: Value to fill in the new array where data is missing (default is -99999.0).
    """
    
    try:
        
        with h5py.File(hdf5_path, 'r+') as base_data:
            bc_group = base_data['Event Conditions']['Unsteady']['Boundary Conditions']['Stage Hydrographs']
            bc_names = list(bc_group.keys())
            print("\n================\nFound Boundary Condition ", bc_names)
    
            adcirc_array = np.array(new_data, dtype=float)
            existing_data = bc_group[bc_names[bc_index]][:]
    
            # Create a new array with the same shape as existing data
            empty_array = np.full(existing_data.shape, fill_value)
            empty_array[:adcirc_array.shape[0], :adcirc_array.shape[1]] = adcirc_array
    
            # Preserve the first column (e.g., time or ID)
            empty_array[:, 0] = existing_data[:, 0]
    
            # Replace the dataset
            bc_group[bc_names[bc_index]][:] = empty_array
    
            # Close file
            base_data.close()
            print(f'Data replaced in hdf\n================')

    except:
        print(f'Error replacing the data in hdf')


#%%

# =============================================================================
# Main
# =============================================================================


# file paths
root_adcirc_csv = r'V:/projects/p00860_coj_2023_cf_jf/02_analysis/HECRAS_250515HECRASV66COMPOUND/CSV/combined_with_time_points_TammyxlRE.xlsx'
model_directory = r'V:\projects\p00860_coj_2023_cf_jf\02_analysis\Test_files_for_RASv66\model_run_files_copy1'
model_name = 'COJCOMPOUNDCOMPUTET.p01.hdf'

# preprocess adcirc data
adcirc_data = process_adcirc_spatial_timeseries(
                    root_adcirc_csv,
                    start_date="2005-10-04T00:00:00", # or simply set this to None
                    end_date="2005-10-05T23:59:59" # or simply set this to None
                )

# replace adcirc data in hdf file
replace_stage_hydrograph_data(
    hdf5_path=rf'{model_directory}/{model_name}',
    new_data=adcirc_data
)
