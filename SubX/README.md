# S2S File Read and Write for ADCIRC

## Overview

This Python script reads and writes S2S (Subseasonal-to-Seasonal) files for ADCIRC simulations. The code focuses on processing wind data and mean sea level pressure (MSLP) from different models and ensembles.

## Prerequisites

- Python 3.x
- Required Python packages: `matplotlib`, `pandas`, `numpy`, `netCDF4`, `mpl_toolkits.basemap`

## Usage

1. **Install Dependencies:**
    ```bash
    pip install matplotlib pandas numpy netCDF4 basemap
    ```

2. **Configure User Inputs:**
    - Modify the `model_sets` list with the desired S2S models.
    - Update the `variable` list with the desired meteorological variables (e.g., `pslmsl`, `uas10m`, `vas10m`).
    - Update the `AnalysisDates` dictionary with the specific analysis dates for each S2S model.
    - Set the `StormNamenYear` variable for the storm name and year.

3. **Run the Script:**
    ```bash
    python script_name.py
    ```

4. **View Plots:**
    - The script generates plots for wind magnitude and mean sea level pressure for each specified S2S model and ensemble.

## Acknowledgments

- The script uses various Python libraries for data processing and visualization.
- It is designed for processing ADCIRC simulations and S2S model output.

## Notes

- Ensure that the required Python packages are installed before running the script.
- This script serves as an example of processing and visualizing S2S model data. Customize it further based on specific requirements.

Feel free to modify and expand the script as needed for your use case. If you encounter any issues or have suggestions for improvement, feel free to [open an issue](https://github.com/your-username/your-repository/issues) or submit a pull request.

## License

This project is licensed under the [MIT License](LICENSE).
