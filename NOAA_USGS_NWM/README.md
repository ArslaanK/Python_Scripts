# NOAA Water Level Data Downloader

This Python script downloads hourly water level data from the NOAA API for specified stations and time ranges. The data is then processed, and the water levels are plotted along with the specified return periods.

## Prerequisites

- Python 3.x
- Required Python packages: `matplotlib`, `fileinput`, `datetime`, `pandas`, `os`, `requests`, `scipy`, `csv`, `numpy`

## Usage

1. **Install Dependencies:**
    ```bash
    pip install matplotlib fileinput pandas requests scipy csv numpy
    ```

2. **Configure User Inputs:**
    - Modify the `rps` list with the desired return periods.
    - Update the `aep_erdc` dictionary with the AEP values for each station.
    - Specify the start and end years for each station in the `strt_end_years_gages` dictionary.

3. **Run the Script:**
    ```bash
    python script_name.py
    ```

4. **View Plots:**
    - The script will generate plots displaying the water levels for the specified NOAA stations.
    - AEP lines for each return period will be included in the plots.

## Acknowledgments

- The script utilizes the NOAA API for retrieving water level data.

## Notes

- Ensure that the required Python packages are installed before running the script.
- The script provides an example of downloading, processing, and plotting NOAA water level data. Customize it further based on specific requirements.

Feel free to modify and expand the script as needed for your use case. If you encounter any issues or have suggestions for improvement, feel free to [open an issue](https://github.com/your-username/your-repository/issues) or submit a pull request.

## License

This project is licensed under the [MIT License](LICENSE).
