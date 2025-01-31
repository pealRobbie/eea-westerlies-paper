"""
calculate_season_average.ipynb

Robert Peal November 2024

- Calculate the seasonal monthly average for some variable

- Generates data compatible with `plot_events_per_year.ipynb' 

"""

import xarray as xr
import numpy as np
import pyproj
pyproj.datadir.set_data_dir("/home/atuin/c104fa/c104fa10/software/conda/envs/atmos_sci/share/proj") ## This line is needed to allow geopandas imports


def main():

    ### User options ################

    # Years to use
    years = np.arange(1980,2023)

    # Path to nc files with westerlies in
    westerlyPath = "/home/atuin/c104fa/c104fa10/data/westerlies/eventDataTCv3/netcdfs/events_"

    # Path to save the seasonal averages
    dataSaveLoc = f"/home/atuin/c104fa/c104fa10/eeaWesterliesPaper/eventsPerYear/seasonAvgs.{years[0]}.{years[-1]}.nc"

    # Seasons
    seasons = [[1,2],[3,4,5],[6,7,8,9],[10,11,12]]
    seasonLabels = ["JF","MAM","JJAS","OND"]

    ### Seasonal average calculation #####

    print(f"Working with event masks from: {westerlyPath}")
    print(f"Working with years: {years}")

    # Empty dict for the data
    seasonAvgs = {}
    nyears = len(years) # Number of years

    # Iterate through the seasons
    for seasonix, season in enumerate(seasons):
        label = seasonLabels[seasonix]

        print(f"Working on season {label} ... ")
        # Load the westerly days and select the westerly mask
        westerlyDays = xr.open_mfdataset([f"{westerlyPath}{year}{month:02}.nc" for year in years for month in season])['westerlyMask'] 

        # Compute the average number of westerly days per month
        nmonths = len(season) * nyears # Total months
        seasonAvgs[label] = (westerlyDays.sum(dim="time").compute()/ nmonths)

    seasonAvgs = xr.Dataset(seasonAvgs)

    print(f"Saving to {dataSaveLoc}")
    seasonAvgs.to_netcdf(dataSaveLoc)

if __name__=="__main__":
    main()