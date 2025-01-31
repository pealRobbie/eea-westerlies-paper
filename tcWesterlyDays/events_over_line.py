# %% [markdown]
# # events_over_line.ipynb
# 
# #### Robert Peal November 2024
# 
# - Timeseries of whether there is an event touching some line

import xarray as xr
import numpy as np
# import pandas as pd
# import sys
# sys.path.append("/home/atuin/c104fa/c104fa10/utils")
# import tctools2 as tct
import numpy as np
# import matplotlib.pyplot as plt
import pyproj
pyproj.datadir.set_data_dir("/home/atuin/c104fa/c104fa10/software/conda/envs/atmos_sci/share/proj") ## This line is needed to allow geopandas imports
# import geopandas as gpd
import argparse

def main():

    parser = argparse.ArgumentParser(description="Create a timeseries showing the days where there were westerlies and TC westerlies over the line")
    
    parser.add_argument(
        '--dataSaveLoc1D', 
        type=str, 
        required=True, 
        help="Location to save the 1-D timeseries to."
    )
    parser.add_argument(
        '--dataSaveLoc2D', 
        type=str, 
        required=True, 
        help="Location to save the 2-D timeseries to"
    )
    parser.add_argument(
        '--westerlyPath', 
        type=str, 
        required=True, 
        help="Path to the westerly event data."
    )

    parser.add_argument(
        '--lon', 
        type=str, 
        required=True, 
        help="Line longitude"
    )

    parser.add_argument(
        '--latMin', 
        type=str, 
        required=True, 
        help="Line min latitude"
    )

    parser.add_argument(
        '--latMax', 
        type=str, 
        required=True, 
        help="Line max latitude"
    )
    # parser.add_argument(
    #     '--swio_stateFile', 
    #     type=str, 
    #     required=True, 
    #     help="Path to the statefile"
    # )
    parser.add_argument(
        '--years', 
        type=int, 
        nargs='+', 
        required=True, 
        help="List of years (e.g., 2001 2002 2003)."
    )

    # parser.add_argument(
    #     '--calculate', 
    #     action='store_true', 
    #     help="Flag to enable calculation functionality."
    # )
    # parser.add_argument(
    #     '--figSaveLoc', 
    #     type=str, 
    #     help="Location to save the figure if plotting is required. If not provided, no plotting occurs."
    # )

    args = parser.parse_args()

    # Path to nc files with westerlies in
    westerlyPath = args.westerlyPath # "/home/atuin/c104fa/c104fa10/data/westerlies/eventDataTCv3/netcdfs/events_"

    # # Path to shapefile
    # region = "EEA"
    # shapefilePath = f"/home/atuin/c104fa/c104fa10/data/shapefiles/{region}.shp"

    # line specs
    lon=30
    latMin=-12
    latMax=5

    years = args.years
    months = np.arange(1,13)

    # dataSaveLoc
    dataSaveLoc1D = args.dataSaveLoc1D #f"/home/atuin/c104fa/c104fa10/eeaWesterliesPaper/tcWesterlyDays/eventsLineSeries.{lon:02.0f}E.{latMin:02.0f}_{latMax:02.0f}N.{years[0]}_{years[-1]}.csv"
    dataSaveLoc2D = args.dataSaveLoc2D #f"/home/atuin/c104fa/c104fa10/eeaWesterliesPaper/tcWesterlyDays/eventsLine.{lon:02.0f}E.{latMin:02.0f}_{latMax:02.0f}N.{years[0]}_{years[-1]}.nc"
    
    # #### Load the data
    print(f"Loading westerly events from {years[0]}-{years[-1]} from {westerlyPath} ... ")
    data = xr.open_mfdataset([f"{westerlyPath}{year}{month:02}.nc" for year in years for month in months])

    # #### Select the line
    print(f"Selecting points along {lon}E from {latMin}-{latMax}N")
    lineData = data.sel(latitude=slice(latMax,latMin))#,method="nearest")
    lineData = lineData.sel(longitude=lon,method="nearest")

    # #### Save days with a westerly on the line
    print("Identifying days with a westerly along the line ...")
    binaryLineData = xr.where(lineData > 0,1,0)
    lineWesterlyBinaryArray = binaryLineData.max(dim="latitude").compute()
    lineWesterlyBinaryDf = lineWesterlyBinaryArray.to_pandas()

    print(f"Saving 2D data to {dataSaveLoc2D}")
    lineData.to_netcdf(dataSaveLoc2D)

    print(f"Saving 1D data to {dataSaveLoc1D}")
    lineWesterlyBinaryDf.to_csv(dataSaveLoc1D)


if __name__=="__main__":

    # # Get the list from the command-line arguments
    # if len(sys.argv) < 2:
    #     print("Usage: python events_over_line.py <years>")
    #     sys.exit(1)

    # try:
    #     # Convert input string into a list of integers

    #     main([int(x) for x in sys.argv[1].split(',')])

    # except ValueError:
    #     print("Error: Please ensure all inputs are integers separated by commas.")
    #     sys.exit(1)
    main()

