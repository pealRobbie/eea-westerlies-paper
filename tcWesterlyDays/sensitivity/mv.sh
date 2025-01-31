#!/usr/bin/bash

dataSavePath1D="/home/atuin/c104fa/c104fa10/eeaWesterliesPaper/tcWesterlyDays/sensitivity/eventsLine1D."

# Define the parameters and their possible values
directionValues=(40 45 50)
magnitudePercentileValues=(60 70 80)
sizeValues=(750 1000 1500)

# Run the code for each value combination
for direction in "${directionValues[@]}"; do
  for percentage in "${magnitudePercentileValues[@]}"; do
    for size in "${sizeValues[@]}"; do

        # Define input arguments
        mv "${dataSavePath1D}dir.${direction}_thresh.${percentage}_size.${size}.nc" "${dataSavePath1D}dir.${direction}_thresh.${percentage}_size.${size}.csv"
    done
  done
done
