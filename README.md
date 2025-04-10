# PM2.5-vertical-distribution-atmospheric-chemistry-
This repository provides a Python script to compute and visualize the  vertical profile of PM2.5 concentrations (in µg/m³) over a chosen location using NASA MERRA-2 3D aerosol reanalysis data.

The script processes NetCDF files for each month in a year , extracts relevant aerosol species and air density data at different pressure levels, and computes the total PM2.5 concentration. 

🌐 Overview
📌 Objective: To analyze how PM2.5 is vertically distributed in the atmosphere across different months.

🌍 Data Source: NASA MERRA-2 Reanalysis (NetCDF format)

🏙️ Region of Interest: User-defined — update latitude and longitude bounds to match specific cities (e.g., Delhi, Kolkata, Patna, Bhubaneswar).

📊 Output:

Line plot showing monthly PM2.5 vertical profiles over pressure levels (hPa)

Saved PNG figure

🔧 How It Works
User sets the latitude and longitude bounds for a target city/region (e.g., Bhubaneswar).

For each month:

The script loads daily MERRA-2 NetCDF files from the vp_2021/ folder

PM2.5 is computed using aerosol components (SO4, BC, OC, DU, SS) and AIRDENS

Daily and monthly averages are computed over the specified grid

curves are plotted together to visualize  variations in PM2.5 concentration across vertical atmospheric layers.

