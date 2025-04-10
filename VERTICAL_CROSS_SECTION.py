# VERTICAL CROSS SECTION PLOT FOR MULTI-YEAR CLIMATIC MEAN OF PM2.5

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
import xarray as xr
from metpy.interpolate import cross_section
import warnings
from tqdm import tqdm

warnings.filterwarnings("ignore")

# === Global Settings and Constants ===
oc_scaling_factor = 1.6
so4_scaling_factor = 1.375

xsection_type = "lon"  # Options: "lon", "lat", "other"

# === Plot Settings ===
custom_cmap = cm.get_cmap("rainbow", 256)
color_array = custom_cmap(np.linspace(0, 1, 256))
white_color = np.array([1.0, 1.0, 1.0, 1])
color_array[:12, :] = white_color
cmap = ListedColormap(color_array)

fontsize = 14
plt.rc("xtick", labelsize=fontsize)
plt.rc("ytick", labelsize=fontsize)
plt.rc("legend", fontsize=fontsize)

# === Location and Domain Setup ===
station_labels = ["NTL"]
station_lats = [29.40]
station_lons = [79.50]

month_labels = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
                'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']

lat_min, lat_max, lon_min, lon_max = 18, 42, 65, 95
domain_extent = [lon_min, lon_max, lat_min, lat_max]

# === File Paths ===
pm25_data_dir = "./KGP/"
wind_data_dir = "./wb_2021/"

# === Initialize Monthly File Lists ===
monthly_pm_files = [[] for _ in range(12)]
monthly_wind_files = [[] for _ in range(12)]

# === Loop to Collect File Names from Multiple Years ===
for month_idx in range(1, 13):
    month_str = f"{month_idx:02d}"
    for yr in range(2020, 2022):
        timestamp = f"{yr}{month_str}"
        pm_files = [f for f in os.listdir(pm25_data_dir) if timestamp in f]
        wind_files = [f for f in os.listdir(wind_data_dir) if timestamp in f]
        monthly_pm_files[month_idx - 1].extend(pm_files)
        monthly_wind_files[month_idx - 1].extend(wind_files)

# === Loop Over Months to Generate Vertical Cross-Sections ===
for month_idx in range(12):

    pressure_month = []
    pm25_month = []
    u_wind_month = []
    v_wind_month = []

    for pm_file, wind_file in tqdm(list(zip(monthly_pm_files[month_idx], monthly_wind_files[month_idx]))):
        # Load data files
        pm_ds = xr.open_dataset(pm25_data_dir + pm_file)
        wind_ds = xr.open_dataset(wind_data_dir + wind_file)

        # Calculate full pressure grid (edge levels to center)
        shape = list(wind_ds.DELP.shape)
        shape[1] += 1
        pressure_full = np.zeros(shape)
        pressure_full[:, 0, :, :] = merra_top_pressure

        for level in range(len(wind_ds.lev)):
            pressure_full[:, level + 1, :, :] = pressure_full[:, level, :, :] + wind_ds.DELP[:, level, :, :]

        center_pressure = (pressure_full[:, :-1] + pressure_full[:, 1:]) / 2
        pm_ds["MER_Pres"] = (("time", "lev", "lat", "lon"), center_pressure)

        # Compute PM2.5 using constituent species
        total_pm25 = (so4_scaling_factor * pm_ds.SO4 + pm_ds.BCPHOBIC + pm_ds.BCPHILIC +
                      oc_scaling_factor * (pm_ds.OCPHILIC + pm_ds.OCPHOBIC) +
                      pm_ds.DU001 + pm_ds.DU002 * 0.38 +
                      pm_ds.SS001 + pm_ds.SS002 + pm_ds.SS003 * 0.83) * \
                     pm_ds.AIRDENS * 1e9
        pm_ds["pm25"] = total_pm25

        # Monthly mean
        pm_avg = pm_ds.resample(time="MS").mean(dim="time").metpy.parse_cf().squeeze()
        wind_avg = wind_ds.resample(time="MS").mean(dim="time").metpy.parse_cf().squeeze()

        # Define cross-section start and end
        if xsection_type == "other":
            start_pt = (24.0, 84.3)
            end_pt = (29.0, 85.7)
        elif xsection_type == "lon":
            start_pt = (20.0, 79.5)
            end_pt = (40.0, 79.5)
        elif xsection_type == "lat":
            start_pt = (29.4, 70.67)
            end_pt = (29.5, 89.5)

        pm_cross = cross_section(pm_avg, start_pt, end_pt).set_coords(("lat", "lon"))
        wind_cross = cross_section(wind_avg, start_pt, end_pt).set_coords(("lat", "lon"))

        # Prepare 2D coordinate arrays for plotting
        lat_grid = np.tile(pm_cross["lat"], (len(pm_avg.lev), 1))
        lon_grid = np.tile(pm_cross["lon"], (len(pm_avg.lev), 1))

        if xsection_type == "lat":
            x_coord_2d = lon_grid
            x_coord_1d = pm_cross["lon"]
            station_ticks = station_lons
        else:
            x_coord_2d = lat_grid
            x_coord_1d = pm_cross["lat"]
            station_ticks = station_lats

        # Collect cross-section data
        pressure_month.append(pm_cross["MER_Pres"] * 1e-2)
        pm25_month.append(pm_cross["pm25"])
        u_wind_month.append(wind_cross["U"])
        v_wind_month.append(wind_cross["V"])

    # === Monthly Averaging Across All Days ===
    pressure_avg = np.mean(np.array(pressure_month), axis=0)
    pm25_avg = np.mean(np.array(pm25_month), axis=0)
    u_wind_avg = np.mean(np.array(u_wind_month), axis=0)
    v_wind_avg = np.mean(np.array(v_wind_month), axis=0)

    # === Plot Section ===
    contour_levels = np.linspace(0, 100, 51)
    plot = plt.contourf(x_coord_2d,
                        pressure_avg,
                        pm25_avg,
                        contour_levels,
                        cmap=cmap,
                        extend="max")

    # Fill ground region
    plt.fill_between(x_coord_1d,
                     1.013e3,
                     pressure_avg[-1, :],
                     facecolor="black")

    plt.colorbar(plot,
                 label="PM2.5 (μg/m³)",
                 ticks=np.arange(0, 120, 20))

    plt.gca().invert_yaxis()
    plt.ylim(1000, 400)

    # Wind vectors
    vertical_idx = list(range(0, 72, 2))
    horizontal_idx = slice(5, 100, 8)

    plt.barbs(x_coord_2d[vertical_idx, horizontal_idx],
              pressure_avg[vertical_idx, horizontal_idx],
              u_wind_avg[vertical_idx, horizontal_idx],
              v_wind_avg[vertical_idx, horizontal_idx],
              color="k",
              length=5,
              pivot="middle",
              zorder=2)

    plt.xlabel('Latitude (°N)' if xsection_type != "lat" else 'Longitude (°E)')
    plt.ylabel('Pressure (hPa)')
    plt.title(f'Vertical Cross-section of Monthly Average PM2.5 - {month_labels[month_idx]}')

    plt.show()
