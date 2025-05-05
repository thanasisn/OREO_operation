#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy    as np
import os
import sys
import glob
import tqdm
import time
import pathlib
import subprocess

##  Load project functions  --------------------------------------------------
sys.path.append("../")
import oreo_mod.utils as Ou
import oreo_mod.calc  as Oc



# Ou.source_code_hash("/home/athan/OREO/operation/DOMOS/Step_01_regrid_ERA5.py")


##  Load configuration profile by host name  ---------------------------------
cnf = Ou.get_configs(
        Ou.parse_arguments(run_profiles_folder = "../run_profiles").profile
    )



from netCDF4 import Dataset
import numpy as np

# Group names and shape settings
groups = ['regionA', 'regionB']
n_alt = 10
n_lat = 20
n_lon = 30

# Coordinate data
altitude = np.linspace(0, 10000, n_alt)
latitudes = np.linspace(-90, 90, n_lat)
longitudes = np.linspace(-180, 180, n_lon)

# First, create root-level lat/lon once
with Dataset('spatial_altitude_data.nc', 'w') as nc_file:
    # Define dimensions in the root group
    nc_file.createDimension('lat', n_lat)
    nc_file.createDimension('lon', n_lon)

    # Create root-level coordinate variables
    lat_var = nc_file.createVariable('latitude', 'f4', ('lat',))
    lon_var = nc_file.createVariable('longitude', 'f4', ('lon',))

    lat_var.units = 'degrees_north'
    lon_var.units = 'degrees_east'

    lat_var[:] = latitudes
    lon_var[:] = longitudes

# Then, open the file iteratively to add groups with only altitude + data
for group_name in groups:
    with Dataset('spatial_altitude_data.nc', 'a') as nc_file:
        if group_name in nc_file.groups:
            print(f"Group '{group_name}' already exists. Skipping.")
            continue

        group = nc_file.createGroup(group_name)

        # Only define altitude locally in each group
        group.createDimension('altitude', n_alt)

        # Create altitude variable in group
        alt_var = group.createVariable('altitude', 'f4', ('altitude',))
        alt_var.units = 'meters'
        alt_var[:] = altitude

        # Link data variable to root dimensions (lat, lon) and local altitude
        data_var = group.createVariable('some_variable', 'f4', ('altitude', 'lat', 'lon'))
        data_var.units = 'arbitrary_units'
        data_var[:, :, :] = np.random.rand(n_alt, n_lat, n_lon)

        print(f"Group '{group_name}' written with data and linked to shared lat/lon.")


for group_name in groups:

    nc_file = Dataset("mytestfile.nc", 'w', format='NETCDF4')

    # Define dimensions in the root group
    nc_file.createDimension('lat', n_lat)
    nc_file.createDimension('lon', n_lon)

    # Create root-level coordinate variables
    lat_var = nc_file.createVariable('latitude', 'f4', ('lat',))
    lon_var = nc_file.createVariable('longitude', 'f4', ('lon',))

    lat_var.units = 'degrees_north'
    lon_var.units = 'degrees_east'

    lat_var[:] = latitudes
    lon_var[:] = longitudes

    if not group_name in nc_file.groups:
        print(f"Group '{group_name}' already exists. Skipping.")

        group = nc_file.createGroup(group_name)

        # Only define altitude locally in each group
        group.createDimension('altitude', n_alt)

        # Create altitude variable in group
        alt_var = group.createVariable('altitude', 'f4', ('altitude',))
        alt_var.units = 'meters'
        alt_var[:] = altitude

        # Link data variable to root dimensions (lat, lon) and local altitude
        data_var = group.createVariable('some_variable', 'f4', ('altitude', 'lat', 'lon'))
        data_var.units = 'arbitrary_units'
        data_var[:, :, :] = np.random.rand(n_alt, n_lat, n_lon)

        print(f"Group '{group_name}' written with data and linked to shared lat/lon.")

    nc_file.close()



# for outer in [10, 20, 30, 40, 50]:
#     print(outer,"\n")
#     for inner in tqdm.tqdm(range(outer), desc="This file", dynamic_ncols = True):
#         print(inner,"\n")
#         time.sleep(1)
# print("done!")
#

# xlon = -62.5
# ylat = 25
# TT = Ou.backscatter_percentiles_lookup(ylat,
#                                     xlon,
#                                    datafile  = cnf.LIVAS.percentiles_fl,
#                                    ylat_step = cnf.D1.LatStep,
#                                    xlon_step = cnf.D1.LonStep)
