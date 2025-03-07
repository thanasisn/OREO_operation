#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 17:29:13 2021
This script uses as input the ERA5 u, v, w, lat, lon, and geopotential and performs two tasks:
    (a) converts geopotential to geometric height a.m.s.l.
    (b) regrids u, v, w wind components into a regular 1x1 deg2 grid.
Here the script is applied for DOMOS domain.
@author: proestakis
"""

import os
import sys
from   datetime import datetime, timedelta
from   dotmap   import DotMap
import yaml
import glob
import netCDF4  as nc
import numpy    as np
import metpy.calc
from   metpy.units import units
from   os          import walk
# import math

SCRIPT_NAME = __file__
tic         = datetime.now()

#  Load configuration by host name  --------------------------------------------
config_file = os.uname()[1] + '.yaml'
with open(config_file, 'r') as config_fl:
    configs = yaml.safe_load(config_fl)
    # Convert dictionary to use dot notation
    cnf = DotMap(configs)
    print("\nRead config file:", config_file, "\n")

#  Check destination folder has been created  ----------------------------------
if not os.path.isdir(cnf.ERA5.path_regrid):
    sys.exit("\nFolder " + cnf.ERA5.path_regrid + " don't exist!\n")




# setting up input and output path of ERA5 files.
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels-monthly-means?tab=overview

# CLIMPACT
mypath      = cnf.ERA5.path_raw
output_path = cnf.ERA5.path_regrid

# checking if input files exist in input dir.
filenames = glob.glob(cnf.ERA5_input.path_raw + '/ERA5_*.nc')
if len(filenames) < 1 :
    sys.exit("\nNo input file found in " + co.ERA5_input.path_raw + " !\n")

# DOMOS
# mypath      = r"D:\PC-backup\Projects\DOMOS\Datasets\2x5_deg_seasonal_mean\ERA5\ERA5_pre_processed"
# output_path = r"D:\PC-backup\Projects\DOMOS\Datasets\2x5_deg_seasonal_mean\ERA5\processed"

# ESA-DOMOS: "... the full coverage of the Atlantic Ocean (including dust emission sources of Africa and S. America, 
# the broader Atlantic Ocean, Caribbean Sea and Gulf of Mexico, confined between latitudes 40°N to 60°S), and of 
# temporal coverage at least between 2010 and 2020".
# Therefore:
# (I)  DOMOS lon: -105E:5:25E
# (II) DOMOS lat:  -65N:5:45N
# (a wider domain is used here to account for (1) all fluxes and (2) the broader domain, and N.Atlandic Dust)

sys.exit("wait")

### CLIMPACT II
lon_array       = np.arange(-10, 45)
lat_array       = np.arange( 30, 50)
DOMOS_lon_array = np.arange(-10, 45, 5)
DOMOS_lat_array = np.arange( 30, 50, 2)

# ## DOMOS ####
# lon_array       = np.arange(-125, 25)
# lat_array       = np.arange( -60, 42)
# DOMOS_lon_array = np.arange(-125, 25, 5) 
# DOMOS_lat_array = np.arange( -62, 42, 2)






for count_filename, filename in enumerate(filenames):

    file        = mypath + '\\' + filename
    dataset     = nc.Dataset(file)
    if int(filename[0:4]) < 2007:
        continue
    print(file)
    

    # https://confluence.ecmwf.int/display/CKB/ERA5%3A+What+is+the+spatial+reference
    # last visit: 03/02/2022.
    # ERA longitude from 0->360 deg to -180->180 deg.    
    dataset_latitude     = dataset['latitude'][:]
    ERA5_idx_lat         = np.where((dataset_latitude >= lat_array[0]) & (dataset_latitude <= lat_array[-1]+1))
    ERA5_idx_lat         = np.ravel(ERA5_idx_lat)
    dataset_latitude     = dataset['latitude'][ERA5_idx_lat]

    dataset_longitude    = dataset['longitude'][:] 
    for lon in dataset_longitude:
        if lon >= 180:
            dataset_longitude[np.where(lon == dataset_longitude)] = lon - 360
    ERA5_idx_lon         = np.where((dataset_longitude >= lon_array[0]) & (dataset_longitude <= lon_array[-1]+1))
    ERA5_idx_lon         = np.ravel(ERA5_idx_lon)
    dataset_longitude    = dataset['longitude'][ERA5_idx_lon] 
    for lon in dataset_longitude:
        if lon >= 180:
            dataset_longitude[np.where(lon == dataset_longitude)] = lon - 360

    dataset_time         = dataset['time'][:]
    dataset_level        = dataset['level'][:]
    dataset_u            = dataset['u'][:,:,ERA5_idx_lat,ERA5_idx_lon]
    dataset_v            = dataset['v'][:,:,ERA5_idx_lat,ERA5_idx_lon]

    # ERA convert Geopotenial to geometric height (a.m.s.l.):
    # https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.geopotential_to_height.html
    # last access: 03/02/2022. 
    dataset_geo_units    = dataset['z'].units
    dataset_geopotential = dataset['z'][:,:,ERA5_idx_lat, ERA5_idx_lon]
    dataset_geopotential = units.Quantity(dataset_geopotential, dataset_geo_units)
    dataset_height       = metpy.calc.geopotential_to_height(dataset_geopotential)
    dataset_height       = np.array(dataset_height)

    # File of previous year - to read December for DJF season
    Year_of_Interest_previous = mypath + '\\' + str(int(filename[0:4])-1)+'.nc'
    dataset_ΙΙ = nc.Dataset(Year_of_Interest_previous)
    dataset_u_ΙΙ            = dataset_ΙΙ['u'][:,:,ERA5_idx_lat,ERA5_idx_lon]
    dataset_v_ΙΙ            = dataset_ΙΙ['v'][:,:,ERA5_idx_lat,ERA5_idx_lon]
    dataset_geopotential_ΙΙ = dataset_ΙΙ['z'][:,:,ERA5_idx_lat,ERA5_idx_lon]
    dataset_geopotential_ΙΙ = units.Quantity(dataset_geopotential_ΙΙ,dataset_geo_units)
    dataset_height_ΙΙ       = metpy.calc.geopotential_to_height(dataset_geopotential_ΙΙ)
    dataset_height_ΙΙ       = np.array(dataset_height_ΙΙ)

    # ERA time:
    # units     = hours since 1900-01-01 00:00:00.0 
    # long_name = time 
    # calendar  = gregorian
    # loop to produce seasonal-mean files

    seasons = ['DJF', 'MAM', 'JJA', 'SON']
    for season_idx, season in enumerate(seasons):

        # initializing: u-mean,SD / v-mean,SD / w-mean,SD / z-mean for 1x1 deg2 grid resolution.
        u_total_mean = np.empty((len(DOMOS_lon_array), len(DOMOS_lat_array), len(dataset_level)))
        u_total_SD   = np.empty((len(DOMOS_lon_array), len(DOMOS_lat_array), len(dataset_level)))
        v_total_mean = np.empty((len(DOMOS_lon_array), len(DOMOS_lat_array), len(dataset_level)))
        v_total_SD   = np.empty((len(DOMOS_lon_array), len(DOMOS_lat_array), len(dataset_level)))
        z_total_mean = np.empty((len(DOMOS_lon_array), len(DOMOS_lat_array), len(dataset_level)))

        # computing and saving: u-mean,SD / v-mean,SD / z-mean for 2x5 deg2 grid resolution.
        for lon in DOMOS_lon_array:

            idx_lon        = np.where((dataset_longitude >= lon) & (dataset_longitude <= lon+5))
            idx_lon        = np.ravel(idx_lon)
            temp_u         = dataset_u[:,:,:,idx_lon]
            temp_v         = dataset_v[:,:,:,idx_lon]
            temp_height    = dataset_height[:,:,:,idx_lon]
            temp_u_II      = dataset_u_ΙΙ[:,:,:,idx_lon]
            temp_v_II      = dataset_v_ΙΙ[:,:,:,idx_lon]
            temp_height_II = dataset_height_ΙΙ[:,:,:,idx_lon]

            for lat in DOMOS_lat_array:

                idx_lat = np.where((dataset_latitude >= lat) & (dataset_latitude <= lat+2))
                idx_lat = np.ravel(idx_lat)

                if season == 'DJF':
                    Months_of_Interest_idx = [0, 1]
                    u_I  = temp_u[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
                    v_I  = temp_v[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
                    z_I  = temp_height[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
                    u_II = temp_u_II[11:12,:,idx_lat,:]
                    v_II = temp_v_II[11:12,:,idx_lat,:]
                    z_II = temp_height_II[11:12,:,idx_lat,:]
                    u    = np.concatenate([u_I,u_II], axis=0)                    
                    v    = np.concatenate([v_I,v_II], axis=0) 
                    z    = np.concatenate([z_I,z_II], axis=0) 
                if season == 'MAM':
                    Months_of_Interest_idx = [3, 5]
                    u = temp_u[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
                    v = temp_v[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
                    z = temp_height[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
                if season == 'JJA':
                    Months_of_Interest_idx = [6, 8]
                    u = temp_u[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
                    v = temp_v[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
                    z = temp_height[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
                if season == 'SON':
                    Months_of_Interest_idx = [9, 11]
                    u = temp_u[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
                    v = temp_v[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
                    z = temp_height[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]

                u_total = [u[:,i,:,:].mean() for i in range(u.shape[1])]
                v_total = [v[:,i,:,:].mean() for i in range(v.shape[1])]
                z_total = [z[:,i,:,:].mean() for i in range(z.shape[1])]
                for idx_level,lev in enumerate(dataset_level):
                    u_total_mean[np.where(lon == DOMOS_lon_array),np.where(lat == DOMOS_lat_array),idx_level] = u_total[idx_level]
                    v_total_mean[np.where(lon == DOMOS_lon_array),np.where(lat == DOMOS_lat_array),idx_level] = v_total[idx_level]
                    z_total_mean[np.where(lon == DOMOS_lon_array),np.where(lat == DOMOS_lat_array),idx_level] = z_total[idx_level]   

                u_total = [u[:,i,:,:].std() for i in range(u.shape[1])]
                v_total = [v[:,i,:,:].std() for i in range(v.shape[1])]
                z_total = [z[:,i,:,:].std() for i in range(z.shape[1])]
                for idx_level,lev in enumerate(dataset_level):
                    u_total_SD[np.where(lon == DOMOS_lon_array),np.where(lat == DOMOS_lat_array),idx_level] = u_total[idx_level]
                    v_total_SD[np.where(lon == DOMOS_lon_array),np.where(lat == DOMOS_lat_array),idx_level] = v_total[idx_level]

        ######################################################################
        #  --- Saving ERA u, v, w, height monthly mean dataset as NetCDF --- #
        ######################################################################

        # creating nc. filename and initializing:
        fn           = output_path + '\\' + filename[0:4] + '_' + season + '.nc'
        ds           = nc.Dataset(fn, 'w', format='NETCDF4')

        # create nc. dimensions:
        longitude    = DOMOS_lon_array + 2.5
        latitude     = DOMOS_lat_array + 1
        lev          = ds.createDimension('lev', len(dataset_level))
        lat          = ds.createDimension('lat', len(latitude)) 
        lon          = ds.createDimension('lon', len(longitude))

        # create nc. variables:
        lats         = ds.createVariable('Latitude', 'f4',   ('lat',),             zlib=True)
        lons         = ds.createVariable('Longitude','f4',   ('lon',),             zlib=True)
        Height       = ds.createVariable('Height',   'f4',   ('lon','lat','lev',), zlib=True)
        U            = ds.createVariable('U',    np.float64, ('lon','lat','lev',), zlib=True)
        U_SD         = ds.createVariable('U_SD', np.float64, ('lon','lat','lev',), zlib=True)
        V            = ds.createVariable('V',    np.float64, ('lon','lat','lev',), zlib=True)
        V_SD         = ds.createVariable('V_SD', np.float64, ('lon','lat','lev',), zlib=True)

        # nc. variables' units
        lats.units   = 'degrees_north'
        lons.units   = 'degrees_east'
        Height.units = 'm'     
        U.units      = 'm s**-1'
        U_SD.units   = 'm s**-1'
        V.units      = 'm s**-1'
        V_SD.units   = 'm s**-1'

        # nc. variables' "long names":
        lats.long_name   = 'Latitude'
        lons.long_name   = 'Longitude'
        Height.long_name = 'Height'
        U.long_name      = 'U component of wind'
        U_SD.long_name   = 'U component of wind SD'
        V.long_name      = 'V component of wind'
        V_SD.long_name   = 'V component of wind SD'

        # nc. variables' "standard names":
        Height.standard_name = 'height'
        U.standard_name      = 'eastward_wind'
        U_SD.standard_name   = 'eastward_wind_SD'
        V.standard_name      = 'northward_wind'
        V_SD.standard_name   = 'northward_wind_SD'

        # nc. saving datasets
        lats[:]      = latitude
        lons[:]      = longitude
        Height[:]    = z_total_mean
        U[:]         = u_total_mean
        U_SD[:]      = u_total_SD
        V[:]         = v_total_mean
        V_SD[:]      = v_total_SD

        ds.close()

        print("End of: " + fn)


# SCRIPT END
out  = datetime.now().strftime("%F %T") + " "
out += os.getlogin() + "@" + os.uname()[1] + " "
out += SCRIPT_NAME + " "
out += str(round((datetime.now() - tic).total_seconds() / 60.0, 2)) + " mins"
print('\n' + out + '\n')
with open(co.LOGs.run, 'a') as runlog:
    runlog.write(out + '\n')
