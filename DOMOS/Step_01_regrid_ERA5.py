#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 17:29:13 2021
This script uses as input the ERA5 u, v, w, lat, lon, and geopotential and
performs two tasks:
    (a) converts geopotential to geometric height a.m.s.l.
    (b) regrids u, v, w wind components into a regular 1x1 deg2 grid.
Here the script is applied for DOMOS domain.
@author: proestakis, thanasisn
"""

import os
import sys
import re
import glob
from   datetime import datetime
import netCDF4  as nc
import numpy    as np
import metpy.calc
import xarray as xr
from   metpy.units import units

# import cartopy.crs as ccrs
# import xarray_regrid


#  Load project functions  ---------------------------------------------------
sys.path.append("../")
import oreo_mod.utils as Ou
import oreo_mod.calc  as Oc
tic = datetime.now()

#  Load configuration profile by host name  ----------------------------------
config_file = f"../run_profiles/{os.uname()[1]}.yaml"
cnf = Ou.get_configs(config_file)

#  Force the reprocess of the inputs
FORCE = cnf.mode.Force

#  Check destination folder exists  ------------------------------------------
if not os.path.isdir(cnf.ERA5.path_regrid):
    sys.exit(f"\nFolder {cnf.ERA5.path_regrid} don't exist !!\n")

#  Choose input files  -------------------------------------------------------
#  Use expanded domain to find input files
fl_North = Oc.border_up(  cnf.D1.North, cnf.D1.LatStep)
fl_South = Oc.border_down(cnf.D1.South, cnf.D1.LatStep)
fl_East  = Oc.border_up(  cnf.D1.East,  cnf.D1.LonStep)
fl_West  = Oc.border_down(cnf.D1.West,  cnf.D1.LonStep)
filenames = glob.glob(f"{cnf.ERA5.path_raw}/ERA5_*_{fl_North}N{fl_South}S{fl_West}W{fl_East}E.nc")
filenames.sort()

if len(filenames) < 1:
    sys.exit(f"\nNo input file found in {cnf.ERA5.path_raw} !!\n")

# ESA-DOMOS: "... the full coverage of the Atlantic Ocean (including dust
# emission sources of Africa and S. America, the broader Atlantic Ocean,
# Caribbean Sea and Gulf of Mexico, confined between latitudes 40°N to 60°S),
# and of temporal coverage at least between 2010 and 2020".
# Therefore:
# (I)  DOMOS lon: -105E:5:25E
# (II) DOMOS lat:  -65N:5:45N
# (a wider domain is used here to account for (1) all fluxes and (2) the
# broader domain, and N.Atlandic Dust)

### DOMOS ####
lat_array       = np.arange(cnf.D1.South, cnf.D1.North)
lon_array       = np.arange(cnf.D1.West,  cnf.D1.East)
DOMOS_lat_array = np.arange(cnf.D1.North, cnf.D1.South, -cnf.D1.LatStep)
DOMOS_lon_array = np.arange(cnf.D1.West,  cnf.D1.East,   cnf.D1.LonStep)

# Process raw ERA5 files  ------------------------------------------------------------
for filein in filenames:
    yyyy = int(re.compile('ERA5_([0-9]*)_.*.nc').search(filein).group(1))

    ##  Limit data time range
    if not cnf.Range.start <= yyyy <= cnf.Range.until:
        continue

    print(f"\nProcessing: {filein}")

    ## load ERA5 main file on a xarray
    DT = xr.open_dataset(filein)

    ##  Get ERA5 spatial step
    lat_res = (np.unique(np.diff(DT.latitude.values))[0])
    lon_res = (np.unique(np.diff(DT.longitude.values))[0])

    ## apply a domain constrains
    # DT = DT.sel(longitude = slice(cnf.ERA5.West,  cnf.ERA5.East ),
    #             latitude  = slice(cnf.ERA5.North, cnf.ERA5.South))

    ## original grid
    # DT.longitude.values
    # DT.latitude.values
    # DT.longitude.values.min()
    # DT.longitude.values.max()
    # len(DT.latitude.values)
    # DT.dims
    # DT.info
    # DT.data_vars
    # DT.coords
    # DT.u
    # DT.z.units
    # DT.u.values
    # DT.u.attrs

    ## test plot
    DT.u.isel(pressure_level = 0, valid_time = 0).plot()
    DT.v.isel(pressure_level = 0, valid_time = 0).plot()

    ##  Compute by season of the year
    seasons = ['Q1_DJF', 'Q2_MAM', 'Q3_JJA', 'Q4_SON']
    for season_idx, season in enumerate(seasons):
        print(f"{yyyy} {season}")

        fileout = os.path.join(
            cnf.ERA5.path_regrid,
            f"ERA5_{yyyy}_{season}_{cnf.D1.North}N{cnf.D1.South}S{cnf.D1.West}W{cnf.D1.East}E.nc"
        )

        ## skip already existing files
        if (not FORCE) and (not Ou.output_needs_update(filein, fileout)):
            continue

        ## select data by season
        if season == 'Q1_DJF':
            # File of previous year to read December for DJF season
            previous_file = list(filter(lambda x:'ERA5_' + str(yyyy - 1) in x, filenames))
            if (len(previous_file)!=1):
                print("SKIP season! No file for the previous year found\n")
                continue

            ## load ERA5 main file on a xarray
            DTpre = xr.open_dataset(previous_file[0])
            ## keep only December of previous year
            DTpre = DTpre.sel(valid_time=f"{yyyy - 1}-12-01")
            ## select main data explicitly
            DTcur = DT.sel(valid_time=[f"{yyyy}-01-01", f"{yyyy}-02-01"])
            ## combine December with current
            DTses = xr.concat([DTpre, DTcur], dim = "valid_time")
            ## create a data stamp
            sesdate = datetime(yyyy, 1, 15)
            del DTpre
            del DTcur

        elif season == 'Q2_MAM':
            ## select main data explicitly
            DTses   = DT.sel(valid_time = slice(f"{yyyy}-03-01", f"{yyyy}-05-01"))
            sesdate = datetime(yyyy, 3, 15)

        elif season == 'Q3_JJA':
            ## select main data explicitly
            DTses   = DT.sel(valid_time = slice(f"{yyyy}-06-01", f"{yyyy}-08-01"))
            sesdate = datetime(yyyy, 7, 15)

        elif season == 'Q4_SON':
            ## select main data explicitly
            DTses   = DT.sel(valid_time = slice(f"{yyyy}-09-01", f"{yyyy}-11-01"))
            sesdate = datetime(yyyy, 10, 15)

        ##  Add geometric height  --------------------------------------------
        DTses = DTses.assign(
            height = xr.DataArray(
                metpy.calc.geopotential_to_height(
                     units.Quantity(DTses.z.values, DTses.z.units)
                ),
                coords = DTses.coords,
            )
        )
        DTses['height'].attrs = {
            'long_name':     'Geometric height',
            'units':         'm',
            'standard_name': 'height'
        }
        DTses.height.values

        # ##  Calculations with xarray  ------------------------------------
        # dd = DTses.coarsen(latitude  = int(-cnf.D1.LatStep / lat_res),
        #                    longitude = int( cnf.D1.LonStep / lon_res),
        #                    boundary  = "trim").mean(skipna = True)
        # ## mean of all months
        # res = dd.mean(dim = ["valid_time"])
        # ## coarsen with 'pad' may get under representation of values due to unequal bins
        # DTses = DTses.assign(valid_time = DTses.valid_time.dt.month)
        #
        # ## this should work
        # DTses.coarsen(latitude  = int(-cnf.D1.LatStep / lat_res),
        #               longitude = int( cnf.D1.LonStep / lon_res),
        #               valid_time = 3,
        #               boundary  = "trim").mean(skipna=True)

        ##  Iterative calculations  ------------------------------------------

        ## init target arrays
        u_total_mean   = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
        u_total_median = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
        u_total_SD     = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
        u_total_N      = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
        v_total_mean   = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
        v_total_median = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
        v_total_SD     = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
        v_total_N      = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
        # z_total_mean   = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
        height_mean    = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))

        ##  Compute stats in each cell
        for ilon, lon in enumerate(DOMOS_lon_array):
            for jlat, lat in enumerate(DOMOS_lat_array):
                for klev, lev in enumerate(DTses.pressure_level):

                    ## This includes each of the limits value two times in order to centre cells
                    cell = DTses.where(
                        (DTses.longitude >= lon) &
                        (DTses.longitude <= lon + cnf.D1.LonStep) &
                        (DTses.latitude  >= lat) &
                        (DTses.latitude  <= lat + cnf.D1.LatStep) &
                        (DTses.pressure_level == lev),
                        drop = True)

                    ## gather all statistics for each cell
                    u_total_mean  [klev, ilon, jlat] = np.mean(                   cell.u.values)
                    u_total_median[klev, ilon, jlat] = np.median(                 cell.u.values)
                    u_total_SD    [klev, ilon, jlat] = np.std(                    cell.u.values)
                    u_total_N     [klev, ilon, jlat] = np.count_nonzero(~np.isnan(cell.u.values))
                    v_total_mean  [klev, ilon, jlat] = np.mean(                   cell.v.values)
                    v_total_median[klev, ilon, jlat] = np.median(                 cell.v.values)
                    v_total_SD    [klev, ilon, jlat] = np.std(                    cell.v.values)
                    v_total_N     [klev, ilon, jlat] = np.count_nonzero(~np.isnan(cell.v.values))
                    # z_total_mean  [klev, ilon, jlat] = np.mean(                   cell.z.values)
                    height_mean   [klev, ilon, jlat] = np.mean(                   cell.height.values)

        # ##  xarray export  ---------------------------------------------------
        # ## add time stamp to the dataset
        # res = res.expand_dims(time = [sesdate])
        # ## store data
        # comp = dict(zlib=True, complevel=5)
        # encoding = {var: comp for var in res.data_vars}
        # res.to_netcdf(fileout, mode = 'w', engine = "netcdf4", encoding = encoding)
        # print(f"Written: {fileout}")

        ##  numpy array export  -----------------------------------------------
        ds           = nc.Dataset(fileout, 'w', format='NETCDF4')

        ##  Define coordinates
        lev          = ds.createDimension('pressure_level', len(DTses.pressure_level))
        latitude     = ds.createDimension('latitude',       len(DOMOS_lat_array))
        longitude    = ds.createDimension('longitude',      len(DOMOS_lon_array))

        DTses.coords
        ##  Create variables data types
        lats         = ds.createVariable('latitude',        'f4', ('latitude', ),                    zlib=True)
        lons         = ds.createVariable('longitude',       'f4', ('longitude',),                    zlib=True)
        height       = ds.createVariable('height',          'f4', ('pressure_level', 'longitude', 'latitude',), zlib=True)
        U_mean       = ds.createVariable('u_mean',    np.float64, ('pressure_level', 'longitude', 'latitude',), zlib=True)
        U_median     = ds.createVariable('u_median',  np.float64, ('pressure_level', 'longitude', 'latitude',), zlib=True)
        U_SD         = ds.createVariable('u_SD',      np.float64, ('pressure_level', 'longitude', 'latitude',), zlib=True)
        U_N          = ds.createVariable('u_N',       np.float64, ('pressure_level', 'longitude', 'latitude',), zlib=True)
        V_mean       = ds.createVariable('v_mean',    np.float64, ('pressure_level', 'longitude', 'latitude',), zlib=True)
        V_median     = ds.createVariable('v_median',  np.float64, ('pressure_level', 'longitude', 'latitude',), zlib=True)
        V_SD         = ds.createVariable('v_SD',      np.float64, ('pressure_level', 'longitude', 'latitude',), zlib=True)
        V_N          = ds.createVariable('v_N',       np.float64, ('pressure_level', 'longitude', 'latitude',), zlib=True)

        ##  Set units attributes
        lats.units     = 'degrees_north'
        lons.units     = 'degrees_east'
        height.units   = 'm'
        U_mean.units   = 'm s**-1'
        U_median.units = 'm s**-1'
        U_SD.units     = 'm s**-1'
        U_N.units      = ''
        V_mean.units   = 'm s**-1'
        V_median.units = 'm s**-1'
        V_SD.units     = 'm s**-1'
        V_N.units      = ''

        ##  Set long name attribute
        lats.long_name     = 'Latitude'
        lons.long_name     = 'Longitude'
        height.long_name   = 'Height'
        U_mean.long_name   = 'U mean component of wind'
        U_median.long_name = 'U median component of wind'
        U_SD.long_name     = 'U SD component of wind'
        U_N.long_name      = 'U component of wind count'
        V_mean.long_name   = 'V mean component of wind'
        V_median.long_name = 'V median component of wind'
        V_SD.long_name     = 'V SD component of wind'
        V_N.long_name      = 'V component of wind count'

        ##  Set standard name attribute
        height.standard_name   = 'height'
        U_mean.standard_name   = 'eastward_wind_mean'
        U_median.standard_name = 'eastward_wind_median'
        U_SD.standard_name     = 'eastward_wind_SD'
        U_N.standard_name      = 'eastward_wind_count'
        V_mean.standard_name   = 'northward_wind_mean'
        V_median.standard_name = 'northward_wind_median'
        V_SD.standard_name     = 'northward_wind_SD'
        V_N.standard_name      = 'northward_wind_count'

        ##  Assign arrays to datasets
        lats[:]     = DOMOS_lat_array + cnf.D1.LatStep / 2
        lons[:]     = DOMOS_lon_array + cnf.D1.LonStep / 2
        height[:]   = height_mean
        U_mean[:]   = u_total_mean
        U_median[:] = u_total_median
        U_SD[:]     = u_total_SD
        U_N[:]      = u_total_N
        V_mean[:]   = v_total_mean
        V_median[:] = v_total_median
        V_SD[:]     = v_total_SD
        V_N[:]      = v_total_N

        ds.close()
        print(f"Written: {fileout}")

    # # https://confluence.ecmwf.int/display/CKB/ERA5%3A+What+is+the+spatial+reference
    # # ERA longitude from 0->360 deg to -180->180 deg.
    # dataset_latitude     = dataset['latitude'][:]
    # ERA5_idx_lat         = np.where((dataset_latitude >= lat_array[0]) & (dataset_latitude <= lat_array[-1]+1))
    # ERA5_idx_lat         = np.ravel(ERA5_idx_lat)
    # dataset_latitude     = dataset['latitude'][ERA5_idx_lat]

    # dataset_longitude    = dataset['longitude'][:]

    # for lon in dataset_longitude:
    #     if lon >= 180:
    #         dataset_longitude[np.where(lon == dataset_longitude)] = lon - 360
    #         print("lon >= 180")

    # ERA5_idx_lon         = np.where((dataset_longitude >= lon_array[0]) & (dataset_longitude <= lon_array[-1]+1))
    # ERA5_idx_lon         = np.ravel(ERA5_idx_lon)
    # dataset_longitude    = dataset['longitude'][ERA5_idx_lon]

    # for lon in dataset_longitude:
    #     if lon >= 180:
    #         dataset_longitude[np.where(lon == dataset_longitude)] = lon - 360
    #         print("lon >= 180")

    # dataset_time         = dataset['valid_time'][:]
    # dataset_level        = dataset['pressure_level'][:]
    # dataset_u            = dataset['u'][:,:,ERA5_idx_lat, ERA5_idx_lon]
    # dataset_v            = dataset['v'][:,:,ERA5_idx_lat, ERA5_idx_lon]


    # # ERA convert Geopotenial to geometric height (a.m.s.l.):
    # # https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.geopotential_to_height.html
    # dataset_geo_units    = dataset['z'].units
    # dataset_geopotential = dataset['z'][:,:,ERA5_idx_lat, ERA5_idx_lon]
    # dataset_geopotential = units.Quantity(dataset_geopotential, dataset_geo_units)
    # dataset_height       = metpy.calc.geopotential_to_height(dataset_geopotential)
    # dataset_height       = np.array(dataset_height)

    # # Year_of_Interest_previous = mypath + '\\' + str(int(filename[0:4])-1)+'.nc'
    # dataset_ΙΙ              = nc.Dataset(previous_file[0])
    # dataset_u_ΙΙ            = dataset_ΙΙ['u'][:,:,ERA5_idx_lat,ERA5_idx_lon]
    # dataset_v_ΙΙ            = dataset_ΙΙ['v'][:,:,ERA5_idx_lat,ERA5_idx_lon]
    # dataset_geopotential_ΙΙ = dataset_ΙΙ['z'][:,:,ERA5_idx_lat,ERA5_idx_lon]
    # dataset_geopotential_ΙΙ = units.Quantity(dataset_geopotential_ΙΙ,dataset_geo_units)
    # dataset_height_ΙΙ       = metpy.calc.geopotential_to_height(dataset_geopotential_ΙΙ)
    # dataset_height_ΙΙ       = np.array(dataset_height_ΙΙ)


    # seasons = ['Q1_DJF', 'Q2_MAM', 'Q3_JJA', 'Q4_SON']
    # for season_idx, season in enumerate(seasons):


    #     # initializing: u-mean,SD / v-mean,SD / w-mean,SD / z-mean for 1x1 deg2 grid resolution.
    #     u_total_mean = np.empty((len(DOMOS_lon_array), len(DOMOS_lat_array), len(dataset_level)))
    #     u_total_SD   = np.empty((len(DOMOS_lon_array), len(DOMOS_lat_array), len(dataset_level)))
    #     v_total_mean = np.empty((len(DOMOS_lon_array), len(DOMOS_lat_array), len(dataset_level)))
    #     v_total_SD   = np.empty((len(DOMOS_lon_array), len(DOMOS_lat_array), len(dataset_level)))
    #     z_total_mean = np.empty((len(DOMOS_lon_array), len(DOMOS_lat_array), len(dataset_level)))

    #     # computing and saving: u-mean,SD / v-mean,SD / z-mean for 2x5 deg2 grid resolution.
    #     for lon in DOMOS_lon_array:

    #         idx_lon        = np.where((dataset_longitude >= lon) & (dataset_longitude <= lon + cnf.ERA5.LonStep))
    #         idx_lon        = np.ravel(idx_lon)
    #         temp_u         = dataset_u[:,:,:,idx_lon]
    #         temp_v         = dataset_v[:,:,:,idx_lon]
    #         temp_height    = dataset_height[:,:,:,idx_lon]
    #         temp_u_II      = dataset_u_ΙΙ[:,:,:,idx_lon]
    #         temp_v_II      = dataset_v_ΙΙ[:,:,:,idx_lon]
    #         temp_height_II = dataset_height_ΙΙ[:,:,:,idx_lon]

    #         for lat in DOMOS_lat_array:

    #             idx_lat = np.where((dataset_latitude >= lat) & (dataset_latitude <= lat + cnf.ERA5.LatStep))
    #             idx_lat = np.ravel(idx_lat)

    #             if season == 'Q1_DJF':
    #                 Months_of_Interest_idx = [0, 1]
    #                 u_I  = temp_u[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
    #                 v_I  = temp_v[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
    #                 z_I  = temp_height[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
    #                 u_II = temp_u_II[11:12,:,idx_lat,:]
    #                 v_II = temp_v_II[11:12,:,idx_lat,:]
    #                 z_II = temp_height_II[11:12,:,idx_lat,:]
    #                 u    = np.concatenate([u_I,u_II], axis=0)
    #                 v    = np.concatenate([v_I,v_II], axis=0)
    #                 z    = np.concatenate([z_I,z_II], axis=0)
    #             if season == 'Q2_MAM':
    #                 Months_of_Interest_idx = [3, 5]
    #                 u = temp_u[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
    #                 v = temp_v[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
    #                 z = temp_height[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
    #             if season == 'Q3_JJA':
    #                 Months_of_Interest_idx = [6, 8]
    #                 u = temp_u[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
    #                 v = temp_v[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
    #                 z = temp_height[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
    #             if season == 'Q4_SON':
    #                 Months_of_Interest_idx = [9, 11]
    #                 u = temp_u[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
    #                 v = temp_v[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
    #                 z = temp_height[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]

    #             u_total = [u[:,i,:,:].mean() for i in range(u.shape[1])]
    #             v_total = [v[:,i,:,:].mean() for i in range(v.shape[1])]
    #             z_total = [z[:,i,:,:].mean() for i in range(z.shape[1])]
    #             for idx_level,lev in enumerate(dataset_level):
    #                 u_total_mean[np.where(lon == DOMOS_lon_array),np.where(lat == DOMOS_lat_array),idx_level] = u_total[idx_level]
    #                 v_total_mean[np.where(lon == DOMOS_lon_array),np.where(lat == DOMOS_lat_array),idx_level] = v_total[idx_level]
    #                 z_total_mean[np.where(lon == DOMOS_lon_array),np.where(lat == DOMOS_lat_array),idx_level] = z_total[idx_level]

    #             u_total = [u[:,i,:,:].std() for i in range(u.shape[1])]
    #             v_total = [v[:,i,:,:].std() for i in range(v.shape[1])]
    #             z_total = [z[:,i,:,:].std() for i in range(z.shape[1])]

    #             for idx_level,lev in enumerate(dataset_level):
    #                 u_total_SD[np.where(lon == DOMOS_lon_array),np.where(lat == DOMOS_lat_array),idx_level] = u_total[idx_level]
    #                 v_total_SD[np.where(lon == DOMOS_lon_array),np.where(lat == DOMOS_lat_array),idx_level] = v_total[idx_level]



#  SCRIPT END  ---------------------------------------------------------------
Ou.goodbye(cnf.LOGs.run, tic=tic, scriptname=__file__)
