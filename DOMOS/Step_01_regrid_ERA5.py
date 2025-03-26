#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 17:29:13 2021

This script uses as input the ERA5 u, v, w, lat, lon, and geopotential and
performs the listed tasks:

  (a) converts geopotential to geometric height a.m.s.l.
  (b) regrids u, v, w wind components into a regular deg2 grid.
  (c) can output both seasonal and monthly values

There is a code duplication between the monthly and seasonal computation, but
due to simplicity of the program there is no need for further optimization at
the moment.

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

##  Load project functions  --------------------------------------------------
sys.path.append("../")
import oreo_mod.utils as Ou
import oreo_mod.calc  as Oc
tic = datetime.now()

##  Load configuration profile by host name  ---------------------------------
config_file = f"../run_profiles/{os.uname()[1]}.yaml"
cnf = Ou.get_configs(config_file)

##  Set switches  ------------------------------------------------------------

## overwrite output files
FORCE = cnf.mode.Force
FORCE = True

## export my each moth
MONTHLY  = False
MONTHLY  = cnf.D1.Monthly

## export be season
SEASONAL = False
# SEASONAL = cnf.D1.Seasonal


##  Check destination folder exists  ------------------------------------------
if not os.path.isdir(cnf.ERA5.path_regrid):
    sys.exit(f"\nFolder {cnf.ERA5.path_regrid} don't exist !!\n")

##  Choose input files  -------------------------------------------------------
##  Use expanded domain to find input files
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

##  Construction of the output grid  ------------------------------------------
lat_array       = np.arange(cnf.D1.South, cnf.D1.North)
lon_array       = np.arange(cnf.D1.West,  cnf.D1.East)
DOMOS_lat_array = np.arange(cnf.D1.North, cnf.D1.South, -cnf.D1.LatStep)
DOMOS_lon_array = np.arange(cnf.D1.West,  cnf.D1.East,   cnf.D1.LonStep)


def Get_Boundaries_of_ERA5_UVW_wind_components(Longitude, Latitude, Height):

    Height_Boundaries    = np.array(np.empty((Height.shape[0], Height.shape[1], Height.shape[2]+1)))
    Height_Boundaries[:] = np.nan
    for count_lon,lon in enumerate(Longitude):
        for count_lat,lat in enumerate(Latitude):
            H       = Height[count_lon,count_lat,:].data
            temp    = np.array(np.empty((H.shape[0]+1)))
            temp[:] = np.nan
            if H[H.shape[0]-1] > 0:
                for count_alt,alt in enumerate(temp):
                    if count_alt == 0:
                        Height_Boundaries[count_lon,count_lat,0] = H[0] + 1000.0
                    if count_alt == temp.shape[0]-1:
                        Height_Boundaries[count_lon,count_lat,count_alt] = 0.0
                    if ((count_alt > 0) & (count_alt < temp.shape[0]-1)):
                        Height_Boundaries[count_lon,count_lat,count_alt] = (H[count_alt-1] + H[count_alt])/2
            else:
                for count_alt,alt in enumerate(temp):
                    if count_alt == 0:
                        Height_Boundaries[count_lon,count_lat,0] = H[0] + 1000.0
                    if count_alt == temp.shape[0]-1:
                        Height_Boundaries[count_lon,count_lat,count_alt] = -100.0
                    if count_alt == temp.shape[0]-2:
                        Height_Boundaries[count_lon,count_lat,count_alt] =    0.0
                    if ((count_alt > 0) & (count_alt < temp.shape[0]-2)):
                        Height_Boundaries[count_lon, count_lat, count_alt] = (H[count_alt-1] + H[count_alt])/2

    return(Height_Boundaries)





##  Process raw ERA5 files  --------------------------------------------------
for filein in filenames:
    yyyy = int(re.compile('ERA5_([0-9]*)_.*.nc').search(filein).group(1))

    ##  Create ouput folders -------------------------------------------------
    seasonl_dir = os.path.join(cnf.ERA5.path_regrid, f"Seasonal_{cnf.D1.LatStep}x{cnf.D1.LonStep}")
    monthly_dir = os.path.join(cnf.ERA5.path_regrid, f"Monthly_{cnf.D1.LatStep}x{cnf.D1.LonStep}")
    os.makedirs(seasonl_dir, exist_ok = True)
    os.makedirs(monthly_dir, exist_ok = True)

    ##  Limit data time range
    if not cnf.Range.start <= yyyy <= cnf.Range.until:
        continue

    print(f"\nProcessing: {filein}")

    ##  load ERA5 main file on a xarray
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


    ##  Export seasonal aggregation  ----------------------------------------
    if SEASONAL:

        ##  Compute by season of the year
        seasons = ['Q1_DJF', 'Q2_MAM', 'Q3_JJA', 'Q4_SON']
        for season_idx, season in enumerate(seasons):
            print(f"  {yyyy} {season}")

            fileout = os.path.join(
                seasonl_dir,
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
                    print(f"SKIP season! No file for the {yyyy - 1} found\n")
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
            # DTses = DTses.assign(
            #     height = xr.DataArray(
            #         metpy.calc.geopotential_to_height(
            #              units.Quantity(DTses.z.values, DTses.z.units)
            #         ),
            #         coords = DTses.coords,
            #     )
            # )
            # DTses['height'].attrs = {
            #     'long_name':     'Geometric height',
            #     'units':         'm',
            #     'standard_name': 'height'
            # }

            DTses = Oc.z_to_height(DTses)
            if DTses.height.min() < 0 :
                sys.exit("\nNegative height found!! Have to resolve this!!\n")


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

            ##  Init target arrays
            u_total_mean   = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
            u_total_median = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
            u_total_SD     = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
            u_total_N      = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
            v_total_mean   = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
            v_total_median = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
            v_total_SD     = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
            v_total_N      = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
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

            temp_h = height[:,0,0]


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


    ##  Export monthly aggregation  ------------------------------------------
    if MONTHLY:

        ##  Compute by month of the year
        for m in range(1, 13):
            print(f"  {yyyy} M {m:02}")

            fileout = os.path.join(
                monthly_dir,
                f"ERA5_{yyyy}_M{m:02}_{cnf.D1.North}N{cnf.D1.South}S{cnf.D1.West}W{cnf.D1.East}E.nc"
            )

            ## skip already existing files
            if (not FORCE) and (not Ou.output_needs_update(filein, fileout)):
                continue

            ## geopotential changes with time
            DT.isel(valid_time = 0, latitude = 0, longitude = 0).z.values
            DT.isel(valid_time = 1, latitude = 0, longitude = 0).z.values


            ## select data explicitly
            ## the date format for the ERA5 is already by month
            DTses   = DT.sel(valid_time = slice(f"{yyyy}-{m:02}-01"))
            sesdate = datetime(yyyy, m, 1)


            ##  Add geometric height  ----------------------------------------
            # DTses = DTses.assign(
            #     height = xr.DataArray(
            #         metpy.calc.geopotential_to_height(
            #              units.Quantity(DTses.z.values, DTses.z.units)
            #         ),
            #         coords = DTses.coords,
            #     )
            # )
            # DTses['height'].attrs = {
            #     'long_name':     'Geometric height',
            #     'units':         'm',
            #     'standard_name': 'height'
            # }

            DTses = Oc.z_to_height(DTses)
            if DTses.height.min() < 0 :
                sys.exit("\nNegative height found!! Have to resolve this!!\n")


            ## ignoring time
            Height_Boundaries = np.empty((DTses.height.longitude.shape[0],
                                          DTses.height.latitude.shape[0],
                                          DTses.height.pressure_level.shape[0] + 1))
            Height_Boundaries.fill(np.nan)

            # Height_Boundaries    = np.array(np.empty((Height.shape[0], Height.shape[1], Height.shape[2]+1)))
            # Height_Boundaries[:] = np.nan
            for count_lon, lon in enumerate(DTses.height.longitude):
                for count_lat, lat in enumerate(DTses.height.latitude):
                    # H       = DTses.height[count_lon, count_lat, :].data
                    # temp    = np.array(np.empty((H.shape[0]+1)))
                    # temp[:] = np.nan

                    H    = DTses.height.isel(valid_time = 0).sel(latitude = lat, longitude = lon).values
                    temp = np.full(H.shape[0] + 1, np.nan)

                    if H[H.shape[0]-1] > 0:
                        print("first")
                        for count_alt, alt in enumerate(temp):
                            if count_alt == 0:
                                Height_Boundaries[count_lon, count_lat, 0] = H[0] + 1000.0
                            if count_alt == temp.shape[0]-1:
                                Height_Boundaries[count_lon, count_lat, count_alt] = 0.0
                            if ((count_alt > 0) & (count_alt < temp.shape[0]-1)):
                                Height_Boundaries[count_lon,count_lat, count_alt] = (H[count_alt-1] + H[count_alt])/2
                    else:
                        print("second")
                        for count_alt, alt in enumerate(temp):
                            if count_alt == 0:
                                Height_Boundaries[count_lon, count_lat, 0] = H[0] + 1000.0
                            if count_alt == temp.shape[0]-1:
                                Height_Boundaries[count_lon,count_lat,count_alt] = -100.0
                            if count_alt == temp.shape[0]-2:
                                Height_Boundaries[count_lon,count_lat,count_alt] =    0.0
                            if ((count_alt > 0) & (count_alt < temp.shape[0]-2)):
                                Height_Boundaries[count_lon, count_lat, count_alt] = (H[count_alt-1] + H[count_alt])/2

                    print(Height_Boundaries[count_lon, count_lat, :])

                    H[H.max() == H][0]
                    HB = np.full(H.shape[0] + 1, np.nan)
                    HB[0] = 0.0
                    HB[H.shape[0]] = H[H.max() == H][0] + 1000.

                    ## assume the same sort
                    if (np.diff(H) > 0).all():
                        print("Ascending heights")
                    else:
                        sys.stop("Descending heights")

                    if not H[H.min() == H].item() > cnf.ERA5.height_at_bottom:
                        sys.stop("Starting point is not lower than cell centre")

                    HB = H + np.diff(np.concat((H, H[H.max() == H] + cnf.ERA5.extra_height_on_top)))/2
                    HB = np.concat(([cnf.ERA5.height_at_bottom], HB))

                    # len(H)
                    # for ih, h in enumerate(H):
                    #     print(ih, h)
                    # H[1:] + np.diff(H)/2

                    # np.diff(H)/2
                    # np.diff(H, append = H[H.max() == H].item())/2
                    # H[H.max() == H].item()

                    for i,j in enumerate(H):
                        print(H[i], HB[i], HB[i+1])


                    sys.exit("HH")

            ##  Iterative calculations  --------------------------------------

            ##  Init target arrays
            u_total_mean   = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
            u_total_median = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
            u_total_SD     = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
            u_total_N      = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
            v_total_mean   = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
            v_total_median = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
            v_total_SD     = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
            v_total_N      = np.empty((len(DTses.pressure_level), len(DOMOS_lon_array), len(DOMOS_lat_array)))
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
                        height_mean   [klev, ilon, jlat] = np.mean(                   cell.height.values)




            ##  numpy array export  ------------------------------------------
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


#  SCRIPT END  ---------------------------------------------------------------
Ou.goodbye(cnf.LOGs.run, tic=tic, scriptname=__file__)
