#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 17:29:13 2021

This script uses as input the ERA5 u, v, w, lat, lon, and geopotential and
performs the listed tasks:

  (a) converts geopotential to geometric height a.m.s.l.
  (b) regrids u, v, w wind components into a regular deg2 grid.
  (c) create height bounds from geometric height
  (d) can output both seasonal and monthly values

There is a code duplication between the monthly and seasonal computation, but
due to simplicity of the program, there is no need for further optimization at
the moment.

@author: proestakis, thanasisn
"""

## TODO remove u_N and v_N from export at production

import os
import sys
import re
import glob
import calendar
from   datetime import datetime
import netCDF4  as nc
import numpy    as np
import xarray   as xr

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

## export by each month
MONTHLY  = False
MONTHLY  = cnf.D1.Monthly

## export by season of the year
SEASONAL = False
SEASONAL = cnf.D1.Seasonal


##  Check destination folder exists  -----------------------------------------
if not os.path.isdir(cnf.ERA5.path_regrid):
    sys.exit(f"\nFolder {cnf.ERA5.path_regrid} don't exist !!\n")

##  Choose input files  ------------------------------------------------------
##  Override random domain with target resolution boundaries
fl_North = Oc.border_up(  cnf.D1.North, cnf.D1.MaxLatStep,   90)
fl_South = Oc.border_down(cnf.D1.South, cnf.D1.MaxLatStep,  -90)
fl_East  = Oc.border_up(  cnf.D1.East,  cnf.D1.MaxLonStep,  180)
fl_West  = Oc.border_down(cnf.D1.West,  cnf.D1.MaxLonStep, -180)

filenames = glob.glob(f"{cnf.ERA5.path_raw}/ERA5_*_lat_{fl_South}_{fl_North}_lon_{fl_West}_{fl_East}.nc")

if len(filenames) < 1:
    sys.exit(f"\nNo input file found in {cnf.ERA5.path_raw} !!\n")


##  Construction of the output grid  ----------------------------------------

##  Make sure these bounds select the whole domain to regrid. Misalignments
##  will not be detected at run time, but v_N and u_N can be used to check.
REGRID_lat_centers = np.arange(fl_North, fl_South + 1, -cnf.D1.LatStep) - cnf.D1.LatStep / 2
REGRID_lon_centers = np.arange(fl_West,  fl_East,       cnf.D1.LonStep) + cnf.D1.LonStep / 2


##  Main regridding function!!  ---------------------------------------------
def regrid():
    """
    Regrid and store ERA5 data.  This is meant to be used only in this script
    to create monthly and seasonal aggregation of ERA5 data with a consistent
    method.  This function defines how ERA5 regridded data are named and stored.

    Returns
    -------
    Writes a new netcdf file with the regridded data.

    """

    # ##  Calculations with xarray  ------------------------------------------
    # ##  have to use a more general approach to do calculation on
    # ##  specified cells as needed with the iterative
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

    ##  Iterative calculations  ----------------------------------------------

    ##  Init target arrays
    height_lower   = np.empty((len(DTses.pressure_level), len(REGRID_lon_centers), len(REGRID_lat_centers)))
    height_mean    = np.empty((len(DTses.pressure_level), len(REGRID_lon_centers), len(REGRID_lat_centers)))
    height_upper   = np.empty((len(DTses.pressure_level), len(REGRID_lon_centers), len(REGRID_lat_centers)))
    u_total_N      = np.empty((len(DTses.pressure_level), len(REGRID_lon_centers), len(REGRID_lat_centers)))
    u_total_SD     = np.empty((len(DTses.pressure_level), len(REGRID_lon_centers), len(REGRID_lat_centers)))
    u_total_mean   = np.empty((len(DTses.pressure_level), len(REGRID_lon_centers), len(REGRID_lat_centers)))
    u_total_median = np.empty((len(DTses.pressure_level), len(REGRID_lon_centers), len(REGRID_lat_centers)))
    v_total_N      = np.empty((len(DTses.pressure_level), len(REGRID_lon_centers), len(REGRID_lat_centers)))
    v_total_SD     = np.empty((len(DTses.pressure_level), len(REGRID_lon_centers), len(REGRID_lat_centers)))
    v_total_mean   = np.empty((len(DTses.pressure_level), len(REGRID_lon_centers), len(REGRID_lat_centers)))
    v_total_median = np.empty((len(DTses.pressure_level), len(REGRID_lon_centers), len(REGRID_lat_centers)))

    ##  Compute stats in each cell  ------------------------------------------
    for ilon, lon in enumerate(REGRID_lon_centers):
        for jlat, lat in enumerate(REGRID_lat_centers):
            for klev, lev in enumerate(DTses.pressure_level):

                ## This includes each of the limit value two times in nearby
                ## cells, in order to pretty centre the cells.
                cell = DTses.where(
                    (DTses.longitude >= lon - cnf.D1.LonStep/2) &
                    (DTses.longitude <= lon + cnf.D1.LonStep/2) &
                    (DTses.latitude  >= lat - cnf.D1.LatStep/2) &
                    (DTses.latitude  <= lat + cnf.D1.LatStep/2) &
                    (DTses.pressure_level == lev),
                    drop = True)
                cell.coords

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

            ##  Create height bounds  ----------------------------------------
            ##  We assume the sort is always correct
            if not (np.diff(height_mean[:, ilon, jlat]) > 0).all():
                sys.exit("Descending heights")

            ##  Assign boundaries to array
            H, height_lower[:, ilon, jlat], height_upper[:, ilon, jlat] = \
                Oc.height_bounds(
                    heights       = height_mean[:, ilon, jlat],
                    remove_bottom = cnf.ERA5.extra_height_at_bottom,
                    add_top       = cnf.ERA5.extra_height_at_top,
                    quiet         = True)

    # ##  xarray export  -----------------------------------------------------
    # ## add time stamp to the dataset
    # res = res.expand_dims(time = [sesdate])
    # ## store data
    # comp = dict(zlib=True, complevel=5)
    # encoding = {var: comp for var in res.data_vars}
    # res.to_netcdf(fileout, mode = 'w', engine = "netcdf4", encoding = encoding)
    # print(f"Written: {fileout}")

    ##  Numpy array export  --------------------------------------------------
    ds = nc.Dataset(fileout, 'w', format='NETCDF4')

    ##  Define coordinates
    ds.createDimension('latitude',       len(REGRID_lat_centers))
    ds.createDimension('pressure_level', len(DTses.pressure_level))
    ds.createDimension('longitude',      len(REGRID_lon_centers))
    ds.createDimension('time',           1)
    ds.createDimension('time_span',      1)

    ##  Create variables data types
    U_N        = ds.createVariable('u_N',       np.float64, ('pressure_level', 'longitude', 'latitude',), zlib=True)  # will be removed
    U_SD       = ds.createVariable('u_SD',      np.float64, ('pressure_level', 'longitude', 'latitude',), zlib=True)
    U_mean     = ds.createVariable('u_mean',    np.float64, ('pressure_level', 'longitude', 'latitude',), zlib=True)
    U_median   = ds.createVariable('u_median',  np.float64, ('pressure_level', 'longitude', 'latitude',), zlib=True)
    V_N        = ds.createVariable('v_N',       np.float64, ('pressure_level', 'longitude', 'latitude',), zlib=True)  # will be removed
    V_SD       = ds.createVariable('v_SD',      np.float64, ('pressure_level', 'longitude', 'latitude',), zlib=True)
    V_mean     = ds.createVariable('v_mean',    np.float64, ('pressure_level', 'longitude', 'latitude',), zlib=True)
    V_median   = ds.createVariable('v_median',  np.float64, ('pressure_level', 'longitude', 'latitude',), zlib=True)
    height     = ds.createVariable('height',          'f4', ('pressure_level', 'longitude', 'latitude',), zlib=True)
    height_low = ds.createVariable('height_low',      'f4', ('pressure_level', 'longitude', 'latitude',), zlib=True)
    height_up  = ds.createVariable('height_up',       'f4', ('pressure_level', 'longitude', 'latitude',), zlib=True)
    lats       = ds.createVariable('latitude',        'f4', ('latitude', ), zlib=True)
    lons       = ds.createVariable('longitude',       'f4', ('longitude',), zlib=True)
    time       = ds.createVariable('time',      np.float64, ('time',))
    time_span  = ds.createVariable('time_span',   np.int32, ('time_span',))

    ##  Set units attributes
    U_SD.units       = 'm s**-1'
    U_mean.units     = 'm s**-1'
    U_median.units   = 'm s**-1'
    V_SD.units       = 'm s**-1'
    V_mean.units     = 'm s**-1'
    V_median.units   = 'm s**-1'
    height.units     = 'km'  ## conversion at the final step
    height_low.units = 'km'  ## conversion at the final step
    height_up.units  = 'km'  ## conversion at the final step
    lats.units       = 'degrees_north'
    lons.units       = 'degrees_east'
    time.calendar    = 'proleptic_gregorian'
    time.units       = 'seconds since 1970-01-01 00:00:00'  ## same as raw ERA5
    time_span.units  = 'month'

    ##  Set long name attribute
    U_SD.long_name       = 'U SD component of wind'
    U_mean.long_name     = 'U mean component of wind'
    U_median.long_name   = 'U median component of wind'
    V_SD.long_name       = 'V SD component of wind'
    V_mean.long_name     = 'V mean component of wind'
    V_median.long_name   = 'V median component of wind'
    height.long_name     = 'Height'
    height_low.long_name = 'Height of lower cell boundary'
    height_up.long_name  = 'Height of upper cell boundary'
    lats.long_name       = 'Latitude'
    lons.long_name       = 'Longitude'
    time.long_name       = 'Time'
    time_span.long_name  = 'The time span of the aggregated data'

    ##  Set standard name attribute
    U_SD.standard_name       = 'eastward_wind_SD'
    U_mean.standard_name     = 'eastward_wind_mean'
    U_median.standard_name   = 'eastward_wind_median'
    V_SD.standard_name       = 'northward_wind_SD'
    V_mean.standard_name     = 'northward_wind_mean'
    V_median.standard_name   = 'northward_wind_median'
    height.standard_name     = 'height'
    height_low.standard_name = 'lower_boundary'
    height_up.standard_name  = 'upper_boundary'
    lats.standard_name       = 'latitude'
    lons.standard_name       = 'longitude'
    time.standard_name       = 'time'
    time_span.standard_name  = 'duration'

    ##  Assign arrays to datasets
    U_N[:]        = u_total_N             # will be removed
    U_SD[:]       = u_total_SD
    U_mean[:]     = u_total_mean
    U_median[:]   = u_total_median
    V_N[:]        = v_total_N             # will be removed
    V_SD[:]       = v_total_SD
    V_mean[:]     = v_total_mean
    V_median[:]   = v_total_median
    height[:]     = height_mean  / 1000.  # to km
    height_low[:] = height_lower / 1000.  # to km
    height_up[:]  = height_upper / 1000.  # to km
    lats[:]       = REGRID_lat_centers    # centre of the cell
    lons[:]       = REGRID_lon_centers    # centre of the cell
    time[:]       = nc.date2num(sesdate, time.units)
    time_span[:]  = stats_duration

    ##  Set global attributes
    my_attrs = dict(title     = "Regridded ERA5 data",
                    type      = data_type,
                    season    = season,
                    details   = stats_message,
                    data_date = str(sesdate),
                    contacts  = cnf.OREO.contact_emails)
    for name, value in my_attrs.items():
        setattr(ds, name, value)

    ##  Do the actual data write
    ds.close()
    print(f"Written: {fileout}")


##  Some variables just to track progress
filesin_N = len(filenames)
filesin_C = 0
seaso_C = 0
month_C = 0
season  = " - "

##  Process raw ERA5 files  --------------------------------------------------
for filein in filenames:
    yyyy = int(re.compile('ERA5_([0-9]*)_.*.nc').search(filein).group(1))

    ##  Create output folders  -----------------------------------------------
    seasonl_dir = os.path.join(cnf.ERA5.path_regrid, f"Seasonal_{cnf.D1.LatStep}x{cnf.D1.LonStep}")
    monthly_dir = os.path.join(cnf.ERA5.path_regrid, f"Monthly_{cnf.D1.LatStep}x{cnf.D1.LonStep}")
    os.makedirs(seasonl_dir, exist_ok = True)
    os.makedirs(monthly_dir, exist_ok = True)

    ##  Limit data time range
    if not cnf.Range.start <= yyyy <= cnf.Range.until:
        continue

    filesin_C += 1  ## process counter
    print(f"\nProcessing {filesin_C}/{filesin_N}: {filein}")

    ##  load ERA5 main file on a xarray
    DT = xr.open_dataset(filein)

    ##  Get ERA5 spatial step
    lat_res = np.unique(np.diff(DT.latitude.values))[0]
    lon_res = np.unique(np.diff(DT.longitude.values))[0]

    ## apply a domain constrains
    # DT = DT.sel(longitude = slice(cnf.ERA5.West,  cnf.ERA5.East ),
    #             latitude  = slice(cnf.ERA5.North, cnf.ERA5.South))

    ## test plot
    # DT.u.isel(pressure_level = 0, valid_time = 0).plot()
    # DT.v.isel(pressure_level = 0, valid_time = 0).plot()

    ## geopotential changes with time!
    # DT.isel(valid_time = 0, latitude = 0, longitude = 0).z.values
    # DT.isel(valid_time = 1, latitude = 0, longitude = 0).z.values

    ##  Export with seasonal aggregation  ------------------------------------
    if SEASONAL:
        stats_message  = "Statistical values (mean, median, SD and N) refer to the seasonal values of 3 calendar months, with the approximate centre the day of 'time' variable."
        stats_duration = 3
        data_type      = "Seasonal"

        ##  Compute by season of the year
        seasons = ['Q1_DJF', 'Q2_MAM', 'Q3_JJA', 'Q4_SON']
        for season_idx, season in enumerate(seasons):
            ##  Output process info  -----------------------------------------
            seaso_C += 1
            complete = 100 * seaso_C / (filesin_N * 4)
            duration = datetime.now() - tic
            total    = duration * 100 / complete
            eta      = total - duration
            eda      = (datetime.now() + eta).replace(microsecond = 0)
            print(f"  {yyyy} {season} {seaso_C}/{filesin_N*4} {complete:6.2f}% D: {format(duration)} R: {format(eta)} F: {eda}")

            ##  Create output file path  -------------------------------------
            fileout = os.path.join(
                seasonl_dir,
                f"ERA5_{yyyy}_{season}_lat_{np.min(REGRID_lat_centers)}_{np.max(REGRID_lat_centers)}_lon_{np.min(REGRID_lon_centers)}_{np.max(REGRID_lon_centers)}.nc"
            )

            ##  Skip already existing files  ---------------------------------
            if (not FORCE) and (not Ou.output_needs_update(filein, fileout)):
                continue

            ##  Select data to use by season  --------------------------------
            if season == 'Q1_DJF':
                # File of previous year to read December for DJF season
                previous_file = list(filter(lambda x:'ERA5_' + str(yyyy - 1) in x, filenames))
                if len(previous_file)!=1:
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

            ##  Add geometric height  ----------------------------------------
            DTses = Oc.z_to_height(DTses)

            ##  Do the regrid and store file  --------------------------------
            regrid()
        ## end season iteration
    ## end seasonal aggregation

    ##  Export monthly aggregation  ------------------------------------------
    if MONTHLY:
        stats_message  = "Statistical values (mean, median, SD and N) refer to the monthly values starting from the first day of 'time' variable."
        stats_duration = 1
        data_type      = "Monthly"

        ##  Iterate all months of a year
        for m in range(1, 13):
            ##  Output process info  -----------------------------------------
            month_C += 1
            complete = 100 * month_C / (filesin_N * 12)
            duration = datetime.now() - tic
            total    = duration * 100 / complete
            eta      = total - duration
            eda      = (datetime.now() + eta).replace(microsecond = 0)
            print(f"  {yyyy} M {m:02} {month_C}/{filesin_N*12} {complete:6.2f}% D: {format(duration)} R: {format(eta)} F: {eda}")
            season   = calendar.month_name[m]

            ##  Create output file path  -------------------------------------
            fileout = os.path.join(
                monthly_dir,
                f"ERA5_{yyyy}_M{m:02}_lat_{np.min(REGRID_lat_centers)}_{np.max(REGRID_lat_centers)}_lon_{np.min(REGRID_lon_centers)}_{np.max(REGRID_lon_centers)}.nc"
            )

            ##  Skip already existing files  ---------------------------------
            if (not FORCE) and (not Ou.output_needs_update(filein, fileout)):
                continue

            ##  Select data to use explicitly  -------------------------------
            ## the date format for the ERA5 is already by month
            DTses   = DT.sel(valid_time = slice(f"{yyyy}-{m:02}-01"))
            sesdate = datetime(yyyy, m, 1)

            ##  Add geometric height  ----------------------------------------
            DTses = Oc.z_to_height(DTses)

            ##  Do the regrid and store file  --------------------------------
            regrid()
        ## end monthly loop
    ## end monthly aggregation
## end raw files iteration

#  SCRIPT END  ---------------------------------------------------------------
Ou.goodbye(cnf.LOGs.run, tic=tic, scriptname=__file__)
