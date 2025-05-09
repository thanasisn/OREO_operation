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
from   datetime import datetime, timezone
import netCDF4  as nc
import numpy    as np
import xarray   as xr
from multiprocessing import Pool
import tqdm

##  Load project functions  --------------------------------------------------
sys.path.append("../")
import oreo_mod.utils as Ou
import oreo_mod.calc  as Oc
tic = datetime.now()

##  Load configuration profile by host name  ---------------------------------
cnf = Ou.get_configs(
        Ou.parse_arguments(run_profiles_folder = "../run_profiles").profile
    )
##  Track the source code status that created each output  -------------------
TRACE = Ou.source_code_hash(__file__)

##  Set switches  ------------------------------------------------------------

## overwrite output files
FORCE = cnf.mode.Force

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



def process_ERA5_file(filein):
    """
    Process an ERA5 data file and produce output files.
    This functions is adapted to run in parallel.
    """

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

        ##  Iterative calculations  ----------------------------------------------

        coords_3d = (len(REGRID_lon_centers), len(REGRID_lat_centers), len(DTses.pressure_level))
        ##  Init target arrays
        height_lower   = np.full(coords_3d, np.nan)
        height_mean    = np.full(coords_3d, np.nan)
        height_upper   = np.full(coords_3d, np.nan)
        u_total_N      = np.full(coords_3d, np.nan)
        u_total_SD     = np.full(coords_3d, np.nan)
        u_total_mean   = np.full(coords_3d, np.nan)
        u_total_median = np.full(coords_3d, np.nan)
        v_total_N      = np.full(coords_3d, np.nan)
        v_total_SD     = np.full(coords_3d, np.nan)
        v_total_mean   = np.full(coords_3d, np.nan)
        v_total_median = np.full(coords_3d, np.nan)

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
                    u_total_mean  [ilon, jlat, klev] = np.mean(                   cell.u.values)
                    u_total_median[ilon, jlat, klev] = np.median(                 cell.u.values)
                    u_total_SD    [ilon, jlat, klev] = np.std(                    cell.u.values)
                    u_total_N     [ilon, jlat, klev] = np.count_nonzero(~np.isnan(cell.u.values))
                    v_total_mean  [ilon, jlat, klev] = np.mean(                   cell.v.values)
                    v_total_median[ilon, jlat, klev] = np.median(                 cell.v.values)
                    v_total_SD    [ilon, jlat, klev] = np.std(                    cell.v.values)
                    v_total_N     [ilon, jlat, klev] = np.count_nonzero(~np.isnan(cell.v.values))
                    height_mean   [ilon, jlat, klev] = np.mean(                   cell.height.values)

                ##  Create height bounds  ----------------------------------------
                ##  We assume the sort is always correct
                if not (np.diff(height_mean[ilon, jlat, :]) > 0).all():
                    sys.exit("Descending heights")

                ##  Assign boundaries to array
                H, height_lower[ilon, jlat, :], height_upper[ilon, jlat, :] = \
                    Oc.height_bounds(
                        heights       = height_mean[ilon, jlat, :],
                        remove_bottom = cnf.ERA5.extra_height_at_bottom,
                        add_top       = cnf.ERA5.extra_height_at_top,
                        quiet         = True)

        ##  Numpy array export  --------------------------------------------------
        ds = nc.Dataset(fileout, 'w', format = 'NETCDF4')

        ##  Define coordinates
        ds.createDimension('latitude',       len(REGRID_lat_centers))
        ds.createDimension('longitude',      len(REGRID_lon_centers))
        ds.createDimension('pressure_level', len(DTses.pressure_level))
        ds.createDimension('time',           1)
        ds.createDimension('time_span',      1)

        ##  Create variables data types
        U_N        = ds.createVariable('u_N',       np.float64, ('longitude', 'latitude', 'pressure_level', ), zlib=True)  # will be removed
        U_SD       = ds.createVariable('u_SD',      np.float64, ('longitude', 'latitude', 'pressure_level', ), zlib=True)
        U_mean     = ds.createVariable('u_mean',    np.float64, ('longitude', 'latitude', 'pressure_level', ), zlib=True)
        U_median   = ds.createVariable('u_median',  np.float64, ('longitude', 'latitude', 'pressure_level', ), zlib=True)
        V_N        = ds.createVariable('v_N',       np.float64, ('longitude', 'latitude', 'pressure_level', ), zlib=True)  # will be removed
        V_SD       = ds.createVariable('v_SD',      np.float64, ('longitude', 'latitude', 'pressure_level', ), zlib=True)
        V_mean     = ds.createVariable('v_mean',    np.float64, ('longitude', 'latitude', 'pressure_level', ), zlib=True)
        V_median   = ds.createVariable('v_median',  np.float64, ('longitude', 'latitude', 'pressure_level', ), zlib=True)
        height     = ds.createVariable('height',          'f4', ('longitude', 'latitude', 'pressure_level', ), zlib=True)
        height_low = ds.createVariable('height_low',      'f4', ('longitude', 'latitude', 'pressure_level', ), zlib=True)
        height_up  = ds.createVariable('height_up',       'f4', ('longitude', 'latitude', 'pressure_level', ), zlib=True)
        lats       = ds.createVariable('latitude',        'f4', ('latitude', ), zlib=True)
        lons       = ds.createVariable('longitude',       'f4', ('longitude',), zlib=True)
        time       = ds.createVariable('time',      np.float64, ('time',))
        time_span  = ds.createVariable('time_span', np.int32, ('time_span',))

        ##  Set units attributes
        U_SD.units       = 'm s**-1'
        U_mean.units     = 'm s**-1'
        U_median.units   = 'm s**-1'
        V_SD.units       = 'm s**-1'
        V_mean.units     = 'm s**-1'
        V_median.units   = 'm s**-1'
        height.units     = 'm'
        height_low.units = 'm'
        height_up.units  = 'm'
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
        height[:]     = height_mean
        height_low[:] = height_lower
        height_up[:]  = height_upper
        lats[:]       = REGRID_lat_centers    # centre of the cell
        lons[:]       = REGRID_lon_centers    # centre of the cell
        time[:]       = nc.date2num(sesdate, time.units)
        time_span[:]  = stats_duration

        ##  Set global attributes
        my_attrs = dict(title          = "Regridded ERA5 data",
                        type           = data_type,
                        season         = season,
                        details        = stats_message,
                        data_date      = str(sesdate),
                        contacts       = cnf.OREO.contact_emails,
                        creation       = format(datetime.now(timezone.utc)),
                        user_host      = os.getlogin() + "@" + os.uname()[1],
                        source_version = TRACE)
        for name, value in my_attrs.items():
            setattr(ds, name, value)

        ##  Do the actual data write
        ds.close()
        print(f"Written: {fileout}")



    ##  Start work for a file   ---------------------------------------------
    yyyy = int(re.compile('ERA5_([0-9]*)_.*.nc').search(filein).group(1))

    ##  Create output folders  -----------------------------------------------
    seasonl_dir = os.path.join(cnf.ERA5.path_regrid, f"Seasonal_{cnf.D1.LatStep}x{cnf.D1.LonStep}")
    monthly_dir = os.path.join(cnf.ERA5.path_regrid, f"Monthly_{cnf.D1.LatStep}x{cnf.D1.LonStep}")
    os.makedirs(seasonl_dir, exist_ok = True)
    os.makedirs(monthly_dir, exist_ok = True)

    ##  Limit data time range
    if not cnf.Range.start <= yyyy <= cnf.Range.until:
        print(f"Year {yyyy} is out of range")
        return "skip year"

    # filesin_C += 1  ## process counter
    # print(f"\nProcessing {filesin_C}/{filesin_N}: {filein}")
    print(f"\nProcessing {filein}")

    ##  load ERA5 main file on a xarray
    DT = xr.open_dataset(filein)

    ##  Get ERA5 spatial step
    lat_res = np.unique(np.diff(DT.latitude.values))[0]
    lon_res = np.unique(np.diff(DT.longitude.values))[0]

    ##  Export with seasonal aggregation  ------------------------------------
    if SEASONAL:
        stats_message  = "Statistical values (mean, median, SD and N) refer to the seasonal values of 3 calendar months, with the approximate centre the day of 'time' variable."
        stats_duration = 3
        data_type      = "Seasonal"

        ##  Compute by season of the year
        seasons = ['Q1_DJF', 'Q2_MAM', 'Q3_JJA', 'Q4_SON']
        for season_idx, season in enumerate(seasons):
            ##  Output process info  -----------------------------------------
            print(f"  {yyyy} {season}")

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
                previous_file = filein.replace(str(yyyy), str(yyyy - 1))
                if not os.path.exists(previous_file):
                    print(f"SKIP season! No file for the {yyyy - 1} found {previous_file}\n")
                    continue

                ## load ERA5 main file on a xarray
                DTpre = xr.open_dataset(previous_file)
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
            print(f"  {yyyy} M {m:02}")

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
## end of file process


##  Process raw ERA5 files  --------------------------------------------------
pool = Pool()
for _ in tqdm.tqdm(pool.imap_unordered(process_ERA5_file, filenames), total=len(filenames)):
    pass

#  SCRIPT END  ---------------------------------------------------------------
Ou.goodbye(cnf.LOGs.run, tic = tic, scriptname = __file__, version = TRACE)
