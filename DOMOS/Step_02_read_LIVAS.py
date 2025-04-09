#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Establishes the required EO-and-ERA5 DOMOS dataset in the same L3 1x1 grid resolution and monthly-mean.

@author: proestakis, thanasisn
"""

import os
import sys
import re

import netCDF4  as nc
import numpy    as np
import pandas   as pd
from   datetime import datetime, timedelta
import glob
# import math
import warnings
import numpy.ma as ma
import xarray as xr
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
VERSION = Ou.source_code_hash(__file__)


##  Set switches  ------------------------------------------------------------

##  Overwrite output files
FORCE = cnf.mode.Force
FORCE = True

##  Export by each month
MONTHLY  = False
MONTHLY  = cnf.D1.Monthly

##  Export by season
SEASONAL = False
SEASONAL = cnf.D1.Seasonal

##  Reduce work for testing
TEST = False
TEST = cnf.mode.Test

##  Allow only one case to run at the time  ----------------------------------
if SEASONAL == MONTHLY:
    print("Seasonal:", SEASONAL)
    print("Monthly: ", MONTHLY)
    sys.exit("Choose only SEASONAL or only MONTHLY")


##  Check destination folder exists  -----------------------------------------
if not os.path.isdir(cnf.OREO.path_output):
    sys.exit(f"\nFolder {cnf.OREO.path_output} don't exist !!\n")

fl_North = Oc.border_up(  cnf.D1.North, cnf.D1.MaxLatStep,   90)
fl_South = Oc.border_down(cnf.D1.South, cnf.D1.MaxLatStep,  -90)
fl_East  = Oc.border_up(  cnf.D1.East,  cnf.D1.MaxLonStep,  180)
fl_West  = Oc.border_down(cnf.D1.West,  cnf.D1.MaxLonStep, -180)
REGRID_lat_centers = np.arange(fl_North, fl_South + 1, -cnf.D1.LatStep) - cnf.D1.LatStep / 2
REGRID_lon_centers = np.arange(fl_West,  fl_East,       cnf.D1.LonStep) + cnf.D1.LonStep / 2


##  Temporal aggregation setup  ----------------------------------------------
if SEASONAL == True:
    print("Work on seasonal data")
    ERA_filenames = glob.glob(
        f"{cnf.ERA5.path_regrid}/Seasonal_{cnf.D1.LatStep}x{cnf.D1.LonStep}/ERA5_*_lat_{np.min(REGRID_lat_centers)}_{np.max(REGRID_lat_centers)}_lon_{np.min(REGRID_lon_centers)}_{np.max(REGRID_lon_centers)}.nc"
    )
elif MONTHLY == True:
    print("Work on monthly data")
    ERA_filenames = glob.glob(
        f"{cnf.ERA5.path_regrid}/Monthly_{cnf.D1.LatStep}x{cnf.D1.LonStep}/ERA5_*_lat_{np.min(REGRID_lat_centers)}_{np.max(REGRID_lat_centers)}_lon_{np.min(REGRID_lon_centers)}_{np.max(REGRID_lon_centers)}.nc"
    )

ERA_filenames.sort()

##  Select ERA5 variables by method
if cnf.ERA5.data == "mean":
    U = "u_mean"
    V = "v_mean"
elif cnf.ERA5.data == "median":
    U = "u_median"
    V = "v_median"

##  Create all LIVAS lats and longs
## TODO this is missing the last/first point
LIVAS_all_lats = np.arange( -89.5,  90)
LIVAS_all_lons = np.arange(-179.5, 180)

##  Directory path of output datasets  ----------------------------------
output_path    = os.path.join(cnf.OREO.path_output, os.path.basename(os.path.dirname(ERA_filenames[0])))
os.makedirs(output_path, exist_ok = True)



## !!!
if TEST: ERA_filenames = [ERA_filenames[0]]

for efid, ERA_file in enumerate(ERA_filenames):
    print(f"\nProcessing: {efid}/{len(ERA_filenames)} {ERA_file}")

    ##  Load ERA5 data  ------------------------------------------------------
    ERA           = xr.open_dataset(ERA_file)
    ERA_Latitude  = ERA.latitude
    ERA_Longitude = ERA.longitude

    ## TODO use to resolve logic
    ERA.season

    # extracting "yyyymm" suffix from ERA5 filename, for finding the satellite-based MM files.
    # ERA_year          = int(re.compile('ERA5_([0-9]*)_.*.nc').search(ERA_file).group(1))
    ERA_year          = pd.DatetimeIndex(ERA.time).year[0]
    ERA_month         = pd.DatetimeIndex(ERA.time).month[0]

    # ERA_season        = re.compile('ERA5_[0-9]*_Q[1-4]_([A-Z]*)_.*.nc').search(ERA_file).group(1)

    if   ERA.season == 'Q1_DJF':
        YoI = [ERA_year - 1, ERA_year, ERA_year]
        MoI = [12, 1, 2]
    elif ERA.season == 'Q2_MAM':
        YoI = [ERA_year, ERA_year, ERA_year]
        MoI = [3, 4, 5]
    elif ERA.season == 'Q3_JJA':
        YoI = [ERA_year, ERA_year, ERA_year]
        MoI = [6, 7, 8]
    elif ERA.season == 'Q4_SON':
        YoI = [ERA_year, ERA_year, ERA_year]
        MoI = [9, 10, 11]
    else:
        sys.exit("Resolve monthly logic")
        ### monthly here?


    fileout = os.path.join(output_path,
                           os.path.basename(ERA_file).replace("ERA5", "ERA5_LIVAS"))

    ## !!!
    if TEST:
        fileout  = fileout.replace(".nc", "_TEST.nc")
        ## do just a part of the ERA file
        ERA_Longitude = ERA_Longitude[0:2]
        ERA_Latitude  = ERA_Latitude [0:3]


    ## skip already existing files
    if (not FORCE) and (not Ou.output_needs_update(ERA_file, fileout)):
        continue

    Percentage_IGBP_1          = np.empty((len(ERA_Longitude), len(ERA_Latitude)))
    Percentage_IGBP_2          = np.empty((len(ERA_Longitude), len(ERA_Latitude)))
    Percentage_IGBP_3          = np.empty((len(ERA_Longitude), len(ERA_Latitude)))
    Percentage_IGBP_4          = np.empty((len(ERA_Longitude), len(ERA_Latitude)))
    Percentage_IGBP_5          = np.empty((len(ERA_Longitude), len(ERA_Latitude)))
    Percentage_IGBP_6          = np.empty((len(ERA_Longitude), len(ERA_Latitude)))
    Percentage_IGBP_7          = np.empty((len(ERA_Longitude), len(ERA_Latitude)))
    Percentage_IGBP_8          = np.empty((len(ERA_Longitude), len(ERA_Latitude)))
    Percentage_IGBP_9          = np.empty((len(ERA_Longitude), len(ERA_Latitude)))
    Percentage_IGBP_10         = np.empty((len(ERA_Longitude), len(ERA_Latitude)))
    Percentage_IGBP_11         = np.empty((len(ERA_Longitude), len(ERA_Latitude)))
    Percentage_IGBP_12         = np.empty((len(ERA_Longitude), len(ERA_Latitude)))
    Percentage_IGBP_13         = np.empty((len(ERA_Longitude), len(ERA_Latitude)))
    Percentage_IGBP_14         = np.empty((len(ERA_Longitude), len(ERA_Latitude)))
    Percentage_IGBP_15         = np.empty((len(ERA_Longitude), len(ERA_Latitude)))
    Percentage_IGBP_16         = np.empty((len(ERA_Longitude), len(ERA_Latitude)))
    Percentage_IGBP_17         = np.empty((len(ERA_Longitude), len(ERA_Latitude)))
    Percentage_IGBP_18         = np.empty((len(ERA_Longitude), len(ERA_Latitude)))
    Percentage_IGBP_1[:]       = np.nan
    Percentage_IGBP_2[:]       = np.nan
    Percentage_IGBP_3[:]       = np.nan
    Percentage_IGBP_4[:]       = np.nan
    Percentage_IGBP_5[:]       = np.nan
    Percentage_IGBP_6[:]       = np.nan
    Percentage_IGBP_7[:]       = np.nan
    Percentage_IGBP_8[:]       = np.nan
    Percentage_IGBP_9[:]       = np.nan
    Percentage_IGBP_10[:]      = np.nan
    Percentage_IGBP_11[:]      = np.nan
    Percentage_IGBP_12[:]      = np.nan
    Percentage_IGBP_13[:]      = np.nan
    Percentage_IGBP_14[:]      = np.nan
    Percentage_IGBP_15[:]      = np.nan
    Percentage_IGBP_16[:]      = np.nan
    Percentage_IGBP_17[:]      = np.nan
    Percentage_IGBP_18[:]      = np.nan

    Final_Number_of_Profiles    = np.empty((len(ERA_Longitude), len(ERA_Latitude)))
    Final_Number_of_L2Profiles  = np.empty((len(ERA_Longitude), len(ERA_Latitude)))
    Final_LIVAS_PD_DOD_532nm    = np.empty((len(ERA_Longitude), len(ERA_Latitude)))
    Final_LIVAS_PD_DOD_532nm_SD = np.empty((len(ERA_Longitude), len(ERA_Latitude)))
    Final_PD_a532nm             = np.empty((len(ERA_Longitude), len(ERA_Latitude), cnf.LIVAS.levels))
    Final_PD_a532nm_SD          = np.empty((len(ERA_Longitude), len(ERA_Latitude), cnf.LIVAS.levels))
    Final_PD_MC                 = np.empty((len(ERA_Longitude), len(ERA_Latitude), cnf.LIVAS.levels))
    Final_PD_MC_SD              = np.empty((len(ERA_Longitude), len(ERA_Latitude), cnf.LIVAS.levels))

    Final_Number_of_Profiles[:]    = np.nan
    Final_Number_of_L2Profiles[:]  = np.nan
    Final_LIVAS_PD_DOD_532nm[:]    = np.nan
    Final_LIVAS_PD_DOD_532nm_SD[:] = np.nan
    Final_PD_a532nm[:]             = np.nan
    Final_PD_a532nm_SD[:]          = np.nan
    Final_PD_MC[:]                 = np.nan
    Final_PD_MC_SD[:]              = np.nan

    Empty_Vertical_array    = np.empty(cnf.LIVAS.levels)
    Empty_Vertical_array[:] = np.nan


    ## TODO use product
    # coords = np.array(np.meshgrid(ERA_Latitude, ERA_Longitude)).T.reshape(-1,2)
    # for lat, lon in coords:
    #     print(lat, lon)


    for lon_id, lon in enumerate(ERA_Longitude.values):
        for lat_id, lat in enumerate(ERA_Latitude.values):

            ##  create list of all LIVAS coords to use for this cell
            Llats = LIVAS_all_lats[np.logical_and(
                LIVAS_all_lats > lat - (cnf.D1.LatStep / 2),
                LIVAS_all_lats < lat + (cnf.D1.LatStep / 2)
                )]

            Llons = LIVAS_all_lons[np.logical_and(
                LIVAS_all_lons > lon - (cnf.D1.LonStep / 2),
                LIVAS_all_lons < lon + (cnf.D1.LonStep / 2)
                )]

            ##  expand all combinations of LIVAS coordinates
            comb = np.array(np.meshgrid(Llats, Llons)).T.reshape(-1, 2)


            ## !!!
            if TEST: comb = comb[0:8]

            ##  Read a LIVAS file  --------------------------------------------
            file_counter = 0
            for LIVAS_lat, LIVAS_lon in comb:
                ec_c = lon_id + lat_id + 2
                ec_t = len(ERA_Longitude) * len(ERA_Latitude)
                li_c = file_counter + 1
                li_t = comb.shape[0]
                print(f"\n{efid}/{len(ERA_filenames)} {ec_c}/{ec_t} {li_c}/{li_t} [{lat} {lon}] <- [{LIVAS_lat} {LIVAS_lon}] ")

                ##  File to read  ---------------------------------------------
                LIVAS_file = os.path.join(
                    cnf.LIVAS.path_lookup,
                    f'LIVAS_CALIPSO_L2_Grid_lon_c_{str(LIVAS_lon)}_lat_c_{str(LIVAS_lat)}.nc')

                ##  Skip missing LIVAS files or issue an error  ---------------
                if not os.path.exists(LIVAS_file):
                    amsg = f"Missing file: {os.path.basename(LIVAS_file)}"
                    if cnf.mode.Test:
                        warnings.warn(amsg)
                        continue
                    else:
                        sys.exit(amsg)

                ##  Load LIVAS file  ------------------------------------------
                LIVAS         = xr.open_datatree(LIVAS_file)

                ##  Create a selection index of profile dates
                id_date_range = (pd.DatetimeIndex(LIVAS.Profile_Time_Parsed).year.isin(  YoI ) &
                                 pd.DatetimeIndex(LIVAS.Profile_Time_Parsed).month.isin( MoI ) )

                print(f"    LIVAS date range: {LIVAS.Profile_Time_Parsed[id_date_range].min().values} -- {LIVAS.Profile_Time_Parsed[id_date_range].max().values}")
                print(f"    Count: {(id_date_range).sum()}")

                ##  Skip files without data to use
                if (id_date_range).sum() == 0:
                    print("No usable data in LIVAS")
                    continue

                ##  Get data selection from LIVAS
                Altitude        = LIVAS.Altitude
                IGBP            = LIVAS.CALIPSO_Flags_and_Auxiliary.Auxiliary.IGBP_Surface_Type

                ##  Will exclude data with negative altitude  -----------------
                id_negative_altitude = Altitude < 0

                ##  Select date range and ignore negative altitude  -----------
                LIVAS_PD_b532nm = LIVAS.LIVAS.Cloud_Free.Pure_Dust_and_Fine_Coarse.Optical_Products.Pure_Dust_Backscatter_Coefficient_532.sel(
                    profile = id_date_range)
                LIVAS_PD_b532nm = xr.where(id_negative_altitude, np.nan, LIVAS_PD_b532nm)

                LIVAS_PD_a532nm = LIVAS.LIVAS.Cloud_Free.Pure_Dust_and_Fine_Coarse.Optical_Products.Pure_Dust_Extinction_Coefficient_532.sel(
                    profile = id_date_range)
                LIVAS_PD_a532nm = xr.where(id_negative_altitude, np.nan, LIVAS_PD_a532nm)

                LIVAS_PD_MC     = LIVAS.LIVAS.Cloud_Free.Pure_Dust_and_Fine_Coarse.Mass_Concentrations.Pure_Dust_Mass_Concentration.sel(
                    profile = id_date_range)
                LIVAS_PD_MC = xr.where(id_negative_altitude, np.nan, LIVAS_PD_MC)

                LIVAS_LR_Dust   = np.unique(LIVAS.LIVAS.Auxiliary.Lidar_Ratio_Assumptions.Lidar_Ratio_Dust.sel(
                    profile = id_date_range))

                ##  Set extinction and mass to zero when no backscater  -------
                idx = np.where(LIVAS_PD_b532nm == 0)
                LIVAS_PD_a532nm[idx] = 0
                LIVAS_PD_MC    [idx] = 0


                ### !!!  should we concerned with that?
                if ma.isMaskedArray(LIVAS_PD_b532nm) == True:
                    print("a maskded array")
                    sys.exit("DDD")
                    LIVAS_PD_a532nm[LIVAS_PD_b532nm.mask == True] = np.nan
                    LIVAS_PD_MC[    LIVAS_PD_b532nm.mask == True] = np.nan


                ##  Gather data for this cell  --------------------------------
                if file_counter == 0:
                    Total_LIVAS_PD_a532nm = LIVAS_PD_a532nm
                    Total_LIVAS_PD_MC     = LIVAS_PD_MC
                    Total_IGBP            = IGBP
                    Total_LIVAS_LR_Dust   = LIVAS_LR_Dust
                else:
                    Total_LIVAS_PD_MC     = xr.concat([LIVAS_PD_MC,         Total_LIVAS_PD_MC]    , "profile")
                    Total_LIVAS_PD_a532nm = xr.concat([LIVAS_PD_a532nm,     Total_LIVAS_PD_a532nm], "profile")
                    Total_IGBP            = xr.concat([Total_IGBP,          IGBP]                 , "profile")
                    Total_LIVAS_LR_Dust   = np.concat([Total_LIVAS_LR_Dust, LIVAS_LR_Dust])

                file_counter += 1
            ##  end iterate LIVAS files for this a cell

            # sys.exit("DDD")

            ##  Prepare selected data  ----------------------------------------



            Number_of_Profiles = np.shape(Total_LIVAS_PD_MC)[0]

            temp = np.copy(Total_LIVAS_PD_MC)
            idx  = ~np.isnan(Total_LIVAS_PD_MC)
            temp[idx] = 0

            idx  = np.isnan(Total_LIVAS_PD_MC)
            temp[idx] = 1
            temp = np.ravel([np.nansum(temp[i,:]) for i in range(np.shape(temp)[0])])

            L2_CF_profiles  = len(temp) - len(np.ravel(np.where(temp == cnf.LIVAS.levels)))


            ##  Ignore any dust hight in the atmospher  ----------------------
            Total_LIVAS_PD_a532nm[Altitude > cnf.LIVAS.height_limit_km] = np.nan
            Total_LIVAS_PD_MC    [Altitude > cnf.LIVAS.height_limit_km] = np.nan

            ##  Create mean vertical profile
            PD_a532nm    = Total_LIVAS_PD_a532nm.mean(dim = "profile", skipna = True)
            PD_MC        = Total_LIVAS_PD_MC    .mean(dim = "profile", skipna = True)

            PD_a532nm_SD = Total_LIVAS_PD_a532nm.std(dim = "profile", skipna = True, ddof = 1)
            PD_MC_SD     = Total_LIVAS_PD_MC    .std(dim = "profile", skipna = True, ddof = 1)

            # PD_a532nm    = np.nanmean(Total_LIVAS_PD_a532nm,axis = 0)
            # PD_MC        = np.nanmean(Total_LIVAS_PD_MC,    axis = 0)
            # PD_MC_SD     = np.nanstd(Total_LIVAS_PD_MC,     axis = 0, ddof = 1)
            # PD_a532nm_SD = np.nanstd(Total_LIVAS_PD_a532nm, axis = 0, ddof = 1)
            # PD_a532nm   [Altitude > cnf.LIVAS.height_limit_km] = np.nan
            # PD_MC       [Altitude > cnf.LIVAS.height_limit_km] = np.nan
            # PD_MC_SD    [Altitude > cnf.LIVAS.height_limit_km] = np.nan
            # PD_a532nm_SD[Altitude > cnf.LIVAS.height_limit_km] = np.nan

            # PD_a532nm.integrate(coord = "ltitude")

            ##  Create DOD with intergration  ---------------------------------


            ## FIXME integration
            arr = np.copy(PD_a532nm)
            arr[np.isnan(arr)] = 0
            DOD_532nm          = np.trapezoid(Altitude, arr)
            # np.unique(np.diff(Altitude))

            sel = ~np.isnan(PD_a532nm)
            PD_a532nm[sel].shape
            Altitude[sel].shape

            if DOD_532nm != np.trapezoid(Altitude[sel], PD_a532nm[sel]):
                sys.exit("ccc")

            arr = np.copy(PD_a532nm)
            arr[np.isnan(arr)] = 0

            arr.shape
            Altitude.shape




            arr = np.copy(PD_a532nm_SD)
            arr[np.isnan(arr)] = 0
            DOD_532nm_SD       = np.trapezoid(Altitude, arr)


            for count_alt in range(cnf.LIVAS.levels):
                Final_PD_MC[   lon_id, lat_id, count_alt]    = PD_MC   [count_alt]
                Final_PD_MC_SD[lon_id, lat_id, count_alt]    = PD_MC_SD[count_alt]
                Final_PD_a532nm                              = PD_a532nm[count_alt]
                Final_PD_a532nm_SD                           = PD_a532nm_SD[count_alt]

            Final_LIVAS_PD_DOD_532nm[lon_id,lat_id]     = DOD_532nm
            Final_LIVAS_PD_DOD_532nm_SD[lon_id,lat_id]  = DOD_532nm_SD

            Percentage_IGBP_1 [lon_id, lat_id] = (np.divide(float(np.count_nonzero(IGBP ==  1)),float(len(IGBP))))*100.0
            Percentage_IGBP_2 [lon_id, lat_id] = (np.divide(float(np.count_nonzero(IGBP ==  2)),float(len(IGBP))))*100.0
            Percentage_IGBP_3 [lon_id, lat_id] = (np.divide(float(np.count_nonzero(IGBP ==  3)),float(len(IGBP))))*100.0
            Percentage_IGBP_4 [lon_id, lat_id] = (np.divide(float(np.count_nonzero(IGBP ==  4)),float(len(IGBP))))*100.0
            Percentage_IGBP_5 [lon_id, lat_id] = (np.divide(float(np.count_nonzero(IGBP ==  5)),float(len(IGBP))))*100.0
            Percentage_IGBP_6 [lon_id, lat_id] = (np.divide(float(np.count_nonzero(IGBP ==  6)),float(len(IGBP))))*100.0
            Percentage_IGBP_7 [lon_id, lat_id] = (np.divide(float(np.count_nonzero(IGBP ==  7)),float(len(IGBP))))*100.0
            Percentage_IGBP_8 [lon_id, lat_id] = (np.divide(float(np.count_nonzero(IGBP ==  8)),float(len(IGBP))))*100.0
            Percentage_IGBP_9 [lon_id, lat_id] = (np.divide(float(np.count_nonzero(IGBP ==  9)),float(len(IGBP))))*100.0
            Percentage_IGBP_10[lon_id, lat_id] = (np.divide(float(np.count_nonzero(IGBP == 10)),float(len(IGBP))))*100.0
            Percentage_IGBP_11[lon_id, lat_id] = (np.divide(float(np.count_nonzero(IGBP == 11)),float(len(IGBP))))*100.0
            Percentage_IGBP_12[lon_id, lat_id] = (np.divide(float(np.count_nonzero(IGBP == 12)),float(len(IGBP))))*100.0
            Percentage_IGBP_13[lon_id, lat_id] = (np.divide(float(np.count_nonzero(IGBP == 13)),float(len(IGBP))))*100.0
            Percentage_IGBP_14[lon_id, lat_id] = (np.divide(float(np.count_nonzero(IGBP == 14)),float(len(IGBP))))*100.0
            Percentage_IGBP_15[lon_id, lat_id] = (np.divide(float(np.count_nonzero(IGBP == 15)),float(len(IGBP))))*100.0
            Percentage_IGBP_16[lon_id, lat_id] = (np.divide(float(np.count_nonzero(IGBP == 16)),float(len(IGBP))))*100.0
            Percentage_IGBP_17[lon_id, lat_id] = (np.divide(float(np.count_nonzero(IGBP == 17)),float(len(IGBP))))*100.0
            Percentage_IGBP_18[lon_id, lat_id] = (np.divide(float(np.count_nonzero(IGBP == 18)),float(len(IGBP))))*100.0

            Final_Number_of_Profiles[  lon_id, lat_id] = Number_of_Profiles
            Final_Number_of_L2Profiles[lon_id, lat_id] = L2_CF_profiles

            # sys.exit("wait")

    ## END iterate domain


    ##  Saving dataset as NetCDF ----------------------------------------------

    ## !!! why to m
    Altitude = Altitude*1000.0

    # creating nc. filename and initiallizing:
    ds           = nc.Dataset(fileout, 'w', format='NETCDF4')

    # create nc. dimensions:
    ds.createDimension('ERA_lev',     ERA.pressure_level.shape[0])
    ds.createDimension('CALIPSO_lev', len(Altitude))
    ds.createDimension('lat',         len(ERA_Latitude))
    ds.createDimension('lon',         len(ERA_Longitude))

    Geolocation_group     = ds.createGroup("Geolocation")
    ERA5_group            = ds.createGroup("ERA5")
    LIVAS_group           = ds.createGroup("LIVAS")
    Land_Ocean_Mask_group = ds.createGroup("Land_Ocean_Mask")

    # create nc. variables:
    lats_id                   = ds.createVariable('Geolocation/Latitude', 'f4', ('lat',))
    lons_id                   = ds.createVariable('Geolocation/Longitude','f4', ('lon',))

    ERA_U_id                  = ds.createVariable('ERA5/U',           np.float64, ('lon','lat','ERA_lev',), zlib=True)
    ERA_U_SD_id               = ds.createVariable('ERA5/U_SD',        np.float64, ('lon','lat','ERA_lev',), zlib=True)
    ERA_V_id                  = ds.createVariable('ERA5/V',           np.float64, ('lon','lat','ERA_lev',), zlib=True)
    ERA_V_SD_id               = ds.createVariable('ERA5/V_SD',        np.float64, ('lon','lat','ERA_lev',), zlib=True)
    ERA_height_low_id         = ds.createVariable('ERA5/height_low',  np.float64, ('lon','lat','ERA_lev',), zlib=True)
    ERA_Height_id             = ds.createVariable('ERA5/height',      np.float64, ('lon','lat','ERA_lev',), zlib=True)
    ERA_height_up_id          = ds.createVariable('ERA5/height_up',   np.float64, ('lon','lat','ERA_lev',), zlib=True)

    LIVAS_Altitude_id         = ds.createVariable('LIVAS/Altitude',                                        'f4',        ('CALIPSO_lev',),zlib=True)
    LIVAS_a532nm_PD_id        = ds.createVariable('LIVAS/Pure_Dust/LIVAS_PD_a532nm',           np.float64,  ('lon','lat','CALIPSO_lev',),zlib=True)
    LIVAS_a532nm_PD_SD_id     = ds.createVariable('LIVAS/Pure_Dust/LIVAS_PD_a532nm_STD',       np.float64,  ('lon','lat','CALIPSO_lev',),zlib=True)
    LIVAS_PD_MC_id            = ds.createVariable('LIVAS/Pure_Dust/LIVAS_PD_MC',               np.float64,  ('lon','lat','CALIPSO_lev',),zlib=True)
    LIVAS_PD_MC_SD_id         = ds.createVariable('LIVAS/Pure_Dust/LIVAS_PD_MC_STD',           np.float64,  ('lon','lat','CALIPSO_lev',),zlib=True)
    LIVAS_N_of_CF_Profiles_id = ds.createVariable('LIVAS/Flags/Number_of_L2_CF_Profiles',      np.float64,  ('lon','lat',)          ,zlib=True)
    LIVAS_N_of_Profiles_id    = ds.createVariable('LIVAS/Flags/Number_of_L2_Profiles',         np.float64,  ('lon','lat',)          ,zlib=True)
    LIVAS_DOD_532nm_mean      = ds.createVariable('LIVAS/Pure_Dust/DOD_532nm_mean',            np.float64,  ('lon','lat',)          ,zlib=True)
    LIVAS_DOD_532nm_SD        = ds.createVariable('LIVAS/Pure_Dust/DOD_532nm_STD',             np.float64,  ('lon','lat',)          ,zlib=True)

    Percentage_IGBP_1_id      = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_1',  np.float64, ('lon','lat',), zlib=True)
    Percentage_IGBP_2_id      = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_2',  np.float64, ('lon','lat',), zlib=True)
    Percentage_IGBP_3_id      = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_3',  np.float64, ('lon','lat',), zlib=True)
    Percentage_IGBP_4_id      = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_4',  np.float64, ('lon','lat',), zlib=True)
    Percentage_IGBP_5_id      = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_5',  np.float64, ('lon','lat',), zlib=True)
    Percentage_IGBP_6_id      = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_6',  np.float64, ('lon','lat',), zlib=True)
    Percentage_IGBP_7_id      = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_7',  np.float64, ('lon','lat',), zlib=True)
    Percentage_IGBP_8_id      = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_8',  np.float64, ('lon','lat',), zlib=True)
    Percentage_IGBP_9_id      = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_9',  np.float64, ('lon','lat',), zlib=True)
    Percentage_IGBP_10_id     = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_10', np.float64, ('lon','lat',), zlib=True)
    Percentage_IGBP_11_id     = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_11', np.float64, ('lon','lat',), zlib=True)
    Percentage_IGBP_12_id     = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_12', np.float64, ('lon','lat',), zlib=True)
    Percentage_IGBP_13_id     = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_13', np.float64, ('lon','lat',), zlib=True)
    Percentage_IGBP_14_id     = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_14', np.float64, ('lon','lat',), zlib=True)
    Percentage_IGBP_15_id     = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_15', np.float64, ('lon','lat',), zlib=True)
    Percentage_IGBP_16_id     = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_16', np.float64, ('lon','lat',), zlib=True)
    Percentage_IGBP_17_id     = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_17', np.float64, ('lon','lat',), zlib=True)
    Percentage_IGBP_18_id     = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_18', np.float64, ('lon','lat',), zlib=True)


    lats_id.units                   = 'degrees_north'
    lons_id.units                   = 'degrees_east'
    ERA_Height_id.units             = 'km'
    ERA_height_low_id.units         = 'km'
    ERA_height_up_id.units          = 'km'
    ERA_U_id.units                  = 'm s**-1'
    ERA_U_SD_id.units               = 'm s**-1'
    ERA_V_id.units                  = 'm s**-1'
    ERA_V_SD_id.units               = 'm s**-1'
    LIVAS_Altitude_id.units         = 'm'
    LIVAS_a532nm_PD_id.units        = 'km-1'
    LIVAS_a532nm_PD_SD_id.units     = 'km-1'
    LIVAS_PD_MC_id.units            = 'micrograms/m^3'
    LIVAS_PD_MC_SD_id.units         = 'micrograms/m^3'
    LIVAS_N_of_CF_Profiles_id.units = 'none'
    LIVAS_N_of_Profiles_id.units    = 'none'
    LIVAS_DOD_532nm_mean.units      = 'none'
    LIVAS_DOD_532nm_SD.units        = 'none'

    Percentage_IGBP_1_id.units      = 'none'
    Percentage_IGBP_2_id.units      = 'none'
    Percentage_IGBP_3_id.units      = 'none'
    Percentage_IGBP_4_id.units      = 'none'
    Percentage_IGBP_5_id.units      = 'none'
    Percentage_IGBP_6_id.units      = 'none'
    Percentage_IGBP_7_id.units      = 'none'
    Percentage_IGBP_8_id.units      = 'none'
    Percentage_IGBP_9_id.units      = 'none'
    Percentage_IGBP_10_id.units     = 'none'
    Percentage_IGBP_11_id.units     = 'none'
    Percentage_IGBP_12_id.units     = 'none'
    Percentage_IGBP_13_id.units     = 'none'
    Percentage_IGBP_14_id.units     = 'none'
    Percentage_IGBP_15_id.units     = 'none'
    Percentage_IGBP_16_id.units     = 'none'
    Percentage_IGBP_17_id.units     = 'none'
    Percentage_IGBP_18_id.units     = 'none'

    lats_id.long_name                   = 'Latitude'
    lons_id.long_name                   = 'Longitude'
    ERA_Height_id.long_name             = 'Height'
    ERA_height_low_id.long_name         = 'Height of lower cell boundary'
    ERA_height_up_id.long_name          = 'Height of upper cell boundary'
    ERA_U_id.long_name                  = 'U component of wind'
    ERA_U_SD_id.long_name               = 'U component of wind SD'
    ERA_V_id.long_name                  = 'V component of wind'
    ERA_V_SD_id.long_name               = 'V component of wind SD'

    LIVAS_Altitude_id.long_name         = 'Height'
    LIVAS_a532nm_PD_id.units            = 'Pure-Dust Extinction Coefficient 532nm'
    LIVAS_a532nm_PD_SD_id.units         = 'Pure-Dust Extinction Coefficient 532nm - SD'
    LIVAS_PD_MC_id.long_name            = 'Pure-Dust Mass Concentration'
    LIVAS_PD_MC_SD_id.long_name         = 'Pure-Dust Mass Concentration - SD'
    LIVAS_N_of_CF_Profiles_id.long_name = 'Number of CALIPSO L2 5km Cloud Free Profiles'
    LIVAS_N_of_Profiles_id.long_name    = 'Number of CALIPSO L2 5km Profiles'
    LIVAS_DOD_532nm_mean.long_name      = 'Dust Optical Depth 532nm - mean'
    LIVAS_DOD_532nm_SD.long_name        = 'Dust Optical Depth 532nm - SD'
    Percentage_IGBP_1_id.long_name      = 'Evergreen-Needleleaf-Forest'
    Percentage_IGBP_2_id.long_name      = 'Evergreen-Broadleaf-Forest'
    Percentage_IGBP_3_id.long_name      = 'Deciduous-Needleleaf-Forest'
    Percentage_IGBP_4_id.long_name      = 'Deciduous-Broadleaf-Forest'
    Percentage_IGBP_5_id.long_name      = 'Mixed-Forest'
    Percentage_IGBP_6_id.long_name      = 'Closed-Shrublands'
    Percentage_IGBP_7_id.long_name      = 'Open-Shrubland (Desert)'
    Percentage_IGBP_8_id.long_name      = 'Woody-Savanna'
    Percentage_IGBP_9_id.long_name      = 'Savanna'
    Percentage_IGBP_10_id.long_name     = 'Grassland'
    Percentage_IGBP_11_id.long_name     = 'Wetland'
    Percentage_IGBP_12_id.long_name     = 'Cropland'
    Percentage_IGBP_13_id.long_name     = 'Urban'
    Percentage_IGBP_14_id.long_name     = 'Crop-Mosaic'
    Percentage_IGBP_15_id.long_name     = 'Permanent-Snow'
    Percentage_IGBP_16_id.long_name     = 'Barren/Desert'
    Percentage_IGBP_17_id.long_name     = 'Water'
    Percentage_IGBP_18_id.long_name     = 'Tundra'

    lats_id.fill_value                   = np.nan
    lons_id.fill_value                   = np.nan

    ERA_Height_id.fill_value             = np.nan
    ERA_height_low_id.fill_value         = np.nan
    ERA_height_up_id.fill_value          = np.nan
    ERA_U_id.fill_value                  = np.nan
    ERA_U_SD_id.fill_value               = np.nan
    ERA_V_id.fill_value                  = np.nan
    ERA_V_SD_id.fill_value               = np.nan

    LIVAS_Altitude_id.fill_value         = np.nan
    LIVAS_a532nm_PD_id.fill_value        = np.nan
    LIVAS_a532nm_PD_SD_id.fill_value     = np.nan
    LIVAS_PD_MC_id.fill_value            = np.nan
    LIVAS_PD_MC_SD_id.fill_value         = np.nan
    LIVAS_N_of_CF_Profiles_id.fill_value = np.nan
    LIVAS_N_of_Profiles_id.fill_value    = np.nan
    LIVAS_DOD_532nm_mean.fill_value      = np.nan
    LIVAS_DOD_532nm_SD.fill_value        = np.nan

    Percentage_IGBP_1_id.fill_value      = np.nan
    Percentage_IGBP_2_id.fill_value      = np.nan
    Percentage_IGBP_3_id.fill_value      = np.nan
    Percentage_IGBP_4_id.fill_value      = np.nan
    Percentage_IGBP_5_id.fill_value      = np.nan
    Percentage_IGBP_6_id.fill_value      = np.nan
    Percentage_IGBP_7_id.fill_value      = np.nan
    Percentage_IGBP_8_id.fill_value      = np.nan
    Percentage_IGBP_9_id.fill_value      = np.nan
    Percentage_IGBP_10_id.fill_value     = np.nan
    Percentage_IGBP_11_id.fill_value     = np.nan
    Percentage_IGBP_12_id.fill_value     = np.nan
    Percentage_IGBP_13_id.fill_value     = np.nan
    Percentage_IGBP_14_id.fill_value     = np.nan
    Percentage_IGBP_15_id.fill_value     = np.nan
    Percentage_IGBP_16_id.fill_value     = np.nan
    Percentage_IGBP_17_id.fill_value     = np.nan
    Percentage_IGBP_18_id.fill_value     = np.nan


    lats_id[:]                   = ERA_Latitude
    lons_id[:]                   = ERA_Longitude

    ## TODO move selection in if TEST for efficiency

    ERA_Height_id[:]             = ERA.height.sel(latitude  = ERA_Latitude,
                                                  longitude = ERA_Longitude)
    ERA_U_id[:]                  = ERA[ U ].sel(latitude  = ERA_Latitude,
                                                longitude = ERA_Longitude)
    ERA_U_SD_id[:]               = ERA.u_SD.sel(latitude  = ERA_Latitude,
                                                longitude = ERA_Longitude)
    ERA_V_id[:]                  = ERA[ V ].sel(latitude  = ERA_Latitude,
                                                longitude = ERA_Longitude)
    ERA_V_SD_id[:]               = ERA.v_SD.sel(latitude  = ERA_Latitude,
                                                longitude = ERA_Longitude)
    ERA_height_low_id[:]          = ERA.height_low.sel(latitude  = ERA_Latitude,
                                                       longitude = ERA_Longitude)
    ERA_height_up_id[:]          = ERA.height_up.sel(latitude  = ERA_Latitude,
                                                      longitude = ERA_Longitude)


    LIVAS_Altitude_id[:]         = Altitude
    LIVAS_a532nm_PD_id[:]        = Final_PD_a532nm
    LIVAS_a532nm_PD_SD_id[:]     = Final_PD_a532nm_SD
    LIVAS_PD_MC_id[:]            = Final_PD_MC
    LIVAS_PD_MC_SD_id[:]         = Final_PD_MC_SD
    LIVAS_N_of_CF_Profiles_id[:] = Final_Number_of_L2Profiles
    LIVAS_N_of_Profiles_id[:]    = Final_Number_of_Profiles
    LIVAS_DOD_532nm_mean[:]      = Final_LIVAS_PD_DOD_532nm
    LIVAS_DOD_532nm_SD[:]        = Final_LIVAS_PD_DOD_532nm_SD

    Percentage_IGBP_1_id[:]      = Percentage_IGBP_1
    Percentage_IGBP_2_id[:]      = Percentage_IGBP_2
    Percentage_IGBP_3_id[:]      = Percentage_IGBP_3
    Percentage_IGBP_4_id[:]      = Percentage_IGBP_4
    Percentage_IGBP_5_id[:]      = Percentage_IGBP_5
    Percentage_IGBP_6_id[:]      = Percentage_IGBP_6
    Percentage_IGBP_7_id[:]      = Percentage_IGBP_7
    Percentage_IGBP_8_id[:]      = Percentage_IGBP_8
    Percentage_IGBP_9_id[:]      = Percentage_IGBP_9
    Percentage_IGBP_10_id[:]     = Percentage_IGBP_10
    Percentage_IGBP_11_id[:]     = Percentage_IGBP_11
    Percentage_IGBP_12_id[:]     = Percentage_IGBP_12
    Percentage_IGBP_13_id[:]     = Percentage_IGBP_13
    Percentage_IGBP_14_id[:]     = Percentage_IGBP_14
    Percentage_IGBP_15_id[:]     = Percentage_IGBP_15
    Percentage_IGBP_16_id[:]     = Percentage_IGBP_16
    Percentage_IGBP_17_id[:]     = Percentage_IGBP_17
    Percentage_IGBP_18_id[:]     = Percentage_IGBP_18

    ##  Set global attributes
    my_attrs = dict(title          = "Regridded ERA5 data with LIVAS lookup data",
                    type           = ERA.type,
                    season         = ERA.season,
                    details        = ERA.details,
                    data_date      = ERA.data_date,
                    contacts       = cnf.OREO.contact_emails,
                    user_host      = os.getlogin() + "@" + os.uname()[1],
                    source_version = VERSION)
    for name, value in my_attrs.items():
        setattr(ds, name, value)

    ds.close()

    print(f"Written: {fileout}")

    # sys.exit("Good")


#  SCRIPT END  ---------------------------------------------------------------
Ou.goodbye(cnf.LOGs.run, tic = tic, scriptname=__file__, version = VERSION)
