#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Intergrade the appropriate LIVAS data with the ERR5 regridded data.

@author: proestakis, thanasisn
"""

import os
import sys
import warnings
import glob
from datetime  import datetime, timezone
from itertools import product

import netCDF4  as nc
import numpy    as np
import pandas   as pd
import numpy.ma as ma
import xarray as xr

from random import shuffle
import matplotlib.pyplot as plt

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
##  Overwrite output files
FORCE = True
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

##  Resolve working domain  --------------------------------------------------
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
## TEST
# shuffle(ERA_filenames)

##  Select ERA5 variables to use by method
if cnf.ERA5.data == "mean":
    U = "u_mean"
    V = "v_mean"
elif cnf.ERA5.data == "median":
    U = "u_median"
    V = "v_median"

##  Create all LIVAS lats and longs  -----------------------------------------
## TODO this is missing the last/first point
LIVAS_all_lats = np.arange( -89.5,  90)
LIVAS_all_lons = np.arange(-179.5, 180)

##  Define percentile test  --------------------------------------------------
percentiles = [ 95, 96, 97, 98, 99, 100 ]

##  Make sure we include all LIVAS point within the domain  ------------------
if fl_East  >= LIVAS_all_lons.max() or \
   fl_West  <= LIVAS_all_lons.min() or \
   fl_North >= LIVAS_all_lats.max() or \
   fl_South <= LIVAS_all_lats.min() :
   sys.exit("Livas files out of bounds, adjust to actual ranges")

##  Directory path of output datasets  ---------------------------------------
output_path = os.path.join(cnf.OREO.path_output, os.path.basename(os.path.dirname(ERA_filenames[0])))
os.makedirs(output_path, exist_ok = True)

# ERA_filenames[1]
## !!! Few files for tests
if TEST: ERA_filenames = ERA_filenames[0:2]

##  Iterate over all regrided ERA files  -------------------------------------
for efid, ERA_file in enumerate(ERA_filenames):
    print(f"\nProcessing: {efid}/{len(ERA_filenames)} {ERA_file}")

    ##  Load ERA5 data  ------------------------------------------------------
    ERA           = xr.open_dataset(ERA_file)
    ERA_Latitude  = ERA.latitude
    ERA_Longitude = ERA.longitude

    ##  Get the ERA file year and month
    ERA_year      = pd.DatetimeIndex(ERA.time).year[0]
    ERA_month     = pd.DatetimeIndex(ERA.time).month[0]

    ##  Resolve the time span of data to use  --------------------------------
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
    elif ERA.type == "Monthly":
        YoI = [ERA_year ]
        MoI = [ERA_month]
    else:
        sys.exit("Can not resolve the time step of the input data!")

    ##  Output file name  ----------------------------------------------------
    fileout = os.path.join(output_path,
                           os.path.basename(ERA_file).replace("ERA5", "ERA5_LIVAS"))

    ## !!! TEST
    if TEST:
        ## doing just a part of the ERA file
        fileout  = fileout.replace(".nc", "_TEST.nc")
        ## create a subdomain for testing
        ERA_Longitude = ERA_Longitude.where(ERA_Longitude.longitude >= -2.5, drop = True)
        ERA_Latitude  = ERA_Latitude.where( ERA_Latitude.latitude   <=  25,  drop = True)
        ERA_Longitude = ERA_Longitude[0:1]
        ERA_Latitude  = ERA_Latitude [0:1]

    ##  Have to remove existing file in order to create multiple groups  -----
    if os.path.exists(fileout):
        os.remove(fileout)
        print(f"Remove existing file {fileout}")

    ##  Skip already existing files  ------
    if (not FORCE) and (not Ou.output_needs_update(ERA_file, fileout)):
        continue

    ##  Get array shape in one call
    coords_2d = (len(ERA_Longitude), len(ERA_Latitude))
    coords_3d = (len(ERA_Longitude), len(ERA_Latitude), cnf.LIVAS.levels)

    ##  Iterate different cutoff of backscatter values  ----------------------
    for i_per, aper in enumerate(percentiles):

        ##  Create structures to gather data for the domain  ------------------
        # Final_Number_of_Profiles    = np.full(coords_2d, np.nan)
        Total_Number_of_Profiles    = np.full(coords_2d, np.nan)
        Final_Number_of_L2Profiles  = np.full(coords_2d, np.nan)
        Final_LIVAS_PD_DOD_532nm    = np.full(coords_2d, np.nan)
        Final_LIVAS_PD_DOD_532nm_SD = np.full(coords_2d, np.nan)

        Final_PD_a532nm             = np.full(coords_3d, np.nan)
        Final_PD_a532nm_SD          = np.full(coords_3d, np.nan)
        Final_PD_b532nm             = np.full(coords_3d, np.nan)
        Final_PD_b532nm_SD          = np.full(coords_3d, np.nan)
        Final_PD_MC                 = np.full(coords_3d, np.nan)
        Final_PD_MC_Fine            = np.full(coords_3d, np.nan)
        Final_PD_MC_Coarse          = np.full(coords_3d, np.nan)
        Final_PD_MC_SD              = np.full(coords_3d, np.nan)
        Final_PD_MC_Fine_SD         = np.full(coords_3d, np.nan)
        Final_PD_MC_Coarse_SD       = np.full(coords_3d, np.nan)

        Final_Surface_Elevation_MEAN_mean = np.full(coords_2d, np.nan)
        Final_Surface_Elevation_MIN_mean  = np.full(coords_2d, np.nan)
        Final_Surface_Elevation_MEAN_SD   = np.full(coords_2d, np.nan)
        Final_Surface_Elevation_MIN_SD    = np.full(coords_2d, np.nan)

        Percentage_IGBP_1           = np.full(coords_2d, np.nan)
        Percentage_IGBP_2           = np.full(coords_2d, np.nan)
        Percentage_IGBP_3           = np.full(coords_2d, np.nan)
        Percentage_IGBP_4           = np.full(coords_2d, np.nan)
        Percentage_IGBP_5           = np.full(coords_2d, np.nan)
        Percentage_IGBP_6           = np.full(coords_2d, np.nan)
        Percentage_IGBP_7           = np.full(coords_2d, np.nan)
        Percentage_IGBP_8           = np.full(coords_2d, np.nan)
        Percentage_IGBP_9           = np.full(coords_2d, np.nan)
        Percentage_IGBP_10          = np.full(coords_2d, np.nan)
        Percentage_IGBP_11          = np.full(coords_2d, np.nan)
        Percentage_IGBP_12          = np.full(coords_2d, np.nan)
        Percentage_IGBP_13          = np.full(coords_2d, np.nan)
        Percentage_IGBP_14          = np.full(coords_2d, np.nan)
        Percentage_IGBP_15          = np.full(coords_2d, np.nan)
        Percentage_IGBP_16          = np.full(coords_2d, np.nan)
        Percentage_IGBP_17          = np.full(coords_2d, np.nan)
        Percentage_IGBP_18          = np.full(coords_2d, np.nan)

        ##  Iterate over all points in ERA5 files  ----------------------------
        ec_c = 0  ## coordinates counter
        for (i_lon, ERA_lon), (j_lat, ERA_lat) in product(enumerate(ERA_Longitude.values), enumerate(ERA_Latitude.values)):

            file_counter    = 0  # this is important for the logic
            profile_counter = 0  # just for statistics

            ##  List of all LIVAS coordinates to use for this cell  -----------
            Llats = LIVAS_all_lats[np.logical_and(
                LIVAS_all_lats > ERA_lat - (cnf.D1.LatStep / 2),
                LIVAS_all_lats < ERA_lat + (cnf.D1.LatStep / 2)
                )]

            Llons = LIVAS_all_lons[np.logical_and(
                LIVAS_all_lons > ERA_lon - (cnf.D1.LonStep / 2),
                LIVAS_all_lons < ERA_lon + (cnf.D1.LonStep / 2)
                )]

            ##  Get percentile limits to use at this location  ----------------
            BLIM = Ou.backscatter_percentiles_lookup(
                ylat      = ERA_lat,
                xlon      = ERA_lon,
                datafile  = cnf.LIVAS.percentiles_fl,
                ylat_step = cnf.D1.LatStep,
                xlon_step = cnf.D1.LonStep)

            ##  expand all combinations of LIVAS coordinates to loop
            comb = np.array(np.meshgrid(Llats, Llons)).T.reshape(-1, 2)

            ## !!! TEST do few of LIVAS files
            # if TEST: comb = comb[0:4]

            ##  Read a LIVAS file  --------------------------------------------
            ec_c += 1                                                      # domain iteration
            ec_t = len(ERA_Longitude) * len(ERA_Latitude)                  # total domain points
            li_t = comb.shape[0]                                           # total cell points
            total_it = ec_t * li_t * len(ERA_filenames) * len(percentiles) # total iterations
            for LIVAS_lat, LIVAS_lon in comb:
                li_c     = file_counter + 1                                # count cell iteration
                curre_it = (efid * ec_t * li_t ) + ((ec_c) * li_t) + li_c  # current iteration
                curre_it = curre_it * (i_per + 1)                          # expand percentile iterations
                complete = 100 * curre_it / total_it                       # percentage of completeness
                duration = datetime.now() - tic                            # current run duration
                if complete == 0:
                    total = duration * 100                                 # dummy duration value
                else:
                    total = duration * 100 / complete                      # estimate of total duration
                eta = total - duration
                eda = (datetime.now() + eta).replace(microsecond = 0)
                print(f"\n> {efid}/{len(ERA_filenames)} {i_per}/{len(percentiles)} {ec_c:>3}/{ec_t} {li_c}/{li_t} [{ERA_lat} {ERA_lon}] <- [{LIVAS_lat} {LIVAS_lon}]")
                print(f"> {ERA_year} {ERA.season} {curre_it}/{total_it} {complete:6.2f}% P: {format(duration)} R: {format(eta)} F: {eda}")

                ##  File to read  ---------------------------------------------
                LIVAS_file = os.path.join(
                    cnf.LIVAS.path_lookup,
                    f"LIVAS_CALIPSO_L2_Grid_lon_c_{str(LIVAS_lon)}_lat_c_{str(LIVAS_lat)}.nc")
                # print(f"> {LIVAS_file}\n")

                ##  Skip missing LIVAS files or issue an error  ---------------
                if not os.path.exists(LIVAS_file):
                    amsg = f"Missing file: {os.path.basename(LIVAS_file)}"
                    if cnf.mode.Test:
                        warnings.warn(amsg)
                        continue
                    else:
                        sys.exit(amsg)

                ##  Load LIVAS file xarray  -----------------------------------
                # LIVAS         = xr.open_datatree(LIVAS_file, engine="netcdf4") ## seg faults?
                LIVAS         = xr.open_datatree(LIVAS_file, engine="h5netcdf")

                ##  Create a selection index of profile dates
                id_date_range = (pd.DatetimeIndex(LIVAS.Profile_Time_Parsed).year.isin(  YoI ) &
                                 pd.DatetimeIndex(LIVAS.Profile_Time_Parsed).month.isin( MoI ) )

                ##  Will keep only profile in these dates
                id_profiles   = LIVAS.profile[id_date_range]

                ##  Skip files without usable data  ---------------------------
                if (id_date_range).sum() == 0:
                    print("No usable data in LIVAS")
                    continue

                ##  Date range of matching LIVAS data
                # print(f"    LIVAS date range: {LIVAS.Profile_Time_Parsed[id_date_range].min().values} -- {LIVAS.Profile_Time_Parsed[id_date_range].max().values}")
                ##  Profiles in each file for this cell
                # print(f"    Count: {id_date_range.sum()}")

                ##  Get data selection from LIVAS  ----------------------------
                ## BEWARE some variables are in km in LIVAS files
                Altitude = LIVAS.Altitude * 1000  ## Convert to meters!!

                ##  Will exclude data with negative altitude
                id_negative_altitude = Altitude < 0

                IGBP = LIVAS.CALIPSO_Flags_and_Auxiliary.Auxiliary.IGBP_Surface_Type

                Surface_elevation = LIVAS.CALIPSO_Flags_and_Auxiliary.Auxiliary.Surface_Elevation_Statistics.sel(
                    profile = id_date_range)

                ## 0: MIN, 1: MAX, 2: MEAN, 3: SD
                Surface_elevation_MIN  = Surface_elevation.sel(statistic = 0)
                Surface_elevation_MEAN = Surface_elevation.sel(statistic = 2)

                ##  Will remove bellow minimum surface elevation of every profile
                Surface_elevation_limit_id = LIVAS.Altitude <= Surface_elevation_MIN

                ##  Select date range and NaN for negative altitude
                LIVAS_PD_b532nm = LIVAS.LIVAS.Cloud_Free.Pure_Dust_and_Fine_Coarse.Optical_Products.Pure_Dust_Backscatter_Coefficient_532.sel(
                    profile = id_profiles)
                LIVAS_PD_b532nm = xr.where(id_negative_altitude, np.nan, LIVAS_PD_b532nm)
                LIVAS_PD_b532nm = LIVAS_PD_b532nm.where(~Surface_elevation_limit_id, np.nan)

                LIVAS_PD_a532nm = LIVAS.LIVAS.Cloud_Free.Pure_Dust_and_Fine_Coarse.Optical_Products.Pure_Dust_Extinction_Coefficient_532.sel(
                    profile = id_profiles)
                LIVAS_PD_a532nm = xr.where(id_negative_altitude, np.nan, LIVAS_PD_a532nm)
                LIVAS_PD_a532nm = LIVAS_PD_a532nm.where( ~Surface_elevation_limit_id, np.nan )

                LIVAS_PD_MC = LIVAS.LIVAS.Cloud_Free.Pure_Dust_and_Fine_Coarse.Mass_Concentrations.Pure_Dust_Mass_Concentration.sel(
                    profile = id_date_range)
                LIVAS_PD_MC = xr.where(id_negative_altitude, np.nan, LIVAS_PD_MC)
                LIVAS_PD_MC = LIVAS_PD_MC.where( ~Surface_elevation_limit_id, np.nan )

                LIVAS_PD_MC_Fine = LIVAS.LIVAS.Cloud_Free.Pure_Dust_and_Fine_Coarse.Mass_Concentrations.Pure_Dust_Fine_Mass_Concentration.sel(
                    profile = id_date_range)
                LIVAS_PD_MC_Fine = xr.where(id_negative_altitude, np.nan, LIVAS_PD_MC_Fine)
                LIVAS_PD_MC_Fine = LIVAS_PD_MC_Fine.where( ~Surface_elevation_limit_id, np.nan )

                LIVAS_PD_MC_Coarse = LIVAS.LIVAS.Cloud_Free.Pure_Dust_and_Fine_Coarse.Mass_Concentrations.Pure_Dust_Coarse_Mass_Concentration.sel(
                    profile = id_date_range)
                LIVAS_PD_MC_Coarse = xr.where(id_negative_altitude, np.nan, LIVAS_PD_MC_Coarse)
                LIVAS_PD_MC_Coarse = LIVAS_PD_MC_Coarse.where( ~Surface_elevation_limit_id, np.nan )

                LIVAS_LR_Dust   = np.unique(LIVAS.LIVAS.Auxiliary.Lidar_Ratio_Assumptions.Lidar_Ratio_Dust.sel(
                    profile = id_date_range))

                ##  Limit extreme values of backscatter  ----------------------
                if aper != 100:
                    ## backscatter > percentile => drop backscatter and all matching data
                    blim = BLIM[f"Percentile_{aper}"]
                    ## make sure we got a scalar
                    if not blim.shape == (1, 1):
                        sys.exit("Dimension error")
                    blim = blim.item()
                    print(f"> Apply backscatter limit of {aper}% or {blim}")

                    ##  Create backscatter limit  -----------------------------
                    id_extreme_backscatter = LIVAS_PD_b532nm > blim

                    ##  Apply backscatter limit  ------------------------------
                    LIVAS_PD_b532nm    = LIVAS_PD_b532nm.where   ( ~id_extreme_backscatter, np.nan )
                    LIVAS_PD_a532nm    = LIVAS_PD_a532nm.where   ( ~id_extreme_backscatter, np.nan )
                    LIVAS_PD_MC        = LIVAS_PD_MC.where       ( ~id_extreme_backscatter, np.nan )
                    LIVAS_PD_MC_Fine   = LIVAS_PD_MC_Fine.where  ( ~id_extreme_backscatter, np.nan )
                    LIVAS_PD_MC_Coarse = LIVAS_PD_MC_Coarse.where( ~id_extreme_backscatter, np.nan )


                ## TODO histogram of backscatter for a region
                ## maps of percentiles

                ##  Set extinction and mass to zero when no backscatter  ----------
                LIVAS_PD_a532nm    = LIVAS_PD_a532nm.   where(LIVAS_PD_b532nm != 0.0, 0)
                LIVAS_PD_MC        = LIVAS_PD_MC.       where(LIVAS_PD_b532nm != 0.0, 0)
                LIVAS_PD_MC_Fine   = LIVAS_PD_MC_Fine.  where(LIVAS_PD_b532nm != 0.0, 0)
                LIVAS_PD_MC_Coarse = LIVAS_PD_MC_Coarse.where(LIVAS_PD_b532nm != 0.0, 0)

                ##  Count total profiles in the cell  -----------------------------
                profile_counter += id_date_range.sum()

                ##  Resolve masked array cases  -----------------------------------
                if ma.isMaskedArray(LIVAS_PD_b532nm) == True:
                    print("A masked array found and applied")
                    LIVAS_PD_a532nm[LIVAS_PD_b532nm.mask == True] = np.nan
                    LIVAS_PD_b532nm[LIVAS_PD_b532nm.mask == True] = np.nan
                    LIVAS_PD_MC[    LIVAS_PD_b532nm.mask == True] = np.nan

                ##  Gather data for this cell  ------------------------------------
                if file_counter == 0:
                    Total_LIVAS_PD_a532nm        = LIVAS_PD_a532nm
                    Total_LIVAS_PD_b532nm        = LIVAS_PD_b532nm
                    Total_LIVAS_PD_MC            = LIVAS_PD_MC
                    Total_LIVAS_PD_MC_Fine       = LIVAS_PD_MC_Fine
                    Total_LIVAS_PD_MC_Coarse     = LIVAS_PD_MC_Coarse
                    Total_IGBP                   = IGBP
                    Total_LIVAS_LR_Dust          = LIVAS_LR_Dust
                    Total_Surface_elevation_MIN  = Surface_elevation_MIN
                    Total_Surface_elevation_MEAN = Surface_elevation_MEAN
                else:
                    Total_LIVAS_PD_a532nm        = xr.concat([Total_LIVAS_PD_a532nm,  LIVAS_PD_a532nm,   ], "profile")
                    Total_LIVAS_PD_b532nm        = xr.concat([Total_LIVAS_PD_b532nm,  LIVAS_PD_b532nm,   ], "profile")
                    Total_LIVAS_PD_MC            = xr.concat([Total_LIVAS_PD_MC,      LIVAS_PD_MC,       ], "profile")
                    Total_LIVAS_PD_MC_Fine       = xr.concat([Total_LIVAS_PD_MC_Fine, LIVAS_PD_MC_Fine,  ], "profile")
                    Total_LIVAS_PD_MC_Coarse     = xr.concat([Total_LIVAS_PD_MC,      LIVAS_PD_MC_Coarse,], "profile")
                    Total_IGBP                   = xr.concat([Total_IGBP,             IGBP               ], "profile")
                    Total_LIVAS_LR_Dust          = np.concat([Total_LIVAS_LR_Dust,    LIVAS_LR_Dust               ])
                    Total_Surface_elevation_MIN  = np.concat([Total_Surface_elevation_MIN,  Surface_elevation_MIN ])
                    Total_Surface_elevation_MEAN = np.concat([Total_Surface_elevation_MEAN, Surface_elevation_MEAN])

                file_counter += 1
            ##  end iterate LIVAS files for this a cell

            ##  Prepare selected data for a cell  =================================

            ##  Count profiles  ---------------------------------------------------
            # Final_Number_of_Profiles[i_lon, j_lat] = np.shape(Total_LIVAS_PD_MC)[1]  ## TODO remove var
            Total_Number_of_Profiles[i_lon, j_lat] = profile_counter

            ##  Count Cloud free profiles -----------------------------------------
            temp = np.isnan(Total_LIVAS_PD_MC).sum(dim = "altitude")
            L2_CF_profiles  = len(temp) - len(temp[np.where(temp == cnf.LIVAS.levels)])
            Final_Number_of_L2Profiles[i_lon, j_lat] = L2_CF_profiles

            ##  Ignore any dust high in the atmosphere  ---------------------------
            Total_LIVAS_PD_a532nm   [Altitude > cnf.LIVAS.height_limit_m] = np.nan
            Total_LIVAS_PD_b532nm   [Altitude > cnf.LIVAS.height_limit_m] = np.nan
            Total_LIVAS_PD_MC       [Altitude > cnf.LIVAS.height_limit_m] = np.nan
            Total_LIVAS_PD_MC_Fine  [Altitude > cnf.LIVAS.height_limit_m] = np.nan
            Total_LIVAS_PD_MC_Coarse[Altitude > cnf.LIVAS.height_limit_m] = np.nan

            ##  Create mean profile on vertical  ----------------------------------
            PD_a532nm    = Total_LIVAS_PD_a532nm   .mean(dim = "profile", skipna = True)
            PD_b532nm    = Total_LIVAS_PD_b532nm   .mean(dim = "profile", skipna = True)
            PD_MC        = Total_LIVAS_PD_MC       .mean(dim = "profile", skipna = True)
            PD_MC_Fine   = Total_LIVAS_PD_MC_Fine  .mean(dim = "profile", skipna = True)
            PD_MC_Coarse = Total_LIVAS_PD_MC_Coarse.mean(dim = "profile", skipna = True)

            PD_a532nm_SD    = Total_LIVAS_PD_a532nm   .std(dim = "profile", skipna = True, ddof = 1)
            PD_b532nm_SD    = Total_LIVAS_PD_b532nm   .std(dim = "profile", skipna = True, ddof = 1)
            PD_MC_SD        = Total_LIVAS_PD_MC       .std(dim = "profile", skipna = True, ddof = 1)
            PD_MC_Fine_SD   = Total_LIVAS_PD_MC_Fine  .std(dim = "profile", skipna = True, ddof = 1)
            PD_MC_Coarse_SD = Total_LIVAS_PD_MC_Coarse.std(dim = "profile", skipna = True, ddof = 1)

            # # Plotting
            # plt.figure(figsize=(6, 8))
            # plt.xlim(0, 0.25)
            # # plt.ylim(0, 10)
            # upper = PD_a532nm + PD_a532nm_SD
            # lower = PD_a532nm - PD_a532nm_SD
            # plt.fill_betweenx(Altitude, lower, upper, color='blue', alpha=0.3)
            # plt.plot(PD_a532nm, Altitude)
            # plt.xlabel("a 532 nm [km-1]")
            # plt.ylabel("Altitude (km)")
            # plt.grid(True)
            # plt.show()

            ##  Create DOD with vertical integration  -----------------------------
            arr = np.copy(PD_a532nm)
            arr[np.isnan(arr)] = 0
            ## Altitude in meters Extinction in kilometres
            DOD_532nm          = np.trapezoid(Altitude / 1000, arr)

            arr = np.copy(PD_a532nm_SD)
            arr[np.isnan(arr)] = 0
            ## Altitude in meters Extinction in kilometres
            DOD_532nm_SD       = np.trapezoid(Altitude / 1000, arr)

            Final_LIVAS_PD_DOD_532nm[   i_lon, j_lat] = DOD_532nm
            Final_LIVAS_PD_DOD_532nm_SD[i_lon, j_lat] = DOD_532nm_SD

            ##  Store vertical distributions  -------------------------------------
            Final_PD_MC[          i_lon, j_lat, :] = PD_MC
            Final_PD_MC_SD[       i_lon, j_lat, :] = PD_MC_SD
            Final_PD_MC_Fine[     i_lon, j_lat, :] = PD_MC_Fine
            Final_PD_MC_Fine_SD[  i_lon, j_lat, :] = PD_MC_Fine_SD
            Final_PD_MC_Coarse[   i_lon, j_lat, :] = PD_MC_Coarse
            Final_PD_MC_Coarse_SD[i_lon, j_lat, :] = PD_MC_Coarse_SD
            Final_PD_a532nm[      i_lon, j_lat, :] = PD_a532nm    / 1000.  # 1/km => 1/m
            Final_PD_a532nm_SD[   i_lon, j_lat, :] = PD_a532nm_SD / 1000.  # 1/km => 1/m
            Final_PD_b532nm[      i_lon, j_lat, :] = PD_b532nm    / 1000.  # 1/(km sr) => 1/(m sr)
            Final_PD_b532nm_SD[   i_lon, j_lat, :] = PD_b532nm_SD / 1000.  # 1/(km sr) => 1/(m sr)

            ##  Store surface elevation stats  ------------------------------------
            Final_Surface_Elevation_MEAN_mean[i_lon, j_lat] = Total_Surface_elevation_MEAN.mean() * 1000.  # km => m
            Final_Surface_Elevation_MIN_mean [i_lon, j_lat] = Total_Surface_elevation_MIN.mean()  * 1000.  # km => m
            Final_Surface_Elevation_MEAN_SD  [i_lon, j_lat] = Total_Surface_elevation_MEAN.std()  * 1000.  # km => m
            Final_Surface_Elevation_MIN_SD   [i_lon, j_lat] = Total_Surface_elevation_MIN.std()   * 1000.  # km => m

            ##  Compute coverage percentage  --------------------------------------
            len_IGBP = float(len(IGBP))
            Percentage_IGBP_1 [i_lon, j_lat] = np.divide(np.count_nonzero(IGBP ==  1), len_IGBP) * 100
            Percentage_IGBP_2 [i_lon, j_lat] = np.divide(np.count_nonzero(IGBP ==  2), len_IGBP) * 100
            Percentage_IGBP_3 [i_lon, j_lat] = np.divide(np.count_nonzero(IGBP ==  3), len_IGBP) * 100
            Percentage_IGBP_4 [i_lon, j_lat] = np.divide(np.count_nonzero(IGBP ==  4), len_IGBP) * 100
            Percentage_IGBP_5 [i_lon, j_lat] = np.divide(np.count_nonzero(IGBP ==  5), len_IGBP) * 100
            Percentage_IGBP_6 [i_lon, j_lat] = np.divide(np.count_nonzero(IGBP ==  6), len_IGBP) * 100
            Percentage_IGBP_7 [i_lon, j_lat] = np.divide(np.count_nonzero(IGBP ==  7), len_IGBP) * 100
            Percentage_IGBP_8 [i_lon, j_lat] = np.divide(np.count_nonzero(IGBP ==  8), len_IGBP) * 100
            Percentage_IGBP_9 [i_lon, j_lat] = np.divide(np.count_nonzero(IGBP ==  9), len_IGBP) * 100
            Percentage_IGBP_10[i_lon, j_lat] = np.divide(np.count_nonzero(IGBP == 10), len_IGBP) * 100
            Percentage_IGBP_11[i_lon, j_lat] = np.divide(np.count_nonzero(IGBP == 11), len_IGBP) * 100
            Percentage_IGBP_12[i_lon, j_lat] = np.divide(np.count_nonzero(IGBP == 12), len_IGBP) * 100
            Percentage_IGBP_13[i_lon, j_lat] = np.divide(np.count_nonzero(IGBP == 13), len_IGBP) * 100
            Percentage_IGBP_14[i_lon, j_lat] = np.divide(np.count_nonzero(IGBP == 14), len_IGBP) * 100
            Percentage_IGBP_15[i_lon, j_lat] = np.divide(np.count_nonzero(IGBP == 15), len_IGBP) * 100
            Percentage_IGBP_16[i_lon, j_lat] = np.divide(np.count_nonzero(IGBP == 16), len_IGBP) * 100
            Percentage_IGBP_17[i_lon, j_lat] = np.divide(np.count_nonzero(IGBP == 17), len_IGBP) * 100
            Percentage_IGBP_18[i_lon, j_lat] = np.divide(np.count_nonzero(IGBP == 18), len_IGBP) * 100

            ##  Verify coverage percentage consistency  ----------------------------
            ## This was tested, results as expected
            # adder = 0.0
            # for i in range(1, 19):
            #     adder += globals()[f"Percentage_IGBP_{i}"][i_lon, j_lat]
            # print(f"Area coverage {adder}")
            # if abs(adder - 100) > 1:
            #     sys.exit("Inconsistent coverage")
        ## END iteration on the whole ERA5 domain

        ##  Saving dataset as NetCDF ----------------------------------------------
        ##  Have to remove existing file for this to work
        ds = nc.Dataset(fileout, 'a', format='NETCDF4')

        ##  Store common data only one time  -----------------
        if not "ERA5" in ds.groups:

            ##  Create global dimensions
            ds.createDimension('ERA_lev',     ERA.pressure_level.shape[0])
            ds.createDimension('CALIPSO_lev', len(Altitude))
            ds.createDimension('lat',         len(ERA_Latitude))
            ds.createDimension('lon',         len(ERA_Longitude))
            ds.createDimension('time',        1)
            ds.createDimension('time_span',   1)

            ##  Create global groups
            ds.createGroup("ERA5")
            ds.createGroup("CALIPSO")
            ds.createGroup("Land_Ocean_Mask")

            ##  Create global variables
            lats_id                   = ds.createVariable('Latitude',  'f4', ('lat',))
            lons_id                   = ds.createVariable('Longitude', 'f4', ('lon',))
            time                      = ds.createVariable('time',       np.float64, ('time',))
            time_span                 = ds.createVariable('time_span',  np.int32,   ('time_span',))

            ERA_U_id                  = ds.createVariable('ERA5/U',          np.float64, ('lon','lat','ERA_lev',), zlib=True)
            ERA_U_SD_id               = ds.createVariable('ERA5/U_SD',       np.float64, ('lon','lat','ERA_lev',), zlib=True)
            ERA_V_id                  = ds.createVariable('ERA5/V',          np.float64, ('lon','lat','ERA_lev',), zlib=True)
            ERA_V_SD_id               = ds.createVariable('ERA5/V_SD',       np.float64, ('lon','lat','ERA_lev',), zlib=True)
            ERA_height_low_id         = ds.createVariable('ERA5/height_low', np.float64, ('lon','lat','ERA_lev',), zlib=True)
            ERA_Height_id             = ds.createVariable('ERA5/height',     np.float64, ('lon','lat','ERA_lev',), zlib=True)
            ERA_height_up_id          = ds.createVariable('ERA5/height_up',  np.float64, ('lon','lat','ERA_lev',), zlib=True)

            CALIPSO_Surface_Elevation_MEAN_mean = ds.createVariable('CALIPSO/Surface_Elevation_MEAN_mean', np.float64, ('lon','lat',), zlib=True)
            CALIPSO_Surface_Elevation_MIN_mean  = ds.createVariable('CALIPSO/Surface_Elevation_MIN_mean',  np.float64, ('lon','lat',), zlib=True)
            CALIPSO_Surface_Elevation_MEAN_SD   = ds.createVariable('CALIPSO/Surface_Elevation_MEAN_SD',   np.float64, ('lon','lat',), zlib=True)
            CALIPSO_Surface_Elevation_MIN_SD    = ds.createVariable('CALIPSO/Surface_Elevation_MIN_SD',    np.float64, ('lon','lat',), zlib=True)

            Percentage_IGBP_1_id  = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_1',  np.float64, ('lon','lat',), zlib=True)
            Percentage_IGBP_2_id  = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_2',  np.float64, ('lon','lat',), zlib=True)
            Percentage_IGBP_3_id  = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_3',  np.float64, ('lon','lat',), zlib=True)
            Percentage_IGBP_4_id  = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_4',  np.float64, ('lon','lat',), zlib=True)
            Percentage_IGBP_5_id  = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_5',  np.float64, ('lon','lat',), zlib=True)
            Percentage_IGBP_6_id  = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_6',  np.float64, ('lon','lat',), zlib=True)
            Percentage_IGBP_7_id  = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_7',  np.float64, ('lon','lat',), zlib=True)
            Percentage_IGBP_8_id  = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_8',  np.float64, ('lon','lat',), zlib=True)
            Percentage_IGBP_9_id  = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_9',  np.float64, ('lon','lat',), zlib=True)
            Percentage_IGBP_10_id = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_10', np.float64, ('lon','lat',), zlib=True)
            Percentage_IGBP_11_id = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_11', np.float64, ('lon','lat',), zlib=True)
            Percentage_IGBP_12_id = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_12', np.float64, ('lon','lat',), zlib=True)
            Percentage_IGBP_13_id = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_13', np.float64, ('lon','lat',), zlib=True)
            Percentage_IGBP_14_id = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_14', np.float64, ('lon','lat',), zlib=True)
            Percentage_IGBP_15_id = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_15', np.float64, ('lon','lat',), zlib=True)
            Percentage_IGBP_16_id = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_16', np.float64, ('lon','lat',), zlib=True)
            Percentage_IGBP_17_id = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_17', np.float64, ('lon','lat',), zlib=True)
            Percentage_IGBP_18_id = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_18', np.float64, ('lon','lat',), zlib=True)

            ##  Set attributes
            lats_id.units                   = 'degrees_north'
            lons_id.units                   = 'degrees_east'
            time.calendar                   = 'proleptic_gregorian'
            time.units                      = 'seconds since 1970-01-01 00:00:00'  ## same as raw ERA5
            time_span.units                 = 'month'

            ERA_Height_id.units             = 'm'
            ERA_height_low_id.units         = 'm'
            ERA_height_up_id.units          = 'm'
            ERA_U_id.units                  = 'm/s'
            ERA_U_SD_id.units               = 'm/s'
            ERA_V_id.units                  = 'm/s'
            ERA_V_SD_id.units               = 'm/s'

            CALIPSO_Surface_Elevation_MEAN_mean.units = 'm'  # changed altitude to meters
            CALIPSO_Surface_Elevation_MIN_mean.units  = 'm'  # changed altitude to meters
            CALIPSO_Surface_Elevation_MEAN_SD.units   = 'm'  # changed altitude to meters
            CALIPSO_Surface_Elevation_MIN_SD.units    = 'm'  # changed altitude to meters

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

            lats_id.long_name               = 'Latitude'
            lons_id.long_name               = 'Longitude'
            ERA_Height_id.long_name         = 'Height'
            ERA_height_low_id.long_name     = 'Height of lower cell boundary'
            ERA_height_up_id.long_name      = 'Height of upper cell boundary'
            ERA_U_id.long_name              = 'U component of wind'
            ERA_U_SD_id.long_name           = 'U component of wind SD'
            ERA_V_id.long_name              = 'V component of wind'
            ERA_V_SD_id.long_name           = 'V component of wind SD'
            time.long_name                  = 'Time'
            time_span.long_name             = 'The time span of the aggregated data'

            CALIPSO_Surface_Elevation_MEAN_mean.long_name = 'Mean of the Surface elevation MEAN'
            CALIPSO_Surface_Elevation_MIN_mean.long_name  = 'Mean of the Surface elevation MIN'
            CALIPSO_Surface_Elevation_MEAN_SD.long_name   = 'SD of the Surface elevation MEAN'
            CALIPSO_Surface_Elevation_MIN_SD.long_name    = 'SD of the Surface elevation MIN'

            Percentage_IGBP_1_id.long_name  = 'Evergreen-Needleleaf-Forest'
            Percentage_IGBP_2_id.long_name  = 'Evergreen-Broadleaf-Forest'
            Percentage_IGBP_3_id.long_name  = 'Deciduous-Needleleaf-Forest'
            Percentage_IGBP_4_id.long_name  = 'Deciduous-Broadleaf-Forest'
            Percentage_IGBP_5_id.long_name  = 'Mixed-Forest'
            Percentage_IGBP_6_id.long_name  = 'Closed-Shrublands'
            Percentage_IGBP_7_id.long_name  = 'Open-Shrubland (Desert)'
            Percentage_IGBP_8_id.long_name  = 'Woody-Savanna'
            Percentage_IGBP_9_id.long_name  = 'Savanna'
            Percentage_IGBP_10_id.long_name = 'Grassland'
            Percentage_IGBP_11_id.long_name = 'Wetland'
            Percentage_IGBP_12_id.long_name = 'Cropland'
            Percentage_IGBP_13_id.long_name = 'Urban'
            Percentage_IGBP_14_id.long_name = 'Crop-Mosaic'
            Percentage_IGBP_15_id.long_name = 'Permanent-Snow'
            Percentage_IGBP_16_id.long_name = 'Barren/Desert'
            Percentage_IGBP_17_id.long_name = 'Water'
            Percentage_IGBP_18_id.long_name = 'Tundra'

            ERA_Height_id.fill_value             = np.nan
            ERA_height_low_id.fill_value         = np.nan
            ERA_height_up_id.fill_value          = np.nan
            ERA_U_id.fill_value                  = np.nan
            ERA_U_SD_id.fill_value               = np.nan
            ERA_V_id.fill_value                  = np.nan
            ERA_V_SD_id.fill_value               = np.nan

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


            ##  Assign data to variables  -------------------------------------
            lats_id[:]                   = ERA_Latitude
            lons_id[:]                   = ERA_Longitude
            time[:]                      = ERA.time
            time_span[:]                 = ERA.time_span

            if TEST:
                ## test in a subdomain
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
                ERA_height_low_id[:]         = ERA.height_low.sel(latitude  = ERA_Latitude,
                                                                  longitude = ERA_Longitude)
                ERA_height_up_id[:]          = ERA.height_up.sel(latitude  = ERA_Latitude,
                                                                  longitude = ERA_Longitude)
            else:
                ## pass ERA5 ncdf as is
                ERA_Height_id[:]             = ERA.height
                ERA_U_id[:]                  = ERA[ U ]
                ERA_U_SD_id[:]               = ERA.u_SD
                ERA_V_id[:]                  = ERA[ V ]
                ERA_V_SD_id[:]               = ERA.v_SD
                ERA_height_low_id[:]         = ERA.height_low
                ERA_height_up_id[:]          = ERA.height_up

            ##  CALIPSO data
            CALIPSO_Surface_Elevation_MEAN_mean[:] = Final_Surface_Elevation_MEAN_mean
            CALIPSO_Surface_Elevation_MIN_mean [:] = Final_Surface_Elevation_MIN_mean
            CALIPSO_Surface_Elevation_MEAN_SD  [:] = Final_Surface_Elevation_MEAN_SD
            CALIPSO_Surface_Elevation_MIN_SD   [:] = Final_Surface_Elevation_MIN_SD

            ## flags coverage
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
                            creation       = format(datetime.now(timezone.utc)),
                            user_host      = os.getlogin() + "@" + os.uname()[1],
                            source_version = TRACE)
            for name, value in my_attrs.items():
                setattr(ds, name, value)

        ##  Create a indipedednt group for each percentile  -------------------
        Lg = ds.createGroup(f"LIVAS_p{aper}")

        ##  need to recreate coords within the group
        Lg.createDimension('CALIPSO_lev', len(Altitude))

        ##  Geo coordinates
        lats_id           = Lg.createDimension('lat', len(ERA_Latitude))
        lons_id           = Lg.createDimension('lon', len(ERA_Longitude))
        lats_id.units     = 'degrees_north'
        lons_id.units     = 'degrees_east'
        lats_id.long_name = 'Latitude'
        lons_id.long_name = 'Longitude'
        lats_id[:]        = ERA_Latitude
        lons_id[:]        = ERA_Longitude


        LIVAS_Altitude_id         = Lg.createVariable('Altitude',                        'f4',                    ('CALIPSO_lev',), zlib=True)
        LIVAS_a532nm_PD_id        = Lg.createVariable('Pure_Dust/LIVAS_PD_a532nm',        np.float64, ('lon','lat','CALIPSO_lev',), zlib=True)
        LIVAS_a532nm_PD_SD_id     = Lg.createVariable('Pure_Dust/LIVAS_PD_a532nm_STD',    np.float64, ('lon','lat','CALIPSO_lev',), zlib=True)
        LIVAS_b532nm_PD_id        = Lg.createVariable('Pure_Dust/LIVAS_PD_b532nm',        np.float64, ('lon','lat','CALIPSO_lev',), zlib=True)
        LIVAS_b532nm_PD_SD_id     = Lg.createVariable('Pure_Dust/LIVAS_PD_b532nm_STD',    np.float64, ('lon','lat','CALIPSO_lev',), zlib=True)
        LIVAS_PD_MC_id            = Lg.createVariable('Pure_Dust/LIVAS_PD_MC',            np.float64, ('lon','lat','CALIPSO_lev',), zlib=True)
        LIVAS_PD_MC_SD_id         = Lg.createVariable('Pure_Dust/LIVAS_PD_MC_STD',        np.float64, ('lon','lat','CALIPSO_lev',), zlib=True)
        LIVAS_PD_MC_Fine_id       = Lg.createVariable('Pure_Dust/LIVAS_PD_MC_Fine',       np.float64, ('lon','lat','CALIPSO_lev',), zlib=True)
        LIVAS_PD_MC_Fine_SD_id    = Lg.createVariable('Pure_Dust/LIVAS_PD_MC_Fine_STD',   np.float64, ('lon','lat','CALIPSO_lev',), zlib=True)
        LIVAS_PD_MC_Coarse_id     = Lg.createVariable('Pure_Dust/LIVAS_PD_MC_Coarse',     np.float64, ('lon','lat','CALIPSO_lev',), zlib=True)
        LIVAS_PD_MC_Coarse_SD_id  = Lg.createVariable('Pure_Dust/LIVAS_PD_MC_Coarse_STD', np.float64, ('lon','lat','CALIPSO_lev',), zlib=True)
        LIVAS_N_of_CF_Profiles_id = Lg.createVariable('Flags/Number_of_L2_CF_Profiles',   np.float64, ('lon','lat',), zlib=True)
        # LIVAS_N_of_Profiles_id    = Lg.createVariable('Flags/Number_of_L2_Profiles',      np.float64, ('lon','lat',), zlib=True)
        LIVAS_total_Profiles_id   = Lg.createVariable('Flags/Total_L2_Profiles_N',        np.float64, ('lon','lat',), zlib=True)
        LIVAS_DOD_532nm_mean      = Lg.createVariable('Pure_Dust/DOD_532nm_mean',         np.float64, ('lon','lat',), zlib=True)
        LIVAS_DOD_532nm_SD        = Lg.createVariable('Pure_Dust/DOD_532nm_STD',          np.float64, ('lon','lat',), zlib=True)

        LIVAS_Altitude_id.units         = 'm'
        LIVAS_a532nm_PD_id.units        = '1/m'       # changed altitude to meters
        LIVAS_a532nm_PD_SD_id.units     = '1/m'       # changed altitude to meters
        LIVAS_b532nm_PD_id.units        = '1/(m sr)'  # changed altitude to meters
        LIVAS_b532nm_PD_SD_id.units     = '1/(m sr)'  # changed altitude to meters
        LIVAS_PD_MC_id.units            = 'micrograms/m^3'
        LIVAS_PD_MC_SD_id.units         = 'micrograms/m^3'
        LIVAS_PD_MC_Fine_id.units       = 'micrograms/m^3'
        LIVAS_PD_MC_Fine_SD_id.units    = 'micrograms/m^3'
        LIVAS_PD_MC_Coarse_id.units     = 'micrograms/m^3'
        LIVAS_PD_MC_Coarse_SD_id.units  = 'micrograms/m^3'
        LIVAS_N_of_CF_Profiles_id.units = 'none'
        # LIVAS_N_of_Profiles_id.units    = 'none'
        LIVAS_DOD_532nm_mean.units      = 'none'
        LIVAS_DOD_532nm_SD.units        = 'none'

        LIVAS_Altitude_id.long_name         = 'Height'
        LIVAS_a532nm_PD_id.long_name        = 'Pure-Dust Extinction Coefficient 532nm'
        LIVAS_a532nm_PD_SD_id.long_name     = 'Pure-Dust Extinction Coefficient 532nm - SD'
        LIVAS_b532nm_PD_id.long_name        = 'Pure-Dust Backscatter Coefficient 532nm'
        LIVAS_b532nm_PD_SD_id.long_name     = 'Pure-Dust Backscatter Coefficient 532nm - SD'
        LIVAS_PD_MC_id.long_name            = 'Pure-Dust Mass Concentration'
        LIVAS_PD_MC_SD_id.long_name         = 'Pure-Dust Mass Concentration - SD'
        LIVAS_PD_MC_Fine_id.long_name       = 'Fine mode Pure-Dust Mass Concentration'
        LIVAS_PD_MC_Fine_SD_id.long_name    = 'Fine mode Pure-Dust Mass Concentration - SD'
        LIVAS_PD_MC_Coarse_id.long_name     = 'Coarse mode Pure-Dust Mass Concentration'
        LIVAS_PD_MC_Coarse_SD_id.long_name  = 'Coarse mode Pure-Dust Mass Concentration - SD'
        LIVAS_DOD_532nm_mean.long_name      = 'Dust Optical Depth 532nm - mean'
        LIVAS_DOD_532nm_SD.long_name        = 'Dust Optical Depth 532nm - SD'
        LIVAS_N_of_CF_Profiles_id.long_name = 'Number of CALIPSO L2 5km Cloud Free Profiles'
        # LIVAS_N_of_Profiles_id.long_name    = 'Number of CALIPSO L2 5km Profiles'
        LIVAS_total_Profiles_id.long_name   = 'Number of profiles in the ERA5 cell'

        LIVAS_Altitude_id.fill_value         = np.nan
        LIVAS_a532nm_PD_id.fill_value        = np.nan
        LIVAS_a532nm_PD_SD_id.fill_value     = np.nan
        LIVAS_b532nm_PD_id.fill_value        = np.nan
        LIVAS_b532nm_PD_SD_id.fill_value     = np.nan
        LIVAS_PD_MC_id.fill_value            = np.nan
        LIVAS_PD_MC_SD_id.fill_value         = np.nan
        LIVAS_N_of_CF_Profiles_id.fill_value = np.nan
        # LIVAS_N_of_Profiles_id.fill_value    = np.nan
        LIVAS_DOD_532nm_mean.fill_value      = np.nan
        LIVAS_DOD_532nm_SD.fill_value        = np.nan

        ##  Assign LIVAS data
        LIVAS_Altitude_id[:]         = Altitude
        LIVAS_a532nm_PD_id[:]        = Final_PD_a532nm
        LIVAS_a532nm_PD_SD_id[:]     = Final_PD_a532nm_SD
        LIVAS_b532nm_PD_id[:]        = Final_PD_b532nm
        LIVAS_b532nm_PD_SD_id[:]     = Final_PD_b532nm_SD
        LIVAS_PD_MC_id[:]            = Final_PD_MC
        LIVAS_PD_MC_SD_id[:]         = Final_PD_MC_SD
        LIVAS_PD_MC_Fine_id[:]       = Final_PD_MC_Fine
        LIVAS_PD_MC_Fine_SD_id[:]    = Final_PD_MC_Fine_SD
        LIVAS_PD_MC_Coarse_id[:]     = Final_PD_MC_Coarse
        LIVAS_PD_MC_Coarse_SD_id[:]  = Final_PD_MC_Coarse_SD
        LIVAS_N_of_CF_Profiles_id[:] = Final_Number_of_L2Profiles
        # LIVAS_N_of_Profiles_id[:]    = Final_Number_of_Profiles
        LIVAS_DOD_532nm_mean[:]      = Final_LIVAS_PD_DOD_532nm
        LIVAS_DOD_532nm_SD[:]        = Final_LIVAS_PD_DOD_532nm_SD
        LIVAS_total_Profiles_id[:]   = Total_Number_of_Profiles

        ds.close()
        print(f"Written: {fileout}")

#  SCRIPT END  ---------------------------------------------------------------
Ou.goodbye(cnf.LOGs.run, tic = tic, scriptname = __file__, version = TRACE)
