#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

@author: proestakis, thanasisn
"""

import os
import sys
import glob
import warnings

import netCDF4  as nc
import numpy    as np
import xarray   as xr
import pandas   as pd

from   datetime    import date, datetime, timedelta
from   pathlib import Path
import numpy.ma as ma
from matplotlib import pyplot as plt
# import geopy.distance


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
# TEST = cnf.mode.Test
# TEST = True

##  Allow only one case to run at the time  ----------------------------------
if SEASONAL == MONTHLY:
    print("Seasonal:", SEASONAL)
    print("Monthly: ", MONTHLY)
    sys.exit("Choose only SEASONAL or only MONTHLY")

##  Check destination folder exists  -----------------------------------------
if not os.path.isdir(cnf.DOMOS.path_output):
    sys.exit(f"\nFolder {cnf.DOMOS.path_output} don't exist !!\n")

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
    DOMOS_filenames = glob.glob(
        f"{cnf.OREO.path_output}/Seasonal_{cnf.D1.LatStep}x{cnf.D1.LonStep}/ERA5_LIVAS_*_lat_{np.min(REGRID_lat_centers)}_{np.max(REGRID_lat_centers)}_lon_{np.min(REGRID_lon_centers)}_{np.max(REGRID_lon_centers)}.nc"
    )
elif MONTHLY == True:
    print("Work on monthly data")
    DOMOS_filenames = glob.glob(
        f"{cnf.OREO.path_output}/Monthly_{cnf.D1.LatStep}x{cnf.D1.LonStep}/ERA5_LIVAS_*_lat_{np.min(REGRID_lat_centers)}_{np.max(REGRID_lat_centers)}_lon_{np.min(REGRID_lon_centers)}_{np.max(REGRID_lon_centers)}.nc"
    )

DOMOS_filenames.sort()
## TEST
# shuffle(ERA_filenames)


# directory paths of input datasets
DOMOS_EO_datasets_path    = r"D:\PC-backup\Projects\DOMOS\Datasets\2x5_deg_seasonal_mean\Step_II_V45"
# directory path of output datasets
output_path               = r"D:\PC-backup\Projects\DOMOS\Datasets\2x5_deg_seasonal_mean\Step_III_V45"
# directory path of output datasets
IGBP_Surface_Type_path    = r"C:\Users\proes\Desktop\PC-backup\Projects\DOMOS\Datasets\2x5_deg_seasonal_mean"


#########     Functions     ########

# def Get_Boundaries_of_ERA5_UVW_wind_components(Longitude,Latitude,Height):

#     Height_Boundaries    = np.array(np.empty((Height.shape[0],Height.shape[1],Height.shape[2]+1)))
#     Height_Boundaries[:] = np.nan
#     for count_lon,lon in enumerate(Longitude):
#         for count_lat,lat in enumerate(Latitude):
#             H       = Height[count_lon,count_lat,:].data
#             temp    = np.array(np.empty((H.shape[0]+1)))
#             temp[:] = np.nan
#             if H[H.shape[0]-1] > 0:
#                 for count_alt,alt in enumerate(temp):
#                     if count_alt == 0:
#                         Height_Boundaries[count_lon,count_lat,0] = H[0] + 1000.0
#                     if count_alt == temp.shape[0]-1:
#                         Height_Boundaries[count_lon,count_lat,count_alt] = 0.0
#                     if ((count_alt > 0) & (count_alt < temp.shape[0]-1)):
#                         Height_Boundaries[count_lon,count_lat,count_alt] = (H[count_alt-1] + H[count_alt])/2
#             else:
#                 for count_alt,alt in enumerate(temp):
#                     if count_alt == 0:
#                         Height_Boundaries[count_lon,count_lat,0] = H[0] + 1000.0
#                     if count_alt == temp.shape[0]-1:
#                         Height_Boundaries[count_lon,count_lat,count_alt] = -100.0
#                     if count_alt == temp.shape[0]-2:
#                         Height_Boundaries[count_lon,count_lat,count_alt] =    0.0
#                     if ((count_alt > 0) & (count_alt < temp.shape[0]-2)):
#                         Height_Boundaries[count_lon,count_lat,count_alt] = (H[count_alt-1] + H[count_alt])/2

#     return(Height_Boundaries)

def Mass_Concentration_along_ERA5_U_and_V_wind_directions(Longitude,Latitude,Total_Mass,U,U_SD,V,V_SD):

    MC_U = Total_Mass.copy()
    MC_V = Total_Mass.copy()

    # in case of masked data, replace with nans #
    if ma.isMaskedArray(Total_Mass) == True:
       MC_U[MC_U.mask  == True ] = 0
       MC_V[MC_V.mask  == True ] = 0

    # Computation of MC_U and MC_V
    for count_lon,lon in enumerate(Longitude):
        for count_lat,lat in enumerate(Latitude):
            for count_u_component,u_component in enumerate(U[count_lon,count_lat,:]):
                if Total_Mass[count_lon,count_lat,count_u_component] <= 0:
                    MC_U[count_lon,count_lat,count_u_component] = 0
                    MC_V[count_lon,count_lat,count_u_component] = 0
                else:
                    u_term = U[count_lon,count_lat,count_u_component]
                    v_term = V[count_lon,count_lat,count_u_component]
                    ### Mass transport along the U direction ###
                    nominator   = np.square(u_term)
                    denominator = np.add(np.square(u_term),np.square(v_term))
                    multiplication_factor = np.divide(nominator,denominator)
                    MC_U[count_lon,count_lat,count_u_component] = np.multiply(multiplication_factor,Total_Mass[count_lon,count_lat,count_u_component])
                    ### Mass transport along the V direction ###
                    nominator   = np.square(v_term)
                    denominator = np.add(np.square(u_term),np.square(v_term))
                    multiplication_factor = np.divide(nominator,denominator)
                    MC_V[count_lon,count_lat,count_u_component] = np.multiply(multiplication_factor,Total_Mass[count_lon,count_lat,count_u_component])


    return(MC_U,MC_V)

def find_distance_between_coordinates(coord1, coord2):
    Distance = geopy.distance.geodesic(coord1, coord2).km
    return(Distance)

def computation_of_Height_Range(Longitude,Latitude,Height,Height_Boundaries):

    Height_Range    = Height.copy()
    Height_Range[:] = np.nan
    for count_x,x in enumerate(Longitude):
        for count_y,y in enumerate(Latitude):
            idx_lat = np.where(y == Latitude)
            idx_lon = np.where(x == Longitude)
            idx_lat = np.squeeze(np.asarray(idx_lat))
            idx_lon = np.squeeze(np.asarray(idx_lon))
            H = Height_Boundaries[idx_lon,idx_lat,:]
            if H[len(H)-1] == -100:
                for count_z,z in enumerate(H):
                    if count_z < len(H)-2:
                        Height_Range[count_x,count_y,count_z+1] = H[count_z]-H[count_z+1]
            if H[len(H)-1] != -100:
                for count_z,z in enumerate(H):
                    if count_z < len(H)-1:
                        Height_Range[count_x,count_y,count_z] = H[count_z]-H[count_z+1]
    return(Height_Range)

def computation_of_dust_deposition(Longitude,Latitude,Percentage_IGBP_17,PD_MC_ERA5res_V,PD_MC_ERA5res_U,Height_Range,LIVAS_DOD_532nm_mean):

    Mask    = np.empty((len(Longitude),len(Latitude)))
    Mask[:] = 0.0

    F_N = PD_MC_ERA5res_V.copy()
    F_N[:] = np.nan
    F_S = PD_MC_ERA5res_V.copy()
    F_S[:] = np.nan
    F_E = PD_MC_ERA5res_U.copy()
    F_E[:] = np.nan
    F_W = PD_MC_ERA5res_U.copy()
    F_W[:] = np.nan

    if np.count_nonzero(~np.isnan(LIVAS_DOD_532nm_mean)) == 0:
       Dust_Deposition    = np.empty((len(Longitude),len(Latitude)))
       Dust_Deposition[:] = np.nan

    if np.count_nonzero(~np.isnan(LIVAS_DOD_532nm_mean)) != 0:
       Dust_Deposition    = np.empty((len(Longitude),len(Latitude)))
       Dust_Deposition[:] = np.nan

       for count_lon,lon in enumerate(Longitude):
           for count_lat,lat in enumerate(Latitude):

               if (count_lon < 1) | (count_lon > len(Longitude)-2):
                   continue
               if (count_lat < 1) | (count_lat > len(Latitude)-2):
                   continue

               idx_lat = count_lat
               idx_lon = count_lon

               if ( LIVAS_DOD_532nm_mean[idx_lon,idx_lat] <= 0.01) | (np.isnan(LIVAS_DOD_532nm_mean[idx_lon,idx_lat]) ):
                   Cell_C   = Percentage_IGBP_17[idx_lon,idx_lat]
                   if (Cell_C < 20):
                       Dust_Deposition[idx_lon,idx_lat] = np.nan
                   else:
                       if (np.isnan(LIVAS_DOD_532nm_mean[idx_lon,idx_lat])):
                           Dust_Deposition[idx_lon,idx_lat] = np.nan
                       else:
                           Dust_Deposition[idx_lon,idx_lat] = 0.0
               else:
                   Cell_C   = Percentage_IGBP_17[idx_lon,idx_lat]
                   if (Cell_C < 20):
                       Dust_Deposition[idx_lon,idx_lat] = np.nan
                   else:
                       N_edge  = find_distance_between_coordinates((lat+1,lon-2.5),(lat+1,lon+2.5))*1000
                       S_edge  = find_distance_between_coordinates((lat-1,lon-2.5),(lat-1,lon+2.5))*1000
                       Horizontal_edge = (N_edge + S_edge)/2
                       E_edge  = find_distance_between_coordinates((lat-1,lon+2.5),(lat+1,lon+2.5))*1000
                       W_edge  = find_distance_between_coordinates((lat-1,lon-2.5),(lat+1,lon-2.5))*1000

                       MC_V_N  = PD_MC_ERA5res_V[idx_lon,idx_lat+1,:]
                       MC_V_S  = PD_MC_ERA5res_V[idx_lon,idx_lat-1,:]
                       MC_U_W  = PD_MC_ERA5res_U[idx_lon-1,idx_lat,:]
                       MC_U_E  = PD_MC_ERA5res_U[idx_lon+1,idx_lat,:]

                       V_N     = V[idx_lon,idx_lat+1,:]
                       V_S     = V[idx_lon,idx_lat-1,:]
                       U_W     = U[idx_lon-1,idx_lat,:]
                       U_E     = U[idx_lon+1,idx_lat,:]

                       Height_Range_N = Height_Range[idx_lon,idx_lat+1,:]
                       Height_Range_S = Height_Range[idx_lon,idx_lat-1,:]
                       Height_Range_W = Height_Range[idx_lon-1,idx_lat,:]
                       Height_Range_E = Height_Range[idx_lon+1,idx_lat,:]

                       ##### Fluxes N  #####
                       if np.count_nonzero(~np.isnan(MC_V_N)) == 0:
                           Flux_N = np.nan
                           F_N[idx_lon,idx_lat,:] = np.nan
                       else:
                           Flux   = np.multiply(np.nansum(np.multiply(np.multiply(np.multiply(MC_V_N,V_N),Height_Range_N),Horizontal_edge)),86.4)
                           temp_N = np.multiply(np.multiply(np.multiply(np.multiply(MC_V_N,V_N),Height_Range_N),Horizontal_edge),86.4)
                           if Flux > 0:
                               Flux_N = Flux*(-1)
                               temp_N = temp_N*(-1)
                           if Flux < 0:
                               Flux_N = Flux*(-1)
                               temp_N = temp_N*(-1)
                           if Flux == 0:
                               Flux_N = 0.0
                           for i in range(len(MC_V_N)):
                               F_N[idx_lon,idx_lat,i] = temp_N[i]

                       ##### Fluxes S  #####
                       if np.count_nonzero(~np.isnan(MC_V_S)) == 0:
                           Flux_S = np.nan
                       else:
                           Flux_S = np.multiply(np.nansum(np.multiply(np.multiply(np.multiply(MC_V_S,V_S),Height_Range_S),Horizontal_edge)),86.4)
                           temp_S = np.multiply(np.multiply(np.multiply(np.multiply(MC_V_S,V_S),Height_Range_S),Horizontal_edge),86.4)
                           for i in range(len(MC_V_S)):
                               F_S[idx_lon,idx_lat,i] = temp_S[i]

                       ##### Fluxes E  #####
                       if np.count_nonzero(~np.isnan(MC_U_E)) == 0:
                           Flux_E = np.nan
                       else:
                           Flux = np.multiply(np.nansum(np.multiply(np.multiply(np.multiply(MC_U_E,U_E),Height_Range_E),E_edge)),86.4)
                           temp_E = np.multiply(np.multiply(np.multiply(np.multiply(MC_U_E,U_E),Height_Range_E),E_edge),86.4)
                           if Flux > 0:
                               Flux_E = Flux*(-1)
                           if Flux < 0:
                               Flux_E = Flux*(-1)
                           if Flux == 0:
                               Flux_E = 0.0
                           for i in range(len(MC_U_E)):
                               F_E[idx_lon,idx_lat,i] = temp_E[i]

                       ##### Fluxes W  #####
                       if np.count_nonzero(~np.isnan(MC_U_W)) == 0:
                           Flux_W = np.nan
                       else:
                           Flux_W = np.multiply(np.nansum(np.multiply(np.multiply(np.multiply(MC_U_W,U_W),Height_Range_W),W_edge)),86.4)
                           temp_W = np.multiply(np.multiply(np.multiply(np.multiply(MC_U_W,U_W),Height_Range_W),W_edge),86.4)
                           for i in range(len(MC_U_W)):
                               F_W[idx_lon,idx_lat,i] = temp_W[i]

                       ### Deposition Area ###
                       Area = (Horizontal_edge*E_edge)

                       ########################
                       ##### Fluxes North #####
                       Flux_North = Flux_N
                       ##### Fluxes South #####
                       Flux_South = Flux_S
                       ##### Fluxes East #####
                       Flux_East  = Flux_E
                       ##### Fluxes West #####
                       Flux_West  = Flux_W

                       Fluxes = np.array([Flux_North, Flux_South, Flux_East, Flux_West])
                       if np.count_nonzero(~np.isnan(Fluxes)) != len(Fluxes):
                           Dust_Deposited = np.nan
                       Dust_Deposited = (Flux_North+Flux_South+Flux_West+Flux_East)/Area
                       if Dust_Deposited < 0:
                           Dust_Deposited = np.nan
                       Dust_Deposition[idx_lon,idx_lat] = Dust_Deposited

    return(F_N,F_S,F_E,F_W,Dust_Deposition)


##########################################################################################
####                                main                                              ####
##########################################################################################




for i_DOMOS, DOMOS_file in enumerate(DOMOS_filenames):



    ### extracting "yyyymm" sufix from DOMOS filename ###
    DOMOS_substrings  = DOMOS_file.split("\\")
    DOMOS_substrings  = DOMOS_substrings[len(DOMOS_substrings)-1]
    DOMOS_substrings  = DOMOS_substrings.split(".")
    DOMOS_substrings  = DOMOS_substrings[len(DOMOS_substrings)-2]
    DOMOS_substrings  = DOMOS_substrings.split("_")
    DOMOS_file_year   = DOMOS_substrings[len(DOMOS_substrings)-2]
    DOMOS_file_season = DOMOS_substrings[len(DOMOS_substrings)-1]
    DOMOS_sufix       = DOMOS_file_year + '_' + DOMOS_file_season

    ### Reading NC variables ###
    DOMOS_dataset = nc.Dataset(DOMOS_file)
    DOMOS         = xr.open_dataset(DOMOS_file)

    DOMOS.data_date

    sys.exit("paths")

    ERA_Latitude  = ERA.latitude
    ERA_Longitude = ERA.longitude

    ##  Get the ERA file year and month
    ERA_year      = pd.DatetimeIndex(DOMOS.data_date).year[0]
    ERA_month     = pd.DatetimeIndex(DOMOS.data_date).month[0]



    Latitude               = DOMOS_dataset['Geolocation/Latitude'][:]
    Longitude              = DOMOS_dataset['Geolocation/Longitude'][:]

    Height                 = DOMOS_dataset['ERA5/ERA5_Height'][:]
    U                      = DOMOS_dataset['ERA5/U'][:]
    U_SD                   = DOMOS_dataset['ERA5/U_SD'][:]
    V                      = DOMOS_dataset['ERA5/V'][:]
    V_SD                   = DOMOS_dataset['ERA5/V_SD'][:]

    LIVAS_Altitude         = DOMOS_dataset['LIVAS/Altitude'][:]
    LIVAS_a532nm_PD        = DOMOS_dataset['LIVAS/Pure_Dust/LIVAS_PD_a532nm'][:]
    LIVAS_a532nm_PD_SD     = DOMOS_dataset['LIVAS/Pure_Dust/LIVAS_PD_a532nm_STD'][:]
    LIVAS_PD_MC            = DOMOS_dataset['LIVAS/Pure_Dust/LIVAS_PD_MC'][:]
    LIVAS_PD_MC_SD         = DOMOS_dataset['LIVAS/Pure_Dust/LIVAS_PD_MC_STD'][:]
    LIVAS_PD_MC_CM         = DOMOS_dataset['LIVAS/Pure_Dust/LIVAS_PD_MC_CM'][:]
    LIVAS_PD_MC_CM_SD      = DOMOS_dataset['LIVAS/Pure_Dust/LIVAS_PD_MC_CM_STD'][:]
    LIVAS_PD_MC_FM         = DOMOS_dataset['LIVAS/Pure_Dust/LIVAS_PD_MC_FM'][:]
    LIVAS_PD_MC_FM_SD      = DOMOS_dataset['LIVAS/Pure_Dust/LIVAS_PD_MC_FM_STD'][:]
    LIVAS_N_of_CF_Profiles = DOMOS_dataset['LIVAS/Flags/Number_of_L2_CF_Profiles'][:]
    LIVAS_N_of_Profiles    = DOMOS_dataset['LIVAS/Flags/Number_of_L2_Profiles'][:]
    LIVAS_DOD_532nm_mean   = DOMOS_dataset['LIVAS/Pure_Dust/DOD_532nm_mean'][:]
    LIVAS_DOD_532nm_SD     = DOMOS_dataset['LIVAS/Pure_Dust/DOD_532nm_STD'][:]

    Percentage_IGBP_1  = DOMOS_dataset['Land_Ocean_Mask/Percentage_IGBP_1'][:]
    Percentage_IGBP_2  = DOMOS_dataset['Land_Ocean_Mask/Percentage_IGBP_2'][:]
    Percentage_IGBP_3  = DOMOS_dataset['Land_Ocean_Mask/Percentage_IGBP_3'][:]
    Percentage_IGBP_4  = DOMOS_dataset['Land_Ocean_Mask/Percentage_IGBP_4'][:]
    Percentage_IGBP_5  = DOMOS_dataset['Land_Ocean_Mask/Percentage_IGBP_5'][:]
    Percentage_IGBP_6  = DOMOS_dataset['Land_Ocean_Mask/Percentage_IGBP_6'][:]
    Percentage_IGBP_7  = DOMOS_dataset['Land_Ocean_Mask/Percentage_IGBP_7'][:]
    Percentage_IGBP_8  = DOMOS_dataset['Land_Ocean_Mask/Percentage_IGBP_8'][:]
    Percentage_IGBP_9  = DOMOS_dataset['Land_Ocean_Mask/Percentage_IGBP_9'][:]
    Percentage_IGBP_10 = DOMOS_dataset['Land_Ocean_Mask/Percentage_IGBP_10'][:]
    Percentage_IGBP_11 = DOMOS_dataset['Land_Ocean_Mask/Percentage_IGBP_11'][:]
    Percentage_IGBP_12 = DOMOS_dataset['Land_Ocean_Mask/Percentage_IGBP_12'][:]
    Percentage_IGBP_13 = DOMOS_dataset['Land_Ocean_Mask/Percentage_IGBP_13'][:]
    Percentage_IGBP_14 = DOMOS_dataset['Land_Ocean_Mask/Percentage_IGBP_14'][:]
    Percentage_IGBP_15 = DOMOS_dataset['Land_Ocean_Mask/Percentage_IGBP_15'][:]
    Percentage_IGBP_16 = DOMOS_dataset['Land_Ocean_Mask/Percentage_IGBP_16'][:]
    Percentage_IGBP_17 = DOMOS_dataset['Land_Ocean_Mask/Percentage_IGBP_17'][:]
    Percentage_IGBP_18 = DOMOS_dataset['Land_Ocean_Mask/Percentage_IGBP_18'][:]

    ### Height Boundaries ###
    Height_Boundaries = Get_Boundaries_of_ERA5_UVW_wind_components(Longitude,Latitude,Height)

    ###########################################################################
    ### LIVAS ###                                                             #
    # Computation_of_LIVAS_based_EO_MC_profiles_in_ERA5_resolution            #
    ###########################################################################

    if ma.isMaskedArray(Height) == True:
       Height[Height.mask  == True ] = np.nan
       Height = Height.data
    if ma.isMaskedArray(LIVAS_Altitude) == True:
          LIVAS_Altitude[LIVAS_Altitude.mask  == True ] = np.nan
          LIVAS_Altitude = LIVAS_Altitude.data

    PD_MC_ERA5res        = Height.copy()
    PD_MC_ERA5res[:]     = 0

    for count_lon,lon in enumerate(Longitude):
        for count_lat,lat in enumerate(Latitude):

            if ma.isMaskedArray(LIVAS_PD_MC) == True:
                MC          = LIVAS_PD_MC[count_lon,count_lat,:].data
            else:
                MC          = LIVAS_PD_MC[count_lon,count_lat,:]
            ERA5_Height = Height[count_lon,count_lat,:]

            for alt_count,alt in enumerate(ERA5_Height):

                if alt_count == 0:
                    PD_MC_ERA5res[count_lon,count_lat,alt_count]     = np.nan
                if alt_count == len(ERA5_Height)-1:
                    H    = (ERA5_Height[alt_count] + ERA5_Height[alt_count-1])/2
                    idx  = np.where(LIVAS_Altitude < H)
                    PD_MC_ERA5res[count_lon,count_lat,alt_count] = np.nanmean(MC[idx])
                if ((alt_count > 0) and (alt_count < len(ERA5_Height)-1)):
                    H1   = (ERA5_Height[alt_count] + ERA5_Height[alt_count-1])/2
                    H2   = (ERA5_Height[alt_count] + ERA5_Height[alt_count+1])/2
                    idx  = np.where((LIVAS_Altitude > H2) & (LIVAS_Altitude < H1))
                    if len(idx) == 0:
                        PD_MC_ERA5res[count_lon,count_lat,alt_count] = 0
                    else:
                        PD_MC_ERA5res[count_lon,count_lat,alt_count] = np.nanmean(MC[idx])

    sys.exit("Work to here")
    ### Height Boundaries and Range ###
    Height_Boundaries = Get_Boundaries_of_ERA5_UVW_wind_components(Longitude,Latitude,Height)
    Height_Range      = computation_of_Height_Range(Longitude,Latitude,Height,Height_Boundaries)

    # Transport of MC along ERA5 U,V wind directions
    # and computation of MC_U_SD and MC_V_SD
    ### LIVAS ###
    [ LIVAS_PD_MC_ERA5res_U, LIVAS_PD_MC_ERA5res_V ] = Mass_Concentration_along_ERA5_U_and_V_wind_directions(Longitude,Latitude,PD_MC_ERA5res,U,U_SD,V,V_SD)

    # Transport of MC along ERA5 U,V wind directions
    # and computation of MC_U_SD and MC_V_SD
    ### LIVAS ###

    ###########################################################
    # "======== Computation of Dust Deposition     =========" #
    ###########################################################

    # [LIVAS_F_N,LIVAS_F_S,LIVAS_F_E,LIVAS_F_W,LIVAS_Dust_Deposition_Rate] = computation_of_dust_deposition(Longitude,Latitude,Percentage_IGBP_17,LIVAS_PD_MC_ERA5res_V,LIVAS_PD_MC_ERA5res_U,Height_Range,LIVAS_DOD_532nm_mean)

    Mask    = np.empty((len(Longitude),len(Latitude)))
    Mask[:] = 0.0

    F_N = LIVAS_PD_MC_ERA5res_V.copy()
    F_N[:] = np.nan
    F_S = LIVAS_PD_MC_ERA5res_V.copy()
    F_S[:] = np.nan
    F_E = LIVAS_PD_MC_ERA5res_U.copy()
    F_E[:] = np.nan
    F_W = LIVAS_PD_MC_ERA5res_U.copy()
    F_W[:] = np.nan

    if np.count_nonzero(~np.isnan(LIVAS_DOD_532nm_mean)) == 0:
       Dust_Deposition    = np.empty((len(Longitude),len(Latitude)))
       Dust_Deposition[:] = np.nan

    if np.count_nonzero(~np.isnan(LIVAS_DOD_532nm_mean)) != 0:
       Dust_Deposition    = np.empty((len(Longitude),len(Latitude)))
       Dust_Deposition[:] = np.nan

       for count_lon,lon in enumerate(Longitude):
           for count_lat,lat in enumerate(Latitude):

               if (count_lon < 1) | (count_lon > len(Longitude)-2):
                   continue
               if (count_lat < 1) | (count_lat > len(Latitude)-2):
                   continue

               idx_lat = count_lat
               idx_lon = count_lon

               if ( LIVAS_DOD_532nm_mean[idx_lon,idx_lat] <= 0.01) | (np.isnan(LIVAS_DOD_532nm_mean[idx_lon,idx_lat]) ):
                   Cell_C_water   = Percentage_IGBP_17[idx_lon,idx_lat]
                   Cell_C_Desert  = Percentage_IGBP_6[idx_lon,idx_lat] + Percentage_IGBP_7[idx_lon,idx_lat] + Percentage_IGBP_16[idx_lon,idx_lat]
                   if (Cell_C_water <= 50) | (Cell_C_Desert >= 10):
                       Dust_Deposition[idx_lon,idx_lat] = np.nan
                   else:
                       if (np.isnan(LIVAS_DOD_532nm_mean[idx_lon,idx_lat])):
                           Dust_Deposition[idx_lon,idx_lat] = np.nan
                       else:
                           Dust_Deposition[idx_lon,idx_lat] = 0.0

               else:
                   Cell_C_water   = Percentage_IGBP_17[idx_lon,idx_lat]
                   Cell_C_Desert  = Percentage_IGBP_6[idx_lon,idx_lat] + Percentage_IGBP_7[idx_lon,idx_lat] + Percentage_IGBP_16[idx_lon,idx_lat]
                   if (Cell_C_water <= 50) | (Cell_C_Desert >= 10):
                       Dust_Deposition[idx_lon,idx_lat] = np.nan
                   else:
                       N_edge  = find_distance_between_coordinates((lat+1,lon-2.5),(lat+1,lon+2.5))*1000
                       S_edge  = find_distance_between_coordinates((lat-1,lon-2.5),(lat-1,lon+2.5))*1000
                       Horizontal_edge = (N_edge + S_edge)/2
                       E_edge  = find_distance_between_coordinates((lat-1,lon+2.5),(lat+1,lon+2.5))*1000
                       W_edge  = find_distance_between_coordinates((lat-1,lon-2.5),(lat+1,lon-2.5))*1000

                       MC_V_N  = LIVAS_PD_MC_ERA5res_V[idx_lon,idx_lat+1,:]
                       MC_V_S  = LIVAS_PD_MC_ERA5res_V[idx_lon,idx_lat-1,:]
                       MC_U_W  = LIVAS_PD_MC_ERA5res_U[idx_lon-1,idx_lat,:]
                       MC_U_E  = LIVAS_PD_MC_ERA5res_U[idx_lon+1,idx_lat,:]

                       V_N     = V[idx_lon,idx_lat+1,:]
                       V_S     = V[idx_lon,idx_lat-1,:]
                       U_W     = U[idx_lon-1,idx_lat,:]
                       U_E     = U[idx_lon+1,idx_lat,:]

                       Height_Range_N = Height_Range[idx_lon,idx_lat+1,:]
                       Height_Range_S = Height_Range[idx_lon,idx_lat-1,:]
                       Height_Range_W = Height_Range[idx_lon-1,idx_lat,:]
                       Height_Range_E = Height_Range[idx_lon+1,idx_lat,:]

                       ##### Fluxes N  #####
                       if np.count_nonzero(~np.isnan(MC_V_N)) == 0:
                           Flux_N = np.nan
                           F_N[idx_lon,idx_lat,:] = np.nan
                       else:
                           Flux   = np.multiply(np.nansum(np.multiply(np.multiply(np.multiply(MC_V_N,V_N),Height_Range_N),Horizontal_edge)),86.4)
                           temp_N = np.multiply(np.multiply(np.multiply(np.multiply(MC_V_N,V_N),Height_Range_N),Horizontal_edge),86.4)
                           if Flux > 0:
                               Flux_N = Flux*(-1)
                               temp_N = temp_N*(-1)
                           if Flux < 0:
                               Flux_N = Flux*(-1)
                               temp_N = temp_N*(-1)
                           if Flux == 0:
                               Flux_N = 0.0
                           for i in range(len(MC_V_N)):
                               F_N[idx_lon,idx_lat,i] = temp_N[i]

                       ##### Fluxes S  #####
                       if np.count_nonzero(~np.isnan(MC_V_S)) == 0:
                           Flux_S = np.nan
                       else:
                           Flux_S = np.multiply(np.nansum(np.multiply(np.multiply(np.multiply(MC_V_S,V_S),Height_Range_S),Horizontal_edge)),86.4)
                           temp_S = np.multiply(np.multiply(np.multiply(np.multiply(MC_V_S,V_S),Height_Range_S),Horizontal_edge),86.4)
                           for i in range(len(MC_V_S)):
                               F_S[idx_lon,idx_lat,i] = temp_S[i]

                       ##### Fluxes E  #####
                       if np.count_nonzero(~np.isnan(MC_U_E)) == 0:
                           Flux_E = np.nan
                       else:
                           Flux = np.multiply(np.nansum(np.multiply(np.multiply(np.multiply(MC_U_E,U_E),Height_Range_E),E_edge)),86.4)
                           temp_E = np.multiply(np.multiply(np.multiply(np.multiply(MC_U_E,U_E),Height_Range_E),E_edge),86.4)
                           if Flux > 0:
                               Flux_E = Flux*(-1)
                           if Flux < 0:
                               Flux_E = Flux*(-1)
                           if Flux == 0:
                               Flux_E = 0.0
                           for i in range(len(MC_U_E)):
                               F_E[idx_lon,idx_lat,i] = temp_E[i]

                       ##### Fluxes W  #####
                       if np.count_nonzero(~np.isnan(MC_U_W)) == 0:
                           Flux_W = np.nan
                       else:
                           Flux_W = np.multiply(np.nansum(np.multiply(np.multiply(np.multiply(MC_U_W,U_W),Height_Range_W),W_edge)),86.4)
                           temp_W = np.multiply(np.multiply(np.multiply(np.multiply(MC_U_W,U_W),Height_Range_W),W_edge),86.4)
                           for i in range(len(MC_U_W)):
                               F_W[idx_lon,idx_lat,i] = temp_W[i]

                       ### Deposition Area ###
                       Area = (Horizontal_edge*E_edge)

                       ########################
                       ##### Fluxes North #####
                       Flux_North = Flux_N
                       ##### Fluxes South #####
                       Flux_South = Flux_S
                       ##### Fluxes East #####
                       Flux_East  = Flux_E
                       ##### Fluxes West #####
                       Flux_West  = Flux_W

                       Fluxes = np.array([Flux_North, Flux_South, Flux_East, Flux_West])
                       if np.count_nonzero(~np.isnan(Fluxes)) != len(Fluxes):
                           Dust_Deposited = np.nan
                       Dust_Deposited = (Flux_North+Flux_South+Flux_West+Flux_East)/Area
                       if Dust_Deposited < 0:
                           Dust_Deposited = np.nan
                       Dust_Deposition[idx_lon,idx_lat] = Dust_Deposited


    DDR = Dust_Deposition.copy()
    Dust_Deposition_Horizontally_smoothed    = Dust_Deposition.copy()
    Dust_Deposition_Horizontally_smoothed[:] = np.nan

    for count_lon,lon in enumerate(Longitude):
        for count_lat,lat in enumerate(Latitude):

            if (count_lon == 0) or (count_lon == len(Longitude)-1):
                continue
            if (count_lat == 0) or (count_lat == len(Latitude)-1):
                continue

            idx_lat = count_lat
            idx_lon = count_lon
            temp    = [ DDR[idx_lon-1,idx_lat],DDR[idx_lon,idx_lat],DDR[idx_lon+1,idx_lat],DDR[idx_lon-1,idx_lat+1],DDR[idx_lon,idx_lat+1],DDR[idx_lon+1,idx_lat+1],DDR[idx_lon-1,idx_lat-1],DDR[idx_lon,idx_lat-1],DDR[idx_lon+1,idx_lat-1]  ]

            Dust_Deposition_Horizontally_smoothed[idx_lon,idx_lat] = np.nanmean(temp)

    IGBP_Surface_Type            = IGBP_Surface_Type_path + '\\' + 'Land_Ocean_Percentage.nc'
    IGBP_Surface_Type_percentage = nc.Dataset(IGBP_Surface_Type)
    Land_Ocean_Percentage        = IGBP_Surface_Type_percentage['/Land_Ocean_Percentage'][:]
    idx  = np.where(Land_Ocean_Percentage == 0.0)
    Dust_Deposition_Horizontally_smoothed[idx] = np.nan

    #####################################
    #  --- Saving dataset as NetCDF --- #
    #####################################

    # creating nc. filename and initiallizing:
    fn           = output_path + '\\' + 'DOMOS_Datasets_DOD-and-DDR-V3-GRID_RESOLUTION_2x5_' + DOMOS_file_season + '_' + DOMOS_file_year + '.nc'
    ds           = nc.Dataset(fn, 'w', format='NETCDF4')

    # create nc. dimensions:
    longitude    = Longitude
    latitude     = Latitude
    ERA5_levels  = ds.createDimension('ERA5_lev', Height.shape[2])
    ERA5_levels_boundaries = ds.createDimension('ERA5lev_boundaries', Height_Boundaries.shape[2])
    CALIPSO_lev  = ds.createDimension('CALIPSO_lev', len(LIVAS_Altitude))
    lat          = ds.createDimension('lat',         len(Latitude))
    lon          = ds.createDimension('lon',         len(Longitude))

    Geolocation_group     = ds.createGroup("Geolocation")
    ERA5_group            = ds.createGroup("ERA5")
    LIVAS_group           = ds.createGroup("LIVAS")
    Land_Ocean_Mask_group = ds.createGroup("Land_Ocean_Mask")
    Dust_Deposition_group = ds.createGroup("Dust_Deposition_Rate")

    # create nc. variables:
    lats_id                   = ds.createVariable('Geolocation/Latitude', 'f4', ('lat',),zlib=True)
    lons_id                   = ds.createVariable('Geolocation/Longitude','f4', ('lon',),zlib=True)

    ERA_Height_id             = ds.createVariable('ERA5/ERA5_Height', 'f4', ('lon','lat','ERA5_lev',),zlib=True)
    Height_Boundaries_id      = ds.createVariable('ERA5/Height_Boundaries', 'f4',       ('lon','lat','ERA5lev_boundaries',),zlib=True)
    ERA_U_id                  = ds.createVariable('ERA5/U',    np.float64, ('lon','lat','ERA5_lev',),zlib=True)
    ERA_U_SD_id               = ds.createVariable('ERA5/U_SD', np.float64, ('lon','lat','ERA5_lev',),zlib=True)
    ERA_V_id                  = ds.createVariable('ERA5/V',    np.float64, ('lon','lat','ERA5_lev',),zlib=True)
    ERA_V_SD_id               = ds.createVariable('ERA5/V_SD', np.float64, ('lon','lat','ERA5_lev',),zlib=True)


    LIVAS_Altitude_id         = ds.createVariable('LIVAS/Altitude',                                        'f4',        ('CALIPSO_lev',),zlib=True)
    LIVAS_a532nm_PD_id        = ds.createVariable('LIVAS/Pure_Dust/LIVAS_PD_a532nm',           np.float64,  ('lon','lat','CALIPSO_lev',),zlib=True)
    LIVAS_a532nm_PD_SD_id     = ds.createVariable('LIVAS/Pure_Dust/LIVAS_PD_a532nm_STD',       np.float64,  ('lon','lat','CALIPSO_lev',),zlib=True)
    LIVAS_PD_MC_id            = ds.createVariable('LIVAS/Pure_Dust/LIVAS_PD_MC',               np.float64,  ('lon','lat','CALIPSO_lev',),zlib=True)
    LIVAS_PD_MC_SD_id         = ds.createVariable('LIVAS/Pure_Dust/LIVAS_PD_MC_STD',           np.float64,  ('lon','lat','CALIPSO_lev',),zlib=True)
    LIVAS_N_of_CF_Profiles_id = ds.createVariable('LIVAS/Flags/Number_of_L2_CF_Profiles',      np.float64,  ('lon','lat',)          ,zlib=True)
    LIVAS_N_of_Profiles_id    = ds.createVariable('LIVAS/Flags/Number_of_L2_Profiles',         np.float64,  ('lon','lat',)          ,zlib=True)
    LIVAS_DOD_532nm_mean_id   = ds.createVariable('LIVAS/Pure_Dust/DOD_532nm_mean',            np.float64,  ('lon','lat',)          ,zlib=True)
    LIVAS_DOD_532nm_SD_id     = ds.createVariable('LIVAS/Pure_Dust/DOD_532nm_STD',             np.float64,  ('lon','lat',)          ,zlib=True)
    LIVAS_PD_MC_ERA5res_id    = ds.createVariable('LIVAS/Pure_Dust/LIVAS_PD_MC_ERA5res',       np.float64, ('lon','lat','ERA5_lev',),zlib=True)
    LIVAS_PD_MC_ERA5res_U_id  = ds.createVariable('LIVAS/Pure_Dust/LIVAS_PD_MC_ERA5res_U',     np.float64, ('lon','lat','ERA5_lev',),zlib=True)
    LIVAS_PD_MC_ERA5res_V_id  = ds.createVariable('LIVAS/Pure_Dust/LIVAS_PD_MC_ERA5res_V',     np.float64, ('lon','lat','ERA5_lev',),zlib=True)

    Land_Ocean_Percentage_id  = ds.createVariable('Land_Ocean_Mask/Land_Ocean_Percentage',     np.float64, ('lon','lat',),zlib=True)

#    Percentage_IGBP_1_id      = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_1',  np.float64, ('lon','lat',), zlib=True)
#    Percentage_IGBP_2_id      = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_2',  np.float64, ('lon','lat',), zlib=True)
#    Percentage_IGBP_3_id      = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_3',  np.float64, ('lon','lat',), zlib=True)
#    Percentage_IGBP_4_id      = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_4',  np.float64, ('lon','lat',), zlib=True)
#    Percentage_IGBP_5_id      = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_5',  np.float64, ('lon','lat',), zlib=True)
#    Percentage_IGBP_6_id      = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_6',  np.float64, ('lon','lat',), zlib=True)
#    Percentage_IGBP_7_id      = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_7',  np.float64, ('lon','lat',), zlib=True)
#    Percentage_IGBP_8_id      = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_8',  np.float64, ('lon','lat',), zlib=True)
#    Percentage_IGBP_9_id      = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_9',  np.float64, ('lon','lat',), zlib=True)
#    Percentage_IGBP_10_id     = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_10', np.float64, ('lon','lat',), zlib=True)
#    Percentage_IGBP_11_id     = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_11', np.float64, ('lon','lat',), zlib=True)
#    Percentage_IGBP_12_id     = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_12', np.float64, ('lon','lat',), zlib=True)
#    Percentage_IGBP_13_id     = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_13', np.float64, ('lon','lat',), zlib=True)
#    Percentage_IGBP_14_id     = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_14', np.float64, ('lon','lat',), zlib=True)
#    Percentage_IGBP_15_id     = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_15', np.float64, ('lon','lat',), zlib=True)
#    Percentage_IGBP_16_id     = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_16', np.float64, ('lon','lat',), zlib=True)
#    Percentage_IGBP_17_id     = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_17', np.float64, ('lon','lat',), zlib=True)
#    Percentage_IGBP_18_id     = ds.createVariable('Land_Ocean_Mask/Percentage_IGBP_18', np.float64, ('lon','lat',), zlib=True)

    LIVAS_Dust_Deposition_id  = ds.createVariable('Dust_Deposition_Rate/LIVAS_Dust_Deposition', np.float64, ('lon','lat',), zlib=True)
    LIVAS_F_N_id              = ds.createVariable('Dust_Deposition_Rate/LIVAS_MFR_N', np.float64, ('lon','lat','ERA5_lev',),zlib=True)
    LIVAS_F_S_id              = ds.createVariable('Dust_Deposition_Rate/LIVAS_MFR_S', np.float64, ('lon','lat','ERA5_lev',),zlib=True)
    LIVAS_F_E_id              = ds.createVariable('Dust_Deposition_Rate/LIVAS_MFR_E', np.float64, ('lon','lat','ERA5_lev',),zlib=True)
    LIVAS_F_W_id              = ds.createVariable('Dust_Deposition_Rate/LIVAS_MFR_W', np.float64, ('lon','lat','ERA5_lev',),zlib=True)

    lats_id.units                   = 'degrees_north'
    lons_id.units                   = 'degrees_east'
    ERA_Height_id.units             = 'm'
    Height_Boundaries_id.units      = 'm'
    ERA_U_id.units                  = 'm s**-1'
    ERA_U_SD_id.units               = 'm s**-1'
    ERA_V_id.units                  = 'm s**-1'
    ERA_V_SD_id.units               = 'm s**-1'
    LIVAS_Altitude_id.units         = 'm'
    LIVAS_a532nm_PD_id.units        = 'km-1'
    LIVAS_a532nm_PD_SD_id.units     = 'km-1'
    LIVAS_PD_MC_id.units            = 'micrograms/m^3'
    LIVAS_PD_MC_SD_id.units         = 'micrograms/m^3'
    LIVAS_PD_MC_ERA5res_id.units    = 'micrograms/m^3'
    LIVAS_PD_MC_ERA5res_U_id.units  = 'micrograms/m^3'
    LIVAS_PD_MC_ERA5res_V_id.units  = 'micrograms/m^3'
    LIVAS_N_of_CF_Profiles_id.units = 'none'
    LIVAS_N_of_Profiles_id.units    = 'none'
    LIVAS_DOD_532nm_mean_id.units   = 'none'
    LIVAS_DOD_532nm_SD_id.units     = 'none'
    LIVAS_Dust_Deposition_id.units  = 'mg/m^2d'
    LIVAS_F_N_id.units              = 'micrograms/s'
    LIVAS_F_S_id.units              = 'micrograms/s'
    LIVAS_F_E_id.units              = 'micrograms/s'
    LIVAS_F_W_id.units              = 'micrograms/s'

    Land_Ocean_Percentage_id.units   = 'none'
#    Percentage_IGBP_1_id.units      = 'none'
#    Percentage_IGBP_2_id.units      = 'none'
#    Percentage_IGBP_3_id.units      = 'none'
#    Percentage_IGBP_4_id.units      = 'none'
#    Percentage_IGBP_5_id.units      = 'none'
#    Percentage_IGBP_6_id.units      = 'none'
#    Percentage_IGBP_7_id.units      = 'none'
#    Percentage_IGBP_8_id.units      = 'none'
#    Percentage_IGBP_9_id.units      = 'none'
#    Percentage_IGBP_10_id.units     = 'none'
#    Percentage_IGBP_11_id.units     = 'none'
#    Percentage_IGBP_12_id.units     = 'none'
#    Percentage_IGBP_13_id.units     = 'none'
#    Percentage_IGBP_14_id.units     = 'none'
#    Percentage_IGBP_15_id.units     = 'none'
#    Percentage_IGBP_16_id.units     = 'none'
#    Percentage_IGBP_17_id.units     = 'none'
#    Percentage_IGBP_18_id.units     = 'none'

    lats_id.long_name                   = 'Latitude'
    lons_id.long_name                   = 'Longitude'
    ERA_Height_id.long_name             = 'Height'
    Height_Boundaries_id.long_name      = 'ERA5 level boundaries'
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
    LIVAS_DOD_532nm_mean_id.long_name   = 'Dust Optical Depth 532nm - mean'
    LIVAS_DOD_532nm_SD_id.long_name     = 'Dust Optical Depth 532nm - SD'
    LIVAS_PD_MC_ERA5res_id.long_name    = 'LIVAS Pure-Dust Mass Concentration profiles - ERA5 resolution'
    LIVAS_PD_MC_ERA5res_U_id.long_name  = 'LIVAS Pure-Dust Mass Concentration profiles - ERA5 resolution - U'
    LIVAS_PD_MC_ERA5res_V_id.long_name  = 'LIVAS Pure-Dust Mass Concentration profiles - ERA5 resolution - V'
    Land_Ocean_Percentage_id.long_name  = 'Land Ocean Mask Percentage'
#    Percentage_IGBP_1_id.long_name      = 'Evergreen-Needleleaf-Forest'
#    Percentage_IGBP_2_id.long_name      = 'Evergreen-Broadleaf-Forest'
#    Percentage_IGBP_3_id.long_name      = 'Deciduous-Needleleaf-Forest'
#    Percentage_IGBP_4_id.long_name      = 'Deciduous-Broadleaf-Forest'
#    Percentage_IGBP_5_id.long_name      = 'Mixed-Forest'
#    Percentage_IGBP_6_id.long_name      = 'Closed-Shrublands'
#    Percentage_IGBP_7_id.long_name      = 'Open-Shrubland (Desert)'
#    Percentage_IGBP_8_id.long_name      = 'Woody-Savanna'
#    Percentage_IGBP_9_id.long_name      = 'Savanna'
#    Percentage_IGBP_10_id.long_name     = 'Grassland'
#    Percentage_IGBP_11_id.long_name     = 'Wetland'
#    Percentage_IGBP_12_id.long_name     = 'Cropland'
#    Percentage_IGBP_13_id.long_name     = 'Urban'
#    Percentage_IGBP_14_id.long_name     = 'Crop-Mosaic'
#    Percentage_IGBP_15_id.long_name     = 'Permanent-Snow'
#    Percentage_IGBP_16_id.long_name     = 'Barren/Desert'
#    Percentage_IGBP_17_id.long_name     = 'Water'
#    Percentage_IGBP_18_id.long_name     = 'Tundra'
    LIVAS_F_N_id.long_name              = 'Mass Flow Rate - North'
    LIVAS_F_S_id.long_name              = 'Mass Flow Rate - South'
    LIVAS_F_E_id.long_name              = 'Mass Flow Rate - East'
    LIVAS_F_W_id.long_name              = 'Mass Flow Rate - West'


    lats_id.fill_value                   = np.nan
    lons_id.fill_value                   = np.nan
    ERA_Height_id.fill_value             = np.nan
    Height_Boundaries_id.fill_value      = np.nan
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
    LIVAS_DOD_532nm_mean_id.fill_value   = np.nan
    LIVAS_DOD_532nm_SD_id.fill_value     = np.nan
    LIVAS_PD_MC_ERA5res_id.fill_value    = np.nan
    LIVAS_PD_MC_ERA5res_U_id.fill_value  = np.nan
    LIVAS_PD_MC_ERA5res_V_id.fill_value  = np.nan
    Land_Ocean_Percentage_id.fill_value  = np.nan
#    Percentage_IGBP_1_id.fill_value      = np.nan
#    Percentage_IGBP_2_id.fill_value      = np.nan
#    Percentage_IGBP_3_id.fill_value      = np.nan
#    Percentage_IGBP_4_id.fill_value      = np.nan
#    Percentage_IGBP_5_id.fill_value      = np.nan
#    Percentage_IGBP_6_id.fill_value      = np.nan
#    Percentage_IGBP_7_id.fill_value      = np.nan
#    Percentage_IGBP_8_id.fill_value      = np.nan
#    Percentage_IGBP_9_id.fill_value      = np.nan
#    Percentage_IGBP_10_id.fill_value     = np.nan
#    Percentage_IGBP_11_id.fill_value     = np.nan
#    Percentage_IGBP_12_id.fill_value     = np.nan
#    Percentage_IGBP_13_id.fill_value     = np.nan
#    Percentage_IGBP_14_id.fill_value     = np.nan
#    Percentage_IGBP_15_id.fill_value     = np.nan
#    Percentage_IGBP_16_id.fill_value     = np.nan
#    Percentage_IGBP_17_id.fill_value     = np.nan
#    Percentage_IGBP_18_id.fill_value     = np.nan
    LIVAS_Dust_Deposition_id.fill_value  = np.nan
    LIVAS_F_N_id.fill_value              = np.nan
    LIVAS_F_S_id.fill_value              = np.nan
    LIVAS_F_E_id.fill_value              = np.nan
    LIVAS_F_W_id.fill_value              = np.nan


    lats_id[:]                   = Latitude
    lons_id[:]                   = Longitude
    ERA_Height_id[:]             = Height
    Height_Boundaries_id[:]      = Height_Boundaries
    ERA_U_id[:]                  = U
    ERA_U_SD_id[:]               = U_SD
    ERA_V_id[:]                  = V
    ERA_V_SD_id[:]               = V_SD
    LIVAS_Altitude_id[:]         = LIVAS_Altitude
    LIVAS_a532nm_PD_id[:]        = LIVAS_a532nm_PD
    LIVAS_a532nm_PD_SD_id[:]     = LIVAS_a532nm_PD
    LIVAS_PD_MC_id[:]            = LIVAS_PD_MC
    LIVAS_PD_MC_SD_id[:]         = LIVAS_PD_MC_SD
    LIVAS_N_of_CF_Profiles_id[:] = LIVAS_N_of_CF_Profiles
    LIVAS_N_of_Profiles_id[:]    = LIVAS_N_of_Profiles
    LIVAS_DOD_532nm_mean_id[:]   = LIVAS_DOD_532nm_mean
    LIVAS_DOD_532nm_SD_id[:]     = LIVAS_DOD_532nm_SD
    LIVAS_PD_MC_ERA5res_id[:]    = PD_MC_ERA5res
    LIVAS_PD_MC_ERA5res_U_id[:]  = LIVAS_PD_MC_ERA5res_U
    LIVAS_PD_MC_ERA5res_V_id[:]  = LIVAS_PD_MC_ERA5res_V
    Land_Ocean_Percentage_id[:]  = Land_Ocean_Percentage
#    Percentage_IGBP_1_id[:]      = Percentage_IGBP_1
#    Percentage_IGBP_2_id[:]      = Percentage_IGBP_2
#    Percentage_IGBP_3_id[:]      = Percentage_IGBP_3
#    Percentage_IGBP_4_id[:]      = Percentage_IGBP_4
#    Percentage_IGBP_5_id[:]      = Percentage_IGBP_5
#    Percentage_IGBP_6_id[:]      = Percentage_IGBP_6
#    Percentage_IGBP_7_id[:]      = Percentage_IGBP_7
#    Percentage_IGBP_8_id[:]      = Percentage_IGBP_8
#    Percentage_IGBP_9_id[:]      = Percentage_IGBP_9
#    Percentage_IGBP_10_id[:]     = Percentage_IGBP_10
#    Percentage_IGBP_11_id[:]     = Percentage_IGBP_11
#    Percentage_IGBP_12_id[:]     = Percentage_IGBP_12
#    Percentage_IGBP_13_id[:]     = Percentage_IGBP_13
#    Percentage_IGBP_14_id[:]     = Percentage_IGBP_14
#    Percentage_IGBP_15_id[:]     = Percentage_IGBP_15
#    Percentage_IGBP_16_id[:]     = Percentage_IGBP_16
#    Percentage_IGBP_17_id[:]     = Percentage_IGBP_17
#    Percentage_IGBP_18_id[:]     = Percentage_IGBP_18
    LIVAS_Dust_Deposition_id[:]  = Dust_Deposition_Horizontally_smoothed

    ds.close()

    print("End of: " + fn)

#  SCRIPT END  ---------------------------------------------------------------
Ou.goodbye(cnf.LOGs.run, tic = tic, scriptname = __file__, version = TRACE)
