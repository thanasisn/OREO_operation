# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 08:54:57 2022
Establishes the required EO-and-ERA5 DOMOS dataset in the same L3 1x1 grid resolution and monthly-mean.

@author: proestakis, thanasisn
"""

import os
import sys
import re

import netCDF4  as nc
import numpy    as np
import os
from   os.path import exists
from   datetime    import datetime, timedelta
from   os          import walk
import glob
import math
import warnings
import metpy
import metpy.calc
from   metpy.units import units
from   pathlib import Path
import numpy.ma as ma
import xarray as xr

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
# MONTHLY  = cnf.D1.Monthly

## export be season
SEASONAL = False
SEASONAL = cnf.D1.Seasonal


##  Allow only one case to run at the time  ----------------------------------
if SEASONAL == MONTHLY:
    print("Seasonal:", SEASONAL)
    print("Monthly: ", MONTHLY)
    sys.exit("Choose only SEASONAL or only MONTHLY")
    

##  Check destination folder exists  -----------------------------------------
if not os.path.isdir(cnf.OREO.path_output):
    sys.exit(f"\nFolder {cnf.ERA5.path_regrid} don't exist !!\n")


##  Temporal aggregation setup  ----------------------------------------------
if SEASONAL == True:
    print("Work on seasonal data")
    
    fERA_ilenames = glob.glob(
        f"{cnf.ERA5.path_regrid}/Seasonal_{cnf.D1.LatStep}x{cnf.D1.LonStep}/ERA5_*_{cnf.D1.North}N{cnf.D1.South}S{cnf.D1.West}W{cnf.D1.East}E.nc"
    )
elif MONTHLY == True:
    print("Work on monthly data")

    ERA_filenames = glob.glob(
        f"{cnf.ERA5.path_regrid}/Monthly_{cnf.D1.LatStep}x{cnf.D1.LonStep}/ERA5_*_{cnf.D1.North}N{cnf.D1.South}S{cnf.D1.West}W{cnf.D1.East}E.nc"
    )

filenames.sort()

if cnf.ERA5.data == "mean":
    U = "u_mean"
    V = "v_mean"
elif cnf.ERA5.data == "median":
    U = "u_median"
    V = "v_median"



##  Directory paths of input datasets
# ERA_path       = r"M:\DOMOS\datasets\2x5_seasonal\ERA5\processed"
LIVAS_path     = cnf.LIVAS.path_lookup
# CAMS_path      = r"M:\DOMOS\datasets\2x5_seasonal\CAMS\processed"

##  Directory path of output datasets
output_path    = os.path.join(cnf.OREO.path_output,
                              os.path.basename(os.path.dirname(filenames[0])))
os.makedirs(output_path, exist_ok = True)


for ERA_file in filenames:

    ##  Load ERA5 data  ------------------------------------------------------
    
    print(f"\nProcessing: {ERA_file}")

    # reading variables of interest.    
    ERA = xr.open_dataset(ERA_file)
    ERA_dataset   = nc.Dataset(ERA_file) 
    
    ERA_Latitude  = ERA_dataset['latitude'][:]
    ERA_Longitude = ERA_dataset['longitude'][:]      
    ERA_Height    = ERA_dataset['height'][:]
    ERA_Latitude  = ERA_dataset['latitude'][:]
    ERA_Longitude = ERA_dataset['longitude'][:]
    ERA_U         = ERA_dataset[ U ][:]
    ERA_U_SD      = ERA_dataset['u_SD'][:]
    ERA_V         = ERA_dataset[ V ][:]
    ERA_V_SD      = ERA_dataset['u_SD'][:]
    
    
    sys.exit("wait")
    
    os.path.basename(ERA_file)
    
    yyyy = int(re.compile('ERA5_([0-9]*)_.*.nc').search(ERA_file).group(1))
    
    
    # extracting "yyyymm" sufix from ERA5 filename, for finding the satellite-based MM files.  
    ERA_year          = int(re.compile('ERA5_([0-9]*)_.*.nc').search(ERA_file).group(1))
    ERA_year_previous = ERA_year - 1
    seas = (re.compile('ERA5_[0-9]*_Q[1-4]_([A-Z]*)_.*.nc').search(ERA_file).group(1))

    ERA_season        = seas
    
    
    if ERA_season == 'DJF':          
        YoI = [ERA_year_previous,ERA_year,ERA_year]
        MoI = [12,1,2]
    if ERA_season == 'MAM':    
        YoI = [ERA_year,ERA_year,ERA_year]
        MoI = [3,4,5]        
    if ERA_season == 'JJA':    
        YoI = [ERA_year,ERA_year,ERA_year]
        MoI = [6,7,8]        
    if ERA_season == 'SON':    
        YoI = [ERA_year,ERA_year,ERA_year]
        MoI = [9,10,11]        

    fn = output_path + '\\' + 'DOMOS_Datasets_' + ERA_substrings[0] + '_' + ERA_season + '.nc' 
    if exists(fn):
        continue

    # ######################
    # ###       CAMS     ###
    # ######################    
    # 
    # CAMS_file      = CAMS_path + '\\' + 'CAMS_' + ERA_substrings[0] + '_' + ERA_season + '.nc' 
    # CAMS_dataset   = nc.Dataset(CAMS_file) 
    # CAMS_latitude  = CAMS_dataset['Latitude'][:]
    # CAMS_longitude = CAMS_dataset['Longitude'][:]      
    # CAMS_Height    = CAMS_dataset['Height'][:]
    # CAMS_Dust_MC   = CAMS_dataset['CAMS_Dust_MC'][:]   

    #######################
    ### LIVAS pure-dust ###
    #######################    

    Percentage_IGBP_1          = np.empty((len(ERA_Longitude),len(ERA_Latitude)))
    Percentage_IGBP_2          = np.empty((len(ERA_Longitude),len(ERA_Latitude)))
    Percentage_IGBP_3          = np.empty((len(ERA_Longitude),len(ERA_Latitude)))
    Percentage_IGBP_4          = np.empty((len(ERA_Longitude),len(ERA_Latitude)))
    Percentage_IGBP_5          = np.empty((len(ERA_Longitude),len(ERA_Latitude)))
    Percentage_IGBP_6          = np.empty((len(ERA_Longitude),len(ERA_Latitude)))
    Percentage_IGBP_7          = np.empty((len(ERA_Longitude),len(ERA_Latitude)))
    Percentage_IGBP_8          = np.empty((len(ERA_Longitude),len(ERA_Latitude)))
    Percentage_IGBP_9          = np.empty((len(ERA_Longitude),len(ERA_Latitude)))
    Percentage_IGBP_10         = np.empty((len(ERA_Longitude),len(ERA_Latitude)))
    Percentage_IGBP_11         = np.empty((len(ERA_Longitude),len(ERA_Latitude)))
    Percentage_IGBP_12         = np.empty((len(ERA_Longitude),len(ERA_Latitude)))
    Percentage_IGBP_13         = np.empty((len(ERA_Longitude),len(ERA_Latitude)))
    Percentage_IGBP_14         = np.empty((len(ERA_Longitude),len(ERA_Latitude)))
    Percentage_IGBP_15         = np.empty((len(ERA_Longitude),len(ERA_Latitude)))
    Percentage_IGBP_16         = np.empty((len(ERA_Longitude),len(ERA_Latitude)))
    Percentage_IGBP_17         = np.empty((len(ERA_Longitude),len(ERA_Latitude)))
    Percentage_IGBP_18         = np.empty((len(ERA_Longitude),len(ERA_Latitude)))
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
       
    Final_Number_of_Profiles      = np.empty(( len(ERA_Longitude), len(ERA_Latitude)))
    Final_Number_of_L2Profiles    = np.empty(( len(ERA_Longitude), len(ERA_Latitude)))  
    Final_LIVAS_PD_DOD_532nm      = np.empty(( len(ERA_Longitude), len(ERA_Latitude)))
    Final_LIVAS_PD_DOD_532nm_SD   = np.empty(( len(ERA_Longitude), len(ERA_Latitude)))      
    Final_PD_a532nm   = np.empty(( len(ERA_Longitude), len(ERA_Latitude), 399))
    Final_PD_a532nm_SD= np.empty(( len(ERA_Longitude), len(ERA_Latitude), 399))
    Final_PD_MC       = np.empty(( len(ERA_Longitude), len(ERA_Latitude), 399))
    Final_PD_MC_FM    = np.empty(( len(ERA_Longitude), len(ERA_Latitude), 399))
    Final_PD_MC_CM    = np.empty(( len(ERA_Longitude), len(ERA_Latitude), 399))
    Final_PD_MC_SD    = np.empty(( len(ERA_Longitude), len(ERA_Latitude), 399))
    Final_PD_MC_FM_SD = np.empty(( len(ERA_Longitude), len(ERA_Latitude), 399))
    Final_PD_MC_CM_SD = np.empty(( len(ERA_Longitude), len(ERA_Latitude), 399))  

    Final_Number_of_Profiles[:]    = np.nan
    Final_Number_of_L2Profiles[:]  = np.nan
    Final_LIVAS_PD_DOD_532nm[:]    = np.nan
    Final_LIVAS_PD_DOD_532nm_SD[:] = np.nan      
    Final_PD_a532nm[:]             = np.nan
    Final_PD_a532nm_SD[:]          = np.nan  
    Final_PD_MC[:]                 = np.nan
  #  Final_PD_MC_FM[:]              = np.nan # fine mode no need
  #  Final_PD_MC_CM[:]              = np.nan
  #  Final_PD_MC_SD[:]              = np.nan
  #  Final_PD_MC_FM_SD[:]           = np.nan
  #  Final_PD_MC_CM_SD[:]           = np.nan
   
    Empty_Vertical_array    = np.empty(399)
    Empty_Vertical_array[:] = np.nan  

    for count_lon,lon in enumerate(ERA_Longitude):        
        for count_lat,lat in enumerate(ERA_Latitude):
 
            LIVAS_lons = [lon-2,lon-1,lon,lon+1,lon+2,lon-2,lon-1,lon,lon+1,lon+2]
            LIVAS_lats = [lat-0.5,lat+0.5,lat-0.5,lat+0.5,lat-0.5,lat+0.5,lat-0.5,lat+0.5,lat-0.5,lat+0.5]
            
            file_counter = 0
            
            idx_lat = count_lat
            idx_lon = count_lon  

            ######################
            ###     LIVAS      ###
            ###################### 
            
            for LIVAS_lons_lats_idx,LIVAS_lons_lats in enumerate(LIVAS_lons):
           
                LIVAS_filename  = 'LIVAS_CALIPSO_L2_Grid_lon_c_' + str(LIVAS_lons[LIVAS_lons_lats_idx]) + '_lat_c_' + str(LIVAS_lats[LIVAS_lons_lats_idx]) + '.nc'
                LIVAS_file      = LIVAS_path + '\\' + LIVAS_filename
                file_existing   = exists(LIVAS_file)
                
                if (not file_existing) | (LIVAS_filename == 'LIVAS_CALIPSO_L2_Grid_lon_c_-106.5_lat_c_-55.5.nc') | (LIVAS_filename == 'LIVAS_CALIPSO_L2_Grid_lon_c_-106.5_lat_c_-48.5.nc') :
                    continue
                else:  
                    
                    LIVAS_dataset       = nc.Dataset(LIVAS_file)
                    Profile_Time_Parsed = LIVAS_dataset['/Profile_Time_Parsed'][:]
                    Months    = np.empty((len(Profile_Time_Parsed)))
                    Years     = np.empty((len(Profile_Time_Parsed)))
                    Months[:] = np.nan
                    Years[:]  = np.nan
                    for count_time,time in enumerate(Profile_Time_Parsed):
                        time = time.split(' ')[0]
                        Months[count_time] = int(time.split('/')[1])
                        Years[count_time]  = int(time.split('/')[0])
                    idx = np.where( (MoI[0] == Months) & (YoI[0] == Years) | (MoI[1] == Months) & (YoI[1] == Years) | (MoI[2] == Months) & (YoI[2] == Years) )
                    idx = np.ravel(idx)
                    if (len(idx) == 0) | (len(idx) == 1):
                        continue
                    Month = Months[idx]
                    Year  = Years[idx]

                    Altitude        = LIVAS_dataset['/Altitude'][:]
                    LIVAS_PD_b532nm = LIVAS_dataset['/LIVAS/Cloud_Free/Pure_Dust_and_Fine_Coarse/Optical_Products/Pure_Dust_Backscatter_Coefficient_532'][idx,:]
                    LIVAS_PD_a532nm = LIVAS_dataset['/LIVAS/Cloud_Free/Pure_Dust_and_Fine_Coarse/Optical_Products/Pure_Dust_Extinction_Coefficient_532'][idx,:] 
                    LIVAS_PD_MC     = LIVAS_dataset['/LIVAS/Cloud_Free/Pure_Dust_and_Fine_Coarse/Mass_Concentrations/Pure_Dust_Mass_Concentration'][idx,:]   
                    LIVAS_PD_MC_FM  = LIVAS_dataset['/LIVAS/Cloud_Free/Pure_Dust_and_Fine_Coarse/Mass_Concentrations/Pure_Dust_Fine_Mass_Concentration'][idx,:]
                    LIVAS_PD_MC_CM  = LIVAS_dataset['/LIVAS/Cloud_Free/Pure_Dust_and_Fine_Coarse/Mass_Concentrations/Pure_Dust_Coarse_Mass_Concentration'][idx,:]

                    IGBP            = LIVAS_dataset['/CALIPSO_Flags_and_Auxiliary/Auxiliary/IGBP_Surface_Type'][:]
                                       
                    idx = np.where(LIVAS_PD_b532nm == 0) ## to check input
                    LIVAS_PD_a532nm[idx] = 0
                    LIVAS_PD_MC[idx]     = 0
                    LIVAS_PD_MC_FM[idx]  = 0
                    LIVAS_PD_MC_CM[idx]  = 0
                    
                    if ma.isMaskedArray(LIVAS_PD_b532nm) == True:
                        LIVAS_PD_a532nm[LIVAS_PD_b532nm.mask == True ] = np.nan
                        LIVAS_PD_MC[LIVAS_PD_b532nm.mask == True ]     = np.nan
                        LIVAS_PD_MC_FM[LIVAS_PD_b532nm.mask == True ]  = np.nan
                        LIVAS_PD_MC_CM[LIVAS_PD_b532nm.mask == True ]  = np.nan
                    
                    if file_counter >  0:
                        Total_LIVAS_PD_a532nm = np.vstack([LIVAS_PD_a532nm,Total_LIVAS_PD_a532nm])
                        Total_LIVAS_PD_MC     = np.vstack([LIVAS_PD_MC,Total_LIVAS_PD_MC])
                        Total_LIVAS_PD_MC_FM  = np.vstack([LIVAS_PD_MC_FM,Total_LIVAS_PD_MC_FM])
                        Total_LIVAS_PD_MC_CM  = np.vstack([LIVAS_PD_MC_CM,Total_LIVAS_PD_MC_CM])
                        Total_IGBP            = np.hstack([Total_IGBP,IGBP])
                        file_counter = file_counter + 1                          
                    if file_counter == 0:
                        Total_LIVAS_PD_a532nm = LIVAS_PD_a532nm
                        Total_LIVAS_PD_MC     = LIVAS_PD_MC
                        Total_LIVAS_PD_MC_FM  = LIVAS_PD_MC_FM
                        Total_LIVAS_PD_MC_CM  = LIVAS_PD_MC_CM  
                        Total_IGBP            = IGBP
                        file_counter = file_counter + 1


                        ## to continut
                        
            Number_of_Profiles = np.shape(Total_LIVAS_PD_MC)[0]
            temp = np.copy(Total_LIVAS_PD_MC)
            idx  = ~np.isnan(Total_LIVAS_PD_MC)
            temp[idx] = 0
            idx  = np.isnan(Total_LIVAS_PD_MC)
            temp[idx] = 1                        
            temp = np.ravel([np.nansum(temp[i,:]) for i in range(np.shape(temp)[0]) ])
            L2_CF_profiles  = len(temp) - len(np.ravel(np.where(temp == 399)))
            
            PD_a532nm    = np.nanmean(Total_LIVAS_PD_a532nm,axis = 0)
            PD_MC        = np.nanmean(Total_LIVAS_PD_MC,    axis = 0)
            PD_MC_FM     = np.nanmean(Total_LIVAS_PD_MC_FM, axis = 0)
            PD_MC_CM     = np.nanmean(Total_LIVAS_PD_MC_CM, axis = 0)            
            PD_MC_SD     = np.nanstd(Total_LIVAS_PD_MC,     axis = 0, ddof = 1)
            PD_MC_FM_SD  = np.nanstd(Total_LIVAS_PD_MC_FM,  axis = 0, ddof = 1)
            PD_MC_CM_SD  = np.nanstd(Total_LIVAS_PD_MC_CM,  axis = 0, ddof = 1)             
            PD_a532nm_SD = np.nanstd(Total_LIVAS_PD_a532nm, axis = 0, ddof = 1)  

            PD_a532nm[Altitude > 10]    = np.nan
            PD_MC[Altitude > 10]        = np.nan
            PD_MC_FM[Altitude > 10]     = np.nan
            PD_MC_CM[Altitude > 10]     = np.nan          
            PD_MC_SD[Altitude > 10]     = np.nan
            PD_MC_FM_SD[Altitude > 10]  = np.nan
            PD_MC_CM_SD[Altitude > 10]  = np.nan            
            PD_a532nm_SD[Altitude > 10] = np.nan

            arr = np.copy(PD_a532nm)
            arr[np.isnan(arr)] = 0                            
            DOD_532nm       = np.trapz(Altitude,arr) 
            arr = np.copy(PD_a532nm_SD)
            arr[np.isnan(arr)] = 0                            
            DOD_532nm_SD    = np.trapz(Altitude,arr) 
            
            for count_alt in range(399):   
                Final_PD_MC[idx_lon,idx_lat,count_alt]       = PD_MC[count_alt]
                Final_PD_MC_FM[idx_lon,idx_lat,count_alt]    = PD_MC_FM[count_alt]
                Final_PD_MC_CM[idx_lon,idx_lat,count_alt]    = PD_MC_CM[count_alt]
                Final_PD_MC_SD[idx_lon,idx_lat,count_alt]    = PD_MC_SD[count_alt]
                Final_PD_MC_FM_SD[idx_lon,idx_lat,count_alt] = PD_MC_FM_SD[count_alt]
                Final_PD_MC_CM_SD[idx_lon,idx_lat,count_alt] = PD_MC_CM_SD[count_alt]
                Final_PD_a532nm                              = PD_a532nm[count_alt]
                Final_PD_a532nm_SD                           = PD_a532nm_SD[count_alt]
            
            Final_LIVAS_PD_DOD_532nm[idx_lon,idx_lat]     = DOD_532nm
            Final_LIVAS_PD_DOD_532nm_SD[idx_lon,idx_lat]  = DOD_532nm_SD
            Percentage_IGBP_1[idx_lon,idx_lat]  = (np.divide(float(np.count_nonzero(IGBP ==  1)),float(len(IGBP))))*100.0
            Percentage_IGBP_2[idx_lon,idx_lat]  = (np.divide(float(np.count_nonzero(IGBP ==  2)),float(len(IGBP))))*100.0
            Percentage_IGBP_3[idx_lon,idx_lat]  = (np.divide(float(np.count_nonzero(IGBP ==  3)),float(len(IGBP))))*100.0
            Percentage_IGBP_4[idx_lon,idx_lat]  = (np.divide(float(np.count_nonzero(IGBP ==  4)),float(len(IGBP))))*100.0
            Percentage_IGBP_5[idx_lon,idx_lat]  = (np.divide(float(np.count_nonzero(IGBP ==  5)),float(len(IGBP))))*100.0
            Percentage_IGBP_6[idx_lon,idx_lat]  = (np.divide(float(np.count_nonzero(IGBP ==  6)),float(len(IGBP))))*100.0
            Percentage_IGBP_7[idx_lon,idx_lat]  = (np.divide(float(np.count_nonzero(IGBP ==  7)),float(len(IGBP))))*100.0
            Percentage_IGBP_8[idx_lon,idx_lat]  = (np.divide(float(np.count_nonzero(IGBP ==  8)),float(len(IGBP))))*100.0
            Percentage_IGBP_9[idx_lon,idx_lat]  = (np.divide(float(np.count_nonzero(IGBP ==  9)),float(len(IGBP))))*100.0
            Percentage_IGBP_10[idx_lon,idx_lat] = (np.divide(float(np.count_nonzero(IGBP == 10)),float(len(IGBP))))*100.0
            Percentage_IGBP_11[idx_lon,idx_lat] = (np.divide(float(np.count_nonzero(IGBP == 11)),float(len(IGBP))))*100.0
            Percentage_IGBP_12[idx_lon,idx_lat] = (np.divide(float(np.count_nonzero(IGBP == 12)),float(len(IGBP))))*100.0
            Percentage_IGBP_13[idx_lon,idx_lat] = (np.divide(float(np.count_nonzero(IGBP == 13)),float(len(IGBP))))*100.0
            Percentage_IGBP_14[idx_lon,idx_lat] = (np.divide(float(np.count_nonzero(IGBP == 14)),float(len(IGBP))))*100.0
            Percentage_IGBP_15[idx_lon,idx_lat] = (np.divide(float(np.count_nonzero(IGBP == 15)),float(len(IGBP))))*100.0
            Percentage_IGBP_16[idx_lon,idx_lat] = (np.divide(float(np.count_nonzero(IGBP == 16)),float(len(IGBP))))*100.0
            Percentage_IGBP_17[idx_lon,idx_lat] = (np.divide(float(np.count_nonzero(IGBP == 17)),float(len(IGBP))))*100.0
            Percentage_IGBP_18[idx_lon,idx_lat] = (np.divide(float(np.count_nonzero(IGBP == 18)),float(len(IGBP))))*100.0

            Final_Number_of_Profiles[idx_lon,idx_lat]   = Number_of_Profiles
            Final_Number_of_L2Profiles[idx_lon,idx_lat] = L2_CF_profiles
            
            print(lon,lat)
            
    #####################################            
    #  --- Saving dataset as NetCDF --- #
    #####################################
    
    Altitude = Altitude*1000.0
    
    # creating nc. filename and initiallizing:                 
    fn           = output_path + '\\' + 'DOMOS_Datasets_' + ERA_substrings[0] + '_' + ERA_season + '.nc' 
    ds           = nc.Dataset(fn, 'w', format='NETCDF4')

    # create nc. dimensions:
    longitude    = ERA_Longitude
    latitude     = ERA_Latitude        
    ERA_lev      = ds.createDimension('ERA_lev',     ERA_Height.shape[2])
    CAMS_lev     = ds.createDimension('CAMS_lev',    CAMS_Height.shape[2])
    CALIPSO_lev  = ds.createDimension('CALIPSO_lev', len(Altitude))
    lat          = ds.createDimension('lat',         len(ERA_Latitude)) 
    lon          = ds.createDimension('lon',         len(ERA_Longitude))
    
    Geolocation_group     = ds.createGroup("Geolocation")
    ERA5_group            = ds.createGroup("ERA5")
    CAMS_group            = ds.createGroup("CAMS")
    LIVAS_group           = ds.createGroup("LIVAS")
    Land_Ocean_Mask_group = ds.createGroup("Land_Ocean_Mask")

    # create nc. variables:
    lats_id                   = ds.createVariable('Geolocation/Latitude', 'f4', ('lat',))
    lons_id                   = ds.createVariable('Geolocation/Longitude','f4', ('lon',))
    
    ERA_Height_id             = ds.createVariable('ERA5/ERA5_Height', 'f4', ('lon','lat','ERA_lev',))
    ERA_U_id                  = ds.createVariable('ERA5/U',    np.float64, ('lon','lat','ERA_lev',))
    ERA_U_SD_id               = ds.createVariable('ERA5/U_SD', np.float64, ('lon','lat','ERA_lev',))
    ERA_V_id                  = ds.createVariable('ERA5/V',    np.float64, ('lon','lat','ERA_lev',))
    ERA_V_SD_id               = ds.createVariable('ERA5/V_SD', np.float64, ('lon','lat','ERA_lev',))      

    CAMS_Height_id            = ds.createVariable('CAMS/CAMS_Height', 'f4',     ('lon','lat','CAMS_lev',))
    CAMS_MC_id                = ds.createVariable('CAMS/CAMS_Dust_MC',    np.float64, ('lon','lat','CAMS_lev',))
 
    LIVAS_Altitude_id         = ds.createVariable('LIVAS/Altitude',                                        'f4',        ('CALIPSO_lev',),zlib=True)
    LIVAS_a532nm_PD_id        = ds.createVariable('LIVAS/Pure_Dust/LIVAS_PD_a532nm',           np.float64,  ('lon','lat','CALIPSO_lev',),zlib=True)
    LIVAS_a532nm_PD_SD_id     = ds.createVariable('LIVAS/Pure_Dust/LIVAS_PD_a532nm_STD',       np.float64,  ('lon','lat','CALIPSO_lev',),zlib=True)        
    LIVAS_PD_MC_id            = ds.createVariable('LIVAS/Pure_Dust/LIVAS_PD_MC',               np.float64,  ('lon','lat','CALIPSO_lev',),zlib=True)
    LIVAS_PD_MC_SD_id         = ds.createVariable('LIVAS/Pure_Dust/LIVAS_PD_MC_STD',           np.float64,  ('lon','lat','CALIPSO_lev',),zlib=True)    
    LIVAS_PD_MC_CM_id         = ds.createVariable('LIVAS/Pure_Dust/LIVAS_PD_MC_CM',            np.float64,  ('lon','lat','CALIPSO_lev',),zlib=True)
    LIVAS_PD_MC_CM_SD_id      = ds.createVariable('LIVAS/Pure_Dust/LIVAS_PD_MC_CM_STD',        np.float64,  ('lon','lat','CALIPSO_lev',),zlib=True)    
    LIVAS_PD_MC_FM_id         = ds.createVariable('LIVAS/Pure_Dust/LIVAS_PD_MC_FM',            np.float64,  ('lon','lat','CALIPSO_lev',),zlib=True)
    LIVAS_PD_MC_FM_SD_id      = ds.createVariable('LIVAS/Pure_Dust/LIVAS_PD_MC_FM_STD',        np.float64,  ('lon','lat','CALIPSO_lev',),zlib=True)        
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
    ERA_Height_id.units             = 'm'
    ERA_U_id.units                  = 'm s**-1' 
    ERA_U_SD_id.units               = 'm s**-1' 
    ERA_V_id.units                  = 'm s**-1' 
    ERA_V_SD_id.units               = 'm s**-1' 
    CAMS_Height_id.units            = 'm'
    CAMS_MC_id.units                = 'micrograms/m^3'
    LIVAS_Altitude_id.units         = 'm'
    LIVAS_a532nm_PD_id.units        = 'km-1'
    LIVAS_a532nm_PD_SD_id.units     = 'km-1'     
    LIVAS_PD_MC_id.units            = 'micrograms/m^3'
    LIVAS_PD_MC_SD_id.units         = 'micrograms/m^3'
    LIVAS_PD_MC_CM_id.units         = 'micrograms/m^3'
    LIVAS_PD_MC_CM_SD_id.units      = 'micrograms/m^3'
    LIVAS_PD_MC_FM_id.units         = 'micrograms/m^3'
    LIVAS_PD_MC_FM_SD_id.units      = 'micrograms/m^3'
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
    ERA_U_id.long_name                  = 'U component of wind'
    ERA_U_SD_id.long_name               = 'U component of wind SD'
    ERA_V_id.long_name                  = 'V component of wind'
    ERA_V_SD_id.long_name               = 'V component of wind SD'
    CAMS_Height_id.long_name            = 'Height'
    CAMS_MC_id.long_name                = 'CAMS Pure-Dust Mass Concentration'
    LIVAS_Altitude_id.long_name         = 'Height'
    LIVAS_a532nm_PD_id.units            = 'Pure-Dust Extinction Coefficient 532nm'
    LIVAS_a532nm_PD_SD_id.units         = 'Pure-Dust Extinction Coefficient 532nm - SD'  
    LIVAS_PD_MC_id.long_name            = 'Pure-Dust Mass Concentration'
    LIVAS_PD_MC_SD_id.long_name         = 'Pure-Dust Mass Concentration - SD'
    LIVAS_PD_MC_CM_id.long_name         = 'Pure-Dust Coarse-Mode Mass Concentration'
    LIVAS_PD_MC_CM_SD_id.long_name      = 'Pure-Dust Coarse-Mode Mass Concentration - SD'
    LIVAS_PD_MC_FM_id.long_name         = 'Pure-Dust Fine-Mode Mass Concentration'
    LIVAS_PD_MC_FM_SD_id.long_name      = 'Pure-Dust Fine-Mode Mass Concentration - SD'
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
    ERA_U_id.fill_value                  = np.nan                  
    ERA_U_SD_id.fill_value               = np.nan               
    ERA_V_id.fill_value                  = np.nan                  
    ERA_V_SD_id.fill_value               = np.nan               
    CAMS_Height_id.fill_value            = np.nan            
    CAMS_MC_id.fill_value                = np.nan                
    LIVAS_Altitude_id.fill_value         = np.nan         
    LIVAS_a532nm_PD_id.fill_value        = np.nan
    LIVAS_a532nm_PD_SD_id.fill_value     = np.nan
    LIVAS_PD_MC_id.fill_value            = np.nan            
    LIVAS_PD_MC_SD_id.fill_value         = np.nan         
    LIVAS_PD_MC_CM_id.fill_value         = np.nan         
    LIVAS_PD_MC_CM_SD_id.fill_value      = np.nan      
    LIVAS_PD_MC_FM_id.fill_value         = np.nan         
    LIVAS_PD_MC_FM_SD_id.fill_value      = np.nan      
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
    ERA_Height_id[:]             = ERA_Height
    ERA_U_id[:]                  = ERA_U    
    ERA_U_SD_id[:]               = ERA_U_SD
    ERA_V_id[:]                  = ERA_V  
    ERA_V_SD_id[:]               = ERA_V_SD
    CAMS_Height_id[:]            = CAMS_Height
    CAMS_MC_id[:]                = CAMS_Dust_MC
    LIVAS_Altitude_id[:]         = Altitude
    LIVAS_a532nm_PD_id[:]        = Final_PD_a532nm
    LIVAS_a532nm_PD_SD_id[:]     = Final_PD_a532nm_SD
    LIVAS_PD_MC_id[:]            = Final_PD_MC
    LIVAS_PD_MC_SD_id[:]         = Final_PD_MC_SD
    LIVAS_PD_MC_CM_id[:]         = Final_PD_MC_CM
    LIVAS_PD_MC_CM_SD_id[:]      = Final_PD_MC_CM_SD
    LIVAS_PD_MC_FM_id[:]         = Final_PD_MC_FM
    LIVAS_PD_MC_FM_SD_id[:]      = Final_PD_MC_FM_SD
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
  
    ds.close()               
    
    print("End of: " + fn)  
    print(datetime.now()) 



#  SCRIPT END  ---------------------------------------------------------------
Ou.goodbye(cnf.LOGs.run, tic=tic, scriptname=__file__)
