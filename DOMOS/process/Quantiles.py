# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 08:54:57 2022
Establishes the required EO-and-ERA5 DOMOS dataset in the same L3 1x1 grid resolution and monthly-mean.
@author: proes
"""

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
import matplotlib.pyplot as plt
from   netCDF4 import Dataset

# Initializing timer:
startTime            = datetime.now()
LIVAS_path           = r"D:\LIVAS\Grid"

longitudes_array_1x1 = np.arange(-179.5,180.5,1)
latitudes_array_1x1  = np.arange(-69.5,70.5, 1)

# longitudes_array_1x1 = np.arange(-19.5,-15.5,1)
# latitudes_array_1x1  = np.arange(10.5,13.5, 1)

longitudes_array     = longitudes_array_1x1
latitudes_array      = latitudes_array_1x1

Percentile_95_1x1    = np.empty((len(longitudes_array),len(latitudes_array)))
Percentile_96_1x1    = np.empty((len(longitudes_array),len(latitudes_array)))
Percentile_97_1x1    = np.empty((len(longitudes_array),len(latitudes_array)))
Percentile_98_1x1    = np.empty((len(longitudes_array),len(latitudes_array)))
Percentile_99_1x1    = np.empty((len(longitudes_array),len(latitudes_array)))
Percentile_95_1x1[:] = np.nan
Percentile_96_1x1[:] = np.nan
Percentile_97_1x1[:] = np.nan
Percentile_98_1x1[:] = np.nan
Percentile_99_1x1[:] = np.nan

for count_lon,lon in enumerate(longitudes_array):        
    for count_lat,lat in enumerate(latitudes_array):
 
        print(lon,lat)        

        LIVAS_filename = 'LIVAS_CALIPSO_L2_Grid_lon_c_' + str(lon) + '_lat_c_' + str(lat) + '.nc'
        LIVAS_file     = LIVAS_path + '\\' + LIVAS_filename

        #if LIVAS_filename == 'LIVAS_CALIPSO_L2_Grid_lon_c_-167.5_lat_c_9.5.nc':
        #    continue

        file_existing   = exists(LIVAS_file)
                
        if (not file_existing):
            continue
        else:                     
            LIVAS_dataset   = nc.Dataset(LIVAS_file)
            Altitude        = LIVAS_dataset['/Altitude'][:]
            LIVAS_PD_b532nm = LIVAS_dataset['/LIVAS/Cloud_Free/Pure_Dust_and_Fine_Coarse/Optical_Products/Pure_Dust_Backscatter_Coefficient_532'][:]

            LIVAS_PD_b532nm[: , Altitude > 10] = np.nan
            values = LIVAS_PD_b532nm[(~np.isnan(LIVAS_PD_b532nm)) & (LIVAS_PD_b532nm != 0)]
            # Now compute the desired percentiles
            percentiles = np.percentile(values, [95, 96, 97, 98, 99])

            Percentile_95_1x1[count_lon,count_lat] = percentiles[0]
            Percentile_96_1x1[count_lon,count_lat] = percentiles[1]
            Percentile_97_1x1[count_lon,count_lat] = percentiles[2]
            Percentile_98_1x1[count_lon,count_lat] = percentiles[3]
            Percentile_99_1x1[count_lon,count_lat] = percentiles[4]



longitudes_array_2x2 = np.arange(-179,181,2)
latitudes_array_2x2  = np.arange(-69,71,2)

# longitudes_array_2x2 = np.arange(-39,-31,2)
# latitudes_array_2x2  = np.arange(9,21,2)

longitudes_array     = longitudes_array_2x2
latitudes_array      = latitudes_array_2x2

Percentile_95_2x2    = np.empty((len(longitudes_array),len(latitudes_array)))
Percentile_96_2x2    = np.empty((len(longitudes_array),len(latitudes_array)))
Percentile_97_2x2    = np.empty((len(longitudes_array),len(latitudes_array)))
Percentile_98_2x2    = np.empty((len(longitudes_array),len(latitudes_array)))
Percentile_99_2x2    = np.empty((len(longitudes_array),len(latitudes_array)))
Percentile_95_2x2[:] = np.nan
Percentile_96_2x2[:] = np.nan
Percentile_97_2x2[:] = np.nan
Percentile_98_2x2[:] = np.nan
Percentile_99_2x2[:] = np.nan

for count_lon,lon in enumerate(longitudes_array):        
    for count_lat,lat in enumerate(latitudes_array):
               
        lons_of_interest = [lon-0.5,lon+0.5,lon-0.5,lon+0.5]
        lats_of_interest = [lat-0.5,lat-0.5,lat+0.5,lat+0.5]        
        file_counter = 0

        print(lon,lat) 
        
        for LIVAS_lons_lats_idx,LIVAS_lons_lats in enumerate(lons_of_interest):
        
            LIVAS_filename  = 'LIVAS_CALIPSO_L2_Grid_lon_c_' + str(lons_of_interest[LIVAS_lons_lats_idx]) + '_lat_c_' + str(lats_of_interest[LIVAS_lons_lats_idx]) + '.nc'
            LIVAS_file      = LIVAS_path + '\\' + LIVAS_filename
            file_existing   = exists(LIVAS_file)

            if (not file_existing):
                continue
            else:                     
                LIVAS_dataset   = nc.Dataset(LIVAS_file)
                Altitude        = LIVAS_dataset['/Altitude'][:]
                LIVAS_PD_b532nm = LIVAS_dataset['/LIVAS/Cloud_Free/Pure_Dust_and_Fine_Coarse/Optical_Products/Pure_Dust_Backscatter_Coefficient_532'][:]

                LIVAS_PD_b532nm[: , Altitude > 10] = np.nan

                if file_counter >  0:
                    Total_LIVAS_PD_b532nm = np.vstack([Total_LIVAS_PD_b532nm,LIVAS_PD_b532nm])
                    file_counter = file_counter + 1                          
                if file_counter == 0:
                    Total_LIVAS_PD_b532nm = LIVAS_PD_b532nm
                    file_counter = file_counter + 1            

        values = LIVAS_PD_b532nm[(~np.isnan(LIVAS_PD_b532nm)) & (LIVAS_PD_b532nm != 0)]
        # Now compute the desired percentiles
        percentiles = np.percentile(values, [95, 96, 97, 98, 99])

        Percentile_95_2x2[count_lon,count_lat] = percentiles[0]
        Percentile_96_2x2[count_lon,count_lat] = percentiles[1]
        Percentile_97_2x2[count_lon,count_lat] = percentiles[2]
        Percentile_98_2x2[count_lon,count_lat] = percentiles[3]
        Percentile_99_2x2[count_lon,count_lat] = percentiles[4]       

longitudes_array_5x2 = np.arange(-177.5,182.5,5)
latitudes_array_5x2  = np.arange(-69,71,2)

# longitudes_array_5x2 = np.arange(-27.5,-12.5,5)
# latitudes_array_5x2  = np.arange(9,21,2)

longitudes_array     = longitudes_array_5x2
latitudes_array      = latitudes_array_5x2

Percentile_95_5x2    = np.empty((len(longitudes_array),len(latitudes_array)))
Percentile_96_5x2    = np.empty((len(longitudes_array),len(latitudes_array)))
Percentile_97_5x2    = np.empty((len(longitudes_array),len(latitudes_array)))
Percentile_98_5x2    = np.empty((len(longitudes_array),len(latitudes_array)))
Percentile_99_5x2    = np.empty((len(longitudes_array),len(latitudes_array)))
Percentile_95_5x2[:] = np.nan
Percentile_96_5x2[:] = np.nan
Percentile_97_5x2[:] = np.nan
Percentile_98_5x2[:] = np.nan
Percentile_99_5x2[:] = np.nan

for count_lon,lon in enumerate(longitudes_array):        
    for count_lat,lat in enumerate(latitudes_array):
        
        print(lon,lat)
        
        lons_of_interest = [lon-2,lon-2,lon-1,lon-1,lon,lon,lon+1,lon+1,lon+2,lon+2]
        lats_of_interest = [lat-0.5,lat+0.5,lat-0.5,lat+0.5,lat-0.5,lat+0.5,lat-0.5,lat+0.5,lat-0.5,lat+0.5]        
        file_counter = 0
        
        for LIVAS_lons_lats_idx,LIVAS_lons_lats in enumerate(lons_of_interest):
        
            LIVAS_filename  = 'LIVAS_CALIPSO_L2_Grid_lon_c_' + str(lons_of_interest[LIVAS_lons_lats_idx]) + '_lat_c_' + str(lats_of_interest[LIVAS_lons_lats_idx]) + '.nc'
            LIVAS_file      = LIVAS_path + '\\' + LIVAS_filename
            file_existing   = exists(LIVAS_file)
            
            if (not file_existing):
                continue
            else:                     
                LIVAS_dataset   = nc.Dataset(LIVAS_file)
                Altitude        = LIVAS_dataset['/Altitude'][:]
                LIVAS_PD_b532nm = LIVAS_dataset['/LIVAS/Cloud_Free/Pure_Dust_and_Fine_Coarse/Optical_Products/Pure_Dust_Backscatter_Coefficient_532'][:]

                LIVAS_PD_b532nm[: , Altitude > 10] = np.nan

                if file_counter >  0:
                    Total_LIVAS_PD_b532nm = np.vstack([Total_LIVAS_PD_b532nm,LIVAS_PD_b532nm])
                    file_counter = file_counter + 1                          
                if file_counter == 0:
                    Total_LIVAS_PD_b532nm = LIVAS_PD_b532nm
                    file_counter = file_counter + 1            

        values = LIVAS_PD_b532nm[(~np.isnan(LIVAS_PD_b532nm)) & (LIVAS_PD_b532nm != 0)]
        # Now compute the desired percentiles
        percentiles = np.percentile(values, [95, 96, 97, 98, 99])

        Percentile_95_5x2[count_lon,count_lat] = percentiles[0]
        Percentile_96_5x2[count_lon,count_lat] = percentiles[1]
        Percentile_97_5x2[count_lon,count_lat] = percentiles[2]
        Percentile_98_5x2[count_lon,count_lat] = percentiles[3]
        Percentile_99_5x2[count_lon,count_lat] = percentiles[4] 

#####################################            
#  --- Saving dataset as NetCDF --- #
#####################################

output_path  = 'D:\PC-backup\Projects\DOMOS\Datasets' 

# creating nc. filename and initiallizing:                 
fn           = output_path + '\\' + 'Quantiles_PD_b532nm.nc' 
ds           = Dataset(fn, mode='w', format='NETCDF4')

# create nc. dimensions: 

#longitudes_array_1x1 = np.arange(-179.5,180.5,1)
#latitudes_array_1x1  = np.arange(-69.5,70.5, 1)
#longitudes_array_2x2 = np.arange(-179,181,2)
#latitudes_array_2x2  = np.arange(-69,71,2)
#longitudes_array_5x2 = np.arange(-177.5,182.5,5)
#latitudes_array_5x2  = np.arange(-69,71,2)

lat_1x1      = ds.createDimension('lat_1x1', len(latitudes_array_1x1 )) 
lon_1x1      = ds.createDimension('lon_1x1', len(longitudes_array_1x1))
lat_2x2      = ds.createDimension('lat_2x2', len(latitudes_array_2x2 )) 
lon_2x2      = ds.createDimension('lon_2x2', len(longitudes_array_2x2))
lat_5x2      = ds.createDimension('lat_5x2', len(latitudes_array_5x2 )) 
lon_5x2      = ds.createDimension('lon_5x2', len(longitudes_array_5x2))

group_1x1    = ds.createGroup("1x1")
group_2x2    = ds.createGroup("2x2")
group_5x2    = ds.createGroup("5x2")

# create nc. variables:

lat_1x1_id              = ds.createVariable('1x1/latitude'     , np.float64, ('lat_1x1',),zlib=True) 
lon_1x1_id              = ds.createVariable('1x1/longitude'    , np.float64, ('lon_1x1',),zlib=True)
lat_2x2_id              = ds.createVariable('2x2/latitude'     , np.float64, ('lat_2x2',),zlib=True) 
lon_2x2_id              = ds.createVariable('2x2/longitude'    , np.float64, ('lon_2x2',),zlib=True)
lat_5x2_id              = ds.createVariable('5x2/latitude'     , np.float64, ('lat_5x2',),zlib=True) 
lon_5x2_id              = ds.createVariable('5x2/longitude'    , np.float64, ('lon_5x2',),zlib=True)
Percentile_95_1x1_id    = ds.createVariable('1x1/Percentile_95', np.float64, ('lon_1x1','lat_1x1',),zlib=True)
Percentile_96_1x1_id    = ds.createVariable('1x1/Percentile_96', np.float64, ('lon_1x1','lat_1x1',),zlib=True)
Percentile_97_1x1_id    = ds.createVariable('1x1/Percentile_97', np.float64, ('lon_1x1','lat_1x1',),zlib=True)
Percentile_98_1x1_id    = ds.createVariable('1x1/Percentile_98', np.float64, ('lon_1x1','lat_1x1',),zlib=True)
Percentile_99_1x1_id    = ds.createVariable('1x1/Percentile_99', np.float64, ('lon_1x1','lat_1x1',),zlib=True)
Percentile_95_2x2_id    = ds.createVariable('2x2/Percentile_95', np.float64, ('lon_2x2','lat_2x2',),zlib=True)
Percentile_96_2x2_id    = ds.createVariable('2x2/Percentile_96', np.float64, ('lon_2x2','lat_2x2',),zlib=True)
Percentile_97_2x2_id    = ds.createVariable('2x2/Percentile_97', np.float64, ('lon_2x2','lat_2x2',),zlib=True)
Percentile_98_2x2_id    = ds.createVariable('2x2/Percentile_98', np.float64, ('lon_2x2','lat_2x2',),zlib=True)
Percentile_99_2x2_id    = ds.createVariable('2x2/Percentile_99', np.float64, ('lon_2x2','lat_2x2',),zlib=True)
Percentile_95_5x2_id    = ds.createVariable('5x2/Percentile_95', np.float64, ('lon_5x2','lat_5x2',),zlib=True)
Percentile_96_5x2_id    = ds.createVariable('5x2/Percentile_96', np.float64, ('lon_5x2','lat_5x2',),zlib=True)
Percentile_97_5x2_id    = ds.createVariable('5x2/Percentile_97', np.float64, ('lon_5x2','lat_5x2',),zlib=True)
Percentile_98_5x2_id    = ds.createVariable('5x2/Percentile_98', np.float64, ('lon_5x2','lat_5x2',),zlib=True)
Percentile_99_5x2_id    = ds.createVariable('5x2/Percentile_99', np.float64, ('lon_5x2','lat_5x2',),zlib=True)

lat_1x1_id.units           = 'degrees_north' 
lon_1x1_id.units           = 'degrees_east'
lat_2x2_id.units           = 'degrees_north' 
lon_2x2_id.units           = 'degrees_east'
lat_5x2_id.units           = 'degrees_north' 
lon_5x2_id.units           = 'degrees_east'
Percentile_95_1x1_id.units = 'km sr**-1' 
Percentile_96_1x1_id.units = 'km sr**-1' 
Percentile_97_1x1_id.units = 'km sr**-1' 
Percentile_98_1x1_id.units = 'km sr**-1' 
Percentile_99_1x1_id.units = 'km sr**-1' 
Percentile_95_2x2_id.units = 'km sr**-1' 
Percentile_96_2x2_id.units = 'km sr**-1' 
Percentile_97_2x2_id.units = 'km sr**-1' 
Percentile_98_2x2_id.units = 'km sr**-1' 
Percentile_99_2x2_id.units = 'km sr**-1' 
Percentile_95_5x2_id.units = 'km sr**-1' 
Percentile_96_5x2_id.units = 'km sr**-1' 
Percentile_97_5x2_id.units = 'km sr**-1' 
Percentile_98_5x2_id.units = 'km sr**-1' 
Percentile_99_5x2_id.units = 'km sr**-1' 
 
lat_1x1_id.long_name           = 'Latitude' 
lon_1x1_id.long_name           = 'Longitude'
lat_2x2_id.long_name           = 'Latitude' 
lon_2x2_id.long_name           = 'Longitude'
lat_5x2_id.long_name           = 'Latitude' 
lon_5x2_id.long_name           = 'Longitude'
Percentile_95_1x1_id.long_name = 'Backscatter Coefficient 532 nm - Percentile 95' 
Percentile_96_1x1_id.long_name = 'Backscatter Coefficient 532 nm - Percentile 96' 
Percentile_97_1x1_id.long_name = 'Backscatter Coefficient 532 nm - Percentile 97' 
Percentile_98_1x1_id.long_name = 'Backscatter Coefficient 532 nm - Percentile 98' 
Percentile_99_1x1_id.long_name = 'Backscatter Coefficient 532 nm - Percentile 99' 
Percentile_95_2x2_id.long_name = 'Backscatter Coefficient 532 nm - Percentile 95' 
Percentile_96_2x2_id.long_name = 'Backscatter Coefficient 532 nm - Percentile 96' 
Percentile_97_2x2_id.long_name = 'Backscatter Coefficient 532 nm - Percentile 97' 
Percentile_98_2x2_id.long_name = 'Backscatter Coefficient 532 nm - Percentile 98' 
Percentile_99_2x2_id.long_name = 'Backscatter Coefficient 532 nm - Percentile 99' 
Percentile_95_5x2_id.long_name = 'Backscatter Coefficient 532 nm - Percentile 95' 
Percentile_96_5x2_id.long_name = 'Backscatter Coefficient 532 nm - Percentile 96' 
Percentile_97_5x2_id.long_name = 'Backscatter Coefficient 532 nm - Percentile 97' 
Percentile_98_5x2_id.long_name = 'Backscatter Coefficient 532 nm - Percentile 98' 
Percentile_99_5x2_id.long_name = 'Backscatter Coefficient 532 nm - Percentile 99'

lat_1x1_id.fill_value           = np.nan
lon_1x1_id.fill_value           = np.nan
lat_2x2_id.fill_value           = np.nan
lon_2x2_id.fill_value           = np.nan
lat_5x2_id.fill_value           = np.nan
lon_5x2_id.fill_value           = np.nan
Percentile_95_1x1_id.fill_value = np.nan 
Percentile_96_1x1_id.fill_value = np.nan 
Percentile_97_1x1_id.fill_value = np.nan 
Percentile_98_1x1_id.fill_value = np.nan 
Percentile_99_1x1_id.fill_value = np.nan 
Percentile_95_2x2_id.fill_value = np.nan 
Percentile_96_2x2_id.fill_value = np.nan 
Percentile_97_2x2_id.fill_value = np.nan 
Percentile_98_2x2_id.fill_value = np.nan 
Percentile_99_2x2_id.fill_value = np.nan 
Percentile_95_5x2_id.fill_value = np.nan 
Percentile_96_5x2_id.fill_value = np.nan 
Percentile_97_5x2_id.fill_value = np.nan 
Percentile_98_5x2_id.fill_value = np.nan 
Percentile_99_5x2_id.fill_value = np.nan 

lat_1x1_id[:]           = latitudes_array_1x1
lon_1x1_id[:]           = longitudes_array_1x1
lat_2x2_id[:]           = latitudes_array_2x2
lon_2x2_id[:]           = longitudes_array_2x2
lat_5x2_id[:]           = latitudes_array_5x2
lon_5x2_id[:]           = longitudes_array_5x2
Percentile_95_1x1_id[:] = Percentile_95_1x1
Percentile_96_1x1_id[:] = Percentile_96_1x1
Percentile_97_1x1_id[:] = Percentile_97_1x1
Percentile_98_1x1_id[:] = Percentile_98_1x1
Percentile_99_1x1_id[:] = Percentile_99_1x1
Percentile_95_2x2_id[:] = Percentile_95_2x2
Percentile_96_2x2_id[:] = Percentile_96_2x2
Percentile_97_2x2_id[:] = Percentile_97_2x2
Percentile_98_2x2_id[:] = Percentile_98_2x2
Percentile_99_2x2_id[:] = Percentile_99_2x2
Percentile_95_5x2_id[:] = Percentile_95_5x2
Percentile_96_5x2_id[:] = Percentile_96_5x2
Percentile_97_5x2_id[:] = Percentile_97_5x2
Percentile_98_5x2_id[:] = Percentile_98_5x2
Percentile_99_5x2_id[:] = Percentile_99_5x2
  
ds.close()               
   
# How much time it took to execute the script in minutes:
print(datetime.now() - startTime)  
            