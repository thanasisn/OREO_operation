#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 17:29:13 2021
This script uses as input the ERA5 u, v, w, lat, lon, and geopotential and performs two tasks:
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
from   metpy.units import units

import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import xarray_regrid


#  Load project functions
sys.path.append("../")
import oreo_mod.utils as Ou
import oreo_mod.calc  as Oc
tic = datetime.now()

#  TEST
# os.chdir("./DOMOS")

#  Force the reprocess of the inputs
FORCE = False
FORCE = True

#  Load configuration profile by host name  ----------------------------------
config_file = "../run_profiles/" + os.uname()[1] + '.yaml'
cnf = Ou.get_configs(config_file)

#  Check destination folder exists  ------------------------------------------
if not os.path.isdir(cnf.ERA5.path_regrid):
    sys.exit("\nFolder " + cnf.ERA5.path_regrid + " don't exist!\n")


#  List input files exist in input dir  --------------------------------------
#  Use expanded domain to find input files
fl_North = Oc.border_up(  cnf.D1.North, cnf.D1.LatStep)
fl_South = Oc.border_down(cnf.D1.South, cnf.D1.LatStep)
fl_East  = Oc.border_up(  cnf.D1.East,  cnf.D1.LonStep)
fl_West  = Oc.border_down(cnf.D1.West,  cnf.D1.LonStep)
filenames = glob.glob(f"{cnf.ERA5.path_raw}/ERA5_*_{fl_North}N{fl_South}S{fl_West}W{fl_East}E.nc")
filenames.sort()

if len(filenames) < 1:
    sys.exit("\nNo input file found in " + cnf.ERA5.path_raw + " !\n")


# ESA-DOMOS: "... the full coverage of the Atlantic Ocean (including dust emission sources of Africa and S. America,
# the broader Atlantic Ocean, Caribbean Sea and Gulf of Mexico, confined between latitudes 40°N to 60°S), and of
# temporal coverage at least between 2010 and 2020".
# Therefore:
# (I)  DOMOS lon: -105E:5:25E
# (II) DOMOS lat:  -65N:5:45N
# (a wider domain is used here to account for (1) all fluxes and (2) the broader domain, and N.Atlandic Dust)


### DOMOS ####
# lon_array       = np.arange(-125, 25)
# lat_array       = np.arange( -60, 42)
# DOMOS_lon_array = np.arange(-125, 25, 5)
# DOMOS_lat_array = np.arange( -62, 42, 2)

lon_array       = np.arange(cnf.D1.South, cnf.D1.North)
lat_array       = np.arange(cnf.D1.West,  cnf.D1.East)
DOMOS_lon_array = np.arange(cnf.D1.South, cnf.D1.North, cnf.D1.LonStep)
DOMOS_lat_array = np.arange(cnf.D1.West,  cnf.D1.East,  cnf.D1.LatStep)


# Process raw ERA5 files  ------------------------------------------------------------
for filein in filenames:
    yyyy = int(re.compile('ERA5_([0-9]*)_.*.nc').search(filein).group(1))
    ## limit data range
    if not cnf.Range.start <= yyyy <= cnf.Range.until:
        continue

    print(f"\nProcessing: {filein}")

    ## load on an numpy array
    dataset = nc.Dataset(filein)

    ## load ERA5 on a xarray
    DT = xr.open_dataset(filein)

    ## TODO find the correct multiple of boundaries to subset!

    DT.longitude.values
    DT.latitude.values

    ##  Ger ERA5 spatial step
    lat_res = (np.unique(np.diff(DT.latitude.values))[0])
    lon_res = (np.unique(np.diff(DT.longitude.values))[0])

    ## apply a domain constrains
    ## remove the one boundary for each dimension
    # DT = DT.sel(longitude = slice(cnf.ERA5.West,  cnf.ERA5.East ),
    #             latitude  = slice(cnf.ERA5.North, cnf.ERA5.South))

    ## original grid
    DT.longitude.values
    DT.latitude.values
   
    
    np.arange(cnf.D1.North, cnf.D1.South, -lon_res)
    np.arange(cnf.D1.West,  cnf.D1.East,  -lat_res)

    np.arange(cnf.D1.North, cnf.D1.South, -cnf.D1.LonStep)
    np.arange(cnf.D1.West,  cnf.D1.East,   cnf.D1.LatStep)

    -cnf.ERA5.LatStep / lat_res

    DT.longitude.values.min()
    DT.longitude.values.max()
    len(DT.longitude.values)
    len(np.arange(cnf.D1.North, cnf.D1.South, -cnf.D1.LonStep))


    DT.latitude.values.min()
    DT.latitude.values.max()
    len(DT.latitude.values)


    sys.exit("stop")
    
    DT.dims
    DT.info
    DT.data_vars
    DT.coords
    DT.u
    DT.z.units


    ## test plot
    DT.u.isel(pressure_level = 0, valid_time = 0).plot()
    DT.v.isel(pressure_level = 0, valid_time = 0).plot()

    ## create height variable
    DT = DT.assign(height = metpy.calc.geopotential_to_height(DT.z))
    DT['height'].attrs = {
        'long_name':     'Height',
        'units':         'm',
        'standard_name': 'height'
    }

    # np.array(geopot)
    # geopot.meters
    # units.Quantity(geopot)
    # DT['heigth'] = geopot

    DT.height.values
    DT.height.attrs
    DT.u.values
    DT.u.attrs
    DT.dims
    DT.data_vars

    # ## test write to nc
    # comp = dict(zlib=True, complevel=5)
    # encoding = {var: comp for var in DT.data_vars}
    # DT.to_netcdf("test.nc", mode = 'w', engine = "netcdf4", encoding = encoding)


    ## grid data

    dd.dims
    dd.data_vars
    dd.coords
    dd.info



    dd = DT.coarsen(latitude  = int(-cnf.D1.LatStep / lat_res), 
                    longitude = int( cnf.D1.LonStep / lon_res),
                    boundary = "trim").mean()
    
    dd.longitude.values    
    np.diff(dd.longitude.values)
    
    dd.latitude.values
    np.diff(dd.latitude.values)
    

    # ## with pad may get under representation of values due to unequal bins
    # bb = DT.coarsen(latitude  = int(-cnf.ERA5.LatStep / lat_res),
    #                 longitude = int( cnf.ERA5.LonStep / lon_res),
    #                 boundary  ='pad').mean()
    # bb.longitude.values    
    # np.diff(bb.longitude.values)


    sys.exit("stop")

    # https://confluence.ecmwf.int/display/CKB/ERA5%3A+What+is+the+spatial+reference
    # ERA longitude from 0->360 deg to -180->180 deg.
    dataset_latitude     = dataset['latitude'][:]
    ERA5_idx_lat         = np.where((dataset_latitude >= lat_array[0]) & (dataset_latitude <= lat_array[-1]+1))
    ERA5_idx_lat         = np.ravel(ERA5_idx_lat)
    dataset_latitude     = dataset['latitude'][ERA5_idx_lat]

    dataset_longitude    = dataset['longitude'][:]

    for lon in dataset_longitude:
        if lon >= 180:
            dataset_longitude[np.where(lon == dataset_longitude)] = lon - 360
            print("lon >= 180")

    ERA5_idx_lon         = np.where((dataset_longitude >= lon_array[0]) & (dataset_longitude <= lon_array[-1]+1))
    ERA5_idx_lon         = np.ravel(ERA5_idx_lon)
    dataset_longitude    = dataset['longitude'][ERA5_idx_lon]

    for lon in dataset_longitude:
        if lon >= 180:
            dataset_longitude[np.where(lon == dataset_longitude)] = lon - 360
            print("lon >= 180")

    dataset_time         = dataset['valid_time'][:]
    dataset_level        = dataset['pressure_level'][:]
    dataset_u            = dataset['u'][:,:,ERA5_idx_lat, ERA5_idx_lon]
    dataset_v            = dataset['v'][:,:,ERA5_idx_lat, ERA5_idx_lon]


    # ERA convert Geopotenial to geometric height (a.m.s.l.):
    # https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.geopotential_to_height.html
    dataset_geo_units    = dataset['z'].units
    dataset_geopotential = dataset['z'][:,:,ERA5_idx_lat, ERA5_idx_lon]
    dataset_geopotential = units.Quantity(dataset_geopotential, dataset_geo_units)
    dataset_height       = metpy.calc.geopotential_to_height(dataset_geopotential)
    dataset_height       = np.array(dataset_height)


    # File of previous year - to read December for DJF season
    previous_file = list(filter(lambda x:'ERA5_'+str(yyyy-1) in x, filenames))
    if (len(previous_file)!=1):
        print("SKIP! No file for previous year exists\n")
        continue

    # Year_of_Interest_previous = mypath + '\\' + str(int(filename[0:4])-1)+'.nc'
    dataset_ΙΙ              = nc.Dataset(previous_file[0])
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
    # sys.exit("stop")

    seasons = ['Q1_DJF', 'Q2_MAM', 'Q3_JJA', 'Q4_SON']
    for season_idx, season in enumerate(seasons):

        fileout = cnf.ERA5.path_regrid + "/ERA5_%s_%s_%sN%sS%sW%sE.nc" % (yyyy, season, cnf.ERA5.North, cnf.ERA5.South, cnf.ERA5.West, cnf.ERA5.East)

        if (not FORCE) and (not Ou.output_needs_update(filein, fileout)):
            continue


        # sys.exit("stop")

        # initializing: u-mean,SD / v-mean,SD / w-mean,SD / z-mean for 1x1 deg2 grid resolution.
        u_total_mean = np.empty((len(DOMOS_lon_array), len(DOMOS_lat_array), len(dataset_level)))
        u_total_SD   = np.empty((len(DOMOS_lon_array), len(DOMOS_lat_array), len(dataset_level)))
        v_total_mean = np.empty((len(DOMOS_lon_array), len(DOMOS_lat_array), len(dataset_level)))
        v_total_SD   = np.empty((len(DOMOS_lon_array), len(DOMOS_lat_array), len(dataset_level)))
        z_total_mean = np.empty((len(DOMOS_lon_array), len(DOMOS_lat_array), len(dataset_level)))

        # computing and saving: u-mean,SD / v-mean,SD / z-mean for 2x5 deg2 grid resolution.
        for lon in DOMOS_lon_array:

            idx_lon        = np.where((dataset_longitude >= lon) & (dataset_longitude <= lon + cnf.ERA5.LonStep))
            idx_lon        = np.ravel(idx_lon)
            temp_u         = dataset_u[:,:,:,idx_lon]
            temp_v         = dataset_v[:,:,:,idx_lon]
            temp_height    = dataset_height[:,:,:,idx_lon]
            temp_u_II      = dataset_u_ΙΙ[:,:,:,idx_lon]
            temp_v_II      = dataset_v_ΙΙ[:,:,:,idx_lon]
            temp_height_II = dataset_height_ΙΙ[:,:,:,idx_lon]

            for lat in DOMOS_lat_array:

                idx_lat = np.where((dataset_latitude >= lat) & (dataset_latitude <= lat + cnf.ERA5.LatStep))
                idx_lat = np.ravel(idx_lat)

                if season == 'Q1_DJF':
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
                if season == 'Q2_MAM':
                    Months_of_Interest_idx = [3, 5]
                    u = temp_u[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
                    v = temp_v[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
                    z = temp_height[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
                if season == 'Q3_JJA':
                    Months_of_Interest_idx = [6, 8]
                    u = temp_u[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
                    v = temp_v[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
                    z = temp_height[Months_of_Interest_idx[0]:Months_of_Interest_idx[1]+1,:,idx_lat,:]
                if season == 'Q4_SON':
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
        ds           = nc.Dataset(fileout, 'w', format='NETCDF4')

        # create nc. dimensions:
        longitude    = DOMOS_lon_array + cnf.ERA5.LonStep / 2
        latitude     = DOMOS_lat_array + cnf.ERA5.LatStep / 2
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
        lats[:]    = latitude
        lons[:]    = longitude
        Height[:]  = z_total_mean
        U[:]       = u_total_mean
        U_SD[:]    = u_total_SD
        V[:]       = v_total_mean
        V_SD[:]    = v_total_SD

        ds.close()

        print(f"\nWritten: {fileout}")


#  SCRIPT END  ---------------------------------------------------------------
Ou.goodbye(cnf.LOGs.run, tic=tic, scriptname=__file__)
