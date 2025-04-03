#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Download ERA5 nc files for the current scheme.

The scheme is assumed by the name of the current host.
The downloaded files are cropped at the CDS server side.

@author: thanasisn
"""

import os
import sys
from   datetime import datetime
import cdsapi

##  Load project functions
sys.path.append("../")
import oreo_mod.utils as Ou
import oreo_mod.calc  as Oc
tic = datetime.now()

##  Load configuration profile by host name  ---------------------------------
config_file = f"../run_profiles/{os.uname()[1]}.yaml"
cnf = Ou.get_configs(config_file)

QUIET = cnf.mode.Quiet

##  Check destination folder exists  -----------------------------------------
if not os.path.isdir(cnf.ERA5.path_raw):
    sys.exit(f"\nFolder {cnf.ERA5.path_raw} don't exist!\n")

print(f"Want ERA5 domain:                      {cnf.D1.North}N {cnf.D1.South}S {cnf.D1.West}W {cnf.D1.East}E")

##  Override random domain with target resolution boundaries  ----------------
cnf.D1.North = Oc.border_up(  cnf.D1.North, cnf.D1.MaxLatStep,   90)
cnf.D1.South = Oc.border_down(cnf.D1.South, cnf.D1.MaxLatStep,  -90)
cnf.D1.East  = Oc.border_up(  cnf.D1.East,  cnf.D1.MaxLonStep,  180)
cnf.D1.West  = Oc.border_down(cnf.D1.West,  cnf.D1.MaxLonStep, -180)

print(f"Domain expanded according to cell resolution, lat: [ {cnf.D1.South}, {cnf.D1.North} ], lon: [ {cnf.D1.West}, {cnf.D1.East} ]")

##  Start cds api client
client = cdsapi.Client(quiet = QUIET)

#  Get each years data  ------------------------------------------------------
for yyyy in range(cnf.Range.start, cnf.Range.until + 1):
    #  Filename of the file to download
    target = os.path.join(
        cnf.ERA5.path_raw,
        f"ERA5_{yyyy}_lat_{cnf.D1.South}_{cnf.D1.North}_lon_{cnf.D1.West}_{cnf.D1.East}.nc"
    )

    #  Skip if we already have the file
    if Ou.true_if_file_exist(target):
        continue

    #  Get ERA5 data
    print("Get new", os.path.basename(target))

    dataset = "reanalysis-era5-pressure-levels-monthly-means"
    request = {
        "product_type": ["monthly_averaged_reanalysis"],
        "variable": [
            "geopotential",
            "u_component_of_wind",
            "v_component_of_wind"
        ],
        "pressure_level": [
              "1",   "2",   "3",
              "5",   "7",  "10",
             "20",  "30",  "50",
             "70", "100", "125",
            "150", "175", "200",
            "225", "250", "300",
            "350", "400", "450",
            "500", "550", "600",
            "650", "700", "750",
            "775", "800", "825",
            "850", "875", "900",
            "925", "950", "975",
            "1000"
        ],
        "year": [yyyy],
        "month": [
            "01", "02", "03",
            "04", "05", "06",
            "07", "08", "09",
            "10", "11", "12"
        ],
        "time": ["00:00"],
        "data_format":     "netcdf",
        "download_format": "unarchived",
        "area": [cnf.D1.North,
                 cnf.D1.West,
                 cnf.D1.South,
                 cnf.D1.East]
    }
    #  Call the api to retrieve data
    client.retrieve(dataset, request, target)

#  SCRIPT END  ---------------------------------------------------------------
Ou.goodbye(cnf.LOGs.run, tic=tic, scriptname=__file__)
