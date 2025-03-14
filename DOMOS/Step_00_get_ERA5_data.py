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

#  Load project functions
sys.path.append("../")
import oreo_mod.utils as Ou
import oreo_mod.calc  as Oc
tic = datetime.now()

#  TEST
# os.chdir("./DOMOS")

#  Load configuration profile by host name  ----------------------------------
config_file = "../run_profiles/" + os.uname()[1] + '.yaml'
cnf = Ou.get_configs(config_file)

#  Check destination folder exists  ------------------------------------------
if not os.path.isdir(cnf.ERA5.path_raw):
    sys.exit("\nFolder " + cnf.ERA5.path_raw + " don't exist!\n")

print(f"Want ERA5 domain: {cnf.ERA5.North}N {cnf.ERA5.South}S {cnf.ERA5.West}W {cnf.ERA5.East}E")

cnf.ERA5.North = Oc.border_up(  cnf.ERA5.North, cnf.ERA5.LonStep)
cnf.ERA5.South = Oc.border_down(cnf.ERA5.South, cnf.ERA5.LonStep)
cnf.ERA5.East  = Oc.border_up(  cnf.ERA5.East,  cnf.ERA5.LatStep)
cnf.ERA5.West  = Oc.border_down(cnf.ERA5.West,  cnf.ERA5.LatStep)

print(f"Expand domain acording to resolution: {cnf.ERA5.North}N {cnf.ERA5.South}S {cnf.ERA5.West}W {cnf.ERA5.East}E")

#  Start cds api client
client = cdsapi.Client(quiet=True)

#  Get each years data  ------------------------------------------------------
for yyyy in range(cnf.Range.start, cnf.Range.until + 1):
    #  Filename of the file to download
    target = os.path.join(
        cnf.ERA5.path_raw,
        f"ERA5_{yyyy}_{cnf.ERA5.North}N{cnf.ERA5.South}S{cnf.ERA5.West}W{cnf.ERA5.East}E.nc"
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
        "area": [cnf.ERA5.North,
                 cnf.ERA5.West,
                 cnf.ERA5.South,
                 cnf.ERA5.East]
    }
    #  Call the api to retrieve data
    client.retrieve(dataset, request, target)

#  SCRIPT END  ---------------------------------------------------------------
Ou.goodbye(cnf.LOGs.run, tic=tic, scriptname=__file__)
