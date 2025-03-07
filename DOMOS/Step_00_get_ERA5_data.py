#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""@package docstring
Download ERA5 nc files for the current scheme.

The scheme is assumed by the name of the current host.

@author: thanasisn
"""

import os
import sys
from   datetime import datetime
from   dotmap   import DotMap
from   domos_functions.size_to_human import size_to_human
import cdsapi
import yaml

# TEST
# os.chdir("./DOMOS")

SCRIPT_NAME = __file__
tic         = datetime.now()

#  Load configuration by host name  --------------------------------------------
config_file = os.uname()[1] + '.yaml'
with open(config_file, 'r') as config_fl:
    configs = yaml.safe_load(config_fl)
    # Convert dictionary to use dot notation
    cnf = DotMap(configs)
    print("\nRead config file:", config_file, "\n")

#  Check destination folder has been created  ----------------------------------
if not os.path.isdir(cnf.ERA5.path_raw):
    sys.exit("\nFolder " + cnf.ERA5.path_raw + " don't exist!\n")

#  Expand domain to download  --------------------------------------------------
print("Get domain:", 
      "%sN %sS %sW %sE" % (cnf.ERA5.North, cnf.ERA5.South, cnf.ERA5.West, cnf.ERA5.East))

# start api
client = cdsapi.Client(quiet=True)

#  Get each year data  ---------------------------------------------------------
for yyyy in range(cnf.Range.start, cnf.Range.until + 1):
    # filename of file to download
    target = os.path.join(
        cnf.ERA5.path_raw, 
        "ERA5_%s_%sN%sS%sW%sE.nc" % 
        (yyyy, cnf.ERA5.North, cnf.ERA5.South, cnf.ERA5.West, cnf.ERA5.East)
    )
    
    # check if we already have the file
    if os.path.isfile(target):
        print("Skip existing file:",
              os.path.basename(target),
              datetime.fromtimestamp((os.path.getmtime(target))).strftime("%F %T"),
              size_to_human(os.path.getsize(target)))
        continue

    # get ERA5 data
    print("Get new", os.path.basename(target))

    dataset = "reanalysis-era5-pressure-levels-monthly-means"
    request = {
        "product_type": ["monthly_averaged_reanalysis"],
        "variable": [
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
    # call the api
    client.retrieve(dataset, request, target)

#  SCRIPT END  -----------------------------------------------------------------
out  = datetime.now().strftime("%F %T") + " "
out += os.getlogin() + "@" + os.uname()[1] + " "
out += SCRIPT_NAME + " "
out += str(round((datetime.now() - tic).total_seconds() / 60.0, 2)) + " mins"
print('\n' + out + '\n')
with open(cnf.LOGs.run, 'a') as runlog:
    runlog.write(out + '\n')
