#!/usr/bin/env python
"""@package docstring
Download ERA5 nc files for the current scheme.

The scheme is assumed by the name of the current host.

@author: thanasisn
"""

import os
import sys
from datetime import datetime
from dotmap import DotMap
from domos_functions.size_to_human import size_to_human
import cdsapi
import yaml

SCRIPT_NAME = __file__
tic         = datetime.now()

#  Load configuration by host name  --------------------------------------------
config_file = os.uname()[1] + '.yaml'
with open(config_file, 'r') as config_fl:
    configs = yaml.safe_load(config_fl)
    # Convert dictionary to use dot notation
    co = DotMap(configs)
    print("\nRead config file:", config_file, '\n')

#  Check destination folder has been created  ---------------------------------
if not os.path.isdir(co.ERA5_input.path_raw):
    sys.exit("\nFolder " + co.ERA5_input.path_raw + " don't exist!\n")

#  Expand domain to download  -------------------------------------------------
NO = co.Domain.North + 1
WE = co.Domain.West  - 1
SO = co.Domain.South - 1
EA = co.Domain.East  + 1
print("Get domain:", "%sN %sS %sW %sE" % (NO, SO, WE, EA))

# start api
client = cdsapi.Client(quiet=True)

#  Get each year data  --------------------------------------------------------
for yyyy in range(co.Range.start, co.Range.until + 1):
    # filename of file to download
    target = os.path.join(co.ERA5_input.path_raw, "ERA5_" + str(yyyy) + ".nc")

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
            "1", "2", "3",
            "5", "7", "10",
            "20", "30", "50",
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
        "area": [NO, WE, SO, EA]
    }

    client.retrieve(dataset, request, target)

# SCRIPT END
out  = datetime.now().strftime("%F %T") + " "
out += os.getlogin() + "@" + os.uname()[1] + " "
out += SCRIPT_NAME + " "
out += str(round((datetime.now() - tic).total_seconds() / 60.0, 2)) + " mins"
print('\n' + out + '\n')
with open(co.LOGs.run, 'a') as runlog:
    runlog.write(out + '\n')
