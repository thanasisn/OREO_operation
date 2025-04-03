#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy    as np
import os
import sys

##  Load project functions  --------------------------------------------------
sys.path.append("../")
import oreo_mod.utils as Ou
import oreo_mod.calc  as Oc

##  Load configuration profile by host name  ---------------------------------
config_file = f"../run_profiles/{os.uname()[1]}.yaml"
cnf = Ou.get_configs(config_file)


cnf.D1.LonStep
cnf.D1.LatStep


lat = 43.
lon = -77.5

print(np.arange( lat - (cnf.D1.LatStep / 2),  lat + (cnf.D1.LatStep / 2) + 1, 1))
print(np.arange( lon - (cnf.D1.LonStep / 2),  lon + (cnf.D1.LonStep / 2) + 1, 1))

np.linspace(1, 3, 1)

lats = np.arange( -89.5,  90)
lons = np.arange(-179.5, 180)

lats[np.logical_and(lats > lat - (cnf.D1.LatStep / 2),
                    lats < lat + (cnf.D1.LatStep / 2))]

lons[np.logical_and(lons > lon - (cnf.D1.LonStep / 2),
                    lons < lon + (cnf.D1.LonStep / 2))]


lats[lats > lon - (cnf.D1.LonStep / 2)]

lats > lon - (cnf.D1.LonStep / 2)
lats < lon + (cnf.D1.LonStep / 2)
lats[np.logical_and( lats > lon - (cnf.D1.LonStep / 2),
                lats < lon + (cnf.D1.LonStep / 2))]



lats.isel(1)

lats < lon + (cnf.D1.LonStep / 2)
