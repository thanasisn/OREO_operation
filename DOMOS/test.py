#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy    as np
import os
import sys
import glob
import tqdm
import time
import pathlib
import subprocess

##  Load project functions  --------------------------------------------------
sys.path.append("../")
import oreo_mod.utils as Ou
import oreo_mod.calc  as Oc





Ou.source_code_hash("/home/athan/OREO/operation/DOMOS/Step_01_regrid_ERA5.py")


##  Load configuration profile by host name  ---------------------------------
cnf = Ou.get_configs(
        Ou.parse_arguments(run_profiles_folder = "../run_profiles").profile
    )

# for outer in [10, 20, 30, 40, 50]:
#     print(outer,"\n")
#     for inner in tqdm.tqdm(range(outer), desc="This file", dynamic_ncols = True):
#         print(inner,"\n")
#         time.sleep(1)
# print("done!")
#
