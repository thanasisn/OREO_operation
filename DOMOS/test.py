#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy    as np
import os
import sys
import glob

##  Load project functions  --------------------------------------------------
sys.path.append("../")
import oreo_mod.utils as Ou
import oreo_mod.calc  as Oc







##  Load configuration profile by host name  ---------------------------------
cnf = Ou.get_configs(
        Ou.parse_arguments(run_profiles_folder = "../run_profiles").profile
    )

print(cnf)





