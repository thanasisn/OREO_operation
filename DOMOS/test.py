#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os

print(os.getcwd())

import sys
sys.path.append("../")

import oreo_mod.utils as Ou


print(Ou.size_to_human(1111))

# TEST
# os.chdir("./DOMOS")
