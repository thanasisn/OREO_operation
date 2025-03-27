#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy    as np
import os

H = np.array([ -100, -20, 10, 30, 100])

Hp = H[H>= 0]
Hn = H[H<  0]

Ht = np.concat(( Hn[Hn.min() == Hn] + -100, Hn, [0],  Hp, Hp[Hp.max() == Hp] + 1000))


  ## create height boundaries
HB = H + np.diff(np.concat((H, H[H.max() == H] + cnf.ERA5.extra_height_on_top))) / 2
HB = np.concat(([cnf.ERA5.height_at_bottom], HB))

np.diff(H) / 2

Hp + np.diff(np.concat((Hp, Hp[Hp.max() == Hp] + 1000)))/2
Hn + np.diff(np.concat((Hn[Hn.min() == Hn] + -100, Hn)))/2

Ht + np.diff(Ht)/2
 

# Hp = np.concat(([0], H[H>0]))
# Hn = np.concat((H[H<0], [0]))

HB = np.concat((
    Hn + np.diff(np.concat((Hn[Hn.min() == Hn] -100, Hn)))/2, 
    [0],
    Hp + np.diff( np.concat((Hp, Hp[Hp.max() == Hp] + 1000)) )/2
))

HBn = np.concat((Hn - np.diff(np.concat((Hn[Hn.min() == Hn] - 100, Hn)))/2, [0]))
HBp = np.concat(([0], Hp + np.diff( np.concat((Hp, Hp[Hp.max() == Hp] + 1000)) )/2))
Hp

for i,y in enumerate(Hp):
    print("L: %9.1f -> [ %9.1f ,  %9.1f ]" % (Hp[i], HBp[i], HBp[i+1]))


for i,y in enumerate(Hn):
    print("L: %9.1f -> [ %9.1f ,  %9.1f ]" % (Hn[i], HBn[i], HBn[i+1]))

HBp[1:]  ## upper
HBp[0:-1] ## lower

HBn[1:]  ## upper
HBn[0:-1] ## lower


upper = np.concat((HBn[1:],   HBp[1:]))
lower = np.concat((HBn[0:-1], HBp[0:-1]))
H
 
for i,y in enumerate(H):
    print("L: %9.1f -> [ %9.1f ,  %9.1f ]" % (H[i], lower[i], upper[i]))
 
 
# TEST
# os.chdir("./DOMOS")
