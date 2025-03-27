#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy    as np
import os
import sys

H = np.array([ -100, -20, 10, 30, 100])
# H = np.array([  10, 30, 100])

Hp = H[H>= 0]
Hn = H[H<  0]

len(Hn)

# Ht = np.concat(( Hn[Hn.min() == Hn] + -100, Hn, [0],  Hp, Hp[Hp.max() == Hp] + 1000))
# 
# 
# ## create height boundaries
# HB = H + np.diff(np.concat((H, H[H.max() == H] + cnf.ERA5.extra_height_on_top))) / 2
# HB = np.concat(([cnf.ERA5.height_at_bottom], HB))
# 
# np.diff(H) / 2
# 
# Hp + np.diff(np.concat((Hp, Hp[Hp.max() == Hp] + 1000)))/2
# Hn + np.diff(np.concat((Hn[Hn.min() == Hn] + -100, Hn)))/2
# 
# Ht + np.diff(Ht)/2
 

# Hp = np.concat(([0], H[H>0]))
# Hn = np.concat((H[H<0], [0]))

# HB = np.concat((
#     Hn + np.diff(np.concat((Hn[Hn.min() == Hn] -100, Hn)))/2, 
#     [0],
#     Hp + np.diff( np.concat((Hp, Hp[Hp.max() == Hp] + 1000)) )/2
# ))

HBn = np.concat((Hn - np.diff(np.concat((Hn[Hn.min() == Hn] - 100, Hn)))/2, [0]))
HBp = np.concat(([0], Hp + np.diff( np.concat((Hp, Hp[Hp.max() == Hp] + 1000)) )/2))

for i,y in enumerate(Hp):
    print("L: %9.1f -> [ %9.1f ,  %9.1f ]" % (Hp[i], HBp[i], HBp[i+1]))


for i,y in enumerate(Hn):
    print("L: %9.1f -> [ %9.1f ,  %9.1f ]" % (Hn[i], HBn[i], HBn[i+1]))


upper = np.concat((HBn[1:],   HBp[1:]))
lower = np.concat((HBn[0:-1], HBp[0:-1]))
H
 
for i,y in enumerate(H):
        print("L: %9.1f -> [ %9.1f ,  %9.1f ]" % (H[i], lower[i], upper[i]))


def height_bounds_through_zero(heights, remove_bottom = 100, add_top = 1000, quiet = True):
    H = heights

    if any(H == 0):
        ## In this case we should change this function to deal with that
        sys.exit("\n\nCan not compute bounds when heights incude zero!!\n")
    
    ## split heights by sign
    Hp = H[H> 0]
    Hn = H[H< 0]
    
    ## create bounds always with zero on the center
    HBp = np.concat(([0], Hp + np.diff( np.concat((Hp, Hp[Hp.max() == Hp] + add_top)) )/2))
    if len(Hn)>0:
        HBn = np.concat((Hn - np.diff(np.concat((Hn[Hn.min() == Hn] - remove_bottom, Hn)))/2, [0]))
        ## saperate bounds for output
        upper = np.concat((HBn[1:],   HBp[1:]))
        lower = np.concat((HBn[0:-1], HBp[0:-1]))
    else:
        ## saperate bounds for output
        upper = HBp[1:]
        lower = HBp[0:-1]

    if not quiet:
        for i,y in enumerate(H):
            print("L: %9.1f -> [ %9.1f ,  %9.1f ]" % (H[i], lower[i], upper[i]))
    
    ## rerurn a simple array with heights and coresponding bounds        
    return(np.array((H, lower, upper)))


fff = height_bounds_through_zero(H)



height_bounds_through_zero(np.array([0, 10, 30, 100]), quiet=False)

any(np.array([0, 10, 30, 100]) == 0)
    
tt = np.array((H, lower, upper))

tt[0,:] 
tt[1,:] 
tt[2,:] 

tt[0]
tt[1]
tt[2]

 
