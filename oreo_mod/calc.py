# -*- coding: utf-8 -*-
"""
Functions for parameters calculation.  Calculation here are essential to the
operation of the project.  They define some aspects of the methodology.

@author: thanasisn
"""

import sys
import warnings
import metpy.calc
import xarray as xr
import numpy  as np
from   metpy.units import units


def z_to_height(DT):
    """
    Return a new xarray where the height is calculated from the geopotential (z
    variable).

    Parameters
    ----------
    DT : xarray
        A xarray with geopotential in the z variable.

    Returns
    -------
    A new xarray with the extra variable of height.
    """
    DT = DT.assign(
        height = xr.DataArray(
            metpy.calc.geopotential_to_height(
                 units.Quantity(DT.z.values, DT.z.units)
            ),
            coords = DT.coords,
        )
    )
    DT['height'].attrs = {
        'long_name':     'Geometric height',
        'units':         'm',
        'standard_name': 'height'
    }
    return DT


def border_up(target, step, limit):
    """
    Return the next coordinate that satisfy the given step, going to the
    positive direction.  Usually for North and East boundaries. Thus the
    limit should be 90 for latitude and 180 for longitude.

    Parameters
    ----------
    target : float
        A target for the coordinate.
    step : TYPE
        A coordinate target with the step/resolution we need.
    limit : float
        The upper coordinate limit.

    Returns
    -------
    A value that includes the target and may be bigger to meet the next step
    value, but always limited by upper 'limit'.

    >>> border_up(-11, 5, -1)
    -10
    >>> border_up(11, 5, 90)
    15
    >>> border_up(91, 5, 90)
    90
    >>> border_up( 89, 5,  90)
    90
    """
    res = target + step - target % step
    if not target * limit > 0:
        warnings.warn("Target coordinate and limit have different signs !!")
    if res > limit:
        res = limit
    return res



def border_down(target, step, limit):
    """
    Return the next coordinate that satisfy the given step, going to the
    negative direction.  Usually for South and West boundaries.  Thus the
    limit should be -90 for latitude and -180 for longitude.

    Parameters
    ----------
    target : float
        A target for the coordinate.
    step : TYPE
        A coordinate target with the step/resolution we need.
    limit : float
        The lower coordinate limit.

    Returns
    -------
    A value that includes the target and may be smaller to meet the previous
    step value, but always limited by lower 'limit'

    >>> border_down(-11, 5, -90)
    -15
    >>> border_down( 11, 5, 1)
    10
    >>> border_down(-91, 5, -90)
    -90
    >>> border_down(-89, 5, -90)
    -90
    """
    res = target - target % step
    if not target * limit > 0:
        warnings.warn("Target coordinate and limit have different signs !!")
    if res < limit:
        res = limit
    return res


def height_bounds(heights, remove_bottom = 100, add_top = 1000, quiet = True):
    """
    Given an array of ascending sorted heights, return an array with the
    heights and their bounds. The lower and upper heights bounds are expanded
    by the given parameters.

    Parameters
    ----------
    heights : float
        A vector of ascending sorted heights.
    remove_bottom : float
        Subtract from bottom height to create cell bounds.
    add_top: float
        Add to top height to create cell bounds.

    Returns
    -------
    An simple array where first column is the original heights, second column
    the lower bounds and third columns the upper bounds.

    >>> height_bounds(np.array([-10, 10, 30, 100]))
    array([[-10.,  10.,  30., 100.],
           [-60.,   0.,  20.,  65.],
           [  0.,  20.,  65., 600.]])
    >>> height_bounds(np.array([-20, -10, 10, 30, 100]))
    array([[-20., -10.,  10.,  30., 100.],
           [-70., -15.,   0.,  20.,  65.],
           [-15.,   0.,  20.,  65., 600.]])
    >>> height_bounds(np.array([-10, 0, 10, 30, 100]))
    array([[-10.,   0.,  10.,  30., 100.],
           [-60.,  -5.,   5.,  20.,  65.],
           [ -5.,   5.,  20.,  65., 600.]])
    >>> height_bounds(np.array([-20, -10, 0, 10, 30, 100]))
    array([[-20., -10.,   0.,  10.,  30., 100.],
           [-70., -15.,  -5.,   5.,  20.,  65.],
           [-15.,  -5.,   5.,  20.,  65., 600.]])
    """
    H = heights

    if any(H == 0):
        """
        Heights include zero, so we use it to create bounds above and bellow
        zero
        """
        ## split heights by sign
        Hp = H[H >= 0]
        Hn = H[H <= 0]

        ## create bounds always with zero on the centre
        HBp = Hp + np.diff(np.concat((Hp, Hp[Hp.max() == Hp] + 1000))) / 2
        if len(Hn) > 0:
            HBn = Hn - np.diff(np.concat((Hn[Hn.min() == Hn] - 100, Hn))) / 2
            ## separate bounds for output
            upper = np.concat((HBn[1:], HBp[:]))
            lower = np.concat((HBn[0:], HBp[0:-1]))
        else:
            ## separate bounds for output
            upper = HBp[:]
            lower = HBp[0:-1]
    else:
        """
        Heights do not include zero, so we add a zero to create bounds above
        and bellow zero
        """
        ## split heights by sign
        Hp = H[H > 0]
        Hn = H[H < 0]

        ## create bounds always with zero on the centre
        HBp = np.concat(([0], Hp + np.diff(np.concat((Hp, Hp[Hp.max() == Hp] + add_top))) / 2))
        if len(Hn) > 0:
            HBn = np.concat((Hn - np.diff(np.concat((Hn[Hn.min() == Hn] - remove_bottom, Hn))) / 2, [0]))
            ## separate bounds for output
            upper = np.concat((HBn[1:],   HBp[1:]))
            lower = np.concat((HBn[0:-1], HBp[0:-1]))
        else:
            ## separate bounds for output
            upper = HBp[1:]
            lower = HBp[0:-1]

    if not quiet:
        for i, y in enumerate(H):
            print("%10.1f  <<<  %10.1f  >>>  %10.1f " % (lower[i], H[i], upper[i]))

    ## return a simple array with heights and corresponding bounds
    return np.array((H, lower, upper))
