# -*- coding: utf-8 -*-
"""
Functions for parameters calculation

@author: thanasisn
"""

import metpy.calc
import xarray as xr
from   metpy.units import units


def z_to_height(DT):
    """
    Return a new xarray where the height is calculated from the
    goeopotential (z variable).

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
    return(DT)


def border_up(limit, step):
    """
    Return the next coordinate that satisfy the given step,
    going to the positive direction.
    Usually for North and East boundaries.

    Parameters
    ----------
    limit : float
        A limit for the coordinate.
    step : TYPE
        A coordinate limit with the step/resolution we need.

    Returns
    -------
    A value that includes the limit and may be bigger to meet the next step value

    >>> border_up(-11, 5)
    -10
    >>> border_up( 11, 5)
    15
    """
    return limit + step - limit % step


def border_down(limit, step):
    """
    Return the next coordinate that satisfy the given step,
    going to the negative direction.
    Usually for South and West boundaries.

    Parameters
    ----------
    limit : float
        A limit for the coordinate.
    step : TYPE
        A coordinate limit with the step/resolution we need.

    Returns
    -------
    A value that includes the limit and may be smaller to meet the previous step value.

    >>> border_down(-11, 5)
    -15
    >>> border_down( 11, 5)
    10
    """
    return limit - limit % step
