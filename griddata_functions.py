
"""
Functions on gridded data files

Created 13 Aug 2019

Author: Rachel White rachel.white@cantab.net
"""

import os, errno
import netCDF4
import numpy as np
import datetime as dt
import pandas as pd
import xarray as xr
import calendar
from rhwhitepackages.CESMconst import *


def ddy(invar):
    # based on https://www.ncl.ucar.edu/Document/Functions/Built-in/uv2dv_cfd.shtml
    # H.B. Bluestein [Synoptic-Dynamic Meteorology in Midlatitudes, 1992, Oxford Univ. Press p113-114]

    # Using this rather ugly approach to allowing different names for the latitude variable
    # Because of the way I have written the xarray code to use variable names, names are hard-coded in places
    try:
        nlats = len(invar['lat'])
        latname = 'lat'
    except KeyError:
        nlats = len(invar['latitude'])
        latname = 'latitude'

    dims = invar.shape
    if len(dims) == 1:
        # treat array values differently if it is a 1D array
        if latname == 'lat':
            ddy = invar.copy(deep=True)
            ddy[:] = np.nan

            ddy.isel(lat=0).values = 0

            for ilat in range(1,nlats-1):

                dy = getlatdist(invar['lat'].isel(lat = ilat+1),
                                    invar['lat'].isel(lat = ilat-1))
                dvar = invar.isel(lat = ilat+1) - invar.isel(lat = ilat-1)
                ddy.isel(lat = ilat).values = dvar/dy - (invar.isel(lat=ilat)/rearth) * np.tan(np.deg2rad(invar.lat.isel(lat=ilat)))

            ddy.isel(lat=nlats-1).values = 0

        elif latname == 'latitude':
            ddy = invar.copy(deep=True)
            ddy[:] = np.nan
            ddy[0] = 0

            for ilat in range(1,nlats-1):
                dy = getlatdist(invar['latitude'].isel(latitude = ilat+1),
                                    invar['latitude'].isel(latitude = ilat-1))
                dvar = invar.isel(latitude = ilat+1) - invar.isel(latitude = ilat-1)
                ddy[ilat] = dvar/dy - (invar.isel(latitude=ilat)/rearth) * np.tan(np.deg2rad(invar.latitude.isel(latitude=ilat)))

            ddy[nlats-1] = 0
    else:
        # more than 1D array:
        
        if latname == 'lat':
            ddy = invar.copy(deep=True)
            ddy.isel(lat=0)[...] = 0

            for ilat in range(1,nlats-1):

                dy = getlatdist(invar['lat'].isel(lat = ilat+1),
                                    invar['lat'].isel(lat = ilat-1))
                dvar = invar.isel(lat = ilat+1) - invar.isel(lat = ilat-1)
                ddy.isel(lat = ilat)[...] = dvar/dy - (invar.isel(lat=ilat)/rearth) * np.tan(np.deg2rad(invar.lat.isel(lat=ilat)))

            ddy.isel(lat=nlats-1)[...] = 0

        elif latname == 'latitude':
            ddy = invar.copy(deep=True)
            ddy[...] = np.nan
            ddy.isel(latitude=0)[...] = 0

            for ilat in range(1,nlats-1):
                dy = getlatdist(invar['latitude'].isel(latitude = ilat+1),
                                    invar['latitude'].isel(latitude = ilat-1))
                dvar = invar.isel(latitude = ilat+1) - invar.isel(latitude = ilat-1)
                ddy.isel(latitude = ilat)[...] = dvar/dy - (invar.isel(latitude=ilat)/rearth) * np.tan(np.deg2rad(invar.latitude.isel(latitude=ilat)))

            ddy.isel(latitude=nlats-1)[...] = 0
    return(ddy)


def ddphi(invar):
    # Using this rather ugly approach to allowing different names for the latitude variable
    # Because of the way I have written the xarray code to use variable names, names are hard-coded in places

    try:
        nlats = len(invar['lat'])
        latname = 'lat'
    except KeyError:
        nlats = len(invar['latitude'])
        latname = 'latitude'


    if latname == 'lat':
        dims_var = invar.dims
        latidx_var = dims_var.index('lat')

        dims_lat = invar.lat.dims
        latidx_lat = dims_lat.index('lat')

        ddphi = invar.copy(deep=True)
        ddphi.isel(lat=0)[...] = 0

        dvar = np.gradient(invar,axis=latidx_var)
        dlat = np.gradient(np.deg2rad(invar.lat),axis=latidx_lat)
        if len(dims_var) > 1:
            dvardphi = (dvar/dlat[:,None]) - (invar/rearth) * np.tan(np.deg2rad(invar.lat))
        else:
            dvardphi = dvar/dlat - (invar/rearth) * np.tan(np.deg2rad(invar.lat))[:,None]

    elif latname == 'latitude':
        dims_var = invar.dims
        latidx_var = dims_var.index('latitude')

        dims_lat = invar.latitude.dims
        latidx_lat = dims_lat.index('latitude')

        ddphi = invar.copy(deep=True)
        ddphi.isel(latitude=0)[...] = 0

        dvar = np.gradient(invar,axis=latidx_var)
        dlat = np.gradient(np.deg2rad(invar.latitude),axis=latidx_lat)

        if len(dims_var) > 1:
            dvardphi = (dvar/dlat[:,None]) - (invar/rearth) * np.tan(np.deg2rad(invar.latitude))
        else:
            dvardphi = (dvar/dlat) - (invar/rearth) * np.tan(np.deg2rad(invar.latitude))


    return(dvardphi)

def ddx(invar):
    # based on https://www.ncl.ucar.edu/Document/Functions/Built-in/uv2dv_cfd.shtml
    # H.B. Bluestein [Synoptic-Dynamic Meteorology in Midlatitudes, 1992, Oxford Univ. Press p113-114]
    nlons = len(invar['lon'])
    nlats = len(invar['lat'])
    ddx = invar.copy(deep=True)
    #ddx.isel(lat=0)[...] = 0

    for ilat in range(1,nlats-1):
        for ilon in range(1,nlons-1):
            dx = getlondist(invar['lon'].isel(lon = ilon+1),
                            invar['lon'].isel(lon = ilon-1),
                            invar['lat'].isel(lat=ilat))
            dvar = invar.isel(lon = ilon+1).isel(lat=ilat) - invar.isel(lon = ilon-1).isel(lat=ilat)
            ddx.isel(lat = ilat)[...] = dvar/dx

        # Now cover ilon = 0 and ilon = nlons-1
        ilon = 0
        dx = getlondist(invar['lon'].isel(lon = ilon+1),
                            invar['lon'].isel(lon = nlons-1),
                            invar['lat'].isel(lat=ilat))
        dvar = invar.isel(lon = ilon+1).isel(lat=ilat) - invar.isel(lon = nlons-1).isel(lat=ilat)
        ddx.isel(lat = ilat)[...] = dvar/dx

        ilon= nlons-1
        dx = getlondist(invar['lon'].isel(lon = 0),
                            invar['lon'].isel(lon = ilon),
                            invar['lat'].isel(lat=ilat))
        dvar = invar.isel(lon = 0).isel(lat=ilat) - invar.isel(lon = ilon-1).isel(lat=ilat)
        ddx.isel(lat = ilat)[...] = dvar/dx

    return(ddx)
