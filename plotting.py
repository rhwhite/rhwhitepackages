"""
Functions to help with plotting in matplotlib

Author: Rachel White, rachelwhite@cantab.net
Creation Date: 12 Oct 2017

"""
import mpl_toolkits
import cartopy.crs as ccrs
import os, errno
import netCDF4
import numpy as np
import datetime as dt
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib

def shiftlons(invar,nlons):
    nlonhalf = nlons/2

    newinvar = np.zeros(invar.shape,np.float)
    newinvar[...,0:nlonhalf] = invar[...,nlonhalf:nlons]
    newinvar[...,nlonhalf:nlons] = invar[...,0:nlonhalf]

    return newinvar

def shiftlonlons(inlon,nlons):
    nlonhalf = nlons/2

    lonsnew = np.zeros(inlon.shape,np.float)
    lonsnew[0:nlonhalf] = inlon[nlonhalf:nlons] - 360.0
    lonsnew[nlonhalf:nlons] = inlon[0:nlonhalf]

    return lonsnew

def drawmap():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_xticks(range(-180,180,60),crs=ccrs.PlateCarree())
    ax.set_yticks(range(-90,90,30),crs=ccrs.PlateCarree())

    return ax

def drawmap_Moll():
    ax = plt.axes(projection=ccrs.Mollweide())
    ax.coastlines()
    ax.set_xticks(range(-180,180,60),crs=ccrs.Mollweide())
    ax.set_yticks(range(-90,90,30),crs=ccrs.Mollweide())

    return ax

def drawmap_sub(ax,proj=ccrs.PlateCarree()):
    ax.coastlines()
    ax.set_xticks(range(-180,180,60),crs=proj)
    ax.set_yticks(range(-90,90,30),crs=proj)

    return ax

def drawpolarmap(pole='N'):
    if pole == 'N':
        central_lat = 90
    else:
        central_lat = -90
    ax = plt.axes(projection=ccrs.Orthographic(central_latitude=central_lat))
    ax.coastlines()
    ax.gridlines()

    return ax


