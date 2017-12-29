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

def drawmap():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_xticks(range(-180,180,60),crs=ccrs.PlateCarree())
    ax.set_yticks(range(-90,90,30),crs=ccrs.PlateCarree())

    return ax

def drawmap_sub(ax):
    ax.coastlines()
    ax.set_xticks(range(-180,180,60),crs=ccrs.PlateCarree())
    ax.set_yticks(range(-90,90,30),crs=ccrs.PlateCarree())

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


