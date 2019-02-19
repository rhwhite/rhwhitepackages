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

# Define new colormap
def get_colormaps():
    # Define new colormap with white centre
    reds = np.array([10,15,30,60,80,130,160,200,230,255,255,255,255,255,255,255,255,255,192,165],np.float)
    greens = np.array([50,75,110,160,180,210,240,250,255,255,255,250,232,192,160,96,50,20,0,0],np.float)
    blues = np.array([120,165,200,240,250,255,255,255,255,255,255,220,120,60,0,0,0,0,0,0],np.float)

    reds = reds/255
    blues = blues/255
    greens = greens/255

    ncolors = len(reds)

    spacing = np.linspace(0,1,ncolors)
    cdict = {'red': ((spacing[0], reds[0], reds[0]),
                     (spacing[1], reds[1], reds[1]),
                     (spacing[2], reds[2], reds[2]),
                     (spacing[3], reds[3], reds[3]),
                     (spacing[4], reds[4], reds[4]),
                     (spacing[5], reds[5], reds[5]),
                     (spacing[6], reds[6], reds[6]),
                     (spacing[7], reds[7], reds[7]),
                     (spacing[8], reds[8], reds[8]),
                     (spacing[9], reds[9], reds[9]),
                     (spacing[10], reds[10], reds[10]),
                     (spacing[11], reds[11], reds[11]),
                     (spacing[12], reds[12], reds[12]),
                     (spacing[13], reds[13], reds[13]),
                     (spacing[14], reds[14], reds[14]),
                     (spacing[15], reds[15], reds[15]),
                     (spacing[16], reds[16], reds[16]),
                     (spacing[17], reds[17], reds[17]),
                     (spacing[18], reds[18], reds[18]),
                     (spacing[19], reds[19], reds[19])
                    ),
             'green': ((spacing[0], greens[0], greens[0]),
                     (spacing[1], greens[1], greens[1]),
                     (spacing[2], greens[2], greens[2]),
                     (spacing[3], greens[3], greens[3]),
                     (spacing[4], greens[4], greens[4]),
                     (spacing[5], greens[5], greens[5]),
                     (spacing[6], greens[6], greens[6]),
                     (spacing[7], greens[7], greens[7]),
                     (spacing[8], greens[8], greens[8]),
                     (spacing[9], greens[9], greens[9]),
                     (spacing[10], greens[10], greens[10]),
                     (spacing[11], greens[11], greens[11]),
                     (spacing[12], greens[12], greens[12]),
                     (spacing[13], greens[13], greens[13]),
                     (spacing[14], greens[14], greens[14]),
                     (spacing[15], greens[15], greens[15]),
                     (spacing[16], greens[16], greens[16]),
                     (spacing[17], greens[17], greens[17]),
                     (spacing[18], greens[18], greens[18]),              
                     (spacing[19], greens[19], greens[19]) 
                      ),
                       
             'blue': ((spacing[0], blues[0], blues[0]),
                     (spacing[1], blues[1], blues[1]),
                     (spacing[2], blues[2], blues[2]),
                     (spacing[3], blues[3], blues[3]),
                     (spacing[4], blues[4], blues[4]),
                     (spacing[5], blues[5], blues[5]),
                     (spacing[6], blues[6], blues[6]),
                     (spacing[7], blues[7], blues[7]),
                     (spacing[8], blues[8], blues[8]),
                     (spacing[9], blues[9], blues[9]),
                     (spacing[10], blues[10], blues[10]),
                     (spacing[11], blues[11], blues[11]),
                     (spacing[12], blues[12], blues[12]),
                     (spacing[13], blues[13], blues[13]),
                     (spacing[14], blues[14], blues[14]),
                     (spacing[15], blues[15], blues[15]),
                     (spacing[16], blues[16], blues[16]),
                     (spacing[17], blues[17], blues[17]),
                     (spacing[18], blues[18], blues[18]),               
                     (spacing[19], blues[19], blues[19])                 
                     )}

    cmap_wc = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)

    reds = np.array([  10, 15, 30, 60, 80,130,160,200,230,245,255,255,255,255,255,255,255,192,165],np.float)
    greens = np.array([50, 75,110,160,180,210,240,250,255,245,250,232,192,160, 96, 50, 20,  0,  0],np.float)
    blues = np.array([120,165,200,240,250,255,255,255,255,245,220,120, 60,  0,  0,  0,  0,  0,  0],np.float)

    reds = reds/255
    blues = blues/255
    greens = greens/255

    ncolors = len(reds)

    spacing = np.linspace(0,1,ncolors)
    cdict = {'red': ((spacing[0], reds[0], reds[0]),
                     (spacing[1], reds[1], reds[1]),
                     (spacing[2], reds[2], reds[2]),
                     (spacing[3], reds[3], reds[3]),
                     (spacing[4], reds[4], reds[4]),
                     (spacing[5], reds[5], reds[5]),
                     (spacing[6], reds[6], reds[6]),
                     (spacing[7], reds[7], reds[7]),
                     (spacing[8], reds[8], reds[8]),
                     (spacing[9], reds[9], reds[9]),
                     (spacing[10], reds[10], reds[10]),
                     (spacing[11], reds[11], reds[11]),
                     (spacing[12], reds[12], reds[12]),
                     (spacing[13], reds[13], reds[13]),
                     (spacing[14], reds[14], reds[14]),
                     (spacing[15], reds[15], reds[15]),
                     (spacing[16], reds[16], reds[16]),
                     (spacing[17], reds[17], reds[17]),
                     (spacing[18], reds[18], reds[18])
                    ),
             'green': ((spacing[0], greens[0], greens[0]),
                     (spacing[1], greens[1], greens[1]),
                     (spacing[2], greens[2], greens[2]),
                     (spacing[3], greens[3], greens[3]),
                     (spacing[4], greens[4], greens[4]),
                     (spacing[5], greens[5], greens[5]),
                     (spacing[6], greens[6], greens[6]),
                     (spacing[7], greens[7], greens[7]),
                     (spacing[8], greens[8], greens[8]),
                     (spacing[9], greens[9], greens[9]),
                     (spacing[10], greens[10], greens[10]),
                     (spacing[11], greens[11], greens[11]),
                     (spacing[12], greens[12], greens[12]),
                     (spacing[13], greens[13], greens[13]),
                     (spacing[14], greens[14], greens[14]),
                     (spacing[15], greens[15], greens[15]),
                     (spacing[16], greens[16], greens[16]),
                     (spacing[17], greens[17], greens[17]),
                     (spacing[18], greens[18], greens[18])
                     ),

            'blue': ((spacing[0], blues[0], blues[0]),
                     (spacing[1], blues[1], blues[1]),
                     (spacing[2], blues[2], blues[2]),
                     (spacing[3], blues[3], blues[3]),
                     (spacing[4], blues[4], blues[4]),
                     (spacing[5], blues[5], blues[5]),
                     (spacing[6], blues[6], blues[6]),
                     (spacing[7], blues[7], blues[7]),
                     (spacing[8], blues[8], blues[8]),
                     (spacing[9], blues[9], blues[9]),
                     (spacing[10], blues[10], blues[10]),
                     (spacing[11], blues[11], blues[11]),
                     (spacing[12], blues[12], blues[12]),
                     (spacing[13], blues[13], blues[13]),
                     (spacing[14], blues[14], blues[14]),
                     (spacing[15], blues[15], blues[15]),
                     (spacing[16], blues[16], blues[16]),
                     (spacing[17], blues[17], blues[17]),
                     (spacing[18], blues[18], blues[18]) 
                     )}

    cmap_def = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)

    # Define new colormap
    reds = np.array([ 245,255,255,255,255,255,255,255,192,165],np.float)
    greens = np.array([245,250,232,192,160,96,50,20,0,0],np.float)
    blues = np.array([245,220,120,60,0,0,0,0,0,0],np.float)

    reds = reds/255
    blues = blues/255
    greens = greens/255

    ncolors = len(reds)

    spacing = np.linspace(0,1,ncolors)
    cdict = {'red': ((spacing[0], reds[0], reds[0]),
                     (spacing[1], reds[1], reds[1]),
                     (spacing[2], reds[2], reds[2]),
                     (spacing[3], reds[3], reds[3]),
                     (spacing[4], reds[4], reds[4]),
                     (spacing[5], reds[5], reds[5]),
                     (spacing[6], reds[6], reds[6]),
                     (spacing[7], reds[7], reds[7]),
                     (spacing[8], reds[8], reds[8]),
                     (spacing[9], reds[9], reds[9])
                    ),
             'green': ((spacing[0], greens[0], greens[0]),
                     (spacing[1], greens[1], greens[1]),
                     (spacing[2], greens[2], greens[2]),
                     (spacing[3], greens[3], greens[3]),
                     (spacing[4], greens[4], greens[4]),
                     (spacing[5], greens[5], greens[5]),
                     (spacing[6], greens[6], greens[6]),
                     (spacing[7], greens[7], greens[7]),
                     (spacing[8], greens[8], greens[8]),
                     (spacing[9], greens[9], greens[9])
                     ),
                       
             'blue': ((spacing[0], blues[0], blues[0]),
                     (spacing[1], blues[1], blues[1]),
                     (spacing[2], blues[2], blues[2]),
                     (spacing[3], blues[3], blues[3]),
                     (spacing[4], blues[4], blues[4]),
                     (spacing[5], blues[5], blues[5]),
                     (spacing[6], blues[6], blues[6]),
                     (spacing[7], blues[7], blues[7]),
                     (spacing[8], blues[8], blues[8]),
                     (spacing[9], blues[9], blues[9]) 
                     )}

    cmap_reds = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return(cmap_def,cmap_wc,cmap_reds)

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


