#-*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os, errno
import numpy as np
import numpy.ma as ma
import netCDF4
from netCDF4 import Dataset
import datetime as dt
import re
import sys
import Ngl
import xarray as xr
import math
import argparse
import resource

# simple function to 'regrid' data by summing conservatively over 
# consecutive gridboxes
def conserveRegrid(indata,latin,lonin,nsumlat,nsumlon):
    lats = indata[latin].values
    lons = indata[lonin].values

    nlons = len(lons)
    nlats = len(lats)

    nlonsnew = math.floor(nlons/nsumlon)
    nlatsnew = math.floor(nlats/nsumlat)

    nlats2 = np.int(nlatsnew * nsumlat)
    nlons2 = np.int(nlonsnew * nsumlon)

    latsnew = np.zeros(nlatsnew,np.float)
    lonsnew = np.zeros(nlonsnew,np.float)

    indatadims = indata.dims
    newdimsize = []
    nnewdims = 0
    countdim=0
    newdims=[]
    for idim in indatadims:
        if idim == latin:
            latdimnum = countdim
            newdims.append('lat')
        elif idim == lonin:
            londimnum = countdim
            newdims.append('lon')
        else:
            nnewdims = nnewdims + 1
            newdimsize.append(len(indata[idim]))
            newdims.append(idim)
        countdim += 1

    nnewdims = len(newdims)
    newdimsize.append(nlatsnew)
    newdimsize.append(nlonsnew)

    # Catch all silly values as these will no longer show up as nans if
    # arrived with normal values

    indata[np.where(indata > 10E30)] = np.nan
    indata[np.where(indata == -9999)] = np.nan
    indata[np.where(indata < -10E30)] = np.nan

    new = np.zeros((newdimsize),np.float)

    #'regrid' by averaging over large boxes (averaging lons and lats)

    inlat = 0
    for ilats in range(0,nlats2,nsumlat):
        inlon = 0
        latsnew[inlat] = np.mean(lats[ilats:ilats+nsumlat])
        for ilons in range(0,nlons2,nsumlon):
            lonsnew[inlon] = np.mean(lons[ilons:ilons+nsumlon])
            new[...,inlat,inlon] = np.mean(
                                indata[...,ilats:ilats+nsumlat,ilons:ilons+nsumlon],
                                axis=(latdimnum,londimnum),dtype=np.float)

            inlon += 1
        inlat += 1

    if nnewdims == 2:
        try:
            newDA = xr.DataArray(new,coords=[('lat',latsnew),
                                               ('lon',lonsnew)],
                                       attrs=[indata.attrs])
        except ValueError:
            newDA = xr.DataArray(new,coords=[('lat',latsnew),
                                               ('lon',lonsnew)])
    elif nnewdims == 3:
        try:
            newDA = xr.DataArray(new,coords=[(newdims[0],indata[newdims[0]]),
                                               ('lat',latsnew),
                                               ('lon',lonsnew)],
                                       attrs=[indata.attrs])
        except ValueError:
            newDA = xr.DataArray(new,coords=[(newdims[0],indata[newdims[0]]),
                                               ('lat',latsnew),
                                               ('lon',lonsnew)])

    return(newDA)

    ncfile.close()

