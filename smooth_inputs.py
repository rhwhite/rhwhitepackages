"""
Functions to interpolate SST and sea-ice fields over land, in case
of any discrepancies in the land-sea mask 

Designed for CESM input data format

Created 20 June 2018

Author: Rachel White rachel.white@cantab.net
"""

import os, errno
import netCDF4
import numpy as np
import datetime as dt
import pandas as pd
import xarray as xr
import math
from scipy import stats
from rhwhitepackages.readwrite import *

# Define functions for interpolating without mask
def findnextland_pop(SST,ilat,ilon,imonth,nlon):
    if ilon >= nlon: return(nlon)
    while (SST[ilon] != -1000):
        ilon += 1
        if ilon >= nlon:
            return(nlon)
    return(ilon)

def findprevland_pop(SST,ilat,ilon,imonth,nlon):
    while (SST[ilon] != -1000):
        ilon -= 1
        if ilon <= -nlon+1:
            break
    return(ilon)

def findnextocean_pop(SST,ilat,ilon,imonth,nlon):
    if ilon >= nlon-1: return(0)
    while (SST[ilon] == -1000):
        ilon += 1
        if ilon >= nlon-1:
            ilon = 0
            break
    return(ilon)

def findprevocean_pop(SST,ilat,ilon,imonth,nlon):
    while (SST[ilon] == -1000):
        ilon -= 1
        if ilon <= -nlon+1:
            break
    return(ilon)

def replacetemps_pop(ocnstart,ocnend,landstart,landend,imonth,ilat,ssts_in):
    starttemp = ssts_in[imonth,ilat,ocnstart].values
    endtemp = ssts_in[imonth,ilat,ocnend].values
    numlons = landend - landstart + 1 # each side matches ocean
    interptemp = np.linspace(starttemp,endtemp,numlons)
    return(interptemp)

def smooth_SSTs(datain,K = True):
    data = datain.copy(deep=True)

    # Maybe need this depending on the type of input data:
    if K:
        data = data - 273.15

    # Mask data with various ways to make sure all invalid data is masked
    ssts.values= np.where(maskin >= 1,-1000.0,ssts.values)

    data.values = np.where(data.values < -100,np.nan,data.values)
    data.values = np.where(data.values < -1.8, -1.8, data.values)

    data.values = np.where(np.isnan(data),-1000,data.values)
    data_new = data.copy(deep=True).load()

    nlat = data.lat.shape[0]
    nlon = data.lon.shape[0]

    for imonth in range(0,12):
        # run simple case: ocean, then land, then ocean
        # deal with latitudes that start with land.
        for ilat in range(0,nlat):
            ilon = 0
            if data_new[imonth,ilat,0] == -1000:
                ocnstart = findprevocean_pop(
                            data_new[imonth,ilat,:],ilat,nlon-1,imonth,nlon)
                ocnend = findnextocean_pop(
                            data_new[imonth,ilat,:],ilat,0,imonth,nlon)
                landstart = ocnstart + 1
                landend = ocnend - 1

                ilon = ocnend

                if ocnend > 0 and ocnend < nlon-1:
                    starttemp = data_new[imonth,ilat,max(0,landstart-1)].values
                    endtemp = data_new[imonth,ilat,min(landend+1,nlon-1)].values
                    nend = nlon-landstart + 1

                    numlons = landend + nend + 1
                    interptemp = np.linspace(starttemp,endtemp,numlons)

                    data_new[imonth,ilat,landstart:nlon] = interptemp[0:nend-1]
                    data_new[imonth,ilat,0:landend+1] = interptemp[nend-1:-1]
                else:
                    data_new[imonth,ilat,:] = -1.8
                    ilon = nlon
            # now continue with first land point after ocean

            while ilon < nlon: # -1
                landstart = findnextland_pop(
                                    data_new[imonth,ilat,:],ilat,ilon,imonth,nlon)
                if landstart != nlon:
                    ocnstart = landstart - 1
                    ocnend = findnextocean_pop(
                                    data_new[imonth,ilat,:],ilat,landstart+1,imonth,nlon)
                    landend = ocnend - 1

                    if landend == -1: landend = nlon-1

                    ilon = ocnend
                    interptemp = replacetemps_pop(ocnstart,ocnend,landstart,landend,imonth,ilat,data_new)
                    data_new[imonth,ilat,landstart:landend+1] = interptemp
                else:

                    # go into next latitude loop
                    break
            #ilon == nlon - 1
            #if data_new[month,ilat,ilon] == -1000.0:
            #    starttemp = data_new[month,ilat,ilon]

    DATA_new = data_new.to_dataset(name='SST_cpl')
    return(DATA_new)

def smooth_CICE(datain,K = True):
    # Smooth ice field
    ice_new = datain

    ice_new.values = np.where(maskin>0,-1000,ice_new.values)
    ice_new.values = np.where(ice_new<-10,-1000,ice_new.values)
    ice_new.values = np.where(np.isnan(ice_new),-1000,ice_new.values)

    nlat = ice_new.lat.shape[0]
    nlon = ice_new.lon.shape[0]


    # Now run interpolation
    for imonth in range(0,12):
        print imonth 
        # run simple case: ocean, then land, then ocean
        # deal with latitudes that start with land.
        for ilat in range(0,nlat):
            ilon = 0

            if ice_new[imonth,ilat,0] == -1000:
                ocnstart = findprevocean_pop(
                            ice_new[imonth,ilat,:],ilat,nlon-1,imonth,nlon)
                ocnend = findnextocean_pop(
                            ice_new[imonth,ilat,:],ilat,0,imonth,nlon)
                landstart = ocnstart + 1
                landend = ocnend - 1

                ilon = ocnend

                if ocnend > 0 and ocnend < nlon-1:
                    starttemp = ice_new[imonth,ilat,max(0,landstart-1)].values
                    endtemp = ice_new[imonth,ilat,min(landend+1,nlon-1)].values
                    nend = nlon-landstart + 1

                    numlons = landend + nend + 1
                    interptemp = np.linspace(starttemp,endtemp,numlons)

                    ice_new[imonth,ilat,landstart:nlon] = interptemp[0:nend-1]
                    ice_new[imonth,ilat,0:landend+1] = interptemp[nend-1:-1]
                else:
                    ice_new[imonth,ilat,:] = 1.0
                    ilon = nlon
            # now continue with first land point after ocean

            while ilon < nlon: #-1:
                landstart = findnextland_pop(
                                    ice_new[imonth,ilat,:],ilat,ilon,imonth,nlon)

                if landstart != nlon:
                    ocnstart = landstart - 1
                    ocnend = findnextocean_pop(
                                    ice_new[imonth,ilat,:],ilat,landstart+1,imonth,nlon)
                    landend = ocnend - 1

                    if landend == -1: landend = nlon-1

                    ilon = ocnend
                    interptemp = replacetemps_pop(ocnstart,ocnend,landstart,landend,imonth,ilat,ice_new)
                    ice_new[imonth,ilat,landstart:landend+1] = interptemp
                else:
                    # go into next latitude loop
                    break

    ICE_new = np.where(maskin == 1.0,ice_new.values, datain.values)
    ICE_new = ice_new.to_dataset(name='ice_cov')
    return(ICE_new)


