"""
Functions to manipulate CESM files

Created 5 Oct 2017

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
from rhwhitepackages.griddata_functions import *

calendar.monthrange(2001,1)[1]
monthstart = [0]
for imonth in range(1,13):
    ndays = calendar.monthrange(2001,imonth)[1]
    monthstart.append(monthstart[imonth-1]+ndays)

def smooth_clim_from_daily(file,ndaysyear,yearin,typein):
    nmonths = 36
    startmonth = 12
    ndays = np.arange(0,ndaysyear)
    Uclim_daily = np.zeros([ndaysyear,
                            Uclim_month.shape[1],
                            Uclim_month.shape[2]])

    monthstart = 0
    for imonth in range(0,12):
        month1 = imonth + 1
        nday1 = cal.monthrange(yearin,month1)[1]

        Uclim_daily[monthstart:monthstart+nday1] = Uclim_month[imonth].values
        monthstart = monthstart + nday1

    # Add on previous year
    for iday in range(0,ndaysyear):
        ndays = np.insert(ndays,0, -1 - iday,axis=0)
    # Add on next years months
    for iday in range(0,ndaysyear):
        ndays = np.append(ndays,ndaysyear + iday + 1)

    Uclim_cont = np.concatenate([Uclim_daily, Uclim_daily, Uclim_daily],0)

    f = interp1d(ndays,Uclim_cont,kind=typein,axis=0)

    days = np.arange(0,ndaysyear)
    return f(days)


def lat_mean(invar):
    latrad = np.deg2rad(invar.lat)
    coslat = np.cos(latrad)
    varnew = invar * coslat
    varmean = varnew.mean(dim='lat')/coslat.mean(dim='lat')
    return varmean

def lat_sum(invar):
    latrad = np.deg2rad(invar.lat)
    coslat = np.cos(latrad)
    varnew = invar * coslat
    varsum = varnew.mean(dim='lat')/coslat.sum(dim='lat')
    return varsum

def getlatdist(inlat1,inlat2):
    circum = 2 * pi * rearth
    inlat1rad = np.deg2rad(inlat1)
    inlat2rad = np.deg2rad(inlat2)
    dist = circum * (inlat1rad - inlat2rad)/(2 * pi)
    return(abs(dist))

def getlondist(inlon1,inlon2,inlat):
    inlatrad = np.deg2rad(inlat)
    inlon1rad = np.deg2rad(inlon1)
    inlon2rad = np.deg2rad(inlon2)
    circum = 2 * pi * rearth * np.cos(inlatrad)
    dist = circum * (inlon1rad - inlon2rad)/(2 * pi)
    return(abs(dist))

def get_latname(invar):
    if lat in nvar.dims:
       latvar = 'lat' 
    elif latitude in invar.dims:
       latvar = 'latitude'
    elif 'Latitude' in nvar.dims:
       latvar = 'Latitude'
    else:
        sys.exit('cannot determine latitude dimension name')
    return(latvar)

def get_circum_area(latstart,latend):
    circumlat = 2.0 * np.pi * rearth * np.cos(np.deg2rad(latstart))
    return(circumlat)

def get_area(latstart,latend,lonstart,lonend):
    dlon = lonend - lonstart
    dlonrad = dlon * (2.0 * np.pi)/360.0
    latstartrad = latstart * (2.0 * pi)/360.0
    latendrad = latend * (2.0 * pi)/360.0
    area = rearth**2 * abs(np.sin(latstartrad) - np.sin(latendrad)) * dlonrad
    return(area)


def seasmean(infile,seas):
    calendar.monthrange(2001,1)[1]
    monthstart = [0]
    for imonth in range(1,13):
        ndays = calendar.monthrange(2001,imonth)[1]
        monthstart.append(monthstart[imonth-1]+ndays)

    if months == 'DJF':
        indices = np.squeeze(np.argwhere((infile.time.values % 365 <= 59) |
            (infile.time.values % 365 > 334)))
        #monthtitle = 'DJF'
    monthselect = infile.isel(time=indices)

    # find first month
    months = (monthselect.time.values % 365)

def selectmonths(infile,months):

    calendar.monthrange(2001,1)[1]
    monthstart = [0]
    for imonth in range(1,13):
        ndays = calendar.monthrange(2001,imonth)[1]
        monthstart.append(monthstart[imonth-1]+ndays)

    if months == 'DJF':
        indices = np.squeeze(np.argwhere(((infile.time.values.astype(int)) % 365 <= 59) |
            (infile.time.values.astype(int) % 365 >= 335)))

    elif months == 'MAM':
        indices = np.squeeze(np.argwhere((infile.time.values.astype(int) % 365 > 59) &
            (infile.time.values.astype(int) % 365 <= 151)))

    elif months == 'JJA':
        indices = np.squeeze(np.argwhere((infile.time.values.astype(int) % 365 > 151) &
            (infile.time.values.astype(int) % 365 <= 243)))

    elif months == 'SON':
        indices = np.squeeze(np.argwhere((infile.time.values.astype(int) % 365 > 243) &
            (infile.time.values.astype(int) % 365 <= 334)))


    elif months == 'NDJFM':
        indices = np.squeeze(np.argwhere((infile.time.values.astype(int) % 365 <= 90) |
            (infile.time.values.astype(int) % 365 > 304)))
        #monthtitle = 'NDJFM'

    elif months =='MJJAS':
        indices = np.squeeze(np.argwhere((infile.time.values.astype(int) % 365 > 120) &
            (infile.time.values.astype(int) % 365 <= 273)))
        #monthtitle = 'MJJAS'

    elif months in range(1,13):
        indices = np.squeeze(np.argwhere((infile.time.values.astype(int) % 365 >=
            monthstart[months-1]) & (infile.time.values.astype(int) % 365 < monthstart[months])))
        #monthtitle = calendar.month_abbr[months]
    else:
        print('month ' + months + 'is not recognized, please try again')
        return ()
    return infile.isel(time=indices)


def get_all_seasons(infile):
    outfile = {}
    for iseas in ['DJF','MAM','JJA','SON']:
        outfile[iseas] = selectmonths(infile,iseas)
    return outfile

def concatmonths(infile,month):
    # check data is in days since Jan 01
    if infile.time.units[0:10] not in ['days since','Days since']:
        sys.exit('didn\'t understand time units, or not in days since...')
    if infile.time.units[16:21] != '01-01':
        print infile.time.units
        sys.exit('didn\'t understand time units, or data don\'t start on Jan 1st')

    # get number of days and first day
    ndays = infile.time.shape[0]
    firstday = infile.time[0]

    # check data is daily
    if infile.time[1]-infile.time[0] != 1:
        sys.exit('not daily data')

    # get first year present in data
    mstart = monthstart[month-1]
    mend = monthstart[month]-1

    while mend < firstday:
        mstart += 365
        mend += 365

    out = Z500['CTL'].sel(time=slice(mstart,mend))

    # get following years
    mstart += 365
    mend += 365
    while mend <= ndays:
        test2 = Z500['CTL'].sel(time=slice(mstart,mend))
        #out = out.combine_first(test2)
        out = xr.concat((test,test2),dim='time')
        mstart += 365
        mend += 365

    return(out)


