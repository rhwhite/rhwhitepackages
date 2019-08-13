"""
Created on Jun 4 2016

@author: rhwhite
"""

import os, errno
import numpy as np
import datetime as dt
import pandas
import xarray as xr
from scipy import stats

dayinmon = [31,28,31,30,31,30,31,31,30,31,30,31]
dayinseas = [90,91,92,92]
dayinyear = 365


def calcann(tstart,tend,var,years):
    ndims = var.shape
    nlons = ndims[-1]
    nlats = ndims[-2]
    ntimes = ndims[0]

    # First check if the data is monthly or daily
    res = getres(years,ntimes,365.0)

    # Check whether first year is a full year:
    yearstart = years[tstart]
    if res == "month":
        yearstart2 = years[tstart+11]
    elif res == "day":
        yearstart2 = years[tstart + 364]
    elif res == "6hour":
        yearstart2 = years[tstart + 364*4]
    elif res == "3hour":
        yearstart2 = years[tstart + 364*8]
    else:
        print "don't have a check for sub-daily data yet, just going to assume 6 hourly"
        yearstart2 = years[tstart + 365*4]  

    if yearstart == yearstart2:
        tstart = tstart
    else:
        print "finding beginning of new year"
        # Find beginning of new year
        while years[tstart] == yearstart:
            tstart = tstart + 1

    # Check whether last year is a full year:
    yearend = years[tend-1]
    if res == "month":
        yearend2 = years[tend-12]
    elif res == "day":
        yearend2 = years[tend - 360]
    elif res == "6hour":
        yearend2 = years[tend - 360*4]
    elif res == "3hour":
        yearend2 = years[tend - 360*8]
    else:
        print "don't have a check for sub-daily data yet, just going to assume 6 hourly"
        yearend2 = years[tend - 360*4]

    if yearend == yearend2:
        tend = tend
    else:
        print "finding beginning of last year"
        # Find beginning of new year
        while years[tend] == yearend:
            tend = tend -1

    startyear = years[tstart]
    endyear = years[tend-1]
    nyears = endyear - startyear + 1
    annvar = np.zeros([nyears,nlats,nlons],np.float)
    annyears = np.zeros(nyears,np.int)
    yearnow = -999
    inyearcount = 0

    for icount in range(tstart,tend):
        if years[icount] == yearnow:
            if res == "month":
                # If monthly need to multiply by dayinmon
                temp = temp + (var[icount,:,:] * dayinmon[inyearcount])
            else:
                temp = temp + var[icount,:,:]

            inyearcount = inyearcount + 1
        else:
            # if not the first time through
            if yearnow != -999:
                writeann(temp,res,yearnow-startyear,inyearcount,annvar)
                #reset inyearcount
                inyearcount = 0
                annyears[yearnow-startyear] = yearnow

            yearnow = years[icount]
            if res == "month":
                temp = var[icount,:,:] * dayinmon[inyearcount]
            else:
                temp = var[icount,:,:]
            inyearcount = inyearcount + 1

        if icount == tend-1:
            print "now reached the end of selection before start of next year, will write anyway if it's a full year"
            writeann(temp,res,yearnow-startyear,inyearcount,annvar) 
            annyears[yearnow-startyear] = yearnow
    print np.sum(annvar)
    return annvar,annyears

def writeann(tempvar,tres,iyear,inyearcount,annout):
    if tres == "month":
        if inyearcount != 12:
            print inyearcount
            print "not 12 months in year"
            print years[icount]
        else:
            annout[iyear,:,:] = tempvar / dayinyear

    elif tres == "day":
        if inyearcount < 364 or inyearcount > 365:
            print inyearcount
            print "not 365 days per year!"
        else:
            annout[iyear,:,:] = tempvar / inyearcount
    elif tres == "6hour":
        if inyearcount < 1456 or inyearcount > 1461:
            print inyearcount
            print "not 365 days per year!"
        else:
            annout[iyear,:,:] = tempvar / inyearcount
    elif tres == "3hour":
        if inyearcount < 2912 or inyearcount > 2928:
            print inyearcount
            print "not 365 days per year!"
        else:
            annout[iyear,:,:] = tempvar / inyearcount

    else:
        annout[iyear,:,:] = tempvar / inyearcount
    #reset inyearcount

def getres(inyears,intimes,mult):
    # if data is at least one year, mult should be 365
    # if data is less than one year, mult should be 28
        if inyears[0] != inyears[13] and inyears[13] != inyears[26]:
            return "month"
        else:
           dayspertime = mult /intimes

        #if inyears[0] != inyears[366] and inyears[366] != inyears[732]:
        if dayspertime < 1.1 and dayspertime > 0.9: 
            return "day"
        elif (dayspertime < 0.4) and (dayspertime > 0.2):
            return "6hour"
        elif (dayspertime < 0.15) and (dayspertime > 0.1):
            return "3hour"
        else:
            print dayspertime
            print "resolution less than 1 day - need to sort this!"
            sys.exit()

def calcseas(tstart,tend,var,years,months,seasonstarts,seasmonths):
    ndims = var.shape
    nlons = ndims[-1]
    nlats = ndims[-2]
    ntimes = ndims[0]

    nyears = years[-1] - years[0] + 1
    ntmonths = (nyears - 1) * 12 + months[-1] - months[0] + 1
    if ntmonths < 4:
        seasvar = np.zeros([1,1,nlats,nlons],np.float)
        seasyears = np.zeros(1,np.int)
    else:
        seasvar = np.zeros([nyears,4,nlats,nlons],np.float)
        seasyears = np.zeros(nyears,np.int)

    # First check if the data is monthly or daily
    res = getres(years,ntimes,ntmonths*365.0/12.0)
    # Check whether first season is a full season:
    monthstart = months[tstart]
    monthstart2 = -100
    seasoncounts = np.zeros(4)
    if monthstart in seasonstarts:
        if res == "month":
            monthstart2 = monthstart
        elif res == "day":
            monthstart2 = months[tstart + 28]
        elif res == "6hour":
            monthstart2 = months[tstart + 28*4]
        elif res == "3hour":
            monthstart2 = months[tstart + 28*8]
        else:
            print "don't have a check for sub-daily data yet, just going to assume 6 hourly"
            monthstart2 = months[tstart + 30*4]

        if monthstart == monthstart2:
            tstart = tstart
        else:
            print "finding beginning of new season"
                # Find beginning of new season
                # if partway through first month of season skip to at least the beginning of the next month
        if monthstart2 != -100: #starting in first season
            if res == "month":
                tstart = tstart + 1
            elif res == "day":
                tstart = tstart + 32
            else:
                tstart = tstart + (32*4)
        while months[tstart] not in seasonstarts:
            tstart = tstart + 1

        seasnow = -999
        inseascount = 0

    icount = tstart
    startyear = years[icount]
    nyears = np.max(years) - np.amin(years)

    if months[icount] == seasonstarts[0]:
        seasnow = 0
    elif months[icount] == seasonstarts[1]:
        seasnow = 1
    elif months[icount] == seasonstarts[2]:
        seasnow = 2
    elif months[icount] == seasonstarts[3]:
        seasnow = 3
    else:
        sys.exit("wrong month to start a season with!")

    if res == "month":
        temp = var[icount,:,:] * dayinmon[months[icount]-1]
    else:
        temp = var[icount,:,:]
    inseascount = inseascount + 1
    minseasoncount = np.min(seasoncounts)

    icount = icount + 1
    if nyears > 1:
        # if running multiple years, do this way to ensure equal number of seasons per year
        while minseasoncount < nyears:
            if months[icount] in seasmonths[seasnow]:
                if res == "month":
                    # If monthly need to multiply by dayinmon
                    temp = temp + (var[icount,:,:] * dayinmon[months[icount]-1])
                else:
                    temp = temp + var[icount,:,:]

                inseascount = inseascount + 1
            else:
                writeseas(temp,res,years[icount]-startyear,seasnow,inseascount,seasvar)
                seasyears[years[icount]-startyear] = years[icount]
                #reset inseascount and seasnow

                inseascount = 0
                seasoncounts[seasnow] = seasoncounts[seasnow] + 1
                minseasoncount = np.min(seasoncounts)
                if months[icount] == seasonstarts[0]:
                    seasnow = 0
                elif months[icount] == seasonstarts[1]:
                    seasnow = 1
                elif months[icount] == seasonstarts[2]:
                    seasnow = 2
                elif months[icount] == seasonstarts[3]:
                    seasnow = 3
                else:
                    sys.error("wrong month to start a season with!")

                if res == "month":
                    temp = var[icount,:,:] * dayinmon[months[icount]-1]
                else:
                    temp = var[icount,:,:]
                inseascount = inseascount + 1
            icount = icount+1

    else:
        # doing one season at a time
        for icount in range(1,ntimes):
            if months[icount] in seasmonths[seasnow]:
                if res == "month":
                    # If monthly need to multiply by dayinmon
                    temp = temp + (var[icount,:,:] * dayinmon[months[icount]-1])
                else:
                    temp = temp + var[icount,:,:]
                inseascount = inseascount + 1
            else:
                print "error: should be doing one season at a time"
                sys.exit()
            if icount == ntimes-1:
                print "reached the end of selection before start of next year, will write anyway if it's a full year"
                writeseas(temp,res,0,0,inseascount,seasvar)
                seasyears[0] = years[icount]

    return seasvar,seasyears

def writeseas(tempvar,tres,iyear,seasnow,inseascount,seasout):
    if tres == "month":
        if inseascount != 3:
            print inseascount
            print "not 3 months in season"
            print iyear,seasnow,inseascount
        else:
            seasout[iyear,seasnow,:,:] = tempvar / dayinseas[seasnow]

    elif tres == "day":
        if inseascount < 90 or inseascount > 93:
            print inseascount
            print "not 91 days per season in daily data!"
        else:
            seasout[iyear,seasnow,:,:] = tempvar / inseascount
    elif tres == "6hour":
        if inseascount < 360 or inseascount > 372:
            print inseascount
            print "not 91 days per season in 6 hourly data!"
        else:
            seasout[iyear,seasnow,:,:] = tempvar / inseascount
    elif tres == "3hour":
        if inseascount < 719 or inseascount > 744:
            print inseascount
            print "not 91 days per season in 3hourly data!"
        else:
            seasout[iyear,seasnow,:,:] = tempvar / inseascount




