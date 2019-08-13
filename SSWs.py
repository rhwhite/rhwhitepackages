# Module to search for and get data on SSWs
# Using the definition of Charlton and Polvani (2007):

# Author rachel.white@cantab.net

# Created July 2017

import numpy as np
import xarray as xr
import math
import sys

def adddays(U,itime,ndays):
    # Find ndays consecutive days with easterlies
    numcons = 0
    torun = True
    while itime < len(U.time):
        if U[itime] > 0:
            numcons += 1
        else:
            numcons = 0
        if numcons >= ndays:
            return(itime,numcons,False)
        itime += 1
    return(itime,numcons,True)

def meanSE(N,in0,in1,in2,in3=0):
    # Calculate mean and standard error of number of SSWs
    # a la Charlton and Polvani (2007)
    p0 = float(in0)/float(N)
    p1 = float(in1)/float(N)
    p2 = float(in2)/float(N)
    p3 = float(in3)/float(N)

    calcmean = p1 + (2 * p2) + (3 * p3)
    calcSE = ((math.sqrt(((0-calcmean)**2 * p0) +
                        ((1-calcmean)**2 * p1) +
                        ((2-calcmean)**2 * p2) + 
                        ((3-calcmean)**2 * p3)))
                /math.sqrt(N))
    return calcmean,calcSE

def findyearSSWs(U,times,count,thresh,lastdate,startdate,toprint,SSWdates):
    # find all SSWs in a single year

    finalwarmingstart = -1
    yearcount = 0
    itime = 0
    # if U starts below 0, iterate until it isn't!
    while U[itime]<0:
        itime +=1
    while itime < len(U.time):
        if U[itime] < 0:
            central,end,itime = findend(U,itime,thresh)
            if end == -1:
                finalwarmingstart = ((times[central]+1) % 365)
            else:
                SSWdates.append(int(times[central]))
                if toprint: print ('SSW, day of year ' +
                                str((times[central]) % 365))
                if lastdate < ((times[central] +1) % 365) < startdate :
                    # it counts as a final warming
                    finalwarmingstart = ((times[central]+1) % 365)
                else:
                    count +=1
                    yearcount +=1

        itime +=1
    return count,yearcount, finalwarmingstart, SSWdates

def findend(U,itime,thresh):
    # find final SSW
    centraltime,endtime = -1,-1
    if U[itime] < 0:
        centraltime = itime

    # Find end date
    while U[itime] < 0:
        itime = itime + 1
        if itime >= len(U.time): return (centraltime,-1,itime)
    endtime = itime

    # Check for final warming: ends after April 30th but started before July
    # Add 10 consective easterly days - must occur before April 30th for event
    # to count 

    newtime,numcons,end = adddays(U,itime,thresh)

    if end:
        return(itime,-1,newtime)
    else:
        # Event counts. Now add 20 consecutive days
        itime,ndays,end = adddays(U,itime,20)
        return(centraltime,endtime,itime)



def findSSWs(U,thresh,Obs=False,startyr = 0):
    # Find SSWs, print the mean number, the standard error, and 
    # return the dates
    # Created for WACCM daily data

    SSWdates = []
    toprint = False
    SSWyears = []
    startdate = 303 # beginning of November
    lastdate = 119 # end of April 
    enddate = 119 # 30th April

    count = 0
    yearcount = 0
    singleyear = 0
    doubleyear = 0
    tripleyear = 0
    final = []
    nyears = len(U.time)//365
    times = U.time

    # Select first year
    if Obs:
        yearU = U.sel(time=slice(str(startyr) + '-01',str(startyr) + '-04'))
        yeartime = times.sel(time=slice(str(startyr) + '-01',
                                            str(startyr) +'-04'))
        yeartime = (yeartime.values - np.datetime64('1980-01-01'))/ np.timedelta64(1, 'D')
    else:
        yearU = U.isel(time=slice(0,120))
        yeartime = times[0:120].values

    count,yearcount,finalW,SSWdates = findyearSSWs(yearU,yeartime,count,thresh,
                                                lastdate,startdate,
                                                toprint, SSWdates)
    if yearcount == 1:
        singleyear +=1
        #if toprint: print('year 0 1 SSW \n')
        SSWyears.append(0)
    elif yearcount ==2:
        doubleyear +=1
        #if toprint: print('year 0 2 SSWs \n')
        SSWyears.append(0)
    elif yearcount ==3:
        tripleyear +=1
        SSWyears.append(0)
    final.append(finalW)

    for iyear in range(0,nyears):
        if Obs:
            yearU = U.sel(time=slice(str(startyr+iyear) +'-11',
                                    str(startyr+iyear+1) + '-04'))
            yeartime = times.sel(time=slice(str(startyr+iyear) + '-11',
                                            str(startyr+iyear+1) +'-04'))
            yeartime = ((yeartime.values - np.datetime64('1980-01-01'))/
                                    np.timedelta64(1, 'D'))

        else:
            yearU = U.isel(time=slice(startdate+(iyear*365),
                                    enddate + ((iyear + 1) * 365)))
            yeartime = (times[startdate+(iyear*365):
                             enddate+((iyear+1)*365)].values)

        count,yearcount,finalW,SSWdates = findyearSSWs(
                               yearU,yeartime,
                               count,thresh,lastdate,startdate,
                               toprint,SSWdates)
        if yearcount == 1:
            singleyear +=1
            SSWyears.append(iyear + 1)
            #if toprint: print('year ' + str(iyear +1) + ' 1 SSW \n')
        elif yearcount ==2:
            doubleyear +=1
            #if toprint: print('year ' + str(iyear +1) + ' 2 SSWs \n')
            SSWyears.append(iyear + 1)
        elif yearcount ==3:
            tripleyear +=1
            SSWyears.append(iyear + 1)
        final.append(finalW)

    if singleyear + 2 * doubleyear +3 * tripleyear != count:
        print(count)
        print(singleyear + 2 * doubleyear +3 * tripleyear)
        sys.exit("problem with counting, maybe a year with more than 3 SSWs?!")

    mean,SE = meanSE(nyears,nyears - singleyear - doubleyear,singleyear,doubleyear)
    print ('mean: ' + str(mean) + ' ; s.e.: ' + str(SE) )

    return(SSWdates)

