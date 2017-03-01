# Module to read and write using xray with useful error messages
# Written by rhwhite rachel.white@cantab.net
from netCDF4 import Dataset
import netCDF4
import numpy as np
import datetime as dt
import pandas
import xray
import sys

def xrayOpen(filenamein,decodetimes=True):
    try:
        if decodetimes:
            filein= xray.open_dataset(filenamein)
        else:
            filein=xray.open_dataset(filenamein,decode_times=False)
    except MemoryError:
        exit("need more memory to read in file " + str(filenamein))
    except RuntimeError:
        print filenamein
        exit("couldn't find file")
    return filein

def xrayMfOpen(filenamein,decodetimes=True):
    try:
        if decodetimes:
            filein= xray.open_mfdataset(filenamein)
        else:
            filein=xray.open_mfdataset(filenamein,decode_times=False)
    except MemoryError:
        exit("need more memory to read in file " + str(filenamein))
    except RuntimeError:
        print filenamein
        exit("couldn't find file")
    return filein


def shiftlons(invar,nlons):
        nlonhalf = nlons/2

        newinvar = np.zeros(invar.shape,np.float)
        if len(invar.shape) == 1:
            newinvar[0:nlonhalf] = invar[nlonhalf:nlons]
            newinvar[nlonhalf:nlons] = invar[0:nlonhalf]
        elif len(invar.shape) == 2:
            newinvar[:,0:nlonhalf] = invar[:,nlonhalf:nlons]
            newinvar[:,nlonhalf:nlons] = invar[:,0:nlonhalf]
        else:
            exit('not set up for arrays with dims > 2 yet')

        return newinvar

def shiftlonlons(inlon,nlons):
        nlonhalf = nlons/2

        lonsnew = np.zeros(inlon.shape,np.float)
        lonsnew[0:nlonhalf] = inlon[nlonhalf:nlons]
        lonsnew[nlonhalf:nlons] = inlon[0:nlonhalf] + 360.0

        return lonsnew


def getunitsdesc(invarname):
        return{
                'gridboxspanSA': ['m2','spatial footprint of event during lifetime multiplied by gridbox surface area'],
                'totalprecipSA': ['m3','total precip of event during lifetime multipled by gridbox surface area'],
                'uniquegridboxspanSA': ['m2','unique gridbox spatial footprint of event during lifetime multipled by gridbox surface area'],
                'gridboxspan': ['# gridboxes','spatial footprint of event during lifetime'],
                'totalprecip': ['mm per event','total precip in event (no spatial weighting)'],
                'uniquegridboxspan': ['# gridboxes','unique gridbox spatial footprint of event during lifetime'],
                'timespan': ['hours','timespan of event'],
                'tstart': ['index','time index of start of event'],
                'tmean': ['index','time index of middle of event'],
                'xcenterstart': ['index','x index of center of event at event start'],
                'xcenterend': ['index','x index of center of event at event end'],
                'ycenterstart': ['index','y index of center of event at event start'],
                'ycenterend': ['index','y index of center of event at event send'],
                'xcentermean': ['index','x index of center of event at event mean time'],
                'ycentermean': ['index','y index of center of event at event mean time'],
                'xmin': ['index','minimum x index over lifetime'],
                'xmax': ['index','maximum x index over lifetime'],
                'ymin': ['index','minimum y index over lifetime'],
                'ymax': ['index','maximum y index over lifetime'],
                'xmaxspeed_4ts': ['m/s','maximum event zonal speed over 4 timesteps'],
                'xmaxspeed_1ts': ['m/s','maximum event zonal speed over 1 timestep'] ,
                }[invarname]


def getdenfilename(mappingi, datai, versioni, fstartyri, fendyri, iboundi, splittypei, uniti, speedtspani, minGBi, tbound1i, tbound2i, monanni,sumlatsi,sumlonsi):
    if splittypei == "maxspeed":
        fileTypeadd = "MaxSpeeds_" + str(speedtspani) + "ts_"
    elif splittypei == "speed":
        fileTypeadd = "Speeds_"
    elif splittypei == "day":
        fileTypeadd = "Sizes_"
    else:
        exit("unexpected splittype")

    if minGBi > 0:
        fileadd = '_min' + str(minGBi) + 'GB'
    else:
        fileadd = ''
    
    if tbound1i[iboundi] < 0:
        tboundtitle = str(tbound1i[iboundi]) + '-' + str(tbound2i[iboundi]) + uniti
    else:
        tboundtitle = str(tbound1i[iboundi]) + '-' + str(tbound2i[iboundi]) + uniti
    
    if sumlatsi > 0:
        addsumlats = '_regrid_' + str(sumlonsi) + 'lons_' + str(sumlatsi) + 'lats'
    else:
        addsumlats = ''

    # get directory
    if splittype == 'day':
        diradd = '/Sizes/'
    elif splittype == 'MaxSpeeds':
        diradd = '/MaxSpeeds/'

    dirIn = ('/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/' + datai + '_output/'
               + versioni + str(fstartyri) + '/proc/' + diradd)

    # concat directory and filename
    fileName = (dirIn + 'DenDirSpd_Map_' + monanni + '_' + fileTypeadd + tboundtitle + '_' + mappingi + '_' 
                    + datai + "_" + str(fstartyri) + '-' + str(fendyri) + '_' + versioni 
                    + fileadd + addsumlats + '.nc')
    return(fileName)

def getPrecipfilename(mappingi, datai, versioni, fstartyri, fendyri, iboundi, splittypei, uniti, speedtspani, minGBi, tbound1i, tbound2i):
    if splittypei == "maxspeed":
        fileTypeadd = "MaxSpeeds_" + str(speedtspani) + "ts_"
    elif splittypei == "speed":
        fileTypeadd = "Speeds_"
    elif splittypei == "day":
        fileTypeadd = "Sizes_"
    else:
        exit("unexpected splittype")

    if tbound1i[iboundi] < 0:
        tboundtitle = str(tbound1i[iboundi]) + '-' + str(tbound2i[iboundi])
    else:
        tboundtitle = str(tbound1i[iboundi]) + '-' + str(tbound2i[iboundi])

    if minGBi > 0:
        fileaddGB = '_min' + str(minGBi) + 'GB'
    else:
        fileaddGB = ''

    if mappingi != 'center':
        fileaddmap = mappingi
    else:
        fileaddmap = ''

    return('Precip_' + fileTypeadd + tboundtitle + uniti + "_" + str(fstartyri) + '-' + str(fendyri) + '_' + versioni + fileaddGB + '.nc')

def getrawPrecipAnn(Data,Version,minlat,maxlat,anstartyr,anendyr):
    base = '/home/disk/eos4/rachel/EventTracking/'
    if Data == "TRMM":
        Fstartyr = 1998
        Fendyr = 2014
        PrecipClimDir = "/home/disk/eos4/rachel/Obs/TRMM/"
        PrecipClimFile = "TRMM_1998-2014_clim_ann_1998-2014.nc"

    elif Data == "TRMMERAIgd":
        Fstartyr = 1998
        Fendyr = 2014
        PrecipClimDir = "/home/disk/eos4/rachel/Obs/TRMM/"
        PrecipClimFile = "regrid2ERAI_TRMM_3B42_1998-2014_annclim.nc"

    elif Data == "ERAI":
        Fstartyr = 1980
        Fendyr = 2014
        PrecipClimDir = '/home/disk/eos4/rachel/Obs/ERAI/Precip_3hrly/'
        PrecipClimFile = 'ERAI_Totalprecip_1980-2015_annmean.nc'

    elif Data == "ERA20C":
        Fstartyr = 1980
        Fendyr = 2011
        PrecipClimDir = '/home/disk/eos4/rachel/Obs/ERA_20C/'
        PrecipClimFile = ('ERA_20C_Ann_Totalprecip_' + str(Fstartyr) + '-' +
                            str(Fendyr) + '.nc')

    elif Data == "CESM":
        Fstartyr = 1990
        Fendyr = 2014

        PrecipClimDir = ('/home/disk/eos4/rachel/EventTracking/Inputs/' +
                            'CESM/f.e13.FAMPIC5.ne120_ne120.1979_2012.001/')
        PrecipClimFile = 'ncra_f.e13.FAMIPC5.ne120_ne120_TotalPrecip_1979-2012.nc'

    else:
        exit('not set up for Data type ' + Data + ' version ' + Version)

    # open up precip file
    FileInPrecip = xrayOpen(PrecipClimDir + PrecipClimFile)

    if Data in ['TRMM']:
        latin = FileInPrecip['latitude']
        invar = 'pcp'

    elif Data in ['TRMMERAIgd']:
        latin = FileInPrecip['lat']
        invar = 'pcp'

    elif Data in ["ERAI","ERA20C"]:
        latin = FileInPrecip['lat']
        invar = 'tpnew'

    elif Data in ["CESM"]:
        latin = FileInPrecip['lat']
        invar = 'PRECT'

    if latin[0] > latin[1]:
        PrecipIn = (FileInPrecip[invar][:,::-1,:].sel(lat=slice(minlat,maxlat))
                                            .sel(time=slice(str(anstartyr),str(anendyr))))
    else:
        PrecipIn = (FileInPrecip[invar].sel(latitude=slice(minlat,maxlat))
                                        .sel(time=slice(str(anstartyr),str(anendyr))))

    try:
        print np.amax(PrecipIn)
        if PrecipIn.units == "mm/hr":
            print 'converting'
            PrecipIn = PrecipIn * 24.0 #convert to mm/day
        elif PrecipIn.units == "mm/day":
            pass
        else:
            error("unexpected unit in Precip file")
    except AttributeError:
        if np.amax(PrecipIn) < 5.0:
            print("guessing we need to convert precip units! You should" + 
                     "check this!")
            PrecipIn = PrecipIn * 24.0 #convert to mm/day
        else:
            print("guessing that we don't need to convert precip units" +
                    " - you should check this!")

    return(PrecipIn)


def getdirectory(splittype):
    if splittype == "maxspeed":
        diradd = "MaxSpeeds"
    elif splittype == "speed":
        diradd = "Speeds"
    elif splittype == "day":
        diradd = "Sizes"
    else:
        exit("unexpected splittype")

    return diradd  
