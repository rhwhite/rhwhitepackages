# Module to read and write using xray with useful error messages
# Written by rhwhite rachel.white@cantab.net
from netCDF4 import Dataset
import netCDF4
import numpy as np
import datetime as dt
import pandas
import xray
import sys
import stats

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

def shiftlons(invar,inlons):
        nlons = inlons.shape[0]
        nlonhalf = nlons/2

        newinvar = np.zeros(invar.shape,np.float)
        newinvar[:,0:nlonhalf] = invar[:,nlonhalf:nlons]
        newinvar[:,nlonhalf:nlons] = invar[:,0:nlonhalf]
        return newinvar

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

