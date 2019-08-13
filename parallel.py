import numpy as np
import numpy.ma as ma
import xarray as xr
import netCDF4
from netCDF4 import Dataset

def calcmissed(procnum,qout,q2out,precipin,threshlow,startin,endin):
#        if lowmem:
#                print "using low memory version"
        temp = np.float32(precipin.isel(time=slice(startin,endin)))
#        else:
#                temp = np.float64(precipnew.isel(time=slice(startin,endin)))
#       if lowmem:
#               print "using low memory version"
#               temp = np.float32(precipnew.sel(time=str(yearin) + '-' + '{:02d}'.format(monthin)))
#       else:
#                temp = np.float64(precipnew.sel(time=str(yearin) + '-' + '{:02d}'.format(monthin)))

        #print "number of bytes:" , temp.nbytes
        temp2 = temp[np.where( temp >= threshlow)]
        sum1 = np.nansum(temp) - np.nansum(temp2) 
        print procnum, "sum1", sum1
        sum2 = np.nansum(temp)
        print procnum, "sum2", sum2
        del(temp,temp2)
	qout.put(sum1)
	q2out.put(sum2)

