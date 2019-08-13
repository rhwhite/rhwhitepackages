# Module to calculate statistics on maps
# Written by rhwhite rachel.white@cantab.net
from netCDF4 import Dataset
import netCDF4
import numpy as np
import datetime as dt
import pandas
import xarray as xr
from scipy import stats

def regressmaps(m,c,r,p,stderr,A,linreg,nlats,nlons):
# Calculate linear trend
        for ilat in range(0,nlats):
                for ilon in range(0,nlons):
                        m[ilat,ilon],c[ilat,ilon],r[ilat,ilon],p[ilat,ilon],stderr[ilat,ilon] = stats.linregress(A,linreg[:,ilat,ilon])

