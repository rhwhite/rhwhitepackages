'''
Custom diagnostics for CESM / CAM model output
This package is built on top of `xr` which provides the underlying
grid-aware data structures.
The method `open_dataset()` wraps the `xr.open_dataset()` method
and attempts to compute a bunch of useful diagnostics in addition to returning
a handle to the raw model output.

Written by R H White, extension and addition to packages by Brian Rose
Downloaded from https://github.com/brian-rose/pyCESM/blob/master/pyCESM/cam_diag.py

Jan 2018
'''

import numpy as np
import xarray as xr
from scipy import integrate
from xr.ufuncs import sin, cos, deg2rad
#from climlab import thermo

C2K = 273.15  # conversion factor from C to K
g2kg = (1.0/1000.0)  # converstion factor from g to kg
cm2m = (1.0/100)   # conversion factor from cm to m
mb_to_Pa = 100.  # conversion factor from mb to Pa

import rhwhitepackages
from rhwhitepackages.physconst import *

def getOHC(indata,heights,depth):
    temp = C2K + indata.TEMP.mean(dim='time') #K
    rho = ((g2kg/(cm2m * cm2m * cm2m)) *
                        indata.RHO.mean(dim='time'))  #g/cm3 -> kg/m3
    HCheight = cpocean * temp * rho * heights[:,None,None] # broadcast height onto lat lon grid 
                                                           # (J/kg/K) * K *
                                                           # kg/m3 * m = J/m2
    totalHC = HCheight.sum(dim='z_t')   # J/m2
    avHC = totalHC/depth  # J/m3

    return(totalHC)


def vertInt(indata)
