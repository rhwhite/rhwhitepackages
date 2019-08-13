# -*- coding: utf-8 -*-
"""
Created on 13 Aug 2019

@author: rachel, rachel.white@cantab.net

Functions to calculation the barotropic refractive index
There are slightly different ways of calculating this in the literature,
potentially due to different ways of approximating to a Mercator projection.
Each method gives similar results

"""

import os, errno
import netCDF4
import numpy as np
import datetime as dt
import pandas as pd
import xarray as xr
#import Ngl
import math

from copy import copy

from scipy import stats
from scipy.interpolate import interp1d
import rhwhitepackages
from rhwhitepackages.readwrite import *
from rhwhitepackages.stats import regressmaps
from rhwhitepackages.griddata_functions import *
from rhwhitepackages.CESMconst import *

def calc_Ks_HK(Uin):

    ## Calculate BetaM
    ## Hoskins and Karoly (see also Vallis (page 551) and Petoukhov et al 2013)

    OMEGA = 7.2921E-5
    a = 6.3781E6
    try:
        lats_r = np.deg2rad(Uin.latitude)
    except AttributeError:
        lats_r = np.deg2rad(Uin.lat)

    coslat = np.cos(lats_r)

    betaM1 = 2.0 * OMEGA * coslat * coslat / a

    Um = Uin / coslat

    cos2Um = Um * coslat * coslat

    # first differentiation
    ddy_1 = ddy(cos2Um)
    # divide by cos2phi
    ddy_1_over_cos2p = ddy_1 * (1.0/(coslat * coslat))

    # second differentiation

    ddy_2 = ddy(ddy_1_over_cos2p)

    betaM = betaM1 - ddy_2
    # Now calculate Ks from BetaM
    # Hoskins and Ambrizzi and Hoskins and Karoly now agree
    Um = Uin / coslat
    Ks2 = a * a * betaM/Um

    Ks = np.sqrt(Ks2)

    return(Ks,Ks2) #,betaM)


def calc_Ks_P(Uin):

    ## Calculate Ks as in Petroukhov et al. 2013

    OMEGA = 7.2921E-5
    a = 6.3781E6
    try:
        lats_r = np.deg2rad(Uin.latitude)
    except AttributeError:
        lats_r = np.deg2rad(Uin.lat)

    coslat = np.cos(lats_r)
    sinlat = np.sin(lats_r)

    dUdPhi = ddphi(Uin)
    d2UdPhi2 = ddphi(dUdPhi)

    Ks1 = 2.0 * OMEGA * coslat * coslat * coslat / (a * Uin)

    Ks2 = (coslat * coslat * d2UdPhi2)/(a * a * Uin)

    Ks3 = (sinlat * coslat * dUdPhi)/(a * a * Uin)

    Ks4 = 1.0/(a * a)

    Ks2 = (a*a) * (Ks1 - Ks2 + Ks3 + Ks4)

    Ks = np.sqrt(Ks2)

    return(Ks,Ks2) #,betaM)

def calc_Ks_HA(Uin):
    print('Use of this function is not recommended as it has not been confirmed ' + 
            'why it differs from the others')
    ## Hoskins and Ambrizzi - THIS DOES NOT MATCH HOSKINS & KAROLY or OTHERS.
    ## I can't see how this would be related to allowing variations in longitude

    OMEGA = 7.2921E-5
    a = 6.3781E6 
    try:
        lats_r = Uin.latitude * 2 * math.pi / 360.0
    except AttributeError:
        lats_r = Uin.lat * 2 * math.pi / 360.0
    coslat = np.cos(lats_r)

    betaM1 = 2.0 * OMEGA

    Um = Uin / (a * coslat)

    cos2Um = Um * coslat * coslat

    # first differentiation

    ddphi_1 = ddphi(cos2Um)

    # divide by cosphi
    ddphi_1_over_cosp = ddphi_1 * (1.0/(coslat))

    # second differentiation

    ddphi_2 = ddphi(ddphi_1_over_cosp)

    # divide by cosphi

    betaM2 = ddphi_2 * (1.0/coslat)

    betaM_HA = (betaM1 - betaM2) * coslat * coslat / a
    # Now calculate Ks from BetaM
    # Hoskins and Ambrizzi and Hoskins and Karoly now agree
    Um = Uin / coslat
    Ks2 = (a * a * betaM_HA/Um)
    Ks = np.sqrt(Ks2)

    return(Ks,Ks2) #,betaM_HA)


