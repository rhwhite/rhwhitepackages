'''
Custom diagnostics for CESM / CAM model output
This package is built on top of `xarray` which provides the underlying
grid-aware data structures.
The method `open_dataset()` wraps the `xr.open_dataset()` method
and attempts to compute a bunch of useful diagnostics in addition to returning
a handle to the raw model output.

Downloaded from https://github.com/brian-rose/pyCESM/blob/master/pyCESM/cam_diag.py

on Oct 4 2017

'''

import numpy as np
import xarray as xr
from scipy import integrate
from xarray.ufuncs import sin, cos, deg2rad
#from climlab import thermo

import physconst
mb_to_Pa = 100.  # conversion factor from mb to Pa


def _pressure_formula(Ak, Bk, P0, PS):
    return Ak*P0 + Bk*PS

def _reorder_pressure(p, run):
    '''Reorder the input field to have the same dimension order as moisture'''
    dims = run.Q.dims
    pdims = ['ilev' if x=='lev' else x for x in dims]
    return p.transpose(*pdims)

def compute_pressure_mid(run):
    '''Convert hybrid sigma-pressure coordinates of mid-levels
    to pressure in Pa.'''
    p = _pressure_formula(Ak=run.hyam, Bk=run.hybm, P0=run.P0, PS=run.PS)
    return p

def compute_pressure_interface(run):
    '''Convert hybrid sigma-pressure coordinates of interface-levels
    to pressure in Pa.'''
    p = _pressure_formula(Ak=run.hyai, Bk=run.hybi, P0=run.P0, PS=run.PS)
    return p

def compute_pressure_intervals(run):
    '''Compute pressure intervals corresponding to each mid-level in Pa.'''
    dP = compute_pressure_interface(run).diff(dim='ilev')
    dP_reordered = _reorder_pressure(dP,run)
    return xr.DataArray(dP_reordered.values, coords=run['Q'].coords)


def open_dataset(filename_or_ob, verbose=True, **kwargs):
    '''Convenience method to open and return an xarray dataset handle,
    with precomputed CAM-specific diagnostic quantities.
    '''
    if verbose:
        print 'Opening dataset ', filename_or_ob
    dataset = xr.open_dataset(filename_or_ob, **kwargs)
    fulldataset = compute_all_diagnostics(dataset, verbose)
    return fulldataset


def compute_all_diagnostics(dataset, verbose=True):
    '''Compute all CAM-specific diagnostic quantities from xarray dataset.'''
    #  xr.Dataset has a property called .T which returns transpose()
    #  but this creates a conflict with data named 'T'
    #  Attempt to rename the T variable to TA
    if ('T' in dataset and 'TA' not in dataset):
        dataset = dataset.rename({'T':'TA'})
    if verbose:
        print 'Variable T renamed to TA.'
    fulldataset = compute_diagnostics(dataset)
    if verbose:
        print 'Gridpoint diagnostics have been computed.'
    try:
        fulldataset = compute_heat_transport_xarray(fulldataset)
        if verbose:
            print 'Heat transport diagnostics have been computed.'
    except:
        if verbose:
            print 'Heat transport computation failed.'
        else:
            pass
    try:
        # overturning mass streamfunction (in 10^9 kg/s or "mass Sverdrup")
        #V = fulldataset.V
        #if 'lon' in V.dims:
        #    #  Take zonal average first
        #    V = V.mean(dim='lon')
        #Psi = overturning(V)
        #Psi = overturning(fulldataset.V)
        Psi = overturning_improved(fulldataset)
        fulldataset = fulldataset.assign(Psi = Psi)
        if verbose:
            print 'Overturning streamfunction has been computed.'
    except:
        if verbose:
            print 'Overturning computation failed.'
        else:
            pass
    #  string data prevents us from doing arithmetic between two datasets
    #  So move any string data to attributes (where it belongs)
    #  Would be better to dynamically search for string data,
    #  but these are the two fields that usually appear in my CAM output
    string_fields = ['time_written', 'date_written']
    for item in string_fields:
        try:
            fulldataset.attrs[item] = fulldataset[item].values
            fulldataset = fulldataset.drop(item)
        except:
            pass

    return fulldataset

def inferred_heat_transport2(energy_in, lat=None, latax=None):
    ''' Similar to inferred_heat_transport, except it calculates
    the integral from S-N and N-S and averages. Requires xarray
    Compute heat transport as integral of local energy imbalance.
    Required input:
        energy_in: energy imbalance in W/m2, positive in to domain
    As either numpy array or xr.DataArray
    If using plain numpy, need to supply these arguments:
        lat: latitude in degrees
        latax: axis number corresponding to latitude in the data
            (axis over which to integrate)
    returns the heat transport in PW.
    Will attempt to return data in xarray.DataArray if possible.

    '''
    if lat is None:
        try: lat = energy_in.lat
        except:
            raise InputError('Must supply xarray with \'lat\' dimension.')
    lat_rad = np.deg2rad(lat)
    coslat = np.cos(lat_rad)
    field = coslat*energy_in
    if latax is None:
        try: latax = field.get_axis_num('lat')
        except:
            raise ValueError('Need to supply axis number for integral over latitude.')
    #  result as plain numpy array

    integral1 = integrate.cumtrapz(field, x=lat_rad, initial=0., axis=latax)
    # Now calculate in opposite direction and then flip back
    field2 = field.sel(lat=slice(None, None, -1))
    tempintegral2 = (integrate.cumtrapz(field2,np.deg2rad(field2.lat),initial=0.,axis=latax))

    # put back into xarray so we can swap the latitudes without making
    # assumptions about which dimension is latitude
    field2.values = tempintegral2
    integral2 = field2.sel(lat=slice(None,None,-1)).values

    integral = 0.5 * (integral1 + integral2)

    result = (1E-15 * 2 * np.math.pi * physconst.rearth**2 * integral)
    if isinstance(field, xr.DataArray):
        result_xarray = field.copy()
        result_xarray.values = result
        return result_xarray
    else:
        return result


def inferred_heat_transport(energy_in, lat=None, latax=None):
    '''Compute heat transport as integral of local energy imbalance.
    Required input:
        energy_in: energy imbalance in W/m2, positive in to domain
    As either numpy array or xr.DataArray
    If using plain numpy, need to supply these arguments:
        lat: latitude in degrees
        latax: axis number corresponding to latitude in the data
            (axis over which to integrate)
    returns the heat transport in PW.
    Will attempt to return data in xarray.DataArray if possible.
    '''
    if lat is None:
        try: lat = energy_in.lat
        except:
            raise InputError('Need to supply latitude array if input data is not self-describing.')
    lat_rad = np.deg2rad(lat)
    coslat = np.cos(lat_rad)
    field = coslat*energy_in
    if latax is None:
        try: latax = field.get_axis_num('lat')
        except:
            raise ValueError('Need to supply axis number for integral over latitude.')
    #  result as plain numpy array

    integral = integrate.cumtrapz(field, x=lat_rad, initial=0., axis=latax)
    result = (1E-15 * 2 * np.math.pi * physconst.rearth**2 * integral)
    if isinstance(field, xr.DataArray):
        result_xarray = field.copy()
        result_xarray.values = result
        return result_xarray
    else:
        return result


#%  overturning mass streamfunction
#Psi_cam4 = 2*pi*a/g.*repmat(cos(lat_cam4/180*pi),[1 26 8]).*permute(cumtrapz(lev_cam4*100,permute(V_cam4,[2 1 3])),[2 1 3])*1E-10;
#Psi_cam4_anom = Psi_cam4 - repmat(Psi_cam4(:,:,1),[1 1 8]);
#%  Compute Stokes mass streamfunction
#VprimeTprime_cam4 = VT_cam4 - V_cam4.*T_cam4;
#VprimeTprime_cam4_anom = VprimeTprime_cam4 - repmat(VprimeTprime_cam4(:,:,1),[1 1 8]);
#VprimeThetaprime_cam4 = VprimeTprime_cam4.*(1000./repmat(lev_cam4',[length(lat_cam4) 1 8])).^(0.286);
#VprimeThetaprime_cam4_anom = VprimeThetaprime_cam4 - repmat(VprimeThetaprime_cam4(:,:,1),[1 1 8]);
#oneoverdthetadp_cam4 = cat(2,1./(diff(THETA_cam4,1,2)./repmat((diff(100*lev_cam4))',[length(lat_cam4) 1 8])),zeros(length(lat_cam4),1,8));
#PsiStokes_cam4 = 2*pi*repmat(cos(lat_cam4*pi/180),[1 length(lev_cam4) 8])*a/g.*VprimeThetaprime_cam4.*oneoverdthetadp_cam4*1E-10;
#PsiStokes_cam4_anom = PsiStokes_cam4 - repmat(PsiStokes_cam4(:,:,1),[1 1 8]);
#% anomalies in TEM overturning:  Psi + Psi_stokes

def _prep_overturning(inputfield, lat=None, lev=None, levax=0):
    if lat is None:
        try: lat = inputfield.lat
        except:
            raise ValueError('Need to supply latitude array if input data is not self-describing.')
    if lev is None:
        try:
            lev = inputfield.lev
        except:
            raise ValueError('Need to supply pressure array if input data is not self-describing.')
    lat_rad = deg2rad(lat)
    coslat = cos(lat_rad)
    field = inputfield * coslat
    try: levax = field.get_axis_num('lev')
    except: pass
    return field, lat, lev, levax

def potential_temperature(TA):
    return TA * (1000./TA.lev)**physconst.cappa

def zonal_average(field):
    if 'lon' in field.dims:
        return field.mean(dim='lon')
    else:
        return field

def overturning_stokes(dataset):
    '''Compute Stokes mass streamfunction.
    Usage:
    Assuming `run` is a xarray.Dataset object containing CAM output:
    psi_stokes = overturning_stokes(run)
    '''
    from collections import OrderedDict

    VprimeTprime = zonal_average(dataset.VT - dataset.V * dataset.TA)
    VprimeThetaprime = zonal_average(potential_temperature(VprimeTprime))
    theta = zonal_average(potential_temperature(dataset.TA))
    field, lat, lev, levax = _prep_overturning(VprimeThetaprime)
    levax_theta = theta.get_axis_num('lev')
    dtheta = theta.diff(dim='lev')
    dp = xr.DataArray(np.diff(dataset.lev*mb_to_Pa),
                        dims='lev', coords={'lev': dtheta.lev})
    one_over_dthetadp = 1./ (dtheta / dp)
    #  Need to pad this with zeros
    #  this is just to get correct dimensions
    shape = list(theta.shape)
    shape[levax_theta] = 1
    newcoords = OrderedDict()
    for dim in theta.dims:
        if dim == 'lev':
            newcoords[dim] = theta.lev[0:1]
        else:
            newcoords[dim] = theta[dim]
    zero_pad = xr.DataArray(np.zeros(shape), coords=newcoords)
    one_over_dthetadp_ext = xr.concat([zero_pad, one_over_dthetadp], dim='lev')

    result = (2*np.pi*physconst.rearth/physconst.gravit *
                VprimeThetaprime * one_over_dthetadp_ext * 1E-9)
    return result

def overturning(V, lat=None, lev=None, levax=0):
    '''compute overturning mass streamfunction (SV)
    Required input:
        V: meridional velocity (m/s) (zonal average)
    As either numpy array or xarray.DataArray
    If using plain numpy, need to supply these arguments:
        lat: latitude in degrees
        lev: pressure levels in mb or hPa
        levax: axis number corresponding to pressure in the data
            (axis over which to integrate)
        levax argument is ignored if V is self-describing xarray.DataArray
    Returns the overturning streamfunction in SV or 10^9 kg/s.
    Will attempt to return data in xarray.DataArray if possible.
    '''
    field, lat, lev, levax = _prep_overturning(V, lat, lev, levax)
    #  Do we need to average zonally first?
    if 'lon' in field.dims:
        #  Take zonal average first
        field = field.mean(dim='lon')
    #  Otherwise we assume V is already zonally averaged
    #  result as plain numpy array
    result = (2*np.pi*physconst.rearth/physconst.gravit *
            integrate.cumtrapz(field, lev*mb_to_Pa, axis=levax,initial=0)*1E-9)
    if isinstance(field, xr.DataArray):
        result_xarray = field.copy()
        result_xarray.values = result
        return result_xarray
    else:
        return result

def overturning_improved(run):
    #  Better method: use the actual layer pressure intervals to weight the
    #  integral
    field = (run.V * run.dP * cos(deg2rad(run.lat)))
    if 'lon' in field.dims:
        field = field.mean(dim='lon')
    factor = 2*np.pi*physconst.rearth/physconst.gravit*1E-9
    #psi = np.cumsum( field, axis=field.get_axis_num('lev'))*factor
    # RHWhite modification - needs to run cumsum on values, not on xarray
    psi = np.cumsum( field.values, axis=field.get_axis_num('lev'))*factor
    return psi


#  SHOULD ALSO SET UNITS FOR EACH FIELD IN METADATA

def compute_diagnostics(run):
    '''Compute a bunch of additional diagnostics from regular CAM ouput.
    Input is an xarray dataset containing a CAM simulation climatology.
    These diagnostics are all computed point by point and so should work
    on any grid.'''

    sinlat = np.sin(np.deg2rad(run.lat))
    coslat = np.cos(np.deg2rad(run.lat))
    dP = compute_pressure_intervals(run)
    SST = run.TS - physconst.tmelt
    #TS_global = global_mean(run.TS,run.lat)
    #SST_global = run.TS_global - const.tempCtoK
    # TOA radiation fluxes
    SWdown_toa = run.SOLIN
    OLR = run.FLNT
    ASR = run.FSNT
    OLRclr = run.FLNTC
    ASRclr = run.FSNTC
    Rtoa = ASR - OLR  # net downwelling radiation
    Rtoaclr = ASRclr - OLRclr
    ASRcld = ASR - ASRclr
    OLRcld = OLR - OLRclr
    Rtoacld = Rtoa - Rtoaclr
    SWup_toa = SWdown_toa - ASR
    ALBtoa = SWup_toa / SWdown_toa

    #  surface energy budget terms, all defined as POSITIVE UP
    #    (from ocean to atmosphere)
    LHF = run.LHFLX
    SHF = run.SHFLX
    LWsfc = run.FLNS
    LWsfc_clr = run.FLNSC
    SWsfc = -run.FSNS
    SWsfc_clr = -run.FSNSC
    SnowFlux =  ((run.PRECSC + run.PRECSL) *
                      physconst.rhoh2o * physconst.latice)
    #  we'll let the SW down term be positive down
    SWdown_sfc = run.FSDS
    SWdown_sfc_clr = run.FSDSC
    SWdown_sfc_cld = SWdown_sfc - SWdown_sfc_clr
    SWup_sfc = SWsfc + SWdown_sfc
    ALBsfc = SWup_sfc / SWdown_sfc
    # all the net surface radiation terms are defined positive up
    LWsfc_cld = LWsfc - LWsfc_clr
    SWsfc_cld = SWsfc - SWsfc_clr
    # net upward radiation from surface
    SurfaceRadiation = LWsfc + SWsfc
    SurfaceRadiation_clr = LWsfc_clr + SWsfc_clr
    SurfaceRadiation_cld = SurfaceRadiation - SurfaceRadiation_clr
    # net upward surface heat flux
    SurfaceHeatFlux = SurfaceRadiation + LHF + SHF + SnowFlux
    # net heat flux into atmosphere
    Fatmin = Rtoa + SurfaceHeatFlux

    #  hydrological cycle, all terms in  kg/m2/s or mm/s
    Evap = run.QFLX
    Precip = (run.PRECC + run.PRECL) * physconst.rhoh2o
    EminusP = Evap - Precip

    newfields = {
        'sinlat': sinlat, 'coslat': coslat, 'dP': dP,
        'SST': SST,
        'SWdown_toa': SWdown_toa, 'SWup_toa': SWup_toa, 'ALBtoa': ALBtoa,
        'OLR': OLR, 'OLRclr': OLRclr, 'OLRcld': OLRcld,
        'ASR': ASR, 'ASRclr': ASRclr, 'ASRcld': ASRcld,
        'Rtoa': Rtoa, 'Rtoaclr': Rtoaclr, 'Rtoacld': Rtoacld,
        'LHF': LHF, 'SHF': SHF,
        'LWsfc': LWsfc, 'LWsfc_clr': LWsfc_clr, 'LWsfc_cld': LWsfc_cld,
        'SWsfc': SWsfc, 'SWsfc_clr': SWsfc_clr, 'SWsfc_cld': SWsfc_cld,
        'SnowFlux': SnowFlux,
        'SWdown_sfc': SWdown_sfc, 'SWdown_sfc_clr': SWdown_sfc_clr,
            'SWdown_sfc_cld': SWdown_sfc_cld,
        'SWup_sfc': SWup_sfc, 'ALBsfc': ALBsfc,
        'SurfaceRadiation': SurfaceRadiation,
        'SurfaceRadiation_clr': SurfaceRadiation_clr,
        'SurfaceRadiation_cld': SurfaceRadiation_cld,
        'SurfaceHeatFlux': SurfaceHeatFlux, 'Fatmin': Fatmin,
        'Evap': Evap, 'Precip': Precip, 'EminusP': EminusP,
    }
    #  use the xr.Dataset.assign() method to create a new dataset
    #  with all the new fields added, and return it.
    return run.assign(**newfields)


def compute_heat_transport_xarray(run):
    '''Compute heat transport diagnostics from regular CAM ouput.
    Input is an xarray dataset containing a CAM simulation climatology.
    These diagnostics involve integration so probably only work on
    regular lat-lon grids (for now).'''

    HT = compute_heat_transport2(run.Rtoa, run.SurfaceHeatFlux, run.EminusP)
    newfields = {
        'HT_total': HT['total'], 'HT_atm': HT['atm'], 'HT_ocean': HT['ocean'],
        'HT_latent': HT['latent'], 'HT_dse': HT['dse'],
    }
    #  use the xr.Dataset.assign() method to create a new dataset
    #  with all the new fields added, and return it.
    return run.assign(**newfields)

def compute_heat_transport2(TOAfluxDown, SurfaceFluxUp, EminusP,
                           lat=None, latax=None):
    # net heat flux into atmosphere
    Fatmin = TOAfluxDown + SurfaceFluxUp

    # heat transport terms
    HT_total = inferred_heat_transport2(TOAfluxDown, lat, latax)
    HT_atm = inferred_heat_transport2(Fatmin, lat, latax)
    HT_ocean = inferred_heat_transport2(-SurfaceFluxUp, lat, latax)
    # atm. latent heat transport from moisture imbal.
    HT_latent = inferred_heat_transport2(EminusP * physconst.latvap, lat, latax)
    # dry static energy transport as residual
    HT_dse = HT_atm - HT_latent
    HT = {
        'total': HT_total, 'atm': HT_atm, 'ocean': HT_ocean,
        'latent': HT_latent, 'dse': HT_dse,
        }
    return HT


def compute_heat_transport(TOAfluxDown, SurfaceFluxUp, EminusP,
                           lat=None, latax=None):
    # net heat flux into atmosphere
    Fatmin = TOAfluxDown + SurfaceFluxUp

    # heat transport terms
    HT_total = inferred_heat_transport(TOAfluxDown, lat, latax)
    HT_atm = inferred_heat_transport(Fatmin, lat, latax)
    HT_ocean = inferred_heat_transport(-SurfaceFluxUp, lat, latax)
    # atm. latent heat transport from moisture imbal.
    HT_latent = inferred_heat_transport(EminusP * physconst.latvap, lat, latax)
    # dry static energy transport as residual
    HT_dse = HT_atm - HT_latent
    HT = {
        'total': HT_total, 'atm': HT_atm, 'ocean': HT_ocean,
        'latent': HT_latent, 'dse': HT_dse,
        }
    return HT

    #  still need to work on the rest!

#    run.dz = -physconst.rair * run.T / physconst.gravit * run.dp / run.p *
#    1.E-3  #  in km
#    run.dTdp_moistadiabat = thermo.pseudoadiabat(run.T,run.p)
#    run.dTdp,ignored = np.gradient(run.T) / run.dp
#    run.dTdp_moistanom = run.dTdp - run.dTdp_moistadiabat
#    run.dTdz = run.dTdp * run.dp / run.dz  # in K / km
#    run.dTdz_moistanom = run.dTdp_moistanom * run.dp / run.dz
#    # convert OMEGA (in Pa/s) to w (in m/s)
#    run.w = -run.omega*const.Rd/const.g*run.T/(run.p*const.mb_to_Pa)
#    # overturning mass streamfunction (in 10^9 kg/s or "mass Sverdrup")
#    run.Psi = overturning(run.V,run.lat,run.lev)
#    #  correct for mass imbalance....
#    run.V_imbal = np.trapz(run.V/run.PS, run.lev*100, axis=0)
#    run.V_bal = run.V - run.V_imbal
#    run.Psi_bal = overturning(run.V_bal,run.lat,run.lev)
#    ind700 = np.nonzero(np.abs(run.lev-700)==np.min(np.abs(run.lev-700)))  #
#    closest vertical level to 700 mb
#    T700 = np.squeeze(run.T[ind700,:])
#    run.EIS = thermo.EIS(run.Ta,T700)
#    run.DSE = const.cp * run.T + const.g * run.Z
#    run.MSE = run.DSE + const.Lhvap * run.Q #  J / kg


#  NOTE would be better to use the xr .rename() method to change variable names
def convert_am2(run):
    '''Translate AM2 model output to the CAM naming conventions, so we can
    use the same diagnostic code.'''
    lev = run.pfull
    TS = run.t_surf
    T = run.temp
    Ta = run.t_ref # near-surface air temperature
    # TOA radiation
    SOLIN = run.swdn_toa
    FLNT = run.olr
    FLNTC = run.olr_clr
    FSNT = run.swdn_toa - run.swup_toa
    FSNTC = run.swdn_toa_clr - run.swup_toa_clr
    #  surface energy budget terms matching CAM sign and unit conventions
    LHFLX = run.evap * physconst.latvap
    SHFLX = run.shflx
    FLNS = -run.lwflx
    FLNSC = run.lwup_sfc_clr - run.lwdn_sfc_clr
    FSDS = run.swdn_sfc
    FSDSC = run.swdn_sfc_clr
    FSNS = run.swdn_sfc - run.swup_sfc
    FSNSC = run.swdn_sfc_clr - run.swup_sfc_clr
    #  snow flux
    PRECSC = run.snow_conv / physconst.rhoh2o
    PRECSL = run.snow_ls / physconst.rhoh2o
    # hydrological cycle
    QFLX = run.evap  # kg/m2/s or mm/s
    PRECC = run.prec_conv / physconst.rhoh2o  # m/s
    PRECL = run.prec_ls / physconst.rhoh2o
    # precipitable water in kg/m2
    TMQ = run.WVP
    # near-surface wind speed
    U10 = run.wind
    # Geopotential height
    Z = run.z_full
    #  moisture and clouds
    RELHUM = run.rh
    Q = run.sphum
    CLOUD = run.cld_amt_dyn
    #  surface pressure
    PS = run.ps
    #  velocity components
    U = run.ucomp
    V = run.vcomp

    newfields = {
    'lev': lev, 'TS': TS, 'T': T, 'Ta': Ta,
    'SOLIN': SOLIN, 'FLNT': FLNT, 'FLNTC': FLNTC, 'FSNT': FSNT, 'FSNTC': FSNTC,
    'LHFLX': LHFLX, 'SHFLX': SHFLX, 'FLNS': FLNS, 'FLNSC': FLNSC,
    'FSDS': FSDS, 'FSDSC': FSDSC, 'FSNS': FSNS, 'FSNSC': FSNSC,
    'PRECSC': PRECSC, 'PRECSL': PRECSL, 'PRECC': PRECC, 'PRECL': PRECL,
    'QFLX': QFLX, 'TMQ': TMQ, 'U10': U10, 'Z': Z,
    'RELHUM': RELHUM, 'Q': Q, 'CLOUD': CLOUD, 'PS': PS,
    'U': U, 'V': V,
    }
    #  use the xr.Dataset.assign() method to create a new dataset
    #  with all the new fields added, and return it.
    return run.assign(**newfields)
