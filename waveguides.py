# Module identify waveguides in refractive index data
# Written by rhwhite rachel.white@cantab.net
import numpy as np
import xarray as xr
import math

# count_waveguides and write out magnitude and latitude
# Only looking for waveguides at 2 specific latitudes
def identify_waveguides_mag(datain,wavenumber):

    wnbrs = wavenumber * wavenumber

    temp = datain.copy(deep=True)

    # Subtropical waveguide
    # Look for a value of wnbrs between 30 and 36
    # Look for at least 6deg (3 consecutive gridboxes) positive values between 34 and 46
    # Look for a value of wnbrs between 44 and 60
    wguideST = np.zeros(len(temp.longitude))
    wguideML = np.zeros(len(temp.longitude))

    wguideSTlat = np.zeros(len(temp.longitude))
    wguideMLlat = np.zeros(len(temp.longitude))
    
    for ilon in range(0,len(temp.longitude)):
        tempsel = temp.isel(longitude=ilon).copy(deep=True)
        # Cut off all values below wnbrs:
        tempsel[:] = np.where(tempsel<wnbrs,wnbrs,tempsel)    
        if np.any(tempsel.sel(latitude=slice(36,30)) == wnbrs):
            for ilat in range(34,46):
                result = np.all(tempsel.sel(latitude=slice(ilat+6,ilat)) > wnbrs)
                if result:
                    if np.any(tempsel.sel(latitude = slice(60,44)) == wnbrs):
                        wguideST[ilon] = np.amax(tempsel.sel(latitude=slice(46,34))) - wnbrs
                        templats = tempsel.sel(latitude=slice(52,34))
                        wguideSTlat[ilon] = templats.latitude.isel(latitude = np.argmax(templats.values))
                        break
                    else:
                        wguideST[ilon] = 0
            else:
                wguideST[ilon] = 0
        else:
            wguideST[ilon] = 0

        if np.any(tempsel.sel(latitude=slice(50,40)) == wnbrs):
            for ilat in range(48,60):
                result = np.all(tempsel.sel(latitude=slice(ilat+6,ilat)) > wnbrs)
                if result:
                    if np.any(tempsel.sel(latitude = slice(70,58)) == wnbrs):
                        wguideML[ilon] = np.amax(tempsel.sel(latitude=slice(60,48))) - wnbrs
                        templats = tempsel.sel(latitude=slice(66,48))
                        wguideMLlat[ilon] = templats.latitude.isel(latitude = np.argmax(templats.values))
                        break
                    else:
                        wguideML[ilon] = 0
            else:
                wguideML[ilon] = 0
        else:
            wguideML[ilon] = 0

    return(wguideST,wguideSTlat,wguideML,wguideMLlat)

# Identify waveguides in data, only looking for waveguides at certain latitudes, can't stray from this!
def waveguide_analysis_fixedlat(Ks2in,wnb):
    maxKs = 200
    Ks2in[...] = np.where(Ks2in < maxKs,Ks2in,maxKs)
    Ks2in[...] = np.where(Ks2in > -maxKs,Ks2in,-maxKs)

    waveguide_freq_ST = np.zeros([len(Ks2in.time),len(Ks2in.longitude)])
    waveguide_freq_ML = np.zeros([len(Ks2in.time),len(Ks2in.longitude)])
    waveguide_lat_ST = np.zeros([len(Ks2in.time),len(Ks2in.longitude)])
    waveguide_lat_ML = np.zeros([len(Ks2in.time),len(Ks2in.longitude)])
    
    for itime in range(0,len(Ks2in.time)):
    #    if itime % 100 == 0: print itime
                   
        waveguide_freq_ST[itime],waveguide_lat_ST[itime],waveguide_freq_ML[itime],waveguide_lat_ML[itime] = identify_waveguides_mag(Ks2in.isel(time=itime),wnb)
    
    DA_ST = xr.DataArray(waveguide_freq_ST, coords=[Ks2in.time, Ks2in.longitude], dims=['time', 'longitude'])
    DA_ML = xr.DataArray(waveguide_freq_ML, coords=[Ks2in.time, Ks2in.longitude], dims=['time', 'longitude'])
    DA_STlat = xr.DataArray(waveguide_lat_ST, coords=[Ks2in.time, Ks2in.longitude], dims=['time', 'longitude'])
    DA_MLlat = xr.DataArray(waveguide_lat_ML, coords=[Ks2in.time, Ks2in.longitude], dims=['time', 'longitude'])

    DS_ST = DA_ST.to_dataset(name='waveguide_freq_ST')
    DS_ML = DA_ML.to_dataset(name='waveguide_freq_ML')
    DS_STlat = DA_STlat.to_dataset(name='waveguide_lat_ST')
    DS_MLlat = DA_MLlat.to_dataset(name='waveguide_lat_ML')
    
    towrite = xr.merge([DS_ST, DS_ML,DS_STlat,DS_MLlat])
    del towrite.time.encoding["contiguous"]
    
    return(towrite) 
