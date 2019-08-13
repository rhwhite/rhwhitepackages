# Module to read and write using xr with useful error messages
# Written by rhwhite rachel.white@cantab.net
from netCDF4 import Dataset
import netCDF4
import numpy as np
import datetime as dt
import pandas
import xarray as xr
import Ngl
from rhwhitepackages.readwrite import *

def setplotrange(resMP,a,b,c):
        resMP.cnMinLevelValF = a
        resMP.cnMaxLevelValF = b
        resMP.cnLevelSpacingF = c

	return(resMP)

def initcontourplot(resMP,plotstartlat,plotstartlon,plotendlat,plotendlon,lats,lons):
        resMP.nglDraw  = False
        resMP.nglFrame = False

        resMP.cnFillOn = True
        resMP.cnMonoFillPattern = True
        resMP.cnMonoFillColor = False

        resMP.cnMissingValFillColor = "white"

        resMP.cnLineLabelsOn = False
        resMP.cnInfoLabelOn         = False    # Turn off informational

        resMP.cnLinesOn = False
        resMP.cnSmoothingOn = False
        resMP.cnLevelSelectionMode = "ManualLevels"

        resMP.pmLabelBarDisplayMode = "Conditional"
        resMP.lbPerimOn = False
        resMP.trYMaxF = plotendlat
        resMP.trYMinF = plotstartlat
        resMP.trXMaxF = plotendlon
        resMP.trXMinF = plotstartlon

        resMP.mpProjection = "CylindricalEquidistant" # Change the map projection.
        resMP.mpCenterLonF = 180.           # Rotate the projection.
        resMP.mpFillOn     = True           # Turn on map fill.
        resMP.mpLimitMode = "LatLon"    # Limit the map view.
        resMP.mpMinLonF = plotstartlon
        resMP.mpMaxLonF = plotendlon
        resMP.mpMinLatF = plotstartlat
        resMP.mpMaxLatF = plotendlat
        resMP.mpOutlineBoundarySets = "National"

        # don't have label bar take up so much space
        resMP.pmLabelBarWidthF = 0.05


        resMP.sfYArray = lats
        resMP.sfXArray = lons

	return(resMP)

def initCScontourplot(resCS,plotstartlat,plotstartdepth,plotendlat,plotenddepth,lats,depths):
        resCS.nglDraw  = False
        resCS.nglFrame = False

        resCS.cnFillOn = True
        resCS.cnMonoFillPattern = True
        resCS.cnMonoFillColor = False
        resCS.cnLineLabelsOn = False
        resCS.cnLinesOn = False
        resCS.cnSmoothingOn = False
        resCS.cnLevelSelectionMode = "ManualLevels"

        resCS.pmLabelBarDisplayMode = "Conditional"
        resCS.lbPerimOn = False
        resCS.trYMaxF = plotenddepth
        resCS.trYMinF = plotstartdepth
        resCS.trXMaxF = plotendlat
        resCS.trXMinF = plotstartlat

	resCS.trYReverse = True

        resCS.sfYArray = depths
        resCS.sfXArray = lats

        return(resCS)

def initLinePlot(resLP,plotstartlat,plotstartlon,plotendlat,plotendlon,lats,lons):
        resLP.nglFrame = False
        resLP.nglDraw = False
        resLP.xyMarkLineMode = "Markers"
        resLP.xyMonoMarkLineMode = True
        resLP.xyMarkerColor = "blue"
        resLP.xyMarkers = 16
        resLP.xyMarkerSizeF = 0.01
        resLP.xyYStyle = "Linear"
        resLP.xyXStyle = "Linear"

        resLP.tiMainFontHeightF = 0.035
        resLP.tiYAxisFontHeightF = 0.03
        resLP.tiXAxisOn = False
        resLP.tmXBLabelFontHeightF = 0.03
        resLP.tmYLLabelFontHeightF = 0.03
        resLP.tmYLMode = "Automatic"
        resLP.tmYLFormat = "@6^g"

        resLP.vpWidthF = 0.9
        resLP.vpHeightF = 0.9

def initvectoroverlay(resOv,plotstartlat,plotstartlon,plotendlat,plotendlon,lats,lons):
	resOv.nglDraw  = False
        resOv.nglFrame = False

        resOv.vcRefLengthF = 0.03

        resOv.vcLevelSelectionMode = "ManualLevels"

        resOv.vcMonoLineArrowColor = True
        resOv.vfYArray = lats
        resOv.vfXArray = lons

        resOv.pmLabelBarDisplayMode = "Never"

        resOv.vcMinMagnitudeF = 2
        resOv.vcMinDistanceF = 0.02       # thin vectors

        resOv.lbLabelFontHeightF = 0.01
        resOv.lbTitleFontHeightF = 0.01

        resOv.tiMainFontHeightF = 0.015

	return(resOv)

def initcontouroverlay(resOv,plotstartlat,plotstartlon,plotendlat,plotendlon,lats,lons):
	resOv.cnInfoLabelOn         = False    # Turn off informational
	resOv.nglDraw  = False
	resOv.nglFrame = False

	resOv.cnLineLabelsOn       = False
	resOv.sfXCStartV = float(lons[0])
	resOv.sfXCEndV = float(lons[len(lons)-1])
	resOv.sfYCStartV = float(lats[0])
	resOv.sfYCEndV = float(lats[len(lats)-1])

	#resOv.mpProjection = "CylindricalEquidistant" # Change the map projection.
	#resOv.mpCenterLonF = 180.           # Rotate the projection.

	resOv.lbLabelBarOn   = False
	#resOv.mpLimitMode = "LatLon"    # Limit the map view.
	#resOv.mpMinLonF = plotslon
	#resOv.mpMaxLonF = plotelon
	#resOv.mpMinLatF = plotslat
	#resOv.mpMaxLatF = plotelat
	#resOv.mpOutlineBoundarySets = "AllBoundaries"

	resOv.cnLevelSelectionMode = "ExplicitLevels" # use explicit levels
	resOv.cnLineDashPattern = 1        # sets negative contours to dash pattern 1

	resOv.cnFillOn = False
	resOv.cnConstFLabelOn = False
	resOv.cnInfoLabelOn       = False        # no info label

	return(resOv)

def contourwithoverlay(wks,resMP,varin, plotoverlay):
	plottemp = Ngl.contour_map(wks,varin,resMP)
        Ngl.overlay(plottemp,plotoverlay)
	return(plottemp)

def gettitle(expin):
	return {
		'CESMtopof19': 'CTL',
		'CESMnoT2f19':'noMT',
		'CESMnoTf19':'noT',
		'CESMnoT4f19':'noM',
		'CESMnoMTf19':'noMT',
		'0':'',
		}[expin]
def setmissing(inarray,missingv):
	inarray[np.where(np.logical_not(np.isfinite(inarray)))] = missingv
	return(inarray)


def plotline(nlines,plotin,yearsin,nsep,title,Datatitle,figtitlein,colormap,xString,yString,unit):
        wkres = Ngl.Resources()
        wkres.wkColorMap = colormap
        wks_type = "eps"
        wks = Ngl.open_wks(wks_type,figtitlein,wkres)

	res = Ngl.Resources()
        plot = []

        res.tiXAxisString = xString
        res.tiYAxisString = yString

	for iline in range(0,nlines):
                yearsplot = yearsin[iline]
                print yearsplot
#               nyearsin = len(yearsplot)
#               yearnums = range(0,nyearsin)
                A = np.array(yearsplot)
                regress = (np.zeros((5,nsep),np.float))

                for ibin in range(0,nsep):
                        if ibin == 0:
                                if (plotdensity):
                                        res.tiYAxisString = "frequency"
                                else:
                                        res.tiYAxisString = "number of events"
                        else:
                                res.tiYAxisString = ""

                        linreg = plotin[iline][ibin][:]
                        print A.shape
                        print linreg.shape
                        regress[:,ibin] = stats.linregress(A,linreg)

                        if ibin == nsep -1:
                                res.tiMainString = '{:^80}'.format('          >' + '{:2.1g}'.format(tbound1[ibin]) + unit + '; p=' + '{:5.3f}'.format(regress[3,ibin]) + "           ")
                        else:
                                res.tiMainString = '{:^80}'.format('{:2.1g}'.format(tbound1[ibin]) + '-' + '{:2.1g}'.format(tbound2[ibin]) + unit + '; p=' + '{:5.3f}'.format(regress[3,ibin]))

                        plot.append(Ngl.xy(wks,yearsplot,plotin[iline][ibin,:],res))

        panelres = Ngl.Resources()
        panelres.nglPanelLabelBar = True
        panelres.nglPanelYWhiteSpacePercent = 8.
        panelres.nglPanelXWhiteSpacePercent = 0.0

        panelres.nglPanelLabelBar                 = False     # Turn on panel labelbar
        panelres.nglPanelTop                      = 1.0
        panelres.nglPanelBottom                      = 0.00
        panelres.nglPanelLeft                   = 0.0
        panelres.nglPanelRight                  = 1.0
        panelres.nglPaperOrientation = "Portrait"
        panelres.nglScale = False
        panelres.nglMaximize = True

        #txres = Ngl.Resources()
        #txres.txFontHeightF = 0.012
        #Ngl.text_ndc(wks,'Annual timeseries of global number of events of various scales from' + Datatitle,0.5,0.94,txres)

        Ngl.panel(wks,plot,[nlines,float(nbounds)],panelres)

        print 'panelled'

def plotmap(plotvars1,plotvars2,
            plotmin1,plotmax1,plotmin2,plotmax2,
            vartitle1,vartitle2,
            title,
            figtitle,
            minlon,maxlon,minlat,maxlat,
            FillValue,panellabels = [],
            labelbarlabels = [],
            labelbarlabels2 = []):

    nplots = plotvars1.shape[0]
    wkres = Ngl.Resources()
    wkres.wkColorMap = "WhiteBlue"
    wks_type = "eps"
    wks = Ngl.open_wks(wks_type,figtitle,wkres)

    # if lons start negative, shift everything over so there isn't a line down
    # the middle of the Pacific

    lons1 = plotvars1.lon.values
    lons2 = plotvars2.lon.values

    lats1 = plotvars1.lat.values
    lats2 = plotvars2.lat.values

    if lons1[0] < 0:
        #nlonhalf1 = len(lons1)/2
        #lonsnew1 = np.zeros(lons1.shape,np.float)
        #lonsnew1[0:nlonhalf1] = lons1[nlonhalf1:nlons1]
        #lonsnew1[nlonhalf1:nlons1] = lons1[0:nlonhalf1] + 360.0

        lonsnew1 = shiftlonlons(lons1,len(lons1))
        lonsnew2 = shiftlonlons(lons2,len(lons2)) 

        for iplot in range(0,nplots):
            plotvars1[iplot] = shiftlons(plotvars1[iplot],len(lons1))
            plotvars2[iplot] = shiftlons(plotvars2[iplot],len(lons2))
    else:
        lonsnew1 = lons1
        lonsnew2 = lons2

    # initialize plotting resources
    res1 = Ngl.Resources()
    res1 = initcontourplot(res1,minlat,minlon,maxlat,maxlon,lats1,lonsnew1)
    res1.sfMissingValueV = FillValue
    res1.lbOrientation   = "Vertical"
    # including some font heights
    res1.lbLabelFontHeightF = 0.01
    res1.lbTitleFontHeightF = 0.01
    res1.tiMainFontHeightF = 0.015
    res1.lbTitlePosition = 'Bottom'
    res1.lbBottomMarginF = 0.0

    # initialize plotting resources
    res2 = Ngl.Resources()
    res2 = initcontourplot(res2,minlat,minlon,maxlat,maxlon,lats2,lonsnew2)
    res2.sfMissingValueV = FillValue
    res2.lbOrientation   = "Vertical"
    # including some font heights
    res2.lbLabelFontHeightF = 0.01
    res2.lbTitleFontHeightF = 0.01
    res2.tiMainFontHeightF = 0.008
    res2.lbTitlePosition = 'Bottom'
    res2.lbBottomMarginF = 0.0


    # turn off grid lines
    res1.mpGridAndLimbOn = False
    res2.mpGridAndLimbOn = False
    # initialize plotting array
    toplot = []
    # fill plotting array
    for iplot in range(0,nplots):
        tempplot = plotvars1[iplot].values
        tempplot[np.where(np.isnan(tempplot))] = FillValue
        # update plot resources with correct lat/lon
        res1.cnMinLevelValF       = plotmin1[iplot]
        res1.cnMaxLevelValF       = plotmax1[iplot]
        res1.cnLevelSpacingF      = ((plotmax1[iplot]-plotmin1[iplot])/10.0)
        res1.tiMainString = (vartitle1[iplot])

        if panellabels != []:
            res1.lbTitleString = labelbarlabels[iplot]
            res1.tiYAxisString  = panellabels[iplot]  # Y axes label.


        toplot.append(Ngl.contour_map(wks,tempplot,res1))

        tempplot = plotvars2[iplot].values
        tempplot[np.where(np.isnan(tempplot))] = FillValue
        res2.cnMinLevelValF       = plotmin2[iplot]          # contour levels.
        res2.cnMaxLevelValF       = plotmax2[iplot]
        res2.cnLevelSpacingF      = ((plotmax2[iplot]-plotmin2[iplot])/10.0)
        res2.tiMainString = vartitle2[iplot]
        if panellabels != []:
            res2.lbTitleString = labelbarlabels2[iplot]
            res2.tiYAxisString  = " "  # so plots are the same
                                                      # size
        toplot.append(Ngl.contour_map(wks,tempplot,res2))

    textres = Ngl.Resources()
    textres.txFontHeightF = 0.015
    Ngl.text_ndc(wks,title,0.5,0.87,textres)

    panelres = Ngl.Resources()
    panelres.nglPanelLabelBar = True
    panelres.nglPanelYWhiteSpacePercent = 0.
    panelres.nglPanelXWhiteSpacePercent = 0.

    panelres.nglPanelLabelBar   = False     # Turn on panel labelbar
    if nplots > 5:
        panelres.nglPanelTop                      = 0.8
        panelres.nglPanelBottom                      = 0.15
    else:
        panelres.nglPanelTop                      = 0.95
        panelres.nglPanelBottom                      = 0.01

    panelres.nglPanelLeft = 0.01
    panelres.nglPanelRight = 0.99

    panelres.nglPanelFigureStrings = (
            ['a.','b.','c.','d.','e.','f.','g.','h.','i.','j.','k.','l.','m.','n.','o.','p.'])
    panelres.nglPanelFigureStringsJust = "TopLeft"
    panelres.nglPanelFigureStringsFontHeightF = 0.008
    panelres.nglPanelFigureStringsParallelPosF = -0.55
    panelres.nglPanelFigureStringsOrthogonalPosF = -0.7
    panelres.nglPanelFigureStringsPerimOn = False   # turn off boxes
    #panelres.amJust = "TopLeft"

    panelres.nglPaperOrientation = "Auto"

    plot = Ngl.panel(wks,toplot,[nplots,2],panelres)

def getFITcolorbars(Datain,minGBin,splittypein,varin):
    cbmin = -1
    cbmax = -1
    if splittypein == 'day':
        if Datain in ["TRMM"]:
            if varin in ['TDensity']:
                cbmin,cbmax = [0.0,98.0,0.0,0.0,0.0,0.0],[300,100.0,1.0,0.2,0.05,1.0]
            elif varin in ['TPrecip']:
                cbmin,cbmax = (
                        [40.0,50.0,0.0,0.0,0.0,0.0],[100,100.0,40.0,40.0,40.0,10.0])
            elif varin in ['LocalDensity']:
                cbmin,cbmax =(
                        [0.0,50.0,0.0,0.0,0.0,0.0],[200,100,30,20.0,10.0,10.0])

        elif Datain in ["TRMMERAIgd"]:
            if minGBin in [4,9]:
                if varin in ['TDensity']:
                    cbmin,cbmax = [0.0,70.0,0.0,0.0,0.0,0.0],[6,100.0,20.0,7.0,4.0,3.0]
                elif varin in ['TPrecip']: 
                    cbmin,cbmax = [0.0,0.0,0.0,0.0,0.0,0.0],[100,80.0,50.0,50.0,50.0,50.0]

            else:
                if varin in ['TDensity']:
                    cbmin,cbmax = [0.0,90.0,0.0,0.0,0.0,0.0],[10,100.0,15.0,5.0,2.0,3.0]
                elif varin in ['TPrecip']: 
                    cbmin,cbmax = [0.0,0.0,0.0,0.0,0.0,0.0],[100,80.0,40.0,50.0,50.0,50.0]

        elif Datain in ["ERAI"]:
            if minGBin in [4,9]:
                if varin in ['TDensity']:
                    cbmin,cbmax = [0.0,70.0,0.0,0.0,0.0,0.0],[6,100.0,20.0,7.0,4.0,3.0]
                elif varin in ['TPrecip']:
                    cbmin,cbmax = [0.0,0.0,0.0,0.0,0.0,0.0],[100,80.0,50.0,50.0,50.0,50.0]

            else:
                if varin in ['TDensity']:
                    cbmin,cbmax = [0.0,90.0,0.0,0.0,0.0,0.0],[10,100.0,15.0,5.0,2.0,3.0]
                elif varin in ['TPrecip']:
                    cbmin,cbmax = [0.0,0.0,0.0,0.0,0.0,0.0],[100,80.0,40.0,50.0,50.0,50.0]

        elif Datain in ["ERA20C"]:
            if varin in ['TDensity']:
                cbmin,cbmax = [0.0,90.0,0.0,0.0,0.0,0.0],[15,100.0,4.0,3.0,2.0,4.0]
            elif varin in ['TPrecip']:
                cbmin,cbmax =[0.0,0.0,0.0,0.0,0.0,0.0],[100,70.0,30.0,50.0,50.0,50.0]

        elif Datain in ["CESM"]:
            if varin in ['TDensity']:
                cbmin,cbmax = [0.0,90.0,0.0,0.0,0.0],[50,100.0,5.0,3.0,0.5]
            elif varin in ['TPrecip']:
                cbmin,cbmax = [0.0,0.0,0.0,0.0,0.0],[100,90.0,40.0,40.0,40.0]  

    elif splittypein == 'maxspeed':
        if Datain == "TRMM":
            if varin in ['TDensity']:
                cbmin,cbmax = ([0.0,0,0,0,95,0,0,0,0],
                               [300,0.5,0.75,5.0,100.0,5.0,0.75,0.5,0.5])
            elif varin in ['TPrecip']:
                cbmin,cbmax = ([0,0,0,0,20,0,0,0,0],
                              [100,30,30,30,100,30,30,30,30])

        elif Datain in ["TRMMERAIgd"]:
            if varin in ['TDensity']:
                cbmin,cbmax  = [0.0,0.0,0.0,90.0,0.0,0.0,0.0],[0.5,1.0,5.0,100.0,5.0,1.0,0.5]
            elif varin in ['TPrecip']:
                cbmin,cbmax = [0.0,0.0,0.0,0.0,0.0,0.0,0.0],[30,30.0,30.0,100.0,30.0,30.0,30.0]

        elif Datain in ["ERAI"]:
            if varin in ['TDensity']:
                cbmin,cbmax = [0.0,0.0,0.0,90.0,0.0,0.0,0.0],[0.5,1.0,5.0,100.0,5.0,1.0,0.5]
            elif varin in ['TPrecip']:
                cbmin,cbmax = [0.0,0.0,0.0,0.0,0.0,0.0,0.0],[30,30.0,30.0,100.0,30.0,30.0,30.0]

        elif Datain == "CESM":
            if varin in ['TDensity']:
                cbmin,cbmax = [0.0,90.0,0.0,0.0,0.0],[50,100.0,5.0,3.0,0.5]
            elif varin in ['TPrecip']:
                cbmin,cbmax = [0.0,0.0,0.0,0.0,0.0],[20,90.0,40.0,40.0,40.0]

        elif Datain in ["ERA20C"]:
            if varin in ['TDensity']:
                cbmin,cbmax = [0.0,90.0,0.0,0.0,0.0,0.0],[15,100.0,4.0,3.0,2.0,4.0]
            elif varin in ['TPrecip']:
                cbmin,cbmax = [0.0,0.0,0.0,0.0,0.0,0.0],[100,70.0,30.0,50.0,50.0,50.0]

    if (cbmin,cbmax) == (-1,-1):
        print(Datain,varin,splittypein,minGBin)
        sys.exit('the combination of input parameters is not yet assigned a ' +
              'colorbar. Please go to rhwhitepackages/plots and add the ' +
              'colorbar values that you want')

    return(cbmin,cbmax)

