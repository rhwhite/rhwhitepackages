# Module to read and write using xray with useful error messages
# Written by rhwhite rachel.white@cantab.net
from netCDF4 import Dataset
import netCDF4
import numpy as np
import datetime as dt
import pandas
import xray
import Ngl

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
        resMP.mpOutlineBoundarySets = "AllBoundaries"

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
            lons,lats,
            minlon,maxlon,minlat,maxlat,
            FillValue):

    nplots = plotvars1.shape[0]
    wkres = Ngl.Resources()
    wkres.wkColorMap = "precip_diff_12lev"
    wks_type = "eps"
    wks = Ngl.open_wks(wks_type,figtitle,wkres)

    # if lons start negative, shift everything over so there isn't a line down
    # the middle of the Pacific
    if lons[0] < 0:
        nlonhalf = nlons/2
        lonsnew = np.zeros(lons.shape,np.float)
        lonsnew[0:nlonhalf] = lons[nlonhalf:nlons]
        lonsnew[nlonhalf:nlons] = lons[0:nlonhalf] + 360.0
        lons = lonsnew

        for iplot in range(0,nplots):
            plotvars1[iplot] = shiftlons(plotvars1[iplot],lons)
            plotvars2[iplot] = shiftlons(plotvars2[iplot],lons)
    else:
        lonsnew = lons

    # initialize plotting resources
    res = Ngl.Resources()
    res = initcontourplot(res,minlat,minlon,maxlat,maxlon,lats,lonsnew)
    res.sfMissingValueV = FillValue

    #    res.pmLabelBarDisplayMode = "Always"

    #res.sfXCStartV = float(lonsnew[0])
    #res.sfXCEndV = float(lonsnew[len(lons)-1])
    #res.sfYCStartV = float(lats[0])
    #res.sfYCEndV = float(lats[len(lats)-1])

    res.lbOrientation   = "Vertical"

    # including some font heights
    res.lbLabelFontHeightF = 0.0125
    res.lbTitleFontHeightF = 0.0125
    res.tiMainFontHeightF = 0.015

    # initialize plotting array
    toplot = []
    # fill plotting array
    for iplot in range(0,nplots):
        tempplot = plotvars1[iplot]
        tempplot[np.where(np.isnan(tempplot))] = FillValue
        res.cnMinLevelValF       = plotmin1[iplot]          # contour levels.
        res.cnMaxLevelValF       = plotmax1[iplot]
        res.cnLevelSpacingF      = ((plotmax1[iplot]-plotmin1[iplot])/10.0)
        res.tiMainString = (vartitle1[iplot])
        toplot.append(Ngl.contour_map(wks,tempplot,res))

        tempplot = plotvars2[iplot]
        tempplot[np.where(np.isnan(tempplot))] = FillValue
        res.cnMinLevelValF       = plotmin2[iplot]          # contour levels.
        res.cnMaxLevelValF       = plotmax2[iplot]
        res.cnLevelSpacingF      = ((plotmax2[iplot]-plotmin2[iplot])/10.0)
        res.tiMainString = vartitle2[iplot]
        toplot.append(Ngl.contour_map(wks,tempplot,res))

    textres = Ngl.Resources()
    textres.txFontHeightF = 0.015
    Ngl.text_ndc(wks,title,0.5,0.87,textres)

    panelres = Ngl.Resources()
    panelres.nglPanelLabelBar = True
    #panelres.nglPanelYWhiteSpacePercent = 5.
    #panelres.nglPanelXWhiteSpacePercent = 5.

    panelres.nglPanelLabelBar   = False     # Turn on panel labelbar
    if nplots > 5:
        panelres.nglPanelTop                      = 0.8
        panelres.nglPanelBottom                      = 0.15
    else:
        panelres.nglPanelTop                      = 0.95
        panelres.nglPanelBottom                      = 0.01

    panelres.nglPanelFigureStrings = ['a','b','c','d','e','f','g','h']
    panelres.nglPanelFigureStringsJust = "TopLeft"
    panelres.nglPanelFigureStringsFontHeightF = 0.008
    panelres.nglPanelFigureStringsParallelPosF = -0.55
    panelres.nglPanelFigureStringsOrthogonalPosF = -0.6
    #panelres.amJust = "TopLeft"

    panelres.nglPaperOrientation = "Auto"

    plot = Ngl.panel(wks,toplot,[nplots,2],panelres)
