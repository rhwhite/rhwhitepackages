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
        resMP.cnLineLabelsOn = False
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


def plotline(nlines,plotin,yearsin,nsep,title,Datatitle,figtitlein.colormap,xString,yString,unit):
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

