# -*- coding: utf-8 -*-
"""
Created on Sat Mar 19 2016

@author: rachel, rachel.white@cantab.net

Note that pressure is non-dimensional, i.e. p = p/1000mb

Quasi-geostrophic formulation. Could calculate quasi-geostrophic streamfunction
from geopotential height. But we've shown that the solution is very close to being geostrophic, 
so it should be fine to use streamfunction calculated from wind components.

"""

import os, errno
from netCDF4 import Dataset
import netCDF4
import numpy as np
import datetime as dt
import Ngl

ideal = 0
addtitle = "Pr2"

timespan = "DJF"
startyr = 2
nyears = 30

timing = ["Stationary","Stationary","Stationary","Stationary"]
#timing = ["Stationary","Transients","Total","Stationary","Transients","Total"]

if ideal == 1:
        Experiments1 = ["CESM_IG54", "CESM_IG44", "CESM_IG29"]
        Experiments2 = ["CESMnotopof19","CESMnotopof19","CESMnotopof19","CESMnotopof19"]
        Diffs = [1,1,1,1]	
	Titles1 = ["No Tibetan plateau","No Mongolian plateau"]
#	Titles1 = ["Idealized Mongolia 63N", "Idealized Mongolia, 53N","Idealized Mongolia 38N"]        
        Titles2 = ["Flat","Flat","Flat","Flat"]
	FigTitle = addtitle + "_U_EPfluxes_zm_Ideal_" + str(nplots) + "_" + timespan + "_" + str(startyr) + "_" + str(nyears) + "_" + str(startlon) + "E-" + str(endlon) + "E"
else:
	Experiments1 = ["CESMtopof19", "CESMtopof19"] #, "CESMtopof19", "CESMtopof19","CESMtopof19", "CESMtopof19"]
	Experiments2 = ["CESMnoTf19","CESMnoT4f19"] #,"CESMnoT2f19","CESMnoT2f19","CESMnoT2f19","CESMnoT2f19"]
	Diffs = [1,1,1,1,1,1]
	Titles1 = ["CTL","CTL","CTL","CTL","CTL","CTL"]
	Titles2 = ["No Tibetan plateau","No Mongolian plateau","NoMT","NoMT","NoMT","NoMT"]
	FigTitle = addtitle + "_U_EPfluxes_zm_Realistic_" + str(nplots) + "_" + timespan + "_" + str(startyr) + "_" + str(nyears) + "_" + str(startlon) + "E-" + str(endlon) + "E"

nexps = len(Experiments1)

DirD = "/home/disk/eos4/rachel/CESM_outfiles/"
DirS = "/home/disk/rachel/CESM_outfiles/" 

slat = 20
elat = 85

FigDir = "/home/disk/eos4/rachel/Figures/Mongolia/"

a = 6.37122e06
PI = 3.14159
P0 = 100000.0	# "Pa"
secs = 86400.0

plot = []
plotdiv = []

def Calc_TakNak(Psi,Psidev,TH,N,p,U,V,lons,lats,z):

# Calculate magnitude of climatological wind
Umag = np.sqrt(np.pow(U,2) + np.pow(V,2))

# Calculate dthetadz
dthetadz = np.zeros(dimsizes(TH),typeof(TH))
for ilat in range(0,len(lats)):
	for ilon in range(0,len(lons)):
		dthetadz(...,ilat,ilon) = np.array(np.gradient(TH(...,ilat,ilon),z(...,ilat,ilon)


for iexp in range(0,nexps):

	for ivar in plotvars:
		print ivar
		
		for ii in range(0,2):
			if ii == 0:
				Experiment = Experiments1[iexp]
			elif ii == 1:
				Experiment = Experiments2[iexp]

			print Experiment
			

			if ivar == "EP": 
				if timing[iexp] != "Stationary":
					filenameD = DirD + Experiment + "/atm/hist/" + "EPfluxes_Daily_" + timespan + str(startyr) + "_" + str(nyears) + "_" + str(startlon) + "E-" + str(endlon) + "E_" + Experiment + ".cam2.h0.nc"
					print filenameD
					fileInD = Dataset(filenameD,'r')
					FphiD = fileInD.variables['Fphi_int'][...]
					FpD = fileInD.variables['Fp_int'][...]
			
					EPdivD = fileInD.variables['EPdiv_int'][...]
			
				filenameS = DirS + Experiment + "/atm/hist/" + "EPfluxesPr_" + timespan + str(startyr) + "_" + str(nyears) + "_" + str(startlon) + "E-" + str(endlon) + "E_" + Experiment + ".cam2.h0.nc"
				print filenameS
				fileInS = Dataset(filenameS,'r')

				lats = fileInS.variables['lat']
				levels = fileInS.variables['level_int']
				lev = levels[...]
				nlats = len(lats)
				nlevels = len(levels)
				if levels[0] > levels[nlevels-1]:
					levground = 0
				else:
					levground = nlevels-1

				print levels

				FphiS = fileInS.variables['Fphi_int'][...]
				FpS = fileInS.variables['Fp_int'][...]

				EPdivS = fileInS.variables['EPdivBEH_int'][...]

				if timing[iexp] == "Transients":
					FphiT = FphiD - FphiS
					FpT = FpD - FpS
					EPdivT = EPdivD - EPdivS

				if timing[iexp] == 'Transients':
					EPdiv = EPdivT
					Fp = FpT
					Fphi = FphiT
				elif timing[iexp] == 'Total':
					EPdiv = EPdivD
					Fp = FpD
					Fphi = FphiD	
				elif timing[iexp] == 'Stationary':
					EPdiv = EPdivS
					Fp = FpS
					Fphi = FphiS

				phi = lats[:] * PI / 180.0
				cphi = np.cos(phi)
				acphi = a * np.cos(phi)
				asphi = a * np.sin(phi)
				#Fphicphi include a factor of cos(phi) for graphical display of arrows, al la Edmon et al 1980

				# Scale according to Edmon et al (and NCL script)
				# If plotting in a log-p display arrows may look divergent when they are not.
				# So try plotting in a p display


				# Scale EPdiv according to Baldwin et al. 1985: plotting Div.F /(acos(phi) e -Z/H) in m2/(s day)
				EPdivscaled = EPdiv * secs / ((acphi) * lev[:,None]/1000.0)  #
				Fpsc = Fp * cphi[None,:]        # using broadcasting to copy cphi into the other dimensions of Fp
				Fphisc = Fphi * cphi[None,:] / a
				
				# Now scale by the relative ranges of the two axes of the plot, Pi/2 radians by 10^5 Pa
				if HScale == 1:
					Fpsc = Fpsc / P0
					Fphisc = Fphisc / (PI/2.0)

				# Now scale by the squareroot of the pressure
				if VScale == 1:
					rhofac = np.sqrt(P0/(levels[:] * 100.0))
					Fpsc = Fpsc * rhofac[:,None]
					Fphisc = Fphisc * rhofac[:,None]
		 
				# Could also scale by a magnification factor above 100hPa, but not so interested in heights above 100hpa

				if ii == 0:
					Fpsc1 = Fpsc
					Fphisc1 = Fphisc
					EPdiv1 = EPdivscaled
				elif ii == 1:
					Fpsc2 = Fpsc
					Fphisc2 = Fphisc
					EPdiv2 = EPdivscaled

			elif ivar == "U":
				for ii in range(0,2):
					if ii == 0:
						Experiment = Experiments1[iexp]
					elif ii == 1:
						Experiment = Experiments2[iexp]

					print Experiment
					filename = DirS + Experiment + "/atm/hist/" + "EPfluxesPr_" + timespan + str(startyr) + "_" + str(nyears) + "_" + str(startlon) + "E-" + str(endlon) + "E_" + Experiment + ".cam2.h0.nc"
					fileIn = Dataset(filename,'r')

					lats = fileIn.variables['lat']
					levels = fileIn.variables['level_int']
					nlats = len(lats)
					nlevels = len(levels)
					if levels[0] > levels[nlevels-1]:
						levground = 0
					else:
						levground = nlevels-1

					U = fileIn.variables['Uzm_int'][...]

#					for ilat in range(0,nlats)
#						for ip in range(0,nlevels)
							

					phi = lats[:] * PI / 180.0

					if ii == 0:
						plot1 = U
					elif ii == 1:
						plot2 = U
		# Now take the differences
		if Diffs[iexp] == 0:
			print "don't take differences"	 
		
			if ivar == "EP":	
				EPdiff = EPdiv1
				Fpscdiff = Fpsc1
				Fphiscdiff = Fphisc1

			elif ivar == "U":
				plotvar = plot1
		else:
			if ivar == "EP":	
				Fpscdiff = Diffs[iexp] * (Fpsc1 - Fpsc2)
				Fphiscdiff = Diffs[iexp] * (Fphisc1 - Fphisc2)

				EPdiff = Diffs[iexp] * (EPdiv1 - EPdiv2)
			elif ivar == "U":
				plotvar = Diffs[iexp] * (plot1 - plot2)

		res = Ngl.Resources()
		res.nglDraw  = False
		res.nglFrame = False
		res.vcRefLengthF = 0.05
		res.tiXAxisString = "Latitude"
		res.tiYAxisString = "pressure (mb)"
		res.trYReverse = True
		res.vcMonoLineArrowColor = True
		#res.vcLevelSelectionMode = "ExplicitLevels"
		res.pmLabelBarDisplayMode = "Never"
		res.lbPerimOn = False
		res.trYMaxF = startp
		res.trYMinF = endp
		res.trXMaxF = elat
		res.trXMinF = slat
		
		res_con = Ngl.Resources()
		res_con.tiXAxisString = "Latitude"
		res_con.tiYAxisString = "pressure (mb)"
		res_con.trYReverse = True	
		res_con.nglDraw = False
		res_con.nglFrame = False
		res_con.cnFillOn = True
		res_con.cnMonoFillPattern = True
		res_con.cnMonoFillColor = False
		res_con.cnLineLabelsOn = False
		res_con.cnLinesOn = False
		res_con.trYReverse = True
		
		#res_con.cnSmoothingOn = True
		res_con.cnLineLabelsOn = False
		res_con.cnLevelSelectionMode = "ManualLevels"

		res_con.trYMaxF = startp
		res_con.trYMinF = endp
		res_con.trXMaxF = elat
		res_con.trXMinF = slat


		levend = 0

		res.vfYArray = levels[:]	#[levend:nlevels]
#		res.vfXArray = lats[:]		#[nlats/2+5:nlats-2]
		
		res_con.sfYArray = levels[:]	#[levend:nlevels]
		res_con.sfXArray = lats[:]		#[nlats/2+5:nlats-2]
		#plotvec = Ngl.vector(wks,Fphiscdiff[levend:nlevels,nlats/2+5:nlats-2],Fpscdiff[levend:nlevels,nlats/2+5:nlats-2],res)	

		if ivar == "EP":
			if ideal == 1:
				res.vcRefMagnitudeF =10
			else:
                                res.vcRefMagnitudeF =20

			# Take out every other latitude
			latsnew = lats[0:nlats:2]
			print latsnew
			Fphiscdiffnew = Fphiscdiff[:,0:nlats:2]
			print Fphiscdiff.shape
			print Fphiscdiffnew.shape
			Fpscdiffnew = Fpscdiff[:,0:nlats:2] 
			res.vfXArray = latsnew[:]  

			#Overplot divergence
			# Hide 1000 mb and top level
			#EPdiff[levground,:] = float('nan')

			if ideal == 1:
				res_con.cnMinLevelValF = -1.0
				res_con.cnMaxLevelValF = 1.0
				res_con.cnLevelSpacingF = 0.2
			else:
                                res_con.cnMinLevelValF = -2.0
                                res_con.cnMaxLevelValF = 2.0
                                res_con.cnLevelSpacingF = 0.4
                        res_con.cnFillOn = True
                        res_con.cnMonoFillPattern = True
                        res_con.cnMonoFillColor = False
                        res_con.cnLineLabelsOn = False

#			res_con.sfXArray = latsnew[:] 

			plotcon = (Ngl.contour(wks,EPdiff[:,:],res_con))
                        plotvec = Ngl.vector(wks,Fphiscdiffnew[:,:],Fpscdiffnew[:,:],res)
	
		#	neg_dash_contours(plotcon)
			Ngl.overlay(plotcon,plotvec)
		
			plot.append(plotcon)

		elif ivar == "U":
			if ideal == 1:
				res_con.cnMinLevelValF = -2.5
				res_con.cnMaxLevelValF = 2.5
				res_con.cnLevelSpacingF = 0.5
			else:
                                res_con.cnMinLevelValF = -5.0
                                res_con.cnMaxLevelValF = 5.0
                                res_con.cnLevelSpacingF = 1.0
			res_con.cnFillOn = True
			res_con.cnMonoFillPattern = True
			res_con.cnMonoFillColor = False
			res_con.cnLineLabelsOn = False

		        plot.append(Ngl.contour(wks,plotvar[:,:],res_con))


print "now here"
textres = Ngl.Resources()
textres.txFontHeightF = 0.012

#Ngl.text_ndc(wks,"Realistic and Ideal experiment diffs, zonal mean: " + str(startlon) + "E to " + str(endlon) + "E",0.5,0.87,textres)
#Ngl.text_ndc(wks,timing,0.5,0.84,textres)


if nplots == 4:
	Ngl.text_ndc(wks,"a.",0.07,0.88,textres)
	Ngl.text_ndc(wks,"b.",0.07,0.50,textres)
	Ngl.text_ndc(wks,"c.",0.57,0.88,textres)
	Ngl.text_ndc(wks,"d.",0.57,0.50,textres)
elif nplots == 3:
        Ngl.text_ndc(wks,"a.",0.34,0.97,textres)
        Ngl.text_ndc(wks,"b.",0.34,0.65,textres)
        Ngl.text_ndc(wks,"c.",0.34,0.33,textres)
elif nplots == 6:
	textres.txFontHeightF = 0.010

        Ngl.text_ndc(wks,"a.",0.24,0.85,textres)
        Ngl.text_ndc(wks,"b.",0.24,0.615,textres)
        Ngl.text_ndc(wks,"c.",0.24,0.385,textres)
        Ngl.text_ndc(wks,"d.",0.54,0.85,textres)
        Ngl.text_ndc(wks,"e.",0.54,0.615,textres)
        Ngl.text_ndc(wks,"f.",0.54,0.385,textres)
panelres = Ngl.Resources()
panelres.nglPanelLabelBar                 = False     # Turn on panel labelbar
panelres.nglPaperOrientation = "Portrait"
panelres.nglPanelYWhiteSpacePercent = 5.
panelres.nglPanelXWhiteSpacePercent = 5.

if nplots == 1:
	textres.txFontHeightF = 0.015

	Ngl.panel(wks,plot,[1,1],panelres)
elif nplots == 2:
        textres.txFontHeightF = 0.008

        Ngl.panel(wks,plot,[1,2],panelres)
elif nplots == 3:
        textres.txFontHeightF = 0.008

        Ngl.panel(wks,plot,[3,1],panelres)

elif nplots == 4:
	textres.txFontHeightF = 0.008

        Ngl.panel(wks,plot,[2,2],panelres)

elif nplots == 6:
        textres.txFontHeightF = 0.008
	panelres.nglPanelTop = 0.85
	panelres.nglPanelBottom = 0.15
        Ngl.panel(wks,plot,[3,2],panelres)


Ngl.end()

