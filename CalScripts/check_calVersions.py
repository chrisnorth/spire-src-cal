#check calibration versions

#output is the min/max of (1 - ver1/ver2) for a range of parameters

#set to also plot some plots
doPlot=True

#set version numbers (type = "file" or "cal")
ver1={'name':"v4",'type':"file"}
ver2={'name':"spire_cal_12_3",'type':"cal"}

#set version numbers for specific files
verBeamProf=[ver1,ver2]
verKBeam=[ver1,ver2]
verKPsrc=[ver1,ver2]
verKExtd=[ver1,ver2]
verApCorr=[ver1,ver2]
verKHfi=[ver1,ver2]

directory=Configuration.getProperty('var.hcss.workdir')

spireBands=["PSW","PMW","PLW"]

#-------------------------------------------------------------------------------
# RadialCorrBeam
if verBeamProf[0]['type']=="file":
    beamProf1=fitsReader('%s//Phot//SCalPhotRadialCorrBeam//SCalPhotRadialCorrBeam_%s.fits'%(directory,verBeamProf[0]['name']))
else:
    beamprof1=spireCal(verBeamProf[0]['name']).getPhot().getRadialCorrBeam()
if verBeamProf[1]['type']=='file':
    beamProf2=fitsReader('%s//Phot//SCalPhotRadialCorrBeam//SCalPhotRadialCorrBeam_%s.fits'%(directory,verBeamProf[1]['name']))
else:
    beamProf2=spireCal(verBeamProf[1]['name']).getPhot().getRadialCorrBeam()

print 'RadialCorrBeam (1-%s/%s):'%(verBeamProf[1]['name'],verBeamProf[0]['name'])
for band in spireBands:
	minLen=MIN([len(beamProf1['core'][band].data),len(beamProf2['core'][band].data)])
	print '  core     (%s): min=%g, max=%g'%\
		(band,min(1.-beamProf2['core'][band].data[0:minLen]/beamProf1['core'][band].data[0:minLen]),\
		max(1.-beamProf2['core'][band].data[0:minLen]/beamProf1['core'][band].data[0:minLen]))
	print '  constant (%s): min=%g, max=%g'%\
		(band,min(1.-beamProf2['constant'][band].data[0:minLen]/beamProf1['constant'][band].data[0:minLen]),\
		max(1.-beamProf2['constant'][band].data[0:minLen]/beamProf1['constant'][band].data[0:minLen]))
	print '  normArea (%s): min=%g, max=%g'%\
		(band,min(1.-beamProf2['normArea'][band].data[0:minLen]/beamProf1['normArea'][band].data[0:minLen]),\
		max(1.-beamProf2['normArea'][band].data[0:minLen]/beamProf1['normArea'][band].data[0:minLen]))

#-------------------------------------------------------------------------------
# ColorCorrBeam
if verKBeam[0]['type']=='file':
    Kbeam1=fitsReader('%s//Phot//SCalPhotColorCorrBeam//SCalPhotColorCorrBeam_%s.fits'%(directory,verKBeam[0]['name']))
else:
    Kbeam1=spireCal(verKBeam[0]['name']).getPhot().getColorCorrBeam()
if verKBeam[1]['type']=='file':
    Kbeam2=fitsReader('%s//Phot//SCalPhotColorCorrBeam//SCalPhotColorCorrBeam_%s.fits'%(directory,verKBeam[1]['name']))
else:
    Kbeam2=spireCal(verKBeam[1]['name']).getPhot().getColorCorrBeam()

print 'ColorCorrBeam (1-%s/%s):'%(verKBeam[1]['name'],verKBeam[0]['name'])
for band in spireBands:
	print '  alpha     (%s): min=%g, max=%g'%\
		(band,min(1.-Kbeam2['alpha'][band].data/Kbeam1['alpha'][band].data),\
		max(1.-Kbeam2['alpha'][band].data/Kbeam1['alpha'][band].data))
	print '  beta_2_00 (%s): min=%g, max=%g'%\
		(band,min(1.-Kbeam2['beta_2_00'][band].data/Kbeam1['beta_2_00'][band].data),\
		max(1.-Kbeam2['beta_2_00'][band].data/Kbeam1['beta_2_00'][band].data))

#-------------------------------------------------------------------------------
# ColorCorrK_point
if verKPsrc[0]['type']=='file':
    KPsrc1=fitsReader('%s//Phot//SCalPhotColorCorrK//SCalPhotColorCorrK_point_%s.fits'%(directory,verKPsrc[0]['name']))
else:
    KPsrc1=spireCal(verKPsrc[0]['name']).getPhot().getColorCorrKList().getProduct('point')
if verKPsrc[1]['type']=='file':
    KPsrc2=fitsReader('%s//Phot//SCalPhotColorCorrK//SCalPhotColorCorrK_point_%s.fits'%(directory,verKPsrc[1]['name']))
else:
    KPsrc2=spireCal(verKPsrc[1]['name']).getPhot().getColorCorrKList().getProduct('point')

print 'ColorCorrK_point (1-%s/%s):'%(verKPsrc[1]['name'],verKPsrc[0]['name'])
for band in spireBands:
	print '  alpha     (%s): min=%g, max=%g'%\
		(band,min(1.-KPsrc2['alpha'][band].data/KPsrc1['alpha'][band].data),\
		max(1.-KPsrc2['alpha'][band].data/KPsrc1['alpha'][band].data))
	print '  beta_2_00 (%s): min=%g, max=%g'%\
		(band,min(1.-KPsrc2['beta_2_00'][band].data/KPsrc1['beta_2_00'][band].data),\
		max(1.-KPsrc2['beta_2_00'][band].data/KPsrc1['beta_2_00'][band].data))

#-------------------------------------------------------------------------------
# ColorCorrK_extended
if verKExtd[0]['type']=='file':
    KExtd1=fitsReader('%s//Phot//SCalPhotColorCorrK//SCalPhotColorCorrK_extended_%s.fits'%(directory,verKExtd[0]['name']))
else:
    KExtd1=spireCal(verKExtd[0]['name']).getPhot().getColorCorrKList().getProduct('extended')
if verKExtd[1]['type']=='file':
    KExtd2=fitsReader('%s//Phot//SCalPhotColorCorrK//SCalPhotColorCorrK_extended_%s.fits'%(directory,verKExtd[1]['name']))
else:
    KExtd2=spireCal(verKExtd[1]['name']).getPhot().getColorCorrKList().getProduct('extended')

print 'ColorCorrK_extended (1-%s/%s):'%(verKExtd[1]['name'],verKExtd[0]['name'])
for band in spireBands:
	print '  alpha     (%s): min=%g, max=%g'%\
		(band,min(1.-KExtd2['alpha'][band].data/KExtd1['alpha'][band].data),\
		max(1.-KExtd2['alpha'][band].data/KExtd1['alpha'][band].data))
	print '  beta_2_00 (%s): min=%g, max=%g'%\
		(band,min(1.-KExtd2['beta_2_00'][band].data/KExtd1['beta_2_00'][band].data),\
		max(1.-KExtd2['beta_2_00'][band].data/KExtd1['beta_2_00'][band].data))

#-------------------------------------------------------------------------------
# ColorCorrAperture_noBG
if verApCorr[0]['type']=='file':
    apCorrNoBG1=fitsReader('%s//Phot//SCalPhotColorCorrAperture//SCalPhotColorCorrAperture_noBG_%s.fits'%(directory,verApCorr[0]['name']))
else:
    apCorrNoBG1=spireCal(verApCorr[0]['name']).getPhot().getColorCorrApertureList().getProduct('noBG')
if verApCorr[1]['type']=='file':
    apCorrNoBG2=fitsReader('%s//Phot//SCalPhotColorCorrAperture//SCalPhotColorCorrAperture_noBG_%s.fits'%(directory,verApCorr[1]['name']))
else:
    apCorrNoBG2=spireCal(verApCorr[1]['name']).getPhot().getColorCorrApertureList().getProduct('noBG')

print 'ColorCorrAperture_noBG (1-%s/%s):'%(verApCorr[1]['name'],verApCorr[0]['name'])
for band in spireBands:
	print '  alpha (%s): min=%g, max=%g'%\
		(band,min(1.-apCorrNoBG2['alpha'][band].data/apCorrNoBG1['alpha'][band].data),\
		max(1.-apCorrNoBG2['alpha'][band].data/apCorrNoBG1['alpha'][band].data))

#-------------------------------------------------------------------------------
# ColorCorrAperture_incBG
if verApCorr[0]['type']=='file':
    apCorrIncBG1=fitsReader('%s//Phot//SCalPhotColorCorrAperture//SCalPhotColorCorrAperture_noBG_%s.fits'%(directory,verApCorr[0]['name']))
else:
    apCorrIncBG1=spireCal(verApCorr[0]['name']).getPhot().getColorCorrApertureList().getProduct('incBG')
if verApCorr[1]['type']=='file':
    apCorrIncBG2=fitsReader('%s//Phot//SCalPhotColorCorrAperture//SCalPhotColorCorrAperture_noBG_%s.fits'%(directory,verApCorr[1]['name']))
else:
    apCorrIncBG2=spireCal(verApCorr[1]['name']).getPhot().getColorCorrApertureList().getProduct('incBG')

print 'ColorCorrAperture_incBG (1-%s/%s):'%(verApCorr[1]['name'],verApCorr[0]['name'])
for band in spireBands:
	print '  alpha (%s): min=%g, max=%g'%\
		(band,min(1.-apCorrIncBG2['alpha'][band].data/apCorrIncBG1['alpha'][band].data),\
		max(1.-apCorrIncBG2['alpha'][band].data/apCorrIncBG1['alpha'][band].data))

#-------------------------------------------------------------------------------
# ColorCorrHfi
if verKHfi[0]['type']=='file':
    KHfi1=fitsReader('%s//Phot//SCalPhotColorCorrHfi//SCalPhotColorCorrHfi_%s.fits'%(directory,verKHfi[0]['name']))
else:
    KHfi1=spireCal(verKHfi[0]['name']).getPhot().getColorCorrHfi()
if verKHfi[1]['type']=='file':
    KHfi2=fitsReader('%s//Phot//SCalPhotColorCorrHfi//SCalPhotColorCorrHfi_%s.fits'%(directory,verKHfi[1]['name']))
else:
    KHfi2=spireCal(verKHfi[1]['name']).getPhot().getColorCorrHfi()

KHfi1i=KHfi1.copy()
KHfi2i=KHfi2.copy()
if len(KHfi2['colorCorr']['Temperature'].data) != len(KHfi1['colorCorr']['Temperature'].data):
	#interpolate ver 1 to match ver2 Temperature grid
	temp1=KHfi1['colorCorr']['Temperature'].data
	temp2=KHfi2['colorCorr']['Temperature'].data
	ix=temp2.where((temp2 >= MIN(temp1)) & (temp2 <= MAX(temp1)))
	KHfi1r_interp=CubicSplineInterpolator(KHfi1['colorCorr']['Temperature'].data,\
		KHfi1['colorCorr']['ratio545_857'].data)
	KHfi1plw_interp=CubicSplineInterpolator(KHfi1['colorCorr']['Temperature'].data,\
		KHfi1['colorCorr']['k545toPLW'].data)
	KHfi1pmw_interp=CubicSplineInterpolator(KHfi1['colorCorr']['Temperature'].data,\
		KHfi1['colorCorr']['k857toPMW'].data)
	KHfi1psw_interp=CubicSplineInterpolator(KHfi1['colorCorr']['Temperature'].data,\
		KHfi1['colorCorr']['k857toPSW'].data)
	KHfi1i['colorCorr']['Temperature'].data=KHfi2['colorCorr']['Temperature'].data[ix]
	KHfi1i['colorCorr']['ratio545_857'].data=KHfi1r_interp(KHfi1i['colorCorr']['Temperature'].data)
	KHfi1i['colorCorr']['k545toPLW'].data=KHfi1plw_interp(KHfi1i['colorCorr']['Temperature'].data)
	KHfi1i['colorCorr']['k857toPMW'].data=KHfi1pmw_interp(KHfi1i['colorCorr']['Temperature'].data)
	KHfi1i['colorCorr']['k857toPSW'].data=KHfi1psw_interp(KHfi1i['colorCorr']['Temperature'].data)

	KHfi2i['colorCorr']['Temperature'].data=KHfi2['colorCorr']['Temperature'].data[ix]
	KHfi2i['colorCorr']['ratio545_857'].data=KHfi2['colorCorr']['ratio545_857'].data[ix]
	KHfi2i['colorCorr']['k545toPLW'].data=KHfi2['colorCorr']['k545toPLW'].data[ix]
	KHfi2i['colorCorr']['k857toPMW'].data=KHfi2['colorCorr']['k857toPMW'].data[ix]
	KHfi2i['colorCorr']['k857toPSW'].data=KHfi2['colorCorr']['k857toPSW'].data[ix]
	print 'ColorCorrHfi (1-%s/%s) [v%s interpolated to match v%s]:'% \
	  (verKHfi[1]['name'],verKHfi[0]['name'],verKHfi[0]['name'],verKHfi[1]['name'])
else:
	print 'ColorCorrHfi (1-%s/%s):'%(verKHfi[1]['name'],verKHfi[0]['name'])

print '  ratio545_857: min=%g, max=%g'%\
	(min(1.-KHfi2i['colorCorr']['ratio545_857'].data/KHfi1i['colorCorr']['ratio545_857'].data),\
	max(1.-KHfi2i['colorCorr']['ratio545_857'].data/KHfi1i['colorCorr']['ratio545_857'].data))
print '  k545toPLW   : min=%g, max=%g'%\
	(min(1.-KHfi2i['colorCorr']['k545toPLW'].data/KHfi1i['colorCorr']['k545toPLW'].data),\
	max(1.-KHfi2i['colorCorr']['k545toPLW'].data/KHfi1i['colorCorr']['k545toPLW'].data))
print '  k857toPMW   : min=%g, max=%g'%\
	(min(1.-KHfi2i['colorCorr']['k857toPMW'].data/KHfi1i['colorCorr']['k857toPMW'].data),\
	max(1.-KHfi2i['colorCorr']['k857toPMW'].data/KHfi1i['colorCorr']['k857toPMW'].data))
print '  k857toPSW   : min=%g, max=%g'%\
	(min(1.-KHfi2i['colorCorr']['k857toPSW'].data/KHfi1i['colorCorr']['k857toPSW'].data),\
	max(1.-KHfi2i['colorCorr']['k857toPSW'].data/KHfi1i['colorCorr']['k857toPSW'].data))


if not doPlot:
    #causes error
    stopHere
    

cols={'PSW':java.awt.Color.BLUE,'PMW':java.awt.Color.GREEN,'PLW':java.awt.Color.RED}
pCore=PlotXY()
for band in spireBands:
	minLen=MIN([len(beamProf1['core'][band].data),len(beamProf2['core'][band].data)])
	pCore.addLayer(LayerXY(beamProf2['core']['radius'].data[:minLen],\
		1.-beamProf2['core'][band].data[:minLen]/beamProf1['core'][band].data[:minLen], \
		color=cols[band],name='%s'%(band)))
	pCore.setTitleText('Beam Profile (core)')
	pCore.yaxis.titleText = '1-%s/%s'%(verBeamProf[1]['name'],verBeamProf[0]['name'])
	pCore.xaxis.titleText = 'radius (arcsec)'
	pCore.setSubtitleText('MIN=%.5g MAX=%.5g'%\
		(min(1.-beamProf2['core'][band].data[:minLen]/beamProf1['core'][band].data[:minLen]),\
		max(1.-beamProf2['core'][band].data[:minLen]/beamProf1['core'][band].data[:minLen])))

pConst=PlotXY()
for band in spireBands:
	minLen=MIN([len(beamProf1['core'][band].data),len(beamProf2['core'][band].data)])
	pConst.addLayer(LayerXY(beamProf2['constant']['radius'].data[:minLen],\
		1.-beamProf2['constant'][band].data[:minLen]/beamProf1['constant'][band].data[:minLen], \
		color=cols[band],name='%s'%(band)))
	pConst.setTitleText('Beam Profile (constant)')
	pConst.yaxis.titleText = '1-%s/%s'%(verBeamProf[1]['name'],verBeamProf[0]['name'])
	pConst.xaxis.titleText = 'radius (arcsec)'
	pConst.setSubtitleText('MIN=%.5g MAX=%.5g'%\
		(min(1.-beamProf2['constant'][band].data[:minLen]/beamProf1['constant'][band].data[:minLen]),\
		max(1.-beamProf2['constant'][band].data[:minLen]/beamProf1['constant'][band].data[:minLen])))

pNorm=PlotXY()
for band in spireBands:
	minLen=MIN([len(beamProf1['core'][band].data),len(beamProf2['core'][band].data)])
	pNorm.addLayer(LayerXY(beamProf2['normArea']['radius'].data[:minLen],\
		1.-beamProf2['normArea'][band].data[:minLen]/beamProf1['normArea'][band].data[:minLen], \
		color=cols[band],name='%s'%(band)))
	pNorm.setTitleText('Beam Profile (normArea)')
	pNorm.yaxis.titleText = '1-%s/%s'%(verBeamProf[1]['name'],verBeamProf[0]['name'])
	pNorm.xaxis.titleText = 'radius (arcsec)'
	pNorm.setSubtitleText('MIN=%.5g MAX=%.5g'%\
		(min(1.-beamProf2['normArea'][band].data[:minLen]/beamProf1['normArea'][band].data[:minLen]),\
		max(1.-beamProf2['normArea'][band].data[:minLen]/beamProf1['normArea'][band].data[:minLen])))

pKBeama=PlotXY()
for band in spireBands:
	pKBeama.addLayer(LayerXY(Kbeam2['alpha']['alpha'].data,\
		1.-Kbeam2['alpha'][band].data/Kbeam1['alpha'][band].data, \
		color=cols[band],name='%s'%(band)))
	pKBeama.setTitleText('Beam Correction (alpha)')
	pKBeama.yaxis.titleText = '1-%s/%s'%(verKBeam[1]['name'],verKBeam[0]['name'])
	pKBeama.xaxis.titleText = 'Spectral index'
	pKBeama.setSubtitleText('MIN=%.5g MAX=%.5g'%\
		(min(1.-Kbeam2['alpha'][band].data/Kbeam1['alpha'][band].data),\
		max(1.-Kbeam2['alpha'][band].data/Kbeam1['alpha'][band].data)))

pKBeamt=PlotXY()
for band in spireBands:
	pKBeamt.addLayer(LayerXY(Kbeam2['beta_2_00']['Temperature'].data,\
		1.-Kbeam2['beta_2_00'][band].data/Kbeam1['beta_2_00'][band].data, \
		color=cols[band],name='%s'%(band)))
	pKBeamt.setTitleText('Beam Correction (beta=2.0)')
	pKBeamt.yaxis.titleText = '1-%s/%s'%(verKBeam[1]['name'],verKBeam[0]['name'])
	pKBeamt.xaxis.titleText = 'Temperature (K)'
	pKBeamt.setSubtitleText('MIN=%.5g MAX=%.5g'%\
		(min(1.-Kbeam2['beta_2_00'][band].data/Kbeam1['beta_2_00'][band].data),\
		max(1.-Kbeam2['beta_2_00'][band].data/Kbeam1['beta_2_00'][band].data)))

pKPsrca=PlotXY()
for band in spireBands:
	pKPsrca.addLayer(LayerXY(KPsrc2['alpha']['alpha'].data,\
		1.-KPsrc2['alpha'][band].data/KPsrc1['alpha'][band].data, \
		color=cols[band],name='%s'%(band)))
	pKPsrca.setTitleText('Point Source Correction (alpha)')
	pKPsrca.yaxis.titleText = '1-%s/%s'%(verKPsrc[1]['name'],verKPsrc[0]['name'])
	pKPsrca.xaxis.titleText = 'Spectral index'
	pKPsrca.setSubtitleText('MIN=%.5g MAX=%.5g'%\
		(min(1.-KPsrc2['alpha'][band].data/KPsrc1['alpha'][band].data),\
		max(1.-KPsrc2['alpha'][band].data/KPsrc1['alpha'][band].data)))

pKPsrct=PlotXY()
for band in spireBands:
	pKPsrct.addLayer(LayerXY(KPsrc2['beta_2_00']['Temperature'].data,\
		1.-KPsrc2['beta_2_00'][band].data/KPsrc1['beta_2_00'][band].data, \
		color=cols[band],name='%s'%(band)))
	pKPsrct.setTitleText('Point Source Correction (beta=2.0)')
	pKPsrct.yaxis.titleText = '1-%s/%s'%(verKPsrc[1]['name'],verKPsrc[0]['name'])
	pKPsrct.xaxis.titleText = 'Temperature (K)'
	pKPsrct.setSubtitleText('MIN=%.5g MAX=%.5g'%\
		(min(1.-KPsrc2['beta_2_00'][band].data/KPsrc1['beta_2_00'][band].data),\
		max(1.-KPsrc2['beta_2_00'][band].data/KPsrc1['beta_2_00'][band].data)))

pKExtda=PlotXY()
for band in spireBands:
	pKExtda.addLayer(LayerXY(KExtd2['alpha']['alpha'].data,\
		1.-KExtd2['alpha'][band].data/KExtd1['alpha'][band].data, \
		color=cols[band],name='%s'%(band)))
	pKExtda.setTitleText('Point Source Correction (alpha)')
	pKExtda.yaxis.titleText = '1-%s/%s'%(verKExtd[1]['name'],verKExtd[0]['name'])
	pKExtda.xaxis.titleText = 'Spectral index'
	pKExtda.setSubtitleText('MIN=%.5g MAX=%.5g'%\
		(min(1.-KExtd2['alpha'][band].data/KExtd1['alpha'][band].data),\
		max(1.-KExtd2['alpha'][band].data/KExtd1['alpha'][band].data)))

pKExtdt=PlotXY()
for band in spireBands:
	pKExtdt.addLayer(LayerXY(KExtd2['beta_2_00']['Temperature'].data,\
		1.-KExtd2['beta_2_00'][band].data/KExtd1['beta_2_00'][band].data, \
		color=cols[band],name='%s'%(band)))
	pKExtdt.setTitleText('Point Source Correction (beta=2.0)')
	pKExtdt.yaxis.titleText = '1-%s/%s'%(verKExtd[1]['name'],verKExtd[0]['name'])
	pKExtdt.xaxis.titleText = 'Temperature (K)'
	pKExtdt.setSubtitleText('MIN=%.5g MAX=%.5g'%\
		(min(1.-KExtd2['beta_2_00'][band].data/KExtd1['beta_2_00'][band].data),\
		max(1.-KExtd2['beta_2_00'][band].data/KExtd1['beta_2_00'][band].data)))

pKHFIr=PlotXY()
for band in spireBands:
	pKHFIr.addLayer(LayerXY(KHfi2i['colorCorr']['Temperature'].data,\
		1.-KHfi2i['colorCorr']['ratio545_857'].data/KHfi1i['colorCorr']['ratio545_857'].data, \
		color=cols[band],name='%s'%(band)))
	pKHFIr.setTitleText('SPIRE-HFI correction (545/857 ratio)')
	pKHFIr.yaxis.titleText = '1-%s/%s'%(verKHfi[1]['name'],verKHfi[0]['name'])
	pKHFIr.xaxis.titleText = 'Temperature (K)'
	pKHFIr.setSubtitleText('MIN=%.5g MAX=%.5g'%\
		(min(1.-KHfi2['colorCorr']['ratio545_857'].data/KHfi1['colorCorr']['ratio545_857'].data),\
		max(1.-KHfi2['colorCorr']['ratio545_857'].data/KHfi1['colorCorr']['ratio545_857'].data)))

pKHFIpsw=PlotXY()
for band in spireBands:
	pKHFIpsw.addLayer(LayerXY(KHfi2i['colorCorr']['Temperature'].data,\
		1.-KHfi2i['colorCorr']['k857toPSW'].data/KHfi1i['colorCorr']['k857toPSW'].data, \
		color=cols[band],name='%s'%(band)))
	pKHFIpsw.setTitleText('SPIRE-HFI correction (857 to PSW)')
	pKHFIpsw.yaxis.titleText = '1-%s/%s'%(verKHfi[1]['name'],verKHfi[0]['name'])
	pKHFIpsw.xaxis.titleText = 'Temperature (K)'
	pKHFIpsw.setSubtitleText('MIN=%.5g MAX=%.5g'%\
		(min(1.-KHfi2['colorCorr']['k857toPSW'].data/KHfi1['colorCorr']['k857toPSW'].data),\
		max(1.-KHfi2['colorCorr']['k857toPSW'].data/KHfi1['colorCorr']['k857toPSW'].data)))

pKHFIpmw=PlotXY()
for band in spireBands:
	pKHFIpmw.addLayer(LayerXY(KHfi2i['colorCorr']['Temperature'].data,\
		1.-KHfi2i['colorCorr']['k857toPMW'].data/KHfi1i['colorCorr']['k857toPMW'].data, \
		color=cols[band],name='%s'%(band)))
	pKHFIpmw.setTitleText('SPIRE-HFI correction (857 to PMW)')
	pKHFIpmw.yaxis.titleText = '1-%s/%s'%(verKHfi[1]['name'],verKHfi[0]['name'])
	pKHFIpmw.xaxis.titleText = 'Temperature (K)'
	pKHFIpmw.setSubtitleText('MIN=%.5g MAX=%.5g'%\
		(min(1.-KHfi2['colorCorr']['k857toPMW'].data/KHfi1['colorCorr']['k857toPMW'].data),\
		max(1.-KHfi2['colorCorr']['k857toPMW'].data/KHfi1['colorCorr']['k857toPMW'].data)))

pKHFIplw=PlotXY()
for band in spireBands:
	pKHFIplw.addLayer(LayerXY(KHfi2i['colorCorr']['Temperature'].data,\
		1.-KHfi2i['colorCorr']['k545toPLW'].data/KHfi1i['colorCorr']['k545toPLW'].data, \
		color=cols[band],name='%s'%(band)))
	pKHFIplw.setTitleText('SPIRE-HFI correction (857 to PMW)')
	pKHFIplw.yaxis.titleText = '1-%s/%s'%(verKHfi[1]['name'],verKHfi[0]['name'])
	pKHFIplw.xaxis.titleText = 'Temperature (K)'
	pKHFIplw.setSubtitleText('MIN=%.5g MAX=%.5g'%\
		(min(1.-KHfi2['colorCorr']['k545toPLW'].data/KHfi1['colorCorr']['k545toPLW'].data),\
		max(1.-KHfi2['colorCorr']['k545toPLW'].data/KHfi1['colorCorr']['k545toPLW'].data)))