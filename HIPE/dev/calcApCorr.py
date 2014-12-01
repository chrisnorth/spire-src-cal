##generate various versions of the aperture corrections
from herschel.spire.ia.scripts.useful.sourceCalBundle import SemiExtendedBundle as semi
from herschel.spire.ia.scripts.useful.sourceCalBundle import SpireHandbookBundle as hb
from herschel.spire.ia.scripts.useful.sourceCalBundle import sources
import os

dataDir=Configuration.getProperty('var.hcss.workdir')

### Set alpha range

#alphaArr=Float1d([-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0])
alphaArr=[-4.0,-1.0,2.0,2.5]
alphaArr=Double1d(alphaArr)
nA=len(alphaArr)

alphaKeys=[]
key2alpha={}
for alpha in alphaArr:
    key='alpha_%g'%alpha
    alphaKeys.append(key)
    key2alpha[key]=alpha

bands=['PSW','PMW','PLW']

apPhotRad={"PSW":22.,"PMW":30.,"PLW":45.}
apPhotBGRad={'in':60.,'out':90.}

################################################################################
####                           Make Beam Profiles                           ####
################################################################################
cal=spireCal(jarFile=os.path.join(dataDir,'spire_cal_13_0_photTest2.jar'))
calPhot=hb.getCal(cal)
beamRad=calPhot.getProduct('RadialCorrBeam').getCoreCorrectionTable()['radius'].data
nRad=len(beamRad)
for alpha in alphaArr:
    for band in bands:
        hb.calcSpireEffBeam(alpha,array=band,verbose=True)
spireEffBeams=hb.spireEffBeams

beamAreaPip={'PSW':calPhot.getProduct('ColorCorrBeam').meta['beamPipelinePswArc'].value,\
  'PMW':calPhot.getProduct('ColorCorrBeam').meta['beamPipelinePmwArc'].value,\
  'PLW':calPhot.getProduct('ColorCorrBeam').meta['beamPipelinePlwArc'].value}

#### Calculate beam areas
beamAreas={'PSW':Double1d(nA),'PMW':Double1d(nA),'PLW':Double1d(nA)}
beamAreasCal={'PSW':Double1d(nA),'PMW':Double1d(nA),'PLW':Double1d(nA)}
kBeams={'PSW':Double1d(nA),'PMW':Double1d(nA),'PLW':Double1d(nA)}
kBeamsCal={'PSW':Double1d(nA),'PMW':Double1d(nA),'PLW':Double1d(nA)}
for band in bands:
    for a in range(nA):
        alpha=alphaArr[a]
        #beamAreas[band][a] = hb.calcOmegaEff(alpha,array=band)
        kBeams[band][a] = hb.calcKBeam(alpha,array=band)
        beamAreas[band][a] = beamAreaPip[band]/kBeams[band][a]
        kBeamsCal[band][a] = calPhot.getProduct('ColorCorrBeam').getAlphaCorrection(alpha,band)
        beamAreasCal[band][a] = beamAreaPip[band]/kBeamsCal[band][a]
        
#### Plot beam areas and Kbeam factors
cols={'PSW':java.awt.Color.BLUE,'PMW':java.awt.Color.GREEN,'PLW':java.awt.Color.RED,\
 'PSW 2':java.awt.Color.ORANGE,'PMW 2':java.awt.Color.MAGENTA,'PLW 2':java.awt.Color.CYAN}
pb=PlotXY()
for band in bands:
    pb.addLayer(LayerXY(alphaArr,beamAreas[band],color=cols[band]))
    pb.addLayer(LayerXY(alphaArr,beamAreasCal[band],color=cols[band+' 2'],line=Style.DASHED))
pb.setXtitle('alpha')
pb.setYtitle('Beam Areas [arcsec^2]')

pbDiff=PlotXY()
for band in bands:
    pbDiff.addLayer(LayerXY(alphaArr,beamAreas[band]/beamAreasCal[band],color=cols[band]))
pbDiff.setXtitle('alpha')
pbDiff.setYtitle('Beam Areas / Cal Tree 13.0')

pk=PlotXY()
for band in bands:
    pk.addLayer(LayerXY(alphaArr,kBeams[band],color=cols[band]))
    pk.addLayer(LayerXY(alphaArr,kBeamsCal[band],color=cols[band+' 2'],line=Style.DASHED))
pk.setXtitle('alpha')
pk.setYtitle('Beam Correction Factor [K_beam]')

pkDiff=PlotXY()
for band in bands:
    pkDiff.addLayer(LayerXY(alphaArr,kBeams[band]/kBeamsCal[band],color=cols[band]))
pkDiff.setXtitle('alpha')
pkDiff.setYtitle('Beam Correction / Cal Tree 13.0')

################################################################################
####                             Make Beam Maps                             ####
################################################################################
beamMaps1arcsec={}
for s in alphaKeys:
    beamMaps1arcsec[s]={}
for band in bands:
    for a in range(nA):
        alpha=alphaArr[a]
        s=alphaKeys[a]
        print 'making map for %s %s'%(s,band)
        beamProf=sources.SourceProfile(radArr=beamRad,profile=spireEffBeams[s][band],key='%s_%s'%(s,band))
        beamMaps1arcsec[s][band]=beamProf.makeImage()
        beamMaps1arcsec[s][band].setWcs(Wcs(crpix1=nRad,crpix2=nRad,crval1=0.,crval2=0.,\
          cdelt1=1./3600,cdelt2=1./3600.,ctype1="RA---TAN",ctype2='DEC--TAN'))


###Make beam maps with nominal pixels
beamMapsNominal={}
pixCtrs={}
nPixNoms={}
for s in alphaKeys:
    beamMapsNominal[s]={}
nominalPix={'PSW':6.,'PMW':10.,'PLW':14.}
for band in spireEffBeams[s]:
    nPixNom=int(nRad/nominalPix[band])
    if nPixNom/2==nPixNom/2.:
        nPixNom-=1
    pixCtr=nPixNom/2. + 0.5
    nPixNoms[band]=nPixNom
    pixCtrs[band]=pixCtr
    print nPixNom,pixCtr
    wcsNominal=Wcs(crpix1=pixCtr,crpix2=pixCtr,crval1=0.,crval2=0.,\
          cdelt1=nominalPix[band]/3600,cdelt2=nominalPix[band]/3600.,\
          ctype1="RA---TAN",ctype2='DEC--TAN',naxis2=nPixNom,naxis1=nPixNom)
    for s in beamMaps1arcsec:
        print 'regridding %s for %s'%(band,s)
        beamMapsNominal[s][band]=regrid(beamMaps1arcsec[s][band],wcs=wcsNominal)

################################################################################
####                     Calculate Aperture Corrections                     ####
################################################################################

###Make apCorr holders
apCorr={}
apCorr['ReCalNoBG'] = TableDataset(description='SPIRE aperture correction product (recalc calibration tree)')
apCorr['ReCalIncBG'] = TableDataset(description='SPIRE aperture correction product (recalc calibration tree)')
apCorr['CalNoBG'] = TableDataset(description='SPIRE aperture correction product (spire_cal_13.0)')
apCorr['CalIncBG'] = TableDataset(description='SPIRE aperture correction product (spire_cal_13.0)')
apCorr['AnaNoBG'] = TableDataset(description='SPIRE aperture correction product (analytical)')
apCorr['AnaIncBG'] = TableDataset(description='SPIRE aperture correction product (analytical)')
apCorr['NomPixNoBG'] = TableDataset(description='SPIRE aperture correction product (nominal pixels)')
apCorr['NomPixIncBG'] = TableDataset(description='SPIRE aperture correction product (nominal pixels)')
apCorr['ReCalFracNoBG'] = TableDataset(description='SPIRE aperture correction product (recalc calibration tree with fractional on)')
apCorr['ReCalFracIncBG'] = TableDataset(description='SPIRE aperture correction product (recalc calibration tree with fractional on)')
for ap in apCorr:
    try:
        apCorr[ap].addColumn('alpha',Column(alphaArr))
        for band in bands:
            apCorr[ap].addColumn(band,Column(Float1d(nA)))
    except:
        pass

###  Calculate aperture corrections
#### Get original calibration values
for band in bands:
    for a in range(nA):
        alpha=alphaArr[a]
        apCorr['CalNoBG'][band].data[a]=\
          calPhot.getProduct('ColorCorrApertureList').getProduct('noBG').getApertColorCorrection(alpha,band)
        apCorr['CalIncBG'][band].data[a]=\
          calPhot.getProduct('ColorCorrApertureList').getProduct('incBG').getApertColorCorrection(alpha,band)
        ####

#### Recreate calibration tree values
for band in bands:
    for a in range(nA):
        alpha=alphaArr[a]
        s=alphaKeys[a]
        print 'calculating CalTree aperture correction for %s %s'%(s,band)
        ####
        apPhotReCal = annularSkyAperturePhotometry(image=beamMaps1arcsec[s][band], \
          fractional=0, centerX=nRad-1, centerY=nRad-1, \
          radiusArcsec=apPhotRad[band], \
          innerArcsec=apPhotBGRad['in'], outerArcsec=apPhotBGRad['out'])
        apCorr['ReCalIncBG'][band].data[a]=beamAreas[band][a]/apPhotReCal.getTargetTotal()
        apCorr['ReCalNoBG'][band].data[a]=beamAreas[band][a]/apPhotReCal.getTargetPlusSkyTotal()
        ####

#### Recreate calibration tree values (fractional=1)
for band in bands:
    for a in range(nA):
        alpha=alphaArr[a]
        s=alphaKeys[a]
        print 'calculating CalTree aperture correction (frac=1)for %s %s'%(s,band)
        apPhotReCalFrac = annularSkyAperturePhotometry(image=beamMaps1arcsec[s][band], \
          fractional=1, centerX=nRad-1, centerY=nRad-1, \
          radiusArcsec=apPhotRad[band], \
          innerArcsec=apPhotBGRad['in'], outerArcsec=apPhotBGRad['out'])
        apCorr['ReCalFracIncBG'][band].data[a]=beamAreas[band][a]/apPhotReCalFrac.getTargetTotal()
        apCorr['ReCalFracNoBG'][band].data[a]=beamAreas[band][a]/apPhotReCalFrac.getTargetPlusSkyTotal()
        ####

#### Calculate analytically
for band in bands:
    for a in range(nA):
        alpha=alphaArr[a]
        print 'calculating analytical aperture correction (frac=1)for %s %s'%(s,band)
        anaApCorr=hb.calcApCorr(alpha,aperture=apPhotRad[band],\
          annulus=[apPhotBGRad['in'],apPhotBGRad['out']],array=band,verbose=False)
        apCorr['AnaIncBG'][band].data[a]=anaApCorr[0]
        apCorr['AnaNoBG'][band].data[a]=anaApCorr[1]

#### Run on nominal pixel sizes
for band in bands:
    for a in range(nA):
        alpha=alphaArr[a]
        s=alphaKeys[a]
        print 'calculating nominal pixel aperture correction (frac=1)for %s %s'%(s,band)
        ####
        apPhotNomPix = annularSkyAperturePhotometry(image=beamMapsNominal[s][band], \
          fractional=0, centerX=pixCtrs[band], centerY=pixCtrs[band], \
          radiusArcsec=apPhotRad[band], \
          innerArcsec=apPhotBGRad['in'], outerArcsec=apPhotBGRad['out'])
        apCorr['NomPixIncBG'][band].data[a]=beamAreas[band][a]/apPhotNomPix.getTargetTotal()
        apCorr['NomPixNoBG'][band].data[a]=beamAreas[band][a]/apPhotNomPix.getTargetPlusSkyTotal()
        ####

types={'Cal':{'col':java.awt.Color.BLACK,'desc':'CalTree (13.0)'},\
  'ReCal':{'col':java.awt.Color.BLUE,'desc':'Recalculated CalTree'},\
  'ReCalFrac':{'col':java.awt.Color.CYAN,'desc':'Recalculated CalTree (Frac=1)'},\
  'Ana':{'col':java.awt.Color.RED,'desc':'Analytical'},\
  'NomPix':{'col':java.awt.Color.GREEN,'desc':'Nominal Pixels'}}
plots={}
for band in bands:
    plots[band]=PlotXY()
    pnInc=band+'IncBG'
    pnNo=band+'NoBG'
    pdInc=band+'IncBGDiff'
    pdNo=band+'NoBGDiff'
    t0='Cal'
    spDone=False
    for t in types:
        lnInc=LayerXY(alphaArr,apCorr[t+'IncBG'][band].data,color=types[t]['col'])
        lnInc.setName(types[t]['desc'])        
        plots[band].addLayer(lnInc,0,0)
        lnNo=LayerXY(alphaArr,apCorr[t+'NoBG'][band].data,color=types[t]['col'])
        plots[band].addLayer(lnNo,0,1)
        ldInc=LayerXY(alphaArr,apCorr[t+'IncBG'][band].data/apCorr[t0+'IncBG'][band].data,color=types[t]['col'])
        plots[band].addLayer(ldInc,1,0)
        ldNo=LayerXY(alphaArr,apCorr[t+'NoBG'][band].data/apCorr[t0+'NoBG'][band].data,color=types[t]['col'])
        plots[band].addLayer(ldNo,1,1)
        if not spDone:
            lnInc.setXtitle('alpha')
            lnInc.setYtitle('ApCorr Factor (inc BG)')
            lnNo.setXtitle('alpha')
            lnNo.setYtitle('ApCorr Factor (no BG)')
            ldInc.setXtitle('alpha')
            ldInc.setYtitle('Diff (inc BG) wrt %s'%types[t0]['desc'])
            ldNo.setXtitle('alpha')
            ldNo.setYtitle('Diff (no BG) wrt %s'%types[t0]['desc'])
            spDone=True
    plots[band].setTitleText('%s ApCorr Factors'%band)
    plots[band].getLegend().setVisible(True)
    #
    #
    for t in types:
        plots[pnNo].addLayer(LayerXY(alphaArr,apCorr[t+'NoBG'][band].data,color=types[t]['col'],name=types[t]['desc']))
    plots[pnNo].setXtitle('alpha')
    plots[pnNo].setYtitle('ApCorr Factor')
    plots[pnNo].setTitleText('%s ApCorr Factor (no BG)'%band)
    plots[pnNo].getLegend().setVisible(True)
    #
    for t in types:
        plots[pdInc].addLayer(LayerXY(alphaArr,apCorr[t+'IncBG'][band].data/apCorr[t0+'IncBG'][band].data,color=types[t]['col'],name='%s/%s'%(types[t]['desc'],types[t0]['desc'])))
    plots[pdInc].setXtitle('alpha')
    plots[pdInc].setYtitle('ApCorr Factor wrt %s'%(types[t0]['desc']))
    plots[pdInc].setTitleText('%s ApCorr Differences (inc BG)'%band)
    plots[pdInc].getLegend().setVisible(True)
    #
    plots[pdNo]=SubPlot()
    for t in types:
        plots[pdNo].addLayer(LayerXY(alphaArr,apCorr[t+'NoBG'][band].data/apCorr[t0+'NoBG'][band].data,color=types[t]['col'],name='%s/%s'%(types[t]['desc'],types[t0]['desc'])))
    plots[pdNo].setXtitle('alpha')
    plots[pdNo].setYtitle('ApCorr Factor wrt %s'%(types[t0]['desc']))
    plots[pdNo].setTitleText('%s ApCorr Differences (no BG)'%band)
    plots[pdNo].getLegend().setVisible(True)

    
for band in bands:
    pacInc.addLayer(LayerXY(alphaArr,apCorr['ReCalIncBG'][band].data,color=cols[band],name='ReCal '+band))
    pacInc.addLayer(LayerXY(alphaArr,apCorr['ReCalFracIncBG'][band].data,color=cols[band+' 2'],name='ReCalFrac '+band))
    pacInc.addLayer(LayerXY(alphaArr,apCorr['CalIncBG'][band].data,color=cols[band],line=Style.DASHED,name='Cal '+band))
    pacInc.addLayer(LayerXY(alphaArr,apCorr['AnaIncBG'][band].data,color=cols[band+' 2'],line=Style.DASHED,name='Ana '+band))
pacInc.setXtitle('alpha')
pacInc.setYtitle('Aperture Correction Factor (inc BG)')
pacInc.getLegend().setVisible(True)

pacNo=PlotXY()
for band in bands:
    pacNo.addLayer(LayerXY(alphaArr,apCorr['ReCalNoBG'][band].data,color=cols[band],name='ReCal '+band))
    pacNo.addLayer(LayerXY(alphaArr,apCorr['ReCalFracNoBG'][band].data,color=cols[band+' 2'],name='ReCalFrac '+band))
    pacNo.addLayer(LayerXY(alphaArr,apCorr['CalNoBG'][band].data,color=cols[band],line=Style.DASHED,name='Cal '+band))
    pacNo.addLayer(LayerXY(alphaArr,apCorr['AnaNoBG'][band].data,color=cols[band+' 2'],line=Style.DASHED,name='Ana '+band))
pacNo.setXtitle('alpha')
pacNo.setYtitle('Aperture Correction Factor (no BG)')
pacNo.getLegend().setVisible(True)

pacDiff=PlotXY()
for band in bands:
    pacDiff.addLayer(LayerXY(alphaArr,apCorr['CalNoBG'][band].data/apCorr['ReCalNoBG'][band].data,color=cols[band],name='Cal/ReCal No '+band))
    pacDiff.addLayer(LayerXY(alphaArr,apCorr['CalIncBG'][band].data/apCorr['ReCalIncBG'][band].data,color=cols[band+' 2'],name='Cal/ReCal Inc '+band))
    pacDiff.addLayer(LayerXY(alphaArr,apCorr['AnaNoBG'][band].data/apCorr['ReCalNoBG'][band].data,color=cols[band],line=Style.DASHED,name='Ana/ReCal No '+band))
    pacDiff.addLayer(LayerXY(alphaArr,apCorr['AnaIncBG'][band].data/apCorr['ReCalIncBG'][band].data,color=cols[band+' 2'],line=Style.DASHED,name='Ana/ReCal Inc '+band))
pacDiff.setXtitle('alpha')
pacDiff.setYtitle('Aperture Correction Factor / Cal Tree 13.0')
pacDiff.getLegend().setVisible(True)

        
        
