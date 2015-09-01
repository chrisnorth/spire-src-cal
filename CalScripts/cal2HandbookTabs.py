#turns cal tables into tex tables for OM

#cal=spireCal(pool='spire_cal_12_2')
#urlHaio ='http://archives.esac.esa.int/hsaint/aio/jsp/'
#archive = HsaReadPool(urlHaio+'metadata.jsp',urlHaio+'product.jsp')
#hsaRead = ProductStorage()
#hsaRead.register(archive)
#lookup = IdLookup ("spire_cal")
#cal = SpireCal.getInstance (hsaRead, lookup) # This should be version spire_cal_13_1
#cal=spireCal(jarFile='/home/chris/hcss/workspace/spire_cal_13_0_photTest2.jar')

from herschel.ia.numeric.toolbox.util.MoreMath import modulo
from herschel.share.unit import *
import os
c = Constant.SPEED_OF_LIGHT.value

################################################################################
##Read calibration tree
cal = spireCal(jarFile='hcss/workspace/spire_cal_14_0_phot.jar')

#set Output path
outPath='/data/Herschel/Calibration/CalProducts/Tables/'

#set Handbook output file
hbFile = 'spire_cal_14_0_phot_HandbookTables.txt'

################################################################################

hbOut=open(os.path.join(outPath,hbFile),'w')
hbOut.write('Cal version: %s\n'%(cal.meta['version'].value))
hbOut.write('Date: %s\n'%(java.util.Date()))
hbOut.write('Handbook Tables\n')

verOut=open(os.path.join(outPath,'version.txt'),'w')
verOut.write('Cal version: %s\n'%(cal.meta['version'].value))
verOut.write('Date: %s\n'%(java.util.Date()))

#set up bands
spireBands=['PSW','PMW','PLW']

#set up wavelengths & frequencies
wl0_um={'PSW':250.0,'PMW':350.,'PLW':500.}
nu0_GHz={}
for band in spireBands:
    nu0_GHz[band] = 1.e-9 * c / (wl0_um[band] * 1.e-6)

#-------------------------------------------------------------------------------
##get Neptune beam areas
beamNepArc= {'PSW':cal.getPhot().getProduct("RadialCorrBeam").meta['beamNeptunePswArc'].double, \
	'PMW':cal.getPhot().getProduct("RadialCorrBeam").meta['beamNeptunePmwArc'].double, \
	'PLW':cal.getPhot().getProduct("RadialCorrBeam").meta['beamNeptunePlwArc'].double}
##get Neptune alpha
alphaNep= {'PSW':cal.getPhot().getProduct("RadialCorrBeam").meta['alphaNeptunePsw'].double, \
	'PMW':cal.getPhot().getProduct("RadialCorrBeam").meta['alphaNeptunePmw'].double, \
	'PLW':cal.getPhot().getProduct("RadialCorrBeam").meta['alphaNeptunePlw'].double}
##get effective frequencies
freqEff = {'PSW':cal.getPhot().getProduct("RadialCorrBeam").meta['freqEffPsw'].double, \
	'PMW':cal.getPhot().getProduct("RadialCorrBeam").meta['freqEffPmw'].double, \
	'PLW':cal.getPhot().getProduct("RadialCorrBeam").meta['freqEffPlw'].double}

##get beam properties
beamMajor = {'PSW':cal.getPhot().getProduct('BeamProfList').getProduct('PSW','fine').meta['FWHM_majorPsw'].double,\
	'PMW':cal.getPhot().getProduct('BeamProfList').getProduct('PMW','fine').meta['FWHM_majorPmw'].double, \
	'PLW':cal.getPhot().getProduct('BeamProfList').getProduct('PLW','fine').meta['FWHM_majorPlw'].double}
beamMinor = {'PSW':cal.getPhot().getProduct('BeamProfList').getProduct('PSW','fine').meta['FWHM_minorPsw'].double,\
	'PMW':cal.getPhot().getProduct('BeamProfList').getProduct('PMW','fine').meta['FWHM_minorPmw'].double,\
	'PLW':cal.getPhot().getProduct('BeamProfList').getProduct('PLW','fine').meta['FWHM_minorPlw'].double}
beamGMean = {'PSW':cal.getPhot().getProduct('BeamProfList').getProduct('PSW','fine').meta['FWHM_gMeanPsw'].double,\
	'PMW':cal.getPhot().getProduct('BeamProfList').getProduct('PMW','fine').meta['FWHM_gMeanPmw'].double,\
	'PLW':cal.getPhot().getProduct('BeamProfList').getProduct('PLW','fine').meta['FWHM_gMeanPlw'].double}
beamFlat = {'PSW':cal.getPhot().getProduct('BeamProfList').getProduct('PSW','fine').meta['FlatteningPsw'].double,\
	'PMW':cal.getPhot().getProduct('BeamProfList').getProduct('PMW','fine').meta['FlatteningPmw'].double,\
	'PLW':cal.getPhot().getProduct('BeamProfList').getProduct('PLW','fine').meta['FlatteningPlw'].double}
#-------------------------------------------------------------------------------
##get beam correction factors
kBeam = cal.getPhot().getProduct("ColorCorrBeam")

##get pipeline beam areas
pipBeamArc = {'PSW':cal.getPhot().getProduct("RadialCorrBeam").meta['beamPipelinePswArc'].double, \
	'PMW':cal.getPhot().getProduct("RadialCorrBeam").meta['beamPipelinePmwArc'].double, \
	'PLW':cal.getPhot().getProduct("RadialCorrBeam").meta['beamPipelinePlwArc'].double}

pipBeamSr = {'PSW':cal.getPhot().getProduct("RadialCorrBeam").meta['beamPipelinePswSr'].double, \
	'PMW':cal.getPhot().getProduct("RadialCorrBeam").meta['beamPipelinePmwSr'].double, \
	'PLW':cal.getPhot().getProduct("RadialCorrBeam").meta['beamPipelinePlwSr'].double}

#calculate effective beam areas
effBeamArc = kBeam.copy()
for n in kBeam.getSets():
	#print n
	for band in spireBands:
		#print band
		effBeamArc[n][band].data = pipBeamArc[band]/kBeam[n][band].data

#-------------------------------------------------------------------------------
##get kColP, kColE tables
kPsrc=cal.getPhot().getProduct("ColorCorrKList")[1]
kExtd=cal.getPhot().getProduct("ColorCorrKList")[0]

#-------------------------------------------------------------------------------
##get k4 parameters
k4P={'PSW':cal.getPhot().getProduct("FluxConvList")[0].meta['k4P_PSW'].double,\
	'PMW':cal.getPhot().getProduct("FluxConvList")[0].meta['k4P_PMW'].double, \
	'PLW':cal.getPhot().getProduct("FluxConvList")[0].meta['k4P_PLW'].double}
k4E={'PSW':cal.getPhot().getProduct("FluxConvList")[0].meta['k4E_PSW'].double,\
	'PMW':cal.getPhot().getProduct("FluxConvList")[0].meta['k4E_PMW'].double, \
	'PLW':cal.getPhot().getProduct("FluxConvList")[0].meta['k4E_PLW'].double}

kMonE={}
kPtoE={}
k4E4P={}
for band in spireBands:
	kMonE[band] = (k4E[band] / pipBeamSr[band]) / 1.e6
	kPtoE[band] = (kMonE[band] / k4P[band])
	k4E4P[band] = k4E[band] / k4P[band]

#-------------------------------------------------------------------------------
##get aperture corrections
apCorr_incBG = cal.getPhot().getProduct('ColorCorrApertureList')[0]
apCorr_noBG = cal.getPhot().getProduct('ColorCorrApertureList')[1]

################################################################################
## Print relevant DRG tables
################################################################################

verOut.write('\n---DRG Tables---\n')

#Table 6.9. SPIRE pipeline conversion factors for point and extended sources
drg6_9='tab_6-9_parameters.csv'
tabOut = open(os.path.join(outPath,drg6_9),'w')
tabOut.write('Parameter,PSW,PMW,PLW\n')
tabOut.write('K4P , %.4f , %.4f , %.4f \n'%\
	(k4P['PSW'],k4P['PMW'],k4P['PLW']))
tabOut.write('KMonE , %.4f , %.4f , %.4f \n'%\
	(kMonE['PSW'],kMonE['PMW'],kMonE['PLW']))
tabOut.write('KPtoE , %.4f , %.4f , %.4f \n'%\
	(kPtoE['PSW'],kPtoE['PMW'],kPtoE['PLW']))
tabOut.write('Beam area (arcsec^2; alpha=-1) , %.4f , %.4f , %.4f \n'%\
	(pipBeamArc['PSW'],pipBeamArc['PMW'],pipBeamArc['PLW']))
tabOut.write('K4E , %.4f , %.4f , %.4f \n'%\
	(k4E['PSW'],k4E['PMW'],k4E['PLW']))
tabOut.write('K4E/K4P (K4EdivK4P) , %.4f , %.4f , %.4f \n'%\
	(k4E4P['PSW'],k4E4P['PMW'],k4E4P['PLW']))
tabOut.close()
print('DRG Table 6.9 written to %s'%(os.path.join(outPath,drg6_9)))
verOut.write('Table 6.9. SPIRE pipeline conversion factors for point and extended sources:\n  %s\n\n'%\
	os.path.join(outPath,drg6_9))

#Table 6.10. Beam Areas assumed by the pipeline (alpha = -1)
drg6_10='tab_6-10_pipeline-beam-areas.csv'
tabOut = open(os.path.join(outPath,drg6_10),'w')
tabOut.write('Spectral Index,Beam Area (arcsec^2),,,Beam Area (/10^-8 sr),,\n')
tabOut.write('F_nu=nu^a,PSW,PMW,PLW,PSW,PMW,PLW\n')
tabOut.write('-1, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f\n'%\
	(pipBeamArc['PSW'],pipBeamArc['PMW'],pipBeamArc['PLW'],pipBeamSr['PSW'],pipBeamSr['PMW'],pipBeamSr['PLW']))
tabOut.close()
print('DRG Table 6.10 written to %s'%(os.path.join(outPath,drg6_10)))
verOut.write('Table 6.10. Beam Areas assumed by the pipeline (alpha = -1):\n  %s\n\n'%\
	os.path.join(outPath,drg6_10))
	
#Table 6.11. Effective Beam Area ratios (beam correction) as function of spectral index (alpha)
drg6_11='tab_6-11_beam-correction.csv'
tabOut = open(os.path.join(outPath,drg6_11),'w')
tabOut.write('Spectral Index,Effective Beam Area Ratio (Omega(-1)/Omega(a)),,\n')
tabOut.write('F_nu=nu^a,PSW,PMW,PLW\n')
kB=kBeam['alpha']
for a in range(len(kB['alpha'].data)):
    tabOut.write('%.1f , %.4f , %.4f , %.4f\n'%\
        (kB['alpha'].data[a],kB['PSW'].data[a],kB['PMW'].data[a],kB['PLW'].data[a]))
tabOut.close()
print('DRG Table 6.11 written to %s'%(os.path.join(outPath,drg6_11)))
verOut.write('Table 6.11. Effective Beam Area ratios (beam correction) as function of spectral index (alpha):\n  %s\n\n'%\
	os.path.join(outPath,drg6_11))
	
#Table 6.12. SPIRE FWHM Parameters for 1 arcsec pixels
drg6_12='tab_6-12_fwhm-parameters.csv'
tabOut = open(os.path.join(outPath,drg6_12),'w')
tabOut.write('Band,FWHM,MeanFWHM,Ellipticity\n')
tabOut.write('(micron),(arcsec),(arcsec),(Ellipticity)\n')
for band in spireBands:
    tabOut.write('%d , %.1fx%.1f , %.1f , %.1f\n'%\
    	(wl0_um[band],beamMajor[band],beamMinor[band],beamGMean[band],beamFlat[band]*100))
tabOut.close()
print('DRG Table 6.12 written to %s'%(os.path.join(outPath,drg6_12)))
verOut.write('Table 6.12. SPIRE FWHM Parameters for 1 arcsec pixels:\n  %s\n\n'%\
	os.path.join(outPath,drg6_12))
	
#Table 6.15. Aperture Corrections for Annular Aperture Photometry
drg6_15='tab_6-15_aperture-corrections.csv'
tabOut = open(os.path.join(outPath,drg6_15),'w')
tabOut.write('Spectral Index,Background included,,,Case of background removed,,\n')
tabOut.write('F_nu=nu^a,PSW,PMW,PLW,PSW,PMW,PLW\n')
aI=apCorr_incBG['alpha']
aN=apCorr_noBG['alpha']
for a in range(len(aI['alpha'].data)):
    tabOut.write('%.1f , %.4f , %.4f , %.4f , %.4f , %.4f , %.4f\n'%\
    	(aI['alpha'].data[a],\
    	aI['PSW'].data[a],aI['PMW'].data[a],aI['PLW'].data[a], \
    	aN['PSW'].data[a],aN['PMW'].data[a],aN['PLW'].data[a]))
tabOut.close()
print('DRG Table 6.15 written to %s'%(os.path.join(outPath,drg6_15)))
verOut.write('Table 6.15. Aperture Corrections for Annular Aperture Photometry:\n  %s\n\n'%\
	os.path.join(outPath,drg6_15))
	
#Table 6.16. Color Corrections for Point Sources and Extended Emission
drg6_16='tab_6-16_colour-corrections.csv'
tabOut = open(os.path.join(outPath,drg6_16),'w')
tabOut.write('Spectral Index,Point Source Colour Correction,,,Extended Emission Colour Correction,,\n')
tabOut.write('F_nu=nu^a,PSW,PMW,PLW,PSW,PMW,PLW\n')
kP=kPsrc['alpha']
kE=kExtd['alpha']
for a in range(len(kP['alpha'].data)):
    tabOut.write('%.1f , %.4f , %.4f , %.4f , %.4f , %.4f , %.4f \n'%\
    	(kP['alpha'].data[a],\
    	kP['PSW'].data[a],kP['PMW'].data[a],kP['PLW'].data[a],\
    	kE['PSW'].data[a],kE['PMW'].data[a],kE['PLW'].data[a]))
tabOut.close()
print('DRG Table 6.16 written to %s'%(os.path.join(outPath,drg6_16)))
verOut.write('Table 6.16. Color Corrections for Point Sources and Extended Emission:\n  %s\n\n'%\
	os.path.join(outPath,drg6_16))
	
################################################################################
## Print relevant Handbook tables
################################################################################

verOut.write('\n---Handbook Tables---\n')

#-------------------------------------------------------------------------------
verOut.write('Table 5.1 (K4 parameters etc.)\n')
hbOut.write('\n----- Table 5.1 (K4 parameters etc.):\n')
hbOut.write('\\begin{tabular}{l|ccc}\n')
hbOut.write('\\hline\\hline\n')
hbOut.write('Band & PSW & PMW & PLW \\\\\n')
hbOut.write('\\hline\n')
hbOut.write('Reference wavelength, $\\lambda_0$ ($\mu$m) & %.1f & %.1f & %.1f \\\\\n'%\
	(wl0_um['PSW'],wl0_um['PMW'],wl0_um['PLW']))
hbOut.write('Reference frequency, $\\nu_0$ (GHz) & %.2f & %.2f & %.2f \\\\\n'%\
	(nu0_GHz['PSW'],nu0_GHz['PMW'],nu0_GHz['PLW']))
hbOut.write('$K_\\mathrm{4P}$ & %.4f & %.4f & %.4f \\\\\n'%\
	(k4P['PSW'],k4P['PMW'],k4P['PLW']))
hbOut.write('$K_\\mathrm{MonE}$ (MJy/sr per Jy/beam) & %.3f & %.3f & %.3f \\\\\n'%\
	(kMonE['PSW'],kMonE['PMW'],kMonE['PLW']))
hbOut.write('$K_\\mathrm{PtoE}$ (MJy/sr per Jy/beam) & %.3f & %.3f & %.3f \\\\\n'%\
	(kPtoE['PSW'],kPtoE['PMW'],kPtoE['PLW']))
hbOut.write('$\\Omega_\\mathrm{pip}$ (arcsec$^2$) & %.2f & %.2f & %.2f \\\\\n'%\
	(pipBeamArc['PSW'],pipBeamArc['PMW'],pipBeamArc['PLW']))
hbOut.write('\\hline\n')
#hbOut.write('$K_\\mathrm{4E}$ & %.4f & %.4f & %.4f \\\\\n'%(k4E['PSW'],k4E['PMW'],k4E['PLW'])
#hbOut.write('$K_\\mathrm{4E}/K_\\mathrm{4P}$ & %.4f & %.4f & %.4f '%(k4E4P['PSW'],k4E4P['PMW'],k4E4P['PLW'])
hbOut.write('\\end{tabular}\n')
print('Handbook Table 5.1 appended to file')

#-------------------------------------------------------------------------------
verOut.write('Table 5.2 (Basic 2-D Gaussian parameters)\n')
hbOut.write('\n----- Table 5.2 (Basic 2-D Gaussian parameters):\n')
hbOut.write('\\begin{tabular}{l|ccc}\n')
hbOut.write('\\hline\\hline\n')
hbOut.write('Band & PSW & PMW & PLW\\\\\n')
hbOut.write('\hline\n')
hbOut.write('$\\alpha_\\mathrm{Nep}$ & %.2f & %.2f & %.2f \\\\\n'%\
	(alphaNep['PSW'],alphaNep['PMW'],alphaNep['PLW']))
hbOut.write('Major$\\times$Minor\\- FWHM (arcsec) & %.1f$\\times$%.1f & %.1f$\\times$%.1f & %.1f$\\times$%.1f \\\\\n'%\
	(beamMinor['PSW'],beamMinor['PSW'],beamMajor['PMW'],beamMinor['PMW'],beamMajor['PLW'],beamMinor['PLW']))
hbOut.write('Geometric mean FWHM ($\\theta_\\mathrm{Nep}$, arcsec) & %.1f & %.1f & %.1f \\\\\n'%\
	(beamGMean['PSW'],beamGMean['PMW'],beamGMean['PLW']))
hbOut.write('Flattening (\\%%)& %.1f & %.1f & %.1f \\\\\n'%\
	(beamFlat['PSW']*100,beamFlat['PMW']*100,beamFlat['PLW']*100))
#hbOut.write('Measured beam solid angle ($\\Omega_\\mathrm{Nep}$, arcsec$^2$) & 450 & 795 & 1665 \\\\\n')
hbOut.write('Measured beam solid angle ($\\Omega_\\mathrm{Nep}$, arcsec$^2$) & %.0f & %.0f & %.0f \\\\\n'%\
	(beamNepArc['PSW'],beamNepArc['PMW'],beamNepArc['PLW']))
hbOut.write('Isophotal frequency$^a$ ($\\nu_\\mathrm{eff}$, GHz) & %.2f & %.2f & %.2f\n'%\
	(freqEff['PSW'],freqEff['PMW'],freqEff['PLW']))
hbOut.write('\\end{tabular}\n')
print('Handbook Table 5.2 appended to file')

#-------------------------------------------------------------------------------
verOut.write('Table 5.3 (Effective beam solid angle, alpha)\n')
hbOut.write('\n----- Table 5.3 (Effective beam solid angle, alpha):\n')
hbOut.write('\\begin{tabular}{r|ccc|ccc}\n')
hbOut.write('\\hline\\hline\n')
hbOut.write('& \\multicolumn{3}{c|}{Effective beam solid angle} & \\multicolumn{3}{c}{Beam correction factor} \\\\\n')
hbOut.write('& \\multicolumn{3}{c|}{$\\Omega_\\mathrm{eff}$ (arcsec$^2$)} & \\multicolumn{3}{c}{$K_\\mathrm{Beam}$} \\\\\n')
hbOut.write('$\\alpha_S$  & PSW & PMW & PLW & PSW & PMW & PLW \\\\\n')
hbOut.write('\\hline\n')
kB=kBeam['alpha']
eB=effBeamArc['alpha']
for a in range(len(eB['alpha'].data)):
	if eB['alpha'].data[a]==-1.:
		hbOut.write('{\\bf %.1f }& {\\bf %.2f} & {\\bf %.2f} & {\\bf %.2f} & {\\bf %.4f} & {\\bf %.4f} & {\\bf %.4f} \\\\\n'%\
			(eB['alpha'].data[a],\
			eB['PSW'].data[a],eB['PMW'].data[a],eB['PLW'].data[a],\
			kB['PSW'].data[a],kB['PMW'].data[a],kB['PLW'].data[a]))
	else:
		hbOut.write('%.1f & %.2f & %.2f & %.2f & %.4f & %.4f & %.4f \\\\\n'%\
			(eB['alpha'].data[a],\
			eB['PSW'].data[a],eB['PMW'].data[a],eB['PLW'].data[a],\
			kB['PSW'].data[a],kB['PMW'].data[a],kB['PLW'].data[a]))
hbOut.write('\\hline\n')
hbOut.write('\\end{tabular}\n')
print('Handbook Table 5.3 appended to file')

#-------------------------------------------------------------------------------
verOut.write('Table 5.4 (Effective beam solid angle, beta=1.5,2.0)\n')
hbOut.write('\n----- Table 5.4 (Effective beam solid angle, beta=1.5,2.0):\n')
hbOut.write('\\begin{tabular}{r|ccc|ccc}\n')
hbOut.write('\\hline\\hline\n')
hbOut.write('& \\multicolumn{3}{c|}{Effective beam solid angle} & \\multicolumn{3}{c}{Beam correction factor} \\\\\n')
hbOut.write('& \\multicolumn{3}{c|}{$\\Omega_\\mathrm{eff}(T,\\beta)$ (arcsec$^2$)} & \\multicolumn{3}{c}{$K_\\mathrm{Beam}(T,\\beta)$} \\\\\n')
hbOut.write('Temp (K) & PSW & PMW & PLW & PSW & PMW & PLW \\\\\n')
kB=kBeam['beta_1_50']
eB=effBeamArc['beta_1_50']
hbOut.write('\\hline\n')
hbOut.write(' & \multicolumn{6}{c}{$\\beta=1.5$}\\\\\n')
hbOut.write('\\hline\n')
for t in range(len(eB['Temperature'].data)):
	if (eB['Temperature'].data[t] < 10) or (eB['Temperature'].data[t] <= 40 and modulo(eB['Temperature'].data[t],5)==0.):
		hbOut.write('%.1f & %.2f & %.2f & %.2f & %.4f & %.4f & %.4f \\\\\n'%\
			(eB['Temperature'].data[t],\
			eB['PSW'].data[t],eB['PMW'].data[t],eB['PLW'].data[t],\
			kB['PSW'].data[t],kB['PMW'].data[t],kB['PLW'].data[t]))
kB=kBeam['beta_2_00']
eB=effBeamArc['beta_2_00']
hbOut.write('\\hline\n')
hbOut.write(' & \multicolumn{6}{c}{$\\beta=1.5$}\\\\\n')
hbOut.write('\\hline\n')
for t in range(len(eB['Temperature'].data)):
	if (eB['Temperature'].data[t] < 10) or (eB['Temperature'].data[t] <= 40 and modulo(eB['Temperature'].data[t],5)==0.):
		hbOut.write('%.1f & %.2f & %.2f & %.2f & %.4f & %.4f & %.4f \\\\\n'%\
			(eB['Temperature'].data[t],\
			eB['PSW'].data[t],eB['PMW'].data[t],eB['PLW'].data[t],\
			kB['PSW'].data[t],kB['PMW'].data[t],kB['PLW'].data[t]))
hbOut.write('\\hline\n')
hbOut.write('\\end{tabular}\n')
print('Handbook Table 5.4 appended to file')

#-------------------------------------------------------------------------------
verOut.write('Table 5.5 (KColP, KColE with alpha)\n')
hbOut.write('\n----- Table 5.5 (KColP, KColE with alpha):\n')
hbOut.write('\\begin{tabular}{c|lll|lll|}\n')
hbOut.write('\\hline\\hline\n')
hbOut.write('& \\multicolumn{3}{c|}{Point Source ($K_\\mathrm{ColP}$)} & \\multicolumn{3}{c|}{Extended source ($K_\\mathrm{ColE}$)} \\\\\n')
hbOut.write('$\\alpha_S$  & PSW & PMW & PLW & PSW & PMW & PLW \\\\\n')
hbOut.write('\\hline\n')
kP=kPsrc['alpha']
kE=kExtd['alpha']
for a in range(len(kP['alpha'].data)):
	if kP['alpha'].data[a]==2.:
		#bold print
		hbOut.write('{\\bf %.1f} & {\\bf %.4f} & {\\bf %.4f} & {\\bf %.4f} & {\\bf %.4f} & {\\bf %.4f} & {\\bf %.4f} \\\\\n'%\
			(kP['alpha'].data[a],\
			kP['PSW'].data[a],kP['PMW'].data[a],kP['PLW'].data[a],\
			kE['PSW'].data[a],kE['PMW'].data[a],kE['PLW'].data[a]))
	else:
		hbOut.write('%.1f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\\n'%\
			(kP['alpha'].data[a],\
			kP['PSW'].data[a],kP['PMW'].data[a],kP['PLW'].data[a],\
			kE['PSW'].data[a],kE['PMW'].data[a],kE['PLW'].data[a]))
hbOut.write('\\end{tabular}\n')
print('Handbook Table 5.5 appended to file')

#-------------------------------------------------------------------------------
verOut.write('Table 5.6 (KColP, kColE with beta=1.5,2.0)\n')
hbOut.write('\n----- Table 5.6 (KColP, kColE with beta=1.5,2.0):\n')
hbOut.write('\\begin{tabular}{r|lll|lll}\n')
hbOut.write('\\hline\\hline\n')
hbOut.write('& \\multicolumn{3}{c|}{Point source ($K_\\mathrm{ColP}$)} &')
hbOut.write('\\multicolumn{3}{c}{Extended source ($K_\\mathrm{ColE}$)} \\\\\n')
hbOut.write('$T$ & PSW & PMW & PLW & PSW & PMW & PLW \\\\\n')
kP=kPsrc['beta_1_50']
kE=kExtd['beta_1_50']
hbOut.write('\\hline\n')
hbOut.write(' & \multicolumn{6}{c}{$\\beta=1.5$}\\\\\n')
hbOut.write('\\hline\n')
for t in range(len(eB['Temperature'].data)):
	if (kP['Temperature'].data[t] < 10) or (kP['Temperature'].data[t] <= 40 and modulo(kP['Temperature'].data[t],5)==0.):
		hbOut.write('%.1f & %.2f & %.2f & %.2f & %.4f & %.4f & %.4f \\\\\n'%\
			(kP['Temperature'].data[t],\
			kP['PSW'].data[t],kP['PMW'].data[t],kP['PLW'].data[t],\
			kE['PSW'].data[t],kE['PMW'].data[t],kE['PLW'].data[t]))
kP=kPsrc['beta_2_00']
kE=kExtd['beta_2_00']
hbOut.write('\\hline\n')
hbOut.write(' & \multicolumn{6}{c}{$\\beta=1.5$}\\\\\n')
hbOut.write('\\hline\n')
for t in range(len(eB['Temperature'].data)):
	if (kP['Temperature'].data[t] < 10) or (kP['Temperature'].data[t] <= 40 and modulo(kP['Temperature'].data[t],5)==0.):
		hbOut.write('%.1f & %.2f & %.2f & %.2f & %.4f & %.4f & %.4f \\\\\n'%\
			(kP['Temperature'].data[t],\
			kP['PSW'].data[t],kP['PMW'].data[t],kP['PLW'].data[t],\
			kE['PSW'].data[t],kE['PMW'].data[t],kE['PLW'].data[t]))
hbOut.write('\\hline\n')
hbOut.write('\\end{tabular}\n')
print('Handbook Table 5.6 appended to file')

#-------------------------------------------------------------------------------
verOut.write('Table 5.7 (Aperture correction)\n')
hbOut.write('\n----- Table 5.7 (Aperture correction):\n')
hbOut.write('\\begin{tabular}{r|ccc|ccc}\n')
hbOut.write('\\hline\\hline\n')
hbOut.write('& \\multicolumn{3}{c|}{Background included} & \\multicolumn{3}{c}{No Background} \\\\\n')
hbOut.write('$\\alpha$ & PSW & PMW & PLW & PSW & PMW & PLW \\\\\n')
hbOut.write('\\hline\n')
aI=apCorr_incBG['alpha']
aN=apCorr_noBG['alpha']
for a in range(len(aI['alpha'].data)):
	hbOut.write('%.1f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\\n'%\
	(aI['alpha'].data[a],\
	aI['PSW'].data[a],aI['PMW'].data[a],aI['PLW'].data[a], \
	aN['PSW'].data[a],aN['PMW'].data[a],aN['PLW'].data[a]))
hbOut.write('\\end{tabular}\n')
print('Handbook Table 5.7 appended to file')

hbOut.close()

print('Handbook tables written to %s'%(os.path.join(outPath,hbFile)))

verOut.write('  %s\n'%(os.path.join(outPath,hbFile)))
verOut.close()
