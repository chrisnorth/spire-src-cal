#
#  This file is part of Herschel Common Science System (HCSS).
#  Copyright 2001-2013 Herschel Science Ground Segment Consortium
#
#  HCSS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
#
#  HCSS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General
#  Public License along with HCSS.
#  If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================
# 
#  Herschel-SPIRE Radial Beam Profile 
# 
#  This routine reads in the Herschel-SPIRE beam profiles and computes the 
#  metadata for the calibration product. It calculates the cumulative beam
#  area as a function of radius, normalised to the final value. Using the
#  measured beam areas and the spectral index of Neptune it calculates effective
#  frequency in each band. This is used to calculate the effective beam area for
#  the pipeline spectrum (alpha=-1).
# 
#  Input:
#    RSRF from SPIRE calibration tree
#    Beam profiles (core and constant) from ascii data files
#    Neptune spectral index and measured area
#    Exponent for scaling of beam FWHM (gamma)
#    Name and version of output file
# 
#  Output:
#    SCalPhotRadialCorrBeam product
# 
#  Calculations:
#   1. Reads in constant and core from ascii data files
#   2. Computes full beam profile as a function of radius, r
#        Beam_full(r) = max(Beam_core(r), Beam_constant(r))
#   3. Calculates the effective frequency, which is the frequency at which the
#      monochromatic beam area equals the broadband beam area measured on
#      Neptune.
#   4. Scales the monochromatic beam with frequency by exponent gamma and
#      calculates the monochromatic beam area
#        Beam_mono(r,nu) = Beam_full(r*(nu/nu_eff)^-gamma)
#           (meaning that the FWHM varies as (nu/nu_eff)^gamma
#        Omega_mono(nu) = int{Beam_mono(r,nu) * 2*pi*r dr}
#   5. By integrating over the RSRF and source spectrum, it calculates the
#      effective beam area for the pipeline spectrum
#        Omega_eff(alpha) = int{Omega_mono(nu) * rsrf(nu) * nu^alpha dnu} / 
#                                 int{rsrf(nu) * nu^alpha dnu}
#  
#   Note that the effective frequency is defined such that
#      Omega_eff(alpha_Nep) = Omega_Nep
#
#  There are some functions defined (common between calibration products)
#   * spireMonoBeam: Calculate monochromatic beam profile & area at given freq
#   * spireMonoAreas: Calculate monochromatic beam areas at range of frequencies
#   * spireEffArea: Calculate the effective beam area for a given spectrum
#   * spireFindEffFreq: Calculate the effective frequency for SPIRE
#
#===============================================================================
# $Id: makeSCalPhotRadialCorrBeam.py,v 1.7 2014/11/14 16:33:28 epoleham Exp $
#
# 
# Edition History
# E. Polehampton   22-10-2013  - First version adapted from Andreas' script - SPCAL-83
# E. Polehampton   01-11-2013  - Add effective frequencies to metadata
# E. Polehampton   21-01-2014  - Add gamma as metadata, and beam areas - SPCAL-95, 96
# E. Polehampton   22-01-2014  - update tables - SPCAL-90
# Chris North      18-02-2014  - 1.0: First version (full calculation)
#                                     Used for spire_cal_12_1
# Chris North      26-02-2014  - 1.1: Added writing of ascii tables and log file
#                                   SPCAL-109  
# Chris North      04-11-2014  - 2.0: SPCAL-126 Used new beams (no constant section)
#                                     Changed Neptune beam areas to B. Schulz values
#                                     "constant" table filled with zeros
#                                     removed constant beam from calculations where appropriate
# E. Polehampton   13-11-2014  - Tidy up metadata (use metadata dictionary)
#===============================================================================

import os
from herschel.share.unit import *
scriptVersionString = "makeSCalPhotRadialCorrBeam.py $Revision: 1.7 $"
metaDict = herschel.spire.ia.util.MetaDataDictionary.getInstance()

#-------------------------------------------------------------------------------
#===============================================================================
#=====                              DIRECTORIES                            =====
#===============================================================================
#-------------------------------------------------------------------------------

#directory = "..//..//..//..//..//..//data//spire//cal//SCal"
#dataDir = "//disks//winchester2//calibration_data//"
#LOCAL VERSION
directory = Configuration.getProperty('var.hcss.workdir')
dataDir = Configuration.getProperty('var.hcss.workdir')

# if inputCalDirTree is True, then calibration products are read from a
#  calibration directory tree
#  e.g. RSRF is read from <directory>/Phot/SCalPhotRsrf/SCalPhotRsrf.fits
# If inputCalDirTree is False, then required calibration products are read from
#  a calibration file (either a local pool or a jarfile)
inputCalDirTree=True
if not inputCalDirTree:
	#read calibration tree from pool
	cal=spireCal(pool='spire_cal_12_3')
	## alternatively, read from jarFile
	#cal=spireCal(pool=os.path.join(dataDir,'spire_cal_12_1_test.jar'))

# if outputCalDirTree is True, then calibration products are written to a
#  calibration directory tree
#  e.g. <directory>/Phot/SCalPhotColorCorrK/SCalPhotRadialCorrBeam_<version>.fits
# If outputCalDirTree is False, then calibration products are written to
#  dataDir e.g. <dataDir>/SCalPhotRadialCorrBeam_<version>.fits
outputCalDirTree=True

writeAscii = False
writeLog = False
#-------------------------------------------------------------------------------
#===============================================================================
#=====                           INPUT PARAMETERS                          =====
#===============================================================================
#-------------------------------------------------------------------------------
# version number
version = "4EP"

#read in core beam
beamCoreIn = TableDataset()
beamCoreIn.addColumn('radius',Column(Double1d.range(1400),unit=Angle.SECONDS_ARC))
beamProfInFileNames = ""
for band in ['PSW','PMW','PLW']:
    beamProfInFileName = '%s_MCore_9.csv'%band
    beamProfInFileNames+", "+beamProfInFileName
    beamProfIn=asciiTableReader(os.path.join(dataDir, beamProfInFileName))['c0'].data
    #set central pixel to 1
    beamProfIn[0]=1.0
    beamCoreIn.addColumn(band,Column(beamProfIn))
beamNewFileCore='beamCore_v%s.csv'%version
asciiTableWriter(beamCoreIn,os.path.join(directory,beamNewFileCore))

#create constant beam filled with zeros
beamConstantIn = TableDataset()
beamConstantIn.addColumn('radius',Column(Double1d.range(1400),unit=Angle.SECONDS_ARC))
for band in ['PSW','PMW','PLW']:
    beamConstantIn.addColumn(band,Column(Double1d(1400)))
beamNewFileConstant='beamConstant_v%s.csv'%version
asciiTableWriter(beamConstantIn,os.path.join(directory,beamNewFileConstant))

# set format version and date format
formatVersion = "1.0"
df  = java.text.SimpleDateFormat("yyyy.MM.dd/HH:mm:ss/z")

# set verbose to print more information during processing
verbose = True

# beam parameters
arcsec2Sr = (Math.PI/(60.*60.*180))**2
# area measured on Neptune
# Taken from Bernhard's presentation at the SDAG on 30 October 2014.
spireAreaEffFreq = {"PSW":454.0*arcsec2Sr, "PMW":803.4*arcsec2Sr, "PLW":1687.0*arcsec2Sr}
spireAreaEffFreqErr = {"PSW":1.0*arcsec2Sr, "PMW":4.0*arcsec2Sr, "PLW":12.0*arcsec2Sr}
# Neptune spectral index (from ESA4 model)
alphaNep={"PSW":1.29, "PMW":1.42, "PLW":1.47}

# Exponent of powerlaw describing FWHM dependence on frequency
# FWHM ~ freq**gamma
gamma = -0.85 


#-------------------------------------------------------------------------------
#===============================================================================
#=====                         READ FILTER PROFILES                        =====
#===============================================================================
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Loading physical constants
h = Constant.H_PLANCK.value
k = Constant.K_BOLTZMANN.value
c = Constant.SPEED_OF_LIGHT.value

# Frequency raster of common frequency grid
deltaNu = 0.1e9		# 0.1 GHz
nuMin   = 300.e9
nuMax   = 1800.e9
nNu     = FIX((nuMax-nuMin)/deltaNu)
freq    = Double1d(range(nNu)) * deltaNu + nuMin

#-------------------------------------------------------------------------------
# Load SPIRE filter functions fron calibration directory structure

# SPIRE band names
spireBands=["PSW","PMW","PLW"]

if inputCalDirTree:
	# SPIRE Photometer RSRF calibration product from cal tree
	rsrfVersion = "3"
	rsrf = fitsReader("%s//Phot//SCalPhotRsrf//SCalPhotRsrf_v%s.fits"%(directory, rsrfVersion))
	if verbose:
		print 'Reading RSRF version %s from calibration directory tree'%(rsrfVersion)
else:
	#read RSRF from calibration tree
	rsrf=cal.getPhot().getProduct('Rsrf')
	rsrfVersion=rsrf.getVersion()
	if verbose:
		print 'Reading RSRF version %s from calibration %s'%(rsrfVersion,cal.getVersion())

spireFreq   = rsrf['rsrf']['frequency'].data*1e9  # Frequency in Hz
#indexes of freq in rsrf
ixR = freq.where((freq>=MIN(spireFreq)) & (freq<=MAX(spireFreq)))

# spire RSRF only
spireFiltOnly={}

#interpolate to freq array
for band in spireBands:
	#create Rsrf interpolation object
	interpRsrf = LinearInterpolator(spireFreq, rsrf.getRsrf(band))
	#make arrays for final object
	spireFiltOnly[band] = Double1d(nNu)
	#interpolate Rsrf to freq array
	spireFiltOnly[band][ixR] = interpRsrf(freq[ixR])

#-------------------------------------------------------------------------------
#===============================================================================
#=====                           DEFINE FUNCTIONS                          =====
#===============================================================================
#-------------------------------------------------------------------------------

# Set up functions to calculate beam profile, effective frequency and effective area
# Function list:
#   * spireMonoBeam: Calculate monochromatic beam profile and area at a given frequency
#   * spireMonoAreas: Calculate monochromatic beam areas over a range of frequencies
#   * spireEffArea: Calculate the effective beam area for a given spectrum
#   * spireFindEffFreq: Calculate the effective frequency for SPIRE

#-------------------------------------------------------------------------------
# Calculate monochromatic beam profile and area at a given frequency
def spireMonoBeam(freqx,beamRad,beamProfs,effFreq,gamma,array):
	"""
	========================================================================
	Implements the full beam model to generate the monochromatic beam profile
	and corresponding monochromatic beam solid angle at a given frequency.

	Inputs:
	  freqx:     (float) frequency [Hz] for which monochromatic beam
	               profile and area should be calculated
	  beamRad:   (array float) radius vector from the input beam profiles
	  beamProfs: (dataset) PhotRadialCorrBeam object
	               [used to retrieve core beam profile]
	  effFreq:   (float) effective frequency [Hz] of array
	  gamma:     (float) Exponent of powerlaw describing FWHM dependence
	               on frequency
	  array:     (string) spire array ('Psw'|'Pmw'|'Plw')

	Outputs:     (list of objects)
	  [0]:       (float) Beam area [arcsec^2] at frequency freqx
	  [1]:       (array float) Monochromatic beam profile at frequency freqx

	Calculation:
          Scales the core beam profile width as (freqx/effFreq)^gamma.
	  Queries the calibration file to generate new core beam profile.
	  Uses constant beam profile where it is larger than core profile.
	  Integrates over radius to caluclate beam area.

	Dependencies:
	  herschel.ia.numeric.toolbox.interp.LinearInterpolator
	  herschel.ia.numeric.toolbox.integr.TrapezoidalIntegrator
	  
	2013/12/19  C. North  initial version
	2014/11/07  C. North  removed beamConst (deprecated) from calculations

	"""

	#calculate the "scaled" radius, as nu^-gamma
	radNew=beamRad*(freqx/effFreq)**-gamma
	maxRad=max(beamRad)
	nRad=len(beamRad)
	#ensure it doesn't go out of range
	radNew[radNew.where(radNew > maxRad)]=maxRad
	#get the corresponding beam profiles
	beamNew=Double1d(nRad)
	for r in range(nRad):
		beamNew[r]=beamProfs.getCoreCorrection(radNew[r],array)
	
	####  DEPRECATED  ####
	#apply the "constant" beam where appropriate
	#beamConst=beamProfs.getConstantCorrectionTable().getColumn(array).data
	#isConst=beamNew.where(beamNew < beamConst)
	#beamNew[isConst]=beamConst[isConst]
	####  /DEPRECATED  ####
	
	#integrate to get solid angle (in arcsec^2)
	
	beamInterp=LinearInterpolator(beamRad,beamNew * 2. * Math.PI * beamRad)
	integrator=TrapezoidalIntegrator(0,maxRad)
	beamMonoArea=integrator.integrate(beamInterp)

	return(beamMonoArea,beamNew)

#-------------------------------------------------------------------------------
# Calculate monochromatic beam areas over a range of frequencies
def spireMonoAreas(freq,beamProfs,effFreq,gamma,array,freqFact=100):

	"""
	========================================================================
	Generates array of monochromatic beam areas over frequency range by
	calculating over a sparser array and interpolating

	Inputs:
	  freq:      (array float) frequency vector [Hz] for which monochromatic
	               beams areas should be calculated
	  beamProfs: (dataset) PhotRadialCorrBeam object from calibration tree
	  effFreq:   (float) effective frequency [Hz] of array
	  gamma:     (float) Exponent of powerlaw describing FWHM dependence
	               on frequency
	  array:     (string) spire array ('Psw'|'Pmw'|'Plw')
	  freqFact:  (int) Factor by which to reduce size of freq.
	               OPTIONAL. Default=100.

	Outputs:     
	             (array float) Monochromatic Beam area [sr] at frequencies
	                corresponding to freq

	Calculation:
          Geneates sparse frequency array of full range
	  Uses spireMonoBeam to calculate monochromatic beam area at sparse freqs
	  Interpolates to full frequency grid

	Dependencies:
	  spireMonoBeam
	  herschel.ia.numeric.toolbox.interp.CubicSplineInterpolator
	  
	2013/12/19  C. North  initial version
	2014/11/07  C. North  removed beamConst (deprecated) from calculations

	"""

	#set up a sparser range of frequencies (otherwise it takes too long)
	nNu=len(freq)
	nNuArea=nNu/freqFact + 1
	#array of indices of full frequency array to use
	iNuArea=Int1d(range(nNuArea))*freqFact
	iNuArea[-1]=nNu-1

	#set up arrays
	beamMonoFreqSparse=Double1d(nNuArea)
	beamMonoAreaSparse=Double1d(nNuArea)

	#get beam radius array from calibration table
	beamRad=beamProfs.getCoreCorrectionTable().getColumn('radius').data
	####  DEPRECATED  ####
	#get constant beam profile from calibration table
	beamConst=beamProfs.getConstantCorrectionTable().getColumn(array).data
	####  /DEPRECATED  ####
	
	# calculate at sparse frequencies
	for fx in range(nNuArea):
		#get corresponding index in full frequency array
		f=iNuArea[fx]
		#populate frequency array
		beamMonoFreqSparse[fx]=freq[f]
		#populate beam area array
		beamMonoAreaSparse[fx]=spireMonoBeam(freq[f],beamRad,beamProfs,effFreq,gamma,array)[0]

	# interpolate to full frequency array and convert to Sr
	beamInterp=CubicSplineInterpolator(beamMonoFreqSparse,beamMonoAreaSparse)
	beamMonoArea=beamInterp(freq)*arcsec2Sr #in sr
	
	return(beamMonoArea)

#-------------------------------------------------------------------------------
# Calculate the effective beam area for a given spectrum
def spireEffArea(freq, transm, monoArea, BB=False, temp=20.0, beta=1.8, alpha=-1.0):
	"""
	========================================================================
	Calculate the effective beam area for a source of a given spectrum

	Inputs:
	  freq:       (array float) frequency vector corresponding to RSRF values [Hz]
	  transm:     (array float) relative spectral response (RSRF) corresponding to freq
	                Note that this should *not* include the aperture efficiency
	  monoArea:   (array float) monochromatic beam solid angle corresponding
	                to frequencies in freq
	  BB:         (boolean) spectral function to use for source spectrum:
	                'True': a modified black body
	                'False' a power-law with exponent alpha=-1
	                OPTIONAL. Default=False
	  temp:       (float) Dust/sky temperature (if BB=True)
	                OPTIONAL. Default=20.0
	  beta:       (float) Dust/sky spectral index (if BB=True)
	                OPTIONAL. Default=1.8
	  alpha:      (float) Exponent of power-law sky background model (if BB=False)
	                OPTIONAL. Default=-1

	Outputs:     
	            (float) Beam area for given spectrum, in same units as monoArea

	Calculation:
	  Calculates the source spectrum (either modifies black body or power law)
          Multiplies the monochromatic beam area by RSRF and source spectrum
	  Integrates over frequency
	  Normalises by integral over frequency of RSRF and source spectrum

	Dependencies:
	  herschel.ia.numeric.toolbox.interp.LinearInterpolator
	  herschel.ia.numeric.toolbox.integr.TrapezoidalIntegrator

	2013/12/19  C. North  initial version

	"""	
	#
	# Calculate sky background model
	#
	if BB == 1:
		#print temp,beta,c,h,k
		fSky  = 2*h * freq**3 / c**2 / (EXP(h*freq/k/temp) - 1.) * freq**beta
	#
	# 2) As a Power-Law
	else:
		fSky  = freq**alpha

	# Integrate monochromatic area over frequency, weighted by rsrf and fSky
	numInterp=LinearInterpolator(freq,transm * fSky * monoArea)
	denomInterp=LinearInterpolator(freq,transm * fSky)
	minFreq=min(freq)
	maxFreq=max(freq)
	integrator=TrapezoidalIntegrator(minFreq,maxFreq)
	numInteg=integrator.integrate(numInterp)
	denomInteg=integrator.integrate(denomInterp)
	effArea = numInteg / denomInteg

	return(effArea)


#-------------------------------------------------------------------------------
# Calculate the effective frequency for SPIRE
def spireFindEffFreq(freq, rsrf, beamProfs, effFreqInit, gamma,
  areaNep, alphaNep, array, simpleBeam=False, freqFact=500,
  initRange=0.001, reqPrec=1.e-6, maxIter=5, verbose=False):
	"""
	========================================================================
	Derive effective frequency for spire bands. This is the frequency at
	which the monochromatic beam area is equal to the area as measured on
	Neptune.

	Inputs:
	  freq:        (array float) frequency vector corresponding to RSRF values [Hz]
	  rsrf:        (array float) relative spectral response (RSRF) corresponding to freq
	                 Note that this should *not* include the aperture efficiency
	  beamProfs:   (dataset) PhotRadialCorrBeam object from calibration tree
	  effFreqInit: (float) Initial estimate of effective frequency [Hz]
	  gamma:       (float) Exponent of powerlaw describing FWHM dependence
	                 on frequency
	  areaNep:     (float) Solid angle of beam measured on Neptune [sr]
	  alphaNep:    (float) Spectral index of Neptune frequency spectrum
	  array:       (string) spire array ('Psw'|'Pmw'|'Plw')
	  simpleBeam:  (boolean) set to use simple beam area model, which
                         scales as freq^(2.gamma)
	  freqFact:    (int) Factor by which to reduce size of freq for calculations.
	                 Only applicable if SimpleBeam=False
	                 OPTIONAL. Default=500.
	  initRange:   (float) Fractional intial range to use in calculations
	                 OPTIONAL. Default=0.01
	  reqPrec:     (float) Required relative precision for convergence.
	                 OPTIONAL. Default=1.e-6
          maxIter:     (int) Maximum interations to try
	                 OPTIONAL. Default=5
	  verbose:     (boolean) set to print more detailed info
                         OPTIONAL. Default=False

	Outputs:
	               (float) Effective frequency [Hz]

	Calculation:
	  Uses effFreqInit +/- initRange to calculate monochromatic Areas
	  Calculates effective beam area for Neptune spectrum (alphaNep)
 	  Compares with measured Neptune beam area (areaNep)
	  Adjusts effFreqInit and iterates until maximum iterations (maxIter) or
	  required precicion (reqPrec) is reached

	Dependencies:
	  spireMonoAreasSimple
	  spireMonoAreas
	  spireEffArea

	2014/01/07  C. North  initial version
	"""
		
	#calculate initial estimates of effective Frequency
	#parameterised be offset from original estimate
	relEff=Double1d([1.-initRange,1.+initRange])
	effFreqs=effFreqInit*relEff
	beamMonoAreaInit=spireMonoAreas(freq,beamProfs,effFreqInit,gamma,array,freqFact=freqFact)
	effAreaInit=spireEffArea(freq,rsrf, beamMonoAreaInit, BB=False, alpha=alphaNep)
	arcsec2Sr = (Math.PI/(60.*60.*180))**2
	print '  Initial %s Effective Frequency: %.2f GHz'%(array,effFreqInit/1.e9)
	print '  Initial %s Effective Area: %.2f arcsec^2'%(array,effAreaInit/arcsec2Sr)
	#calculate effective beam area for initial estimates of effFreq
	if simpleBeam:
		beamMonoArea0=spireMonoAreasSimple(freq,effFreqs[0],areaNep,gamma)
		beamMonoArea1=spireMonoAreasSimple(freq,effFreqs[1],areaNep,gamma)
	else:
		beamMonoArea0=spireMonoAreas(freq,beamProfs,effFreqs[0],gamma,array,freqFact=freqFact)
		beamMonoArea1=spireMonoAreas(freq,beamProfs,effFreqs[1],gamma,array,freqFact=freqFact)
	beamMonoDiff=beamMonoArea1-beamMonoArea0
	effAreas=Double1d(2)
	effAreas[0]=spireEffArea(freq,rsrf, beamMonoArea0, BB=False, alpha=alphaNep)
	effAreas[1]=spireEffArea(freq,rsrf, beamMonoArea1, BB=False, alpha=alphaNep)
	if verbose:
		print 'Calculating effective frequency for %s'%array
		print '  Initial %s Effective Frequencies: [%.2f : %.2f] GHz'%(array,effFreqs[0]/1.e9,effFreqs[1]/1.e9)
		print '  Initial %s Effective Areas: [%.2f : %.2f] arcsec^2'%(array,effAreas[0]/arcsec2Sr,effAreas[1]/arcsec2Sr)
	iter=0
	done=False
	while ((done==False) and (iter <= maxIter)):
		iter=iter+1
		#difference from measured beam area
		diffAreas=effAreas-areaNep
		relAreas=diffAreas/areaNep
		#calculate new esitmate of rel
		grad=(diffAreas[1]-diffAreas[0])/(relEff[1]-relEff[0])
		relEffNew=relEff[1] - diffAreas[1]/grad

		#move values in arrays
		relEff[0]=relEff[1]
		effAreas[0]=effAreas[1]

		#calculate new effective beam area
		relEff[1]=relEffNew
		effFreqs=effFreqInit*relEff

		if simpleBeam:
			beamMonoNew=spireMonoAreasSimple(freq,effFreqs[1],areaNep,gamma)
		else:
			beamMonoNew=spireMonoAreas(freq,beamProfs,effFreqs[1],gamma,array,freqFact=freqFact)
		effAreas[1]=spireEffArea(freq,rsrf, beamMonoNew, BB=False, alpha=alphaNep)

		diffAreas=effAreas-areaNep
		relAreas=diffAreas/areaNep
		if verbose:
			print '    iter %d: %.4f [effFreq %.2f GHz], Area=%.2f, RelDiff=%.4g'%(iter,relEff[1],effFreqs[1]/1.e9,effAreas[1]/arcsec2Sr,relAreas[1])
		if (Math.abs(relAreas[1]) < reqPrec):
			done=True
	if ((iter > maxIter) and (done==False)):
		print "  Warning: maximum iterations [%d] exceeded without conversion [%g]"%(maxIter,reqPrec)

	if verbose:
		print '  Final %s effFreq: %.4f'%(array,effFreqs[1]/1.e9)
		print '  Resulting %s Neptune area: %.2f [rel Diff: %.3g]'%(array,effAreas[1]/arcsec2Sr,relAreas[1])
	return(effFreqs[1])

#-------------------------------------------------------------------------------
#===============================================================================
#=====                           END OF FUNCTIONS                          =====
#===============================================================================
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#===============================================================================
#=====                        CALCULATE BEAM PROFILES                      =====
#===============================================================================
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

#create product and set basic metadata
# define the start and end dates for the product
# starting at beginning of PFM1
startDate = df.parse("2005.02.22/00:00:00/GMT")
endDate   = df.parse("2020.01.01/00:00:00/GMT")

beamProfs = herschel.spire.ia.dataset.PhotRadialCorrBeam()
beamProfs.meta["modelName"].value = "FM"
beamProfs.meta["creator"].value = scriptVersionString
beamProfs.meta["creationDate"].value = FineTime(java.util.Date())
beamProfs.meta["startDate"].value = FineTime(startDate)
beamProfs.meta["endDate"].value = FineTime(endDate)
beamProfs.setVersion(version)
beamProfs.setFormatVersion(formatVersion)
beamProfs.meta['dataOrigin'] = metaDict.newParameter('dataOrigin', "%s; RSRF v%s"%(beamProfInFileNames[2:],rsrfVersion))
beamProfs.meta["author"]  = metaDict.newParameter("author", "Chris North")

#read in beams from file
beamProfs["constant"] = beamConstantIn
beamProfs["core"] = beamCoreIn
#update labels
beamProfs["constant"].setDescription("[DEPRECATED, set to zero] Frequency-independent constant part of the radial beam profile")
beamProfs["core"].setDescription("Frequency-dependent core part of the radial beam profile")


# Re-compute normalised beam areas for new beams
# get beam radius
beamRad=beamProfs.getCoreCorrectionTable().getColumn('radius').data
nRad=len(beamRad)

# Create normArea table
beamNormArea=TableDataset(description="Area as a function of radius, normalised by final value")
# add radius column
beamNormArea.addColumn("radius",Column(Float1d(beamRad),unit=Angle.SECONDS_ARC))

print 'Calculating normalised beam area...'
for band in spireBands:
	print band
	# add column to table
	beamNormArea.addColumn(band,Column(Float1d(nRad)))
	# get core beam
	beamComb=beamProfs.getCoreCorrectionTable().getColumn(band).data.copy()
	####  DEPRECATED  ####
	# get const beam
	#beamConst=beamProfs.getConstantCorrectionTable().getColumn(band).data
	# work out where constant beam applies
	#isConst = beamComb.where(beamComb < beamConst)
	# apply constant beam where applicable
	#beamComb[isConst] = beamConst[isConst]
	####  /DEPRECATED  ####
	# interpolate and integrate
	beamInterp=CubicSplineInterpolator(Double1d(beamRad),beamComb*2.*Math.PI*beamRad)
	integrator=TrapezoidalIntegrator(0,max(beamRad))
	beamTotArea=integrator.integrate(beamInterp)
	for r in range(nRad):
		integrator=TrapezoidalIntegrator(0,beamRad[r])
		beamNormArea[band].data[r]=integrator.integrate(beamInterp)/beamTotArea
# write table to beamProfs object
beamProfs["normArea"]=beamNormArea

#-------------------------------------------------------------------------------
# Add metadata
# add gamma to metadata
beamProfs.meta['gamma'] = metaDict.newParameter('gamma', gamma)
# add Neptune beam area to metaData (square arseconds)
beamProfs.meta['beamNeptunePswArc']= metaDict.newParameter('beamNeptunePswArc', spireAreaEffFreq['PSW']/arcsec2Sr)
beamProfs.meta['beamNeptunePmwArc']= metaDict.newParameter('beamNeptunePmwArc', spireAreaEffFreq['PMW']/arcsec2Sr)
beamProfs.meta['beamNeptunePlwArc']= metaDict.newParameter('beamNeptunePlwArc', spireAreaEffFreq['PLW']/arcsec2Sr)
# add Neptune beam area error to metaData (square arseconds)
beamProfs.meta['beamNeptunePswArcErr']= metaDict.newParameter('beamNeptunePswArcErr', spireAreaEffFreqErr['PSW']/arcsec2Sr)
beamProfs.meta['beamNeptunePmwArcErr']= metaDict.newParameter('beamNeptunePmwArcErr', spireAreaEffFreqErr['PMW']/arcsec2Sr)
beamProfs.meta['beamNeptunePlwArcErr']= metaDict.newParameter('beamNeptunePlwArcErr', spireAreaEffFreqErr['PLW']/arcsec2Sr)
# add Neptune beam area to metaData (steradians)
beamProfs.meta['beamNeptunePswSr']= metaDict.newParameter('beamNeptunePswSr', spireAreaEffFreq['PSW'])
beamProfs.meta['beamNeptunePmwSr']= metaDict.newParameter('beamNeptunePmwSr', spireAreaEffFreq['PMW'])
beamProfs.meta['beamNeptunePlwSr']= metaDict.newParameter('beamNeptunePlwSr', spireAreaEffFreq['PLW'])
#add Neptune spectral index to metadata (steradians)
beamProfs.meta['alphaNeptunePsw']= metaDict.newParameter('alphaNeptunePsw', alphaNep['PSW'])
beamProfs.meta['alphaNeptunePmw']= metaDict.newParameter('alphaNeptunePmw', alphaNep['PMW'])
beamProfs.meta['alphaNeptunePlw']= metaDict.newParameter('alphaNeptunePlw', alphaNep['PLW'])

#-------------------------------------------------------------------------------
#===============================================================================
#=====                    CALCULATE EFFECTIVE FREQUENCIES                  =====
#===============================================================================
#-------------------------------------------------------------------------------
#uses existing numbers as initial guess
spireEffFreq = {"PSW":1224.06*1.e9,\
	"PMW":878.07*1.e9,\
	"PLW":609.86*1.e9}
print '\nGenerating new effective frequencies...'

for band in spireBands:
	spireEffFreq[band] = spireFindEffFreq(freq, spireFiltOnly[band],
	  beamProfs,spireEffFreq[band], gamma, spireAreaEffFreq[band],
	  alphaNep[band], band, verbose=verbose)

# update RadialBeamCorr metadata (in GHz)
beamProfs.meta['freqEffPsw'] = metaDict.newParameter('freqEffPsw', spireEffFreq['PSW']/1.e9)
beamProfs.meta['freqEffPmw'] = metaDict.newParameter('freqEffPmw', spireEffFreq['PMW']/1.e9)
beamProfs.meta['freqEffPlw'] = metaDict.newParameter('freqEffPlw', spireEffFreq['PLW']/1.e9)

#-------------------------------------------------------------------------------
#===============================================================================
#=====                  CALCULATE MONOCHROMATIC BEAM AREAS                 =====
#===============================================================================
#-------------------------------------------------------------------------------

#calculate monochromatic beam areas using full or simple beam treatment
beamMonoArea = {}
beamAreaPipSr = {}
beamAreaPipArc = {}
print '\nCalculating monochromatic beam areas...'
for band in spireBands:
	print band
	#monochromatic beam areas
	beamMonoArea[band] = spireMonoAreas(freq, beamProfs, 
	  spireEffFreq[band], gamma, band)
	#pipeline beam areas
	beamAreaPipSr[band]=spireEffArea(freq, spireFiltOnly[band], \
	  beamMonoArea[band], BB=False, alpha=-1)
	beamAreaPipArc[band]=beamAreaPipSr[band]/arcsec2Sr

beamProfs.meta['beamPipelinePswSr'] = metaDict.newParameter('beamPipelinePswSr', beamAreaPipSr['PSW'])
beamProfs.meta['beamPipelinePmwSr'] = metaDict.newParameter('beamPipelinePmwSr', beamAreaPipSr['PMW'])
beamProfs.meta['beamPipelinePlwSr'] = metaDict.newParameter('beamPipelinePlwSr', beamAreaPipSr['PLW'])
beamProfs.meta['beamPipelinePswArc']= metaDict.newParameter('beamPipelinePswArc', beamAreaPipArc['PSW'])
beamProfs.meta['beamPipelinePmwArc']= metaDict.newParameter('beamPipelinePmwArc', beamAreaPipArc['PMW'])
beamProfs.meta['beamPipelinePlwArc']= metaDict.newParameter('beamPipelinePlwArc', beamAreaPipArc['PLW'])

#-------------------------------------------------------------------------------
# Output to FITS files
if outputCalDirTree:
	filename = java.io.File(r"%s//Phot//SCalPhotRadialCorrBeam//SCalPhotRadialCorrBeam_v%s.fits"%(directory, version))
else:
	filename = java.io.File(r"%s//SCalPhotRadialCorrBeam_v%s.fits"%(dataDir, version))

beamProfs.meta['fileName'] = StringParameter(value=filename.name, description="Name of file when exported")
fitsWriter = FitsArchive()
fitsWriter.rules.append(metaDict.getFitsDictionary())
fitsWriter.save(filename.toString(), beamProfs)

#-------------------------------------------------------------------------------
# Write tables to ascii files 
if writeAscii:
   asciiFileCore='RadialCorrBeam_core_v%s.csv'%version
   asciiFileConstant='RadialCorrBeam_constant_v%s.csv'%version
   asciiFileNormArea='RadialCorrBeam_normArea_v%s.csv'%version
   if outputCalDirTree:
      dirRadialCorrBeam='%s//Phot//SCalPhotRadialCorrBeam'%directory
   else:
      dirRadialCorrBeam=dataDir
   asciiTableWriter(table=beamProfs['core'],file=os.path.join(dirRadialCorrBeam,asciiFileCore))
   asciiTableWriter(table=beamProfs['constant'],file=os.path.join(dirRadialCorrBeam,asciiFileConstant))
   asciiTableWriter(table=beamProfs['normArea'],file=os.path.join(dirRadialCorrBeam,asciiFileNormArea))
   if verbose:
      print 'RadialCorrBeam written to ascii files'

#-------------------------------------------------------------------------------
# Write log file
# Define functions to write syntax
def logMeta(name, meta, type='double'):
	"""
	Produce a metadata string for the log file
	
	Inputs:
	  name: (string) name of table in obj
	  meta: (Mtadata) metadata object containing key [name]
	  type: (string) type of metadate [double|string|int]
	        OPTIONAL. Default is double.

	Outputs:
	  	(string) String containing name, value, unit and description

	2014/01/20  C. North  initial version

	"""

	if type=='double':
		#check if there's a unit
		if herschel.share.util.StringUtil.asString(meta[name].unit)=='null':
			unit=''
		else:
			unit=meta[name].unit.getName()
		metaStr='%s: %.9g [%s] (Type=Double, Desc="%s")'%(name, meta[name].double, unit, meta[name].description)

	elif type=='int':
		#check if there's a unit
		if herschel.share.util.StringUtil.asString(meta[name].unit)=='null':
			unit=''
		else:
			unit=meta[name].unit.getName()
		metaStr='%s: %d [%s] (Type=Integer, Desc="%s")'%(name, meta[name].int, unit, meta[name].description)

	elif type=='string':
		metaStr='%s: "%s" (Type=String, Desc="%s")'%(name, meta[name].string, meta[name].description)

	return(metaStr)

def logTable(name, file, obj):
	"""
	Produce a table string for the log file
	
	Inputs:
	  name: (string) name of table in obj
	  file: (string) filename of ascii file containing table
	  obj:  (Product) calibration product containing table [name]

	Outputs:
	  	(string) String containing name, file and description

	2014/01/20  C. North  initial version

	"""
	tableStr='%s: %s (Desc="%s")'%(name, file, obj[name].getDescription())
	return(tableStr)

if writeLog:
   beamLogFile = 'RadialCorrBeam_log_v%s.dat'%version
   beamLog=open(os.path.join(dirRadialCorrBeam,beamLogFile),'w')
   beamLog.write('Beam Radial Profiles Version %s\n'%version)
   beamLog.write('Creation Date: %s\n'%java.util.Date())
   beamLog.write('\nTables:\n')
   beamLog.write('  %s\n'%(logTable('core',asciiFileCore,beamProfs)))
   beamLog.write('  %s\n'%(logTable('constant',asciiFileConstant,beamProfs)))
   beamLog.write('  %s\n'%(logTable('normArea',asciiFileNormArea,beamProfs)))
   beamLog.write('\nMetadata:\n')
   beamLog.write('  %s\n'%(logMeta('gamma',beamProfs.meta,'double')))
   beamLog.write('  %s\n'%(logMeta('beamNeptunePswArc',beamProfs.meta,'double')))
   beamLog.write('  %s\n'%(logMeta('beamNeptunePmwArc',beamProfs.meta,'double')))
   beamLog.write('  %s\n'%(logMeta('beamNeptunePlwArc',beamProfs.meta,'double')))
   beamLog.write('  %s\n'%(logMeta('beamNeptunePswSr',beamProfs.meta,'double')))
   beamLog.write('  %s\n'%(logMeta('beamNeptunePmwSr',beamProfs.meta,'double')))
   beamLog.write('  %s\n'%(logMeta('beamNeptunePlwSr',beamProfs.meta,'double')))
   beamLog.write('  %s\n'%(logMeta('alphaNeptunePsw',beamProfs.meta,'double')))
   beamLog.write('  %s\n'%(logMeta('alphaNeptunePmw',beamProfs.meta,'double')))
   beamLog.write('  %s\n'%(logMeta('alphaNeptunePlw',beamProfs.meta,'double')))
   beamLog.write('  %s\n'%(logMeta('freqEffPsw',beamProfs.meta,'double')))
   beamLog.write('  %s\n'%(logMeta('freqEffPmw',beamProfs.meta,'double')))
   beamLog.write('  %s\n'%(logMeta('freqEffPlw',beamProfs.meta,'double')))
   beamLog.write('  %s\n'%(logMeta('beamPipelinePswSr',beamProfs.meta,'double')))
   beamLog.write('  %s\n'%(logMeta('beamPipelinePmwSr',beamProfs.meta,'double')))
   beamLog.write('  %s\n'%(logMeta('beamPipelinePlwSr',beamProfs.meta,'double')))
   beamLog.write('  %s\n'%(logMeta('beamPipelinePswArc',beamProfs.meta,'double')))
   beamLog.write('  %s\n'%(logMeta('beamPipelinePmwArc',beamProfs.meta,'double')))
   beamLog.write('  %s\n'%(logMeta('beamPipelinePlwArc',beamProfs.meta,'double')))
   beamLog.close()
