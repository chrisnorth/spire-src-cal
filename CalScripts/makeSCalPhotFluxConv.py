# 
#  This file is part of Herschel Common Science System (HCSS).
#  Copyright 2001-2011 Herschel Science Ground Segment Consortium
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
# $Id: makeSCalPhotFluxConv.py,v 1.17 2014/11/14 16:08:03 epoleham Exp $
# Copyright (c) 2008 CCLRC
#
# Script to produce a FITS file for the photometer flux conversion calibration
# product
#
# Author: E. Polehampton following script by Arnie Schwartz & Nanyao Lu
#                      (last revised by Nanyao Lu on July 16, 2008)
#
# History:
# E.Polehampton  25-09-2008       - First Version
# E.Polehampton  22-10-2008       - modified for new delivery from Nanyao which is a FITS file with
#                                   correct structure - just need to add standard metadata
# E.Polehampton  20-03-2009       - update to latest file from Nanyao (v1.6) and to version 1
# E.Polehampton  04-06-2009       - Update to read Nanyao's input data v1.8
#                                 - include biasMode and time editions
# E.Polehampton  05-06-2009       - Update creator string
# E.Polehampton  12-09-2009       - Update description (SPIRE-1490), add editions for flight (SPIRE-1952).
# E.Polehampton  17-09-2009       - correction in latest edition for nominal mode (2.1->2.2) (SPIRE-1972)
# E.Polehampton  25-11-2009       - update with latest files (SPIRE-2191)
# E.Polehampton  07-05-2010       - update with numbers from George (test version)
# E.Polehampton  09-09-2010       - update with final numbers from George (SPIRE-2801)
# E.Polehampton  25-01-2011       - update numbers from George (SPIRE-3150)
# E.Polehampton  31-03-2011       - fix respControlStamp (SPCAL-11)
# E.Polehampton  03-05-2013       - include George's new numbers (SPCAL-76)
# E.Polehampton  21-01-2014       - put K4P and K4E into main metadata, and update their values (SPCAL-91)
# E.Polehampton  26-02-2014       - update K4E  (SPCAL-110)
# E.Polehampton  14-11-2014       - read K4 params from text file written by ColorCorrK script (SPCAL-126)
#
# -----------------------------------------------------------------------------
#
# This isn't used at the moment.
def hpXcalKcorr(freq0, freq, transm, BB=True, temp=20.0, beta=1.8, alpha=-1.0, gamma=0.0):
	"""
	================================================================================
	Calculation of the K-correction factor from isophotal flux to a monochromatic 
	flux-density at a given reference frequency (data to be multiplied!)
	This routine is needed by hpXcalColorCorr.py
	
	Inputs:
	  freq0:    (float) waveband reference frequency [Hz] for which monochromatic
	            flux-density is given
	  freq:     (array float) frequency vector corresponding to RSRF values [Hz]
	  transm:   (array float) relative spectral response (RSRF) corresponding to freq
	  BB:       (boolean) spectral function to use for source spectrum:
	                'True': a modified black body
	                'False' a power-law with exponent alpha=-1
	  temp:     (float) Dust/sky temperature [optional; default is 20K; 
	                only for modified black body]
	  beta:     (float) Dust/sky spectral index [optional; default is 1.8; 
	                only for modified black body]
	  alpha:    (float) Exponent of power-law sky background model
	  gamma:    (float) Exponent of powerlaw describing FWHM dependence on frequency
	
	Outputs:
	 (list)     1st item: K-correction factor
	            2nd item: Sky emission at reference fequency (fSky0)
	
	Calculation:
	  Depending on the state of the input parameter BB, either the spectrum of a
	  Planck function multiplied by frequency to the power of beta, or a power-law
	  spectrum with spectral index alpha is calculated for all values in the vector 
	  frequ. In addition the same value is calculated at the discrete frequency 
	  freq0. Then the product of this value and the integral of the RSRF over all
	  frequencies, divided by the integral over all products of frequency and RSRF	
	  is calculated. Note that the integrals are coded as simple sums as the 
	  HIPE numeric integral doesn't allow too may discrete points and the RSRF
	  is sampled to quite some detail.
	
	
	2012/04/18  L. Conversi  initial version in file hp_xcal_Total_v1.6.py.txt
	2012/08/22  B. Schulz    added comments, reformatted, renamed kCorr to hpXcalKcorr
	2012/08/24  B. Schulz    fixed defaults for temp and beta, changed inputs to frequency,
	                         updated header and comments, renamed inputFreq to freq
	                         and inputFilt to transm added powerlaw sky background
	2012/08/27  B. Schulz    brush -up on header and comments
	2013/03/28  B. Schulz    removed implicit limitation to fixed frequency interval
	2013/04/05  B. Schulz    implemented proper tabulated integration function
	2013/06/18  L. Conversi  implemented gamma parameter in the denominator
	                         (previsuly applie directly to RSRFs, i.e. in the numerator too)
	
	================================================================================
	"""
	#
	# Calculate sky background model
	#
	# 1) As a modified Blackbody
	if BB == 1:
		fSky  = 2*h * freq**3 / c**2 / (EXP(h*freq/k/temp) - 1.) * freq**beta
		fSky0 = 2*h * freq0**3 / c**2 / (EXP(h*freq0/k/temp) - 1.) * freq0**beta
	#
	# 2) As a Power-Law
	else:
		fSky  = freq**alpha
		fSky0 = freq0**alpha
	#
	# Using the K-correction as defined in th OM (eq. 5.5, 5.15 & 5.16)
	kWave = fSky0 * IntTabulated(freq)(transm) / IntTabulated(freq)(transm * fSky * (freq/freq0)**(2*gamma))
	#
	# Return the result as a 2-element list of K-correction and flux at freq0
	return (kWave, fSky0)
# ====================================================================================
#
import java
import herschel
from herschel.share.fltdyn.time import FineTime
metaDict = herschel.spire.ia.util.MetaDataDictionary.getInstance()
fitsWriter = herschel.ia.io.fits.FitsArchive()
fitsWriter.rules.append(metaDict.getFitsDictionary())
from herschel.share.unit import *
c = Constant.SPEED_OF_LIGHT.value

scriptVersionString = "makeSCalPhotFluxConv.py $Revision: 1.17 $"

#directory = "..//..//..//..//..//..//data//spire//cal//SCal"
#dataDir = "//disks//winchester2//calibration_data//"
#LOCAL VERSION
directory = Configuration.getProperty('var.hcss.workdir')
dataDir = Configuration.getProperty('var.hcss.workdir')

# define the meta data
df        = java.text.SimpleDateFormat("yyyy.MM.dd/HH:mm:ss/z")

startDates = {}
endDates   = {}
versions   = {}
fileOrigin = {}
inputFiles = {}
### NOMINAL ###
#                        start of PFM1 campaign               OD6-98                               OD98 - 
startDates['nominal'] = [df.parse("2005.02.22/00:00:00/GMT"), df.parse("2009.05.19/10:59:00/GMT"), df.parse("2009.08.19/14:43:16/GMT")]
endDates['nominal']   = [df.parse("2009.05.19/10:59:00/GMT"), df.parse("2009.08.19/14:43:16/GMT"), df.parse("2020.01.01/00:00:00/GMT")]
versions['nominal']   = ["6"                                , "5"                                , "10"]
inputFiles['nominal'] = ["SCalPhotFluxConv_v1.8.fits", "ScalPhotFluxConv_Nominal_V2.0.fits", "ScalPhotFluxConv_Nominal_V2.3.fits"]
fileOrigin['nominal'] = ["v1.8"                      , "v2.0"                              , "GJB_26Feb13"]
    
### BRIGHT ###
#                       start of PFM1 campaign -116          OD116 - 
startDates['bright'] = [df.parse("2005.02.22/00:00:00/GMT"), df.parse("2009.09.06/13:35:35/GMT")]
endDates['bright']   = [df.parse("2009.09.06/13:35:35/GMT"), df.parse("2020.01.01/00:00:00/GMT")]
versions['bright']   = ["5"                                , "6"  ]
inputFiles['bright']  = ["SCalPhotFluxConv_v1.8_bsmode.fits", "ScalPhotFluxConv_BrightSource_V2.3.fits"]
fileOrigin['bright']  = ["v1.8"                             , "GJB_26Feb13"]

photRespControlStamp = open(dataDir+"PhotRespControlStamp.dat", "r")
respControlStamp = FineTime(Long(photRespControlStamp.read()))
photRespControlStamp.close()

################################################################################
# Calculate the K parameters for the metadata (following Luca & Bernhard's script)
spireRefFreq = c/Double1d([250.,350.,500.])*1e6 
gamma = -0.85 #
deltaNu = 0.1e9		# 0.1 GHz
nuMin   = 150.e9#
nuMax   = 1500.e9
nNu     = FIX((nuMax-nuMin)/deltaNu)
freq    = Double1d(range(nNu)) * deltaNu + nuMin
# Photometer RSRF
rsrfVersion = "3" 
rsrf = fitsReader("%s//Phot//SCalPhotRsrf//SCalPhotRsrf_v%s.fits"%(directory, rsrfVersion))
spireFreq   = rsrf['rsrf']['frequency'].data*1e9  # Frequency in Hz
ix = freq.where((freq>=MIN(spireFreq)) & (freq<=MAX(spireFreq)))
# Get RSRF for normal point sources
# Interpolate to common frequency grid
interpPLW = CubicSplineInterpolator(spireFreq, rsrf['rsrf']['plw'].data)
interpPMW = CubicSplineInterpolator(spireFreq, rsrf['rsrf']['pmw'].data)
interpPSW = CubicSplineInterpolator(spireFreq, rsrf['rsrf']['psw'].data)
spireFiltPLW = Double1d(nNu)
spireFiltPMW = Double1d(nNu)
spireFiltPSW = Double1d(nNu)
spireFiltPLW[ix] = interpPLW(freq[ix])
spireFiltPMW[ix] = interpPMW(freq[ix])
spireFiltPSW[ix] = interpPSW(freq[ix])
# Aperture efficiency table
apertureEfficiencyVersion = "1"
apertureEfficiency = fitsReader("%s//Phot//SCalPhotApertureEfficiency//SCalPhotApertureEfficiency_v%s.fits"%(directory, apertureEfficiencyVersion))
spireApEffFreq = apertureEfficiency['frequency']['frequency'].data * 1e9 #comes in [GHz]
spireApEffPsw  = apertureEfficiency['frequency']["PSW"].data
spireApEffPmw  = apertureEfficiency['frequency']["PMW"].data
spireApEffPlw  = apertureEfficiency['frequency']["PLW"].data
# Fold in and interpolate aperture efficiency
ix = freq.where((freq>=MIN(spireApEffFreq)) & (freq<=MAX(spireApEffFreq)))
interpPLW = CubicSplineInterpolator(spireApEffFreq, spireApEffPlw)
interpPMW = CubicSplineInterpolator(spireApEffFreq, spireApEffPmw)
interpPSW = CubicSplineInterpolator(spireApEffFreq, spireApEffPsw)
# Also for point source filter profiles
spireFiltPLW[ix] = interpPLW(freq[ix]) * spireFiltPLW[ix]
spireFiltPMW[ix] = interpPMW(freq[ix]) * spireFiltPMW[ix]
spireFiltPSW[ix] = interpPSW(freq[ix]) * spireFiltPSW[ix]
# Calculate K-correction factors for extended source assuming alpha=-1
#k4E = {"PSW": hpXcalKcorr(spireRefFreq[0], freq, spireFiltPSW, False, gamma=gamma)[0],\
#       "PMW": hpXcalKcorr(spireRefFreq[1], freq, spireFiltPMW, False, gamma=gamma)[0],\
#       "PLW": hpXcalKcorr(spireRefFreq[2], freq, spireFiltPLW, False, gamma=gamma)[0]}
# Calculate K-correction factors for point source assuming alpha=-1
#k4P = {"PSW": hpXcalKcorr(spireRefFreq[0], freq, spireFiltPSW, False)[0],\
#       "PMW": hpXcalKcorr(spireRefFreq[1], freq, spireFiltPMW, False)[0],\
#       "PLW": hpXcalKcorr(spireRefFreq[2], freq, spireFiltPLW, False)[0]}

#k4E = {"PSW":1.0087219700386618,\
#       "PMW":1.01099639794548,\
#       "PLW":1.0049005210207904}
#k4P = {"PSW":1.01015837,\
#       "PMW":1.00947124,\
#       "PLW":1.005577}

# Read in the K4 parameters generated by the colorCorrK script.
colorCorrKversion = "v4"
file = open("%s/colourCorrectionK4Parameters.txt"%(dataDir),"r")
lines = file.readlines()
print lines
colorCorrKread = (lines[0].split()[3])
print colorCorrKversion, colorCorrKread
if colorCorrKversion == colorCorrKread:
   k4E = {"PSW":float(lines[2].split()[1]),\
       "PMW":float(lines[2].split()[2]),\
       "PLW":float(lines[2].split()[3])}
   k4P = {"PSW":float(lines[3].split()[1]),\
       "PMW":float(lines[3].split()[2]),\
       "PLW":float(lines[3].split()[3])}
else:
   raise Exception("\n\nCould not read in K4 parameters... version mismatch!!!\n")

file.close()

print "Extended K4 parameters = %s"%k4E
print "Point source K4 parameters = %s"%k4P
################################################################################


# ------------ Read in the file -------------------------------------------------------------
biasEditions = ['nominal', 'bright']

for biasEdition in biasEditions:
    for i in range(len(versions[biasEdition])):
        fluxConvTable = fitsWriter.load(dataDir+inputFiles[biasEdition][i])

        # -------------------------------------------------------------------------------------------

        fluxConvTable.meta["creator"].value   = scriptVersionString
        # add metadata descriptions where they do not exist:
        fluxConvTable.meta["creator"].description = "Generator of this product"
        fluxConvTable.meta["modelName"].value = "FM"
        fluxConvTable.meta["modelName"].description = "Model name attached to this product"
        fluxConvTable.meta["startDate"].value = FineTime(startDates[biasEdition][i])
        fluxConvTable.meta["startDate"].description = "Start date of this product"
        fluxConvTable.meta["endDate"].value   = FineTime(endDates[biasEdition][i])
        fluxConvTable.meta["endDate"].description = "End date of this product"
        # set the creation date (as modifying original product).
        fluxConvTable.meta["creationDate"].value   = FineTime(java.util.Date())
        fluxConvTable.meta["description"].value = "Photometer Flux Conversion Calibration Table"
        fluxConvTable.meta["dependency"]= herschel.ia.dataset.StringParameter(value="biasMode, time",description="Keywords on which product depends")
        fluxConvTable.meta["biasMode"]  = herschel.ia.dataset.StringParameter(value=biasEdition, description="Nominal/bright source mode")
        if fileOrigin[biasEdition][i] == "GJB_26Feb13":
            fluxConvTable.meta["author"]  = herschel.ia.dataset.StringParameter(value="George Bendo", description="Author of the data")
            fluxConvTable.meta["sourceTable"].value  = "fluxcal_feb2013.txt"
            #fluxConvTable.meta["respControlStamp"].value  = FineTime(java.util.Date())
            # manually fix the control stamp to correct bug and match the tempDrift product (SPCAL-11)
            #tempDriftCorr = fitsReader(directory+"//Phot//SCalPhotTempDriftCorr//SCalPhotTempDriftCorr_nominal_20090819_v5.fits")
            fluxConvTable.meta["respControlStamp"].value  = respControlStamp
        else:
            fluxConvTable.meta["author"]  = herschel.ia.dataset.StringParameter(value="Nanyao Lu", description="Author of the data")        
        fluxConvTable.meta["dataOrigin"]  = herschel.ia.dataset.StringParameter(value="%s"%fileOrigin[biasEdition][i], description="Origin of the data")
        fluxConvTable.setVersion(versions[biasEdition][i])

        # fill the file with George's in flight values
        #if fileOrigin[biasEdition][i] == "GJB_16Jul10":
        #if fileOrigin[biasEdition][i] == "GJB_25Jan11":
        if fileOrigin[biasEdition][i] == "GJB_26Feb13":
            for array in ["PSW", "PMW", "PLW"]:
                #file  = open("PhotFluxConv_results_nominal_16July2010.txt","r")
                #file  = open("PhotFluxConv_results_nominal_25Jan11.txt","r")
                if biasEdition == "nominal":
                    file  = asciiTableReader(file=dataDir+"fluxcal_%s_feb2013.csv"%array.lower(), tableType='CSV')
                if biasEdition == "bright":
                    file  = asciiTableReader(file=dataDir+"fluxcal_%sbright_feb2013.csv"%array.lower(), tableType='CSV')
                for j in range(len(file['names'].data)):
                    names = fluxConvTable[array]['names'].data
                    sel = names.where(names.eq(file['names'].data[j]))
                    fluxConvTable[array]['v0'].data[sel]      = file['v0'].data[j]
                    fluxConvTable[array]['k1'].data[sel]      = file['k1'].data[j]
                    fluxConvTable[array]['k1Error'].data[sel] = ABS(file['k1Error'].data[j])
                    fluxConvTable[array]['k2'].data[sel]      = file['k2'].data[j]
                    fluxConvTable[array]['k2Error'].data[sel] = ABS(file['k2Error'].data[j])
                    fluxConvTable[array]['k3'].data[sel]      = file['k3'].data[j]
                    fluxConvTable[array]['k3Error'].data[sel] = ABS(file['k3Error'].data[j])
                    fluxConvTable[array]['vMin'].data[sel]    = file['vMin'].data[j]
                    fluxConvTable[array]['vMax'].data[sel]    = file['vMax'].data[j]
                    #print values[0], sel.length(), fluxConvTable[values[0][:3]]['v0'].data[sel]
        # Add K4 parameters to metadata for all editions    
        for array in ["PSW", "PMW", "PLW"]:
            fluxConvTable[array].meta["k4E_%s"%array] = DoubleParameter(k4E[array], "%s colour correction for standard nu*F_nu=const. reference SED of an extended source."%array)
            fluxConvTable[array].meta["k4P_%s"%array] = DoubleParameter(k4P[array], "%s colour correction for standard nu*F_nu=const. reference SED of a point source."%array)
            # And also add to main metadata
            fluxConvTable.meta["k4E_%s"%array] = DoubleParameter(k4E[array], "%s colour correction for standard nu*F_nu=const. reference SED of an extended source."%array)
            fluxConvTable.meta["k4P_%s"%array] = DoubleParameter(k4P[array], "%s colour correction for standard nu*F_nu=const. reference SED of a point source."%array)



        # ******* write to FITS *******
        df2 = java.text.SimpleDateFormat("yyyyMMdd")
        filename = java.io.File(r"%s//Phot//SCalPhotFluxConv//SCalPhotFluxConv_%s_%s_v%s.fits"%(directory, biasEdition, df2.format(startDates[biasEdition][i]),versions[biasEdition][i]))
        print
        print "writing to FITS...."
        fitsWriter.save(filename.toString(), fluxConvTable)
        print
        print "written: %s"%filename.toString()
        print


