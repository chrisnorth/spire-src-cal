#
#  This file is part of Herschel Common Science System (HCSS).
#  Copyright 2001-2014 Herschel Science Ground Segment Consortium
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

##Example script for HIPE source calibration scripts

#-------------------------------------------------------------------------------
#===============================================================================
#=====                   POINT AND FULLY-EXTENDED SOURCES                  =====
#===============================================================================
#-------------------------------------------------------------------------------

#
##Need to make sure SpireHandbookBundle_dev.py is in HIPE path
##Import it as a module
import SpireHandbookBundle_dev as hb
## An alternative is to import all the functions (and global variables directly)
#from SpireHandbookBundle_dev import *
## or just opened the script and execute in HIPE (or use to execfile())
## However, the examples below assume it has been import as hb

## If required, specify a calibration version (default is to use spire_cal_12_2)
## 1) Can be read from pool
#hb.getCal(calPool='spire_cal_12_2')
## 2) or from HSA
#hb.getCal(calTree='spire_cal_12_2')
## 3) or from jarfile
#hb.getCal(calFile='spire_cal_12_2')
## 4) or read a calibration and check
#cal=spireCal()
#hp.getCal(cal=cal)

## If no calibration specified, default is used. Equivalent to
#hb.getCal()

## The calibration is then used for all the processes in the module

## By default, the module uses the full beam treatment to calculate the 
## monochromatic beam areas, and since this is relatively slow, only calculates
## them once. This is done as required, but can be achieved manually by
hb.calcBeamMonoArea()
## To use a simple beam treatment, which scales the Neptune areas by nu^2*gamma,
## set beamType='Simple'
#hb.calcBeamMonoArea(beamType='Simple')

## The monochromatic beam areas can be recalculated by re-running the
## calcBeamMonoArea() function with beamType='Simple' or 'Full'

## 

## Calculate pipeline parameters
print 'Testing pipeline parameters'
print 'K4P=',hb.calcK4P()
print 'KMonE=',hb.calcKMonE()
print 'K4E=',hb.calcK4E()
print 'KPtoE=',hb.calcKPtoE()
print 'Omega(alpha=-1)=',hb.calcOmegaEff(-1.0)

## Set range of power laws to use in thie example
alphaArr=Float1d(range(-4,5)) #range of alphas to use.

## Calculate effective beams
print '\nTesting hb.calcOmegaEff'
# 1) single power law spectrum
print 'OmegaEff(-2)=',hb.calcOmegaEff(-2.0)
# 2) multiple power law spectrum
omegaEff = hb.calcOmegaEff(alphaArr)
# 3) single modified black body spectrum
print 'OmegaEff_BB(1.75,20)=',hb.calcOmegaEff_BB(1.75,20.0)
# 4) multiple modified black body spectrum
omegaEff_BB = hb.calcOmegaEff_BB(1.75,[20.0,30.,40.])

## Calculate beam corrections
print '\nTesting hb.calcKBeam'
# 1) single power law spectrum
print 'KBeam(-2)=',hb.calcKBeam(-2.0)
# 2) multiple power law spectrum
KBeam = hb.calcKBeam(alphaArr)
# 3) single modified black body spectrum
print 'KBeam(1.75,20)=',hb.calcKBeam_BB(1.75,20.0)
# 4) multiple modified black body spectrum
KBeam_BB = hb.calcKBeam_BB(1.75,[20.0,30.,40.])

## Calculate point source colour corrections
print '\nTesting hb.calcKColP'
# 1) single power law spectrum
print 'KColP(-2)=',hb.calcKColP(-2.0)
# 2) multiple power law spectrum
KColP = hb.calcKColP(alphaArr)
# 3) single modified black body spectrum
print 'KColP(1.75,20)=',hb.calcKColP_BB(1.75,20.0)
# 4) multiple modified black body spectrum
KColP_BB = hb.calcKColP_BB(1.75,[20.0,30.,40.])

## Calculate fully-extended source colour corrections
print '\nTesting hb.calcKColE'
# 1) single power law spectrum
print 'KColE(-2)=',hb.calcKColE(-2.0)
# 2) multiple power law spectrum
KColE = hb.calcKColE(alphaArr)
# 3) single modified black body spectrum
print 'KColE(1.75,20)=',hb.calcKColE_BB(1.75,20.0)
# 4) multiple modified black body spectrum
KColE_BB = hb.calcKColE_BB(1.75,[20.0,30.,40.])

## Calculate aperture corrections (can be slow)
print '\nTesting aperture corrections'
# 1) single power law spectrum
apCorr=hb.calcApCorr(-2.0,verbose=True)
print 'Ap Corr [incBG] (-2)=',apCorr[0]
print 'Ap Corr [noBG] (-2)=',apCorr[1]
# 2) multiple power law spectrum
#apCorr = hb.calcApCorr(alphaArr,verbose=True)
#apCorrIncBG=apCorr[0]
#apCorrNoBG=apCorr[1]
# 3) single modified black body spectrum
#apCorr_BB = hb.calcApCorr_BB(1.75,20.0,verbose=True)
#print 'Ap Corr [incBG] (1.75,20)=',apCorr_BB[0]
#print 'Ap Corr [noBG] (1.75,20)=',apCorr_BB[1]
# 4) multiple modified black body spectrum
#apCorr_BB = hb.calcApCorr_BB(1.75,[20.0,30.,40.],verbose=True)
#apCorrIncBG_BB=apCorr_BB[0]
#apCorrNoBG_BB=apCorr_BB[1]

#-------------------------------------------------------------------------------
#===============================================================================
#=====                          INTERNAL VARIABLES                         =====
#===============================================================================
#-------------------------------------------------------------------------------
## As part of its processing the module uses some global variables
## If functions are imported directly (i.e. using "from <module> import *"
## then these are directly available. Otherwise they are accessed similarly to
## accessing contents of a class

## Using the access methods below, the variables will be returned, or if not
## present, calculated from scratch

## If they need to be recalculated, the global variable must be deleted, and the
## function re-run.

## The photometer calibration tree is stored in the spireCalPhot global variable
calPhot=hb.spireCalPhot
## it can be accessed/calculated via the getCal() function
calPhot=hb.getCal()

## The frequency raster is stored in spireFreq global variable
freq=hb.spireFreq
## it can be accessed/calculated via the getSpireFreq() function, which returns
## a single array
freq=hb.getSpireFreq()

## The reference frequencies for the SPIRE bands are stored in spireRefFreq
refFreq=hb.spireRefFreq
## they are accessed/calculated via the getSpireRefFreq() function, which returns
## a dict with keys "PSW", "PMW", "PLW"
refFreq=hb.getSpireRefFreq()

## The effective frequencies for the SPIRE beams are stored in spireEffFreq
effFreq=hb.spireEffFreq
## they are accessed/calculated via the getSpireEffFreq() function, which returns
## a dict with keys "PSW", "PMW", "PLW"
effFreq=hb.getSpireEffFreq()

## The RSRF filter profiles for the SPIRE bands are stored in spireFiltOnly
rsrf=hb.spireFiltOnly
## The same profiles multiplied by the aperture efficiency are spireFilt
## they are accessed/calculated via the getSpireFilt() function, which returns
## a dict with keys "PSW", "PMW", "PLW"
rsrf=hb.getSpireFilt(rsrfOnly=True)
rsrf_apEff=hb.getSpireFilt(rsrfOnly=False)

## The monochromatic beam areas (which take time to calculate) are stored in
## the beamMonoArea global variable
beamMonoArea=hb.beamMonoArea
## They area accessed/calculated via the calcBeamMonoArea() function, which returns
## a dict with keys "PSW", "PMW", "PLW", and a "beamType" key with value "Simple"
## or "Full"
beamMonoArea=hb.calcBeamMonoArea()
## By default, this will return the previously calculated values, or if none
## exist will use the full beam treatment.
## Setting the beamType parameter will force the return (and if necessary, the
## calculation of the 'Simple' or 'Full' monochromatic beam areas
beamMonoAreaSimple=hb.calcBeamMonoArea(beamType='Simple')
beamMonoAreaSimple=hb.calcBeamMonoArea(beamType='Full')
## Note that doing so will change the values used by any subsequent calculations
## which use the monochromatic beam areas