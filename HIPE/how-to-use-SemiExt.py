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
##Need to make sure SpireHandbookBundle.py is in HIPE path
#specify directory where script are stored and add to path
directory  = Configuration.getProperty('var.hcss.workdir')
import sys
sys.path.append(directory)
##Import Handbook module
import SemiExtendedBundle as semi
import SpireHandbookBundle as hb
import sources

#-------------------------------------------------------------------------------
#===============================================================================
#=====                  SPECIFY CALIBRATION TREE (OPTIONAL)                =====
#===============================================================================
#-------------------------------------------------------------------------------

## OPTIONAL: Specify a specific calibration version or calibration tree.
## Default is to use the version obtained by running spireCal()
## 1) Can be read from pool
#semio.getCal(calPool='spire_cal_12_3')

## 2) or from HSA
#semi.getCal(calTree='spire_cal_12_3')

## 3) or from jar file
#semi.getCal(calFile='spire_cal_12_3')

## 4) or read a calibration, make modifications, use that
#cal_mod=spireCal()
#<make modifications>
#semi.getCal(cal=cal_mod)

## The calibration is then used for all the processes in the module

#-------------------------------------------------------------------------------
#===============================================================================
#=====                          IMPORT BEAM PROFILES                       =====
#===============================================================================
#-------------------------------------------------------------------------------

#get calibration tree
calPhot=semi.getCal()
beams=calPhot.getProduct('RadialCorrBeam').getCoreCorrectionTable()
beamProfs={}
for band in ['PSW','PMW','PLW']:
    beamProfs[band]=sources.SourceProfile(radArr=beams['radius'].data,\
      profile=beams[band].data)

#make maps of beams
beamMaps={}
for band in ['PSW','PMW','PLW']:
    beamMaps[band]=beamProfs[band].makeImage()

#-------------------------------------------------------------------------------
#===============================================================================
#=====                        GENERATE SOURCE PROFILE                      =====
#===============================================================================
#-------------------------------------------------------------------------------

#List source types
sources.Source().listSrcTypes()

#create Gaussian source with width (1sigma) of 50. arcsec)
srcGauss=sources.Source('GAUSSIAN',50.)
#print FWHM
print 'Gaussian source FWHM = %f'%srcGauss.calcFwhm()

#create TopHat source with width of 50. arcsec)
srcTophat=sources.Source('TOPHAT',50.)
#print FWHM
print 'Tophat source FWHM = %f'%srcTophat.calcFwhm()

#generate SourceProfile object (using beam radArr)
srcGaussProf=srcGauss.calcProfile(radArr=beamProfs['PSW'].radArr)
srcTophatProf=srcTophat.calcProfile(radArr=beamProfs['PSW'].radArr)

#-------------------------------------------------------------------------------
#===============================================================================
#=====                      CALCULATE COLOUR CORRECTIONS                   =====
#===============================================================================
#-------------------------------------------------------------------------------
#calculate full-extended correction
alpha=2.0
kColEFull=hb.calcKColE(alpha,verbose=False)
print 'KColE (Full)=',kColEFull
KColEPartialGaussian=semi.calcKColESrc(alpha,srcGauss,verbose=False)
print 'KColE (Gaussian)=',KColEPartialGaussian
KColEPartialTophat=semi.calcKColESrc(alpha,srcTophatProf,verbose=False)
print 'KColE (Top Hat)=',KColEPartialTophat

#-------------------------------------------------------------------------------
#===============================================================================
#=====                  CONVOLVE BEAM SND SOURCE PROFILE                   =====
#===============================================================================
#-------------------------------------------------------------------------------

#convolve beam with Gaussian source
srcGaussConv={}
for band in ['PSW','PMW','PLW']:
    srcGaussConv[band]=sources.convolveProfiles(beamProfs[band],srcGaussProf,verbose=True)

#convolve beam with Tophat source
srcTophatConv={}
for band in ['PSW','PMW','PLW']:
    srcTophatConv[band]=sources.convolveProfiles(beamProfs[band],srcTophatProf,verbose=True)

#make maps of convolved sources
srcGaussConvMaps={}
srcTophatConvMaps={}
for band in ['PSW','PMW','PLW']:
    srcGaussConvMaps[band]=srcGaussConv[band].makeImage()
    srcTophatConvMaps[band]=srcTophatConv[band].makeImage()
    
#-------------------------------------------------------------------------------
#===============================================================================
#=====                  CONVOLVE BEAM SND SOURCE PROFILE                   =====
#===============================================================================
#-------------------------------------------------------------------------------

