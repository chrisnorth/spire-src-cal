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
#===============================================================================
# 
#  Herschel-SPIRE Semi-Extended source Colour Correction factors 
# 
#  This routine calculates colour correction tables for partially extented sources
#  to convert from the standard pipeline flux densities (which are quoted for a 
#  nu*F_nu=const. spectrum) into flux densities for a range of source spectra.
#  
#  These produce monochromatic peak surface brightnesses at the SPIRE reference 
#  wavelengths of 250, 350 and 500 micron. The source spectra include power law 
#  spectra, and modified black body spectra with a range of emissivities.
#
#  The script relies on functions in SpireHandbookBundle_dev.py
#
#  Input:
#    RSRF and Aperture efficiency profiles from SPIRE calibration tree
#    Beam profiles and beam model parameters from SPIRE calibration tree
#    Flux conversion parameters from SPIRE calibration tree
#    Name and version of output file
#
#  Output:
#    Functions which can be used to calculate colour correction parameters
#    Some of them also generate global variables
#
#  Usage:
#    A script is provided at the end for usage of the functions
#
#===============================================================================
# 
#  Herschel-SPIRE SemiExtended Source Colour Corrections
# 
#  This module provides functions to calculate the colour corrections for a
#  range of source spectra, for partially extended sources.
# 
#
#===============================================================================
# 
#  Edition History
#   Chris North      15-04-2014  - Initial version
#
#===============================================================================

#-------------------------------------------------------------------------------
#===============================================================================
#=====                      IMPORT HIPE & JAVA MODULES                     =====
#===============================================================================
#-------------------------------------------------------------------------------
import os
import herschel
from herschel.ia.numeric import Double1d,Float1d,Int1d
from herschel.spire.ia.cal import SpireCalTask
spireCal = SpireCalTask()
from herschel.ia.numeric.toolbox.interp import LinearInterpolator,CubicSplineInterpolator
from herschel.ia.numeric.toolbox.integr import TrapezoidalIntegrator
from java.lang.Math import PI
from java.lang import Double,Float,String
from herschel.share.unit import Constant
from herschel.ia.gui.plot import *
import java.awt.Color
from herschel.ia.numeric.toolbox.basic import Floor,Min,Max,Exp
FLOOR=herschel.ia.numeric.toolbox.basic.Floor.PROCEDURE
EXP=herschel.ia.numeric.toolbox.basic.Exp.PROCEDURE
MAX=herschel.ia.numeric.toolbox.basic.Max.FOLDR
MIN=herschel.ia.numeric.toolbox.basic.Min.FOLDR

from SpireHandbookBundle_dev import *
from beamfunctions import spireMonoBeamSrc,spireMonoSrcAreas
import sources_dev as srcMod

#scriptVersionString = "SemiExtendedBundle.py $Revision: 1.0 $"
    
#-------------------------------------------------------------------------------
#===============================================================================
#=====                  CALCULATE MONOCHROMATIC BEAM AREAS                 =====
#===============================================================================
#-------------------------------------------------------------------------------


def calcBeamSrcMonoArea(srcKey,verbose=False):
    #calculate monochromatic beam areas using full or simple beam treatment
    #print '\nCalculating monochromatic beam areas...'

    global arcsec2Sr
    global beamMonoSrcArea

    try:
        arcsec2Sr
    except:
        #define arcsec2Sr
        arcsec2Sr = (PI/(60.*60.*180))**2
    
    try:
        beamMonoSrcArea[srcKey]
        if (verbose):print'Using existing monochromatic areas (key %s)'%srcKey
        #exists, don't recalculate
    except:
        #doesn't exist, recalculate
        #check if main structure exists
        try:
            beamMonoSrcArea
            #exists, don't create
        except:
            #doesn't exist, must create
            beamMonoSrcArea={}
        beamProfs = spireCalPhot.getProduct("RadialCorrBeam")
        beamRad=beamProfs.getCoreCorrectionTable().getColumn('radius').data
        #define srcProf to be same length as beamRad
    
        srcProf=srcMod.src2Prof(srcKey,beamRad,verbose=True)
        
        gamma = beamProfs.meta['gamma'].double
        beamMonoSrcArea[srcKey] = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
        for band in spireBands:
            if (verbose):print 'Calculating monochromatic beam areas for %s source for %s band'%(srcKey,band)
            #monochromatic beam areas
            beamMonoSrcArea[srcKey][band] = spireMonoSrcAreas(getSpireFreq(), beamProfs, 
              getSpireEffFreq()[band], gamma, srcProf, band)
        if (verbose):print 'All beam areas calculated for %s source'%srcKey
    return beamMonoSrcArea[srcKey]

def calcKColESrc(alphaK,srcKey,verbose=False,table=False):
    #-----------------------------------------------------------------------
    #print '\nCalculating extended source colour correction parameters over alpha...'
    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #freq=getSpireFreq()
    #spireFilt=getSpireFilt(cal)

    k4E_Tot=calcKMonE()

    try:
        na=len(alphaK)
        aList=True
    except:
        aList=False

    if not aList:
        # alphaK is scalar
        kColESrc = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
        #kBeamK = calcKBeam(alphaK)
        for band in spireBands:
            k4EaTot_x=calcSpireKcorr(getSpireRefFreq()[band], getSpireFreq(),\
             getSpireFilt()[band], BB=False, alpha=alphaK,\
             ext=True, monoArea=calcBeamSrcMonoArea(srcKey,verbose=verbose)[band])[0]/1.e6
            kColESrc[band] = k4EaTot_x / k4E_Tot[band]
        if (verbose): print 'Calculated KColE for alpha=%f and source Width: '%alphaK,kColESrc
    else:
        # alphaK is list
        kColESrc = {'PSW': Double1d(na,Double.NaN), 'PMW': Double1d(na,Double.NaN), 'PLW': Double1d(na,Double.NaN)}
        for a in range(na):
            for band in spireBands:
                k4EaTot_x=calcSpireKcorr(getSpireRefFreq()[band], getSpireFreq(),\
                 getSpireFilt()[band], BB=False, alpha=alphaK[a],\
                 ext=True, monoArea=calcBeamSrcMonoArea(srcKey,verbose=verbose)[band])[0]/1.e6
                kColESrc[band][a] = k4EaTot_x / k4E_Tot[band]
            if (verbose): print 'Calculated KColE for alpha=%f: '%alphaK[a],kColESrc["PSW"][a],kColESrc["PMW"][a],kColESrc["PLW"][a]

    if not table:
        #return as is
        return kColESrc
    else:
        #create and returnTableDataset
        kColE_table=TableDataset()
        kColE_table.setDescription("Partially Extended Source Colour Correction (Spectral Index)")
        kColE_table.addColumn("alpha",Column(Double1d(alphaK)))
        for band in spireBands:
            kColE_table.addColumn(band,Column(kColESrc[band]))
        return kColE_table

def semiExtendedTest():

    #make Gaussian source
    alphaK=2.0
    srcWidths=[1.,10.,100.,1000]
    KColE_full=calcKColE(alphaK,verbose=True)
    KColE_partial=[]
    for wid in srcWidths:
        srcKey=srcMod.setSrc('Gaussian',[wid])
        KColE_partial.append(calcKColESrc(alphaK,srcKey,verbose=True,table=False))