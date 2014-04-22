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
from java.lang import Double,Float,String,Integer
from herschel.share.unit import Constant
from herschel.ia.gui.plot import *
import java.awt.Color
from herschel.ia.numeric.toolbox.basic import Floor,Min,Max,Exp
FLOOR=herschel.ia.numeric.toolbox.basic.Floor.PROCEDURE
EXP=herschel.ia.numeric.toolbox.basic.Exp.PROCEDURE
MAX=herschel.ia.numeric.toolbox.basic.Max.FOLDR
MIN=herschel.ia.numeric.toolbox.basic.Min.FOLDR

import SpireHandbookBundle_dev as sh
import beamfunctions as beam
import sources_dev as srcMod
from SpireHandbookBundle_dev import spireBands
#scriptVersionString = "SemiExtendedBundle.py $Revision: 1.0 $"
    
#-------------------------------------------------------------------------------
#===============================================================================
#=====                  CALCULATE MONOCHROMATIC BEAM AREAS                 =====
#===============================================================================
#-------------------------------------------------------------------------------


def calcBeamSrcMonoArea(src,verbose=False):
    """
    ========================================================================
    calcBeamSrcMonoArea(src,verbose=False):
        Calculate monochromatic beam areas over frequency raster, convolved with
          source profile
        Stored as entry in global variable beamMonoSrcArea, so only calculated
          once per source

    Inputs:
      src:     (Source) Source object to multiply beam profile by
      verbose: (boolean) print more information. Default=False

    Outputs:     
      (dict) Monochromatic beam areas over frequency raster.
        [one float array per band, plus beamType='full'|'simple']
        
    Calculation:
      + If entry for Source alread exists in global variable beamMonoArea,
          returns previosly calculated version
      + If not:
          - Retrieves photometer calibration product (using getCal())
          - Extracts radial beam profile
          - Calculates profile for Source object
          + If scaleWidth of source is <1 step of beam profile radius array
            - ignores beam profile, just calculates source area
          + If not:
            + If scaleWidth of source is <10 steps of beam profile radius array
              - increase sampling for beam radius array
            - Calculates area of beam profile multipled by source profile

    Dependencies:
        sources_dev.Source - Source object class
        sources_dev.SourceProfile - Source profile object class
        SpireHandbookBundle_dev.spireBands() - get list of spire bands
        SpireHandbookBundle_dev.getCal() - retrieve photometer calibration product
        SpireHandbookBundle_dev.getSpireFreq() - retrieve frequency raster
        SpireHandbookBundle_dev.getSpireEffFreq() - retrieve effective frequencies
        fineBeam() - increase sampling of radius of RadialCorrBeam calibration product
        beamfunctions.spireMonoSrcAreas() - calculate monochromatic beam areas multiplied by source
    """
    #calculate monochromatic beam areas using full or simple beam treatment
    #print '\nCalculating monochromatic beam areas...'

    global beamMonoSrcArea

    try:
        beamMonoSrcArea[src.key]
        if (verbose):print'Using existing monochromatic areas (key %s)'%src.key
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
                
        beamProfs = sh.getCal().getProduct("RadialCorrBeam")
        beamRad=beamProfs.getCoreCorrectionTable().getColumn('radius').data
        #define srcProf to be same length as beamRad
        srcProf=src.calcProfile(beamRad)
        if src.scaleWidth<beamRad[10]:
            #quite a small source (< 10 steps in radius array)
            beamProfsFine=fineBeam(beamProfs,src.scaleWidth)
            beamRadFine=beamProfsFine.getCoreCorrectionTable().getColumn('radius').data
            #define srcProfFine to be same length as beamRadFine
            srcProfFine=src.calcProfile(beamRadFine)
            srcAreaFine=srcProfFine.calcArea()*sh.arcsec2Sr()
            beamMonoSrcArea[src.key] = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
            if src.scaleWidth<beamRad[1]:
                #VERY small (scale width < 1 step in radius array)
                if (verbose):print 'Very small source [%g arcsec]. Using source area only (ignoring beam)'%src.scaleWidth
                for band in spireBands():
                    spireFreq=sh.getSpireFreq()
                    beamMonoSrcArea[src.key][band]=Double1d(len(spireFreq),srcAreaFine)
            else:
                #quote small source (scale width 1-10 steps in radius array
                for band in spireBands():
                    if (verbose):print 'Calculating monochromatic beam areas for small [%g arcsec] %s source for %s band'%(src.scaleWidth,src.key,band)
                    #monochromatic beam areas
                    gamma = beamProfsFine.meta['gamma'].double
                    beamMonoSrcArea[src.key][band] = beam.spireMonoSrcAreas(sh.getSpireFreq(), beamProfsFine, 
                      sh.getSpireEffFreq()[band], gamma, srcProfFine, band)
                
        else:
            #calculate area of beam and source
            gamma = beamProfs.meta['gamma'].double
            beamMonoSrcArea[src.key] = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
            for band in spireBands():
                if (verbose):print 'Calculating monochromatic beam areas for %s source for %s band'%(src.key,band)
                #monochromatic beam areas
                beamMonoSrcArea[src.key][band] = beam.spireMonoSrcAreas(sh.getSpireFreq(), beamProfs, 
                  sh.getSpireEffFreq()[band], gamma, srcProf, band)
            
        
        if (verbose):print 'All beam areas calculated for %s source'%src.key        
    return beamMonoSrcArea[src.key]

def fineBeam(beamProfsIn,scaleWidth,verbose=False):
    """
    ========================================================================
    fineBeam(beamProfsIn,scaleWidth,verbose=False):
        Makes copy of RadialBeamCorr product and interpolates to finer radius
          sampling for inner region

    Inputs:
      beamProfsIn: (PhotRadialCorrBeam) beamProfsIn
      scaleWidth:  (double) width to
      verbose:     (boolean) print more information. Default=False

    Outputs:     
      (dict) Monochromatic beam areas over frequency raster.
        [one float array per band, plus beamType='full'|'simple']
        
    Calculation:
        + Copies PhotRadialCorrBeam product
        + If scaleWidth of source is <=10 steps of beam profile radius array
          - Resamples inner 10 steps to scaleWidth/2
          - Interpolates all columns to finer radius array


    Dependencies:
        herschel.ia.numeric.toolbox.interp.CubicSplineInterpolator
    """
    #adjust sampling of fineBeam according to scaleWidth
    beamRadIn=beamProfsIn.getCoreCorrectionTable().getColumn('radius').data

    iFine=11
    if scaleWidth < beamRadIn[iFine]:
        print 'fine-sampling beam profiles'
        #source scale Width is smaller than 10 beamRad steps
        #add finely-sampled region of profile
        #get coarse part of array
        radArrCoarse=beamRadIn[iFine+1:]
        #generate fine part of array
        maxRadFine=beamRadIn[iFine]
        radFineStep=scaleWidth/2
        nRadFine=int(maxRadFine/radFineStep)
        radArrFine=Double1d().range(nRadFine)
        radArrFine=radArrFine*radFineStep
        
        #concatenate arrays
        beamRad=radArrFine
        beamRad.append(radArrCoarse)

        #interpolate beamProfs
        beamProfs=beamProfsIn.copy()
        beamProfs.getCoreCorrectionTable().getColumn('radius').data=beamRad
        beamProfs.getConstantCorrectionTable().getColumn('radius').data=beamRad
        beamProfs.getNormAreaCorrectionTable().getColumn('radius').data=beamRad
        for band in spireBands():
            #get input arrays
            beamCoreIn=beamProfsIn.getCoreCorrectionTable().getColumn(band).data
            beamConstIn=beamProfsIn.getConstantCorrectionTable().getColumn(band).data
            beamNormIn=beamProfsIn.getNormAreaCorrectionTable().getColumn(band).data
            #make interpolation objects
            beamCoreInterp=CubicSplineInterpolator(beamRadIn,beamCoreIn)
            beamConstInterp=CubicSplineInterpolator(beamRadIn,beamConstIn)
            beamNormInterp=CubicSplineInterpolator(beamRadIn,Double1d(beamNormIn))
            #do interpolation
            beamCoreNew=beamCoreInterp(beamRad)
            beamConstNew=beamConstInterp(beamRad)
            beamNormNew=beamNormInterp(beamRad)
            #replace profile arrays
            beamProfs.getCoreCorrectionTable().getColumn(band).data=beamCoreNew
            beamProfs.getConstantCorrectionTable().getColumn(band).data=beamConstNew
            beamProfs.getNormAreaCorrectionTable().getColumn(band).data=beamNormNew
    else:
        #get original RadialBeamCorr product
        beamProfs=beamProfsIn.copy()
        
    return(beamProfs)

def calcKColPSrc(alphaK,src,verbose=False,table=False):
    """
    ========================================================================
    calcKColPSrc(alphaK,src,verbose=False,table=False):
        Calculates flux density colour correction for source with power law
          spectrum and a given (frequency-independent) source profile

    Inputs:
        alphaK:  (scalar/float array) Spectral index(es) to use
        src:     (Source object class) Source object to use
        verbose: (boolean) Set to print more information. Default=True.
        table:   (boolean) Set to output TableDataset. Default=False.

    Outputs:
        If table==False:
            (dict) colour correction factor (one scalar/list per band,
              depending on whether alphaK is scalar/list)
        If table==True:
            (TableDataset) colour correction factor (one column per band,
              one row per alphaK value)
        
    Calculation:
        [See Handbook Eq 5.40-43]
        Calculates K4P parameter for point source and pipeline spectrum (alpha=-1)
        Calculates source area (analytically if possible, if not then by integrating profile)
        Calculates flux conversion to get peak surface brightness for
          partially-resolved source with power law spectrum and profile defined
          by Source object class
        Multiplies peak surface brightness by source area to get total flux
          density of source
        Divides by the k4P Parameter to get colour correction
        If table is set, returns a TableDataset
        
        KColPSrc converts from point source monochromatic pipeline flux density
          to total flux density of partially-resolved source.

    Dependencies:
        sources_dev.Source - Source object class
        sources_dev.SourceProfile - Source profile object class
        SpireHandbookBundle_dev.spireBands() - get list of spire bands
        SpireHandbookBundle_dev.getCal() - retrieve photometer calibration product
        SpireHandbookBundle_dev.getSpireFreq() - retrieve frequency raster
        SpireHandbookBundle_dev.getSpireRefFreq() - retrieve reference frequencies
        calcBeamSrcMonoArea() - calculate area of source multiplied by beam
        herschel.ia.dataset.TableDataset()
    """
    #-----------------------------------------------------------------------
    #print '\nCalculating extended source colour correction parameters over alpha...'
    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #freq=getSpireFreq()
    #spireFilt=getSpireFilt(cal)

    k4P=sh.calcK4P()
    print 'k4P:',k4P
    try:
        na=len(alphaK)
        aList=True
    except:
        aList=False

    #calculate source area analytically
    srcArea=src.calcArea()
    if Double.isNaN(srcArea):
        #make profile and integrate
        beamRad=sh.getCal().getProduct("RadialCorrBeam").getCoreCorrectionTable().getColumn('radius').data
        srcArea=src.calcProfile(radArr).calcArea()

    if not aList:
        # alphaK is scalar
        kColPSrc = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
        #kBeamK = calcKBeam(alphaK)
        for band in spireBands():
            #calculate peak surface brightness
            k4EaTot_x=sh.calcSpireKcorr(sh.getSpireRefFreq()[band], sh.getSpireFreq(),\
              sh.getSpireFilt()[band], BB=False, alpha=alphaK,\
              ext=True, monoArea=calcBeamSrcMonoArea(src,verbose=verbose)[band])[0]
            #multiply by source area to get total flux density of source
            k4Pa_x = k4EaTot_x * srcArea*sh.arcsec2Sr()
            #calculate colour correction
            kColPSrc[band] = k4Pa_x / k4P[band]
        if (verbose): print 'Calculated KColP for alpha=%f and source %s: '%(alphaK,src.key),kColPSrc
    else:
        # alphaK is list
        kColPSrc = {'PSW': Double1d(na,Double.NaN), 'PMW': Double1d(na,Double.NaN), 'PLW': Double1d(na,Double.NaN)}
        for a in range(na):
            for band in spireBands():
                k4EaTot_x=sh.calcSpireKcorr(sh.getSpireRefFreq()[band], sh.getSpireFreq(),\
                  sh.getSpireFilt()[band], BB=False, alpha=alphaK[a],\
                  ext=True, monoArea=calcBeamSrcMonoArea(src,verbose=verbose)[band])[0]
                #multiply by source area to get total flux density of source
                k4Pa_x = k4EaTot_x * srcArea*sh.arcsec2Sr()
                #calculate colour correction
                kColPSrc[band][a] = k4Pa_x / k4P[band]
            if (verbose): print 'Calculated KColPSrc for alpha=%f and source %s: '%(alphaK[a],src.key),kColPSrc["PSW"][a],kColPSrc["PMW"][a],kColPSrc["PLW"][a]

    if not table:
        #return as is
        return kColPSrc
    else:
        #create and returnTableDataset
        kColPSrc_table=TableDataset()
        kColPSrc_table.setDescription("Partially Extended Flux Density Colour Correction (Spectral Index)")
        kColPSrc_table.addColumn("alpha",Column(Double1d(alphaK)))
        for band in spireBands():
            kColPSrc_table.addColumn(band,Column(kColPSrc[band]))
        return kColP_table

def calcKColPSrc_BB(betaK,tempK,src,verbose=False,table=False):
    """
    ========================================================================
    calcKColPSrc_BB(betaK,tempK,src,verbose=False,table=False):
        Calculates flux density colour correction for source with modified
          black body spectrum and a given (frequency-independent) source profile

    Inputs:
        betaK:   (float) Emmisivity spectral index to use
        tempK:   (scalar/float array) Temperature(s) to use
        src:     (Source object class) Source object to use
        verbose: (boolean) Set to print more information. Default=True.
        table:   (boolean) Set to output TableDataset. Default=False.

    Outputs:
        If table==False:
            (dict) colour correction factor (one scalar/list per band,
              depending on whether tempK is scalar/list)
        If table==True:
            (TableDataset) colour correction factor (one column per band,
              one row per tempK value)
        
    Calculation:
        [See Handbook Eq 5.40-43]
        Calculates K4P parameter for point source and pipeline spectrum (alpha=-1)
        Calculates source area (analytically if possible, if not then by integrating profile)
        Calculates flux conversion to get peak surface brightness for
          partially-resolved source with power law spectrum and profile defined
          by Source object class
        Multiplies peak surface brightness by source area to get total flux
          density of source
        Divides by the k4P Parameter to get colour correction
        If table is set, returns a TableDataset
        
        KColPSrc converts from point source monochromatic pipeline flux density
          to total flux density of partially-resolved source.

    Dependencies:
        sources_dev.Source - Source object class
        sources_dev.SourceProfile - Source profile object class
        SpireHandbookBundle_dev.spireBands() - get list of spire bands
        SpireHandbookBundle_dev.getCal() - retrieve photometer calibration product
        SpireHandbookBundle_dev.getSpireFreq() - retrieve frequency raster
        SpireHandbookBundle_dev.getSpireRefFreq() - retrieve reference frequencies
        calcBeamSrcMonoArea() - calculate area of source multiplied by beam
        herschel.ia.dataset.TableDataset()
    """
    #-----------------------------------------------------------------------
    #print '\nCalculating extended source colour correction parameters over alpha...'
    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #freq=getSpireFreq()
    #spireFilt=getSpireFilt(cal)

    k4P=sh.calcK4P()
    print 'k4P:',k4P
    try:
        nt=len(tempK)
        tList=True
    except:
        tList=False

    #calculate source area analytically
    srcArea=src.calcArea()
    if Double.isNaN(srcArea):
        #make profile and integrate
        beamRad=sh.getCal().getProduct("RadialCorrBeam").getCoreCorrectionTable().getColumn('radius').data
        srcArea=src.calcProfile(radArr).calcArea()

    if not tList:
        # tempK is scalar
        kColPSrc_BB = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
        #kBeamK = calcKBeam(alphaK)
        for band in spireBands():
            #calculate peak surface brightness
            k4EaTot_x=sh.calcSpireKcorr(sh.getSpireRefFreq()[band], sh.getSpireFreq(),\
              sh.getSpireFilt()[band], BB=True, beta=betaK, temp=tempK,\
              ext=True, monoArea=calcBeamSrcMonoArea(src,verbose=verbose)[band])[0]
            #multiply by source area to get total flux density of source
            k4Pa_x = k4EaTot_x * srcArea*sh.arcsec2Sr()
            #calculate colour correction
            kColPSrc_BB[band] = k4Pa_x / k4P[band]
        if (verbose): print 'Calculated KColPSrc for modBB with T=%f K, beta=%f and source %s : '%(tempK,betaK,src.key),kColPSrc_BB
    else:
        # tempK is list
        kColPSrc_BB = {'PSW': Double1d(nt,Double.NaN), 'PMW': Double1d(nt,Double.NaN), 'PLW': Double1d(nt,Double.NaN)}
        for t in range(nt):
            for band in spireBands():
                k4EaTot_x=sh.calcSpireKcorr(sh.getSpireRefFreq()[band], sh.getSpireFreq(),\
                  sh.getSpireFilt()[band], BB=True, beta=betaK, temp=tempK[t],\
                  ext=True, monoArea=calcBeamSrcMonoArea(src,verbose=verbose)[band])[0]
                #multiply by source area to get total flux density of source
                k4Pa_x = k4EaTot_x * srcArea*sh.arcsec2Sr()
                #calculate colour correction
                kColPSrc_BB[band][t] = k4Pa_x / k4P[band]
            if (verbose): print 'Calculated KColP for modBB with T=%f K, beta=%f and source %s: '%(tempK[t],betaK,src.key),kColPSrc_BB["PSW"][t],kColPSrc_BB["PMW"][t],kColPSrc_BB["PLW"][t]

    if not table:
        #return as is
        return kColPSrc_BB
    else:
        #create and returnTableDataset
        kColPSrc_BB_table=TableDataset()
        kColPSrc_BB_table.setDescription("Partially Extended Flux Density Colour Correction (Modified Black Body, beta=%.2f)"%betaK)
        kColPSrc_BB_table.meta['beta']=DoubleParameter(betaK,"Emissivity spectral index")
        kColPSrc_BB_table.addColumn("Temperature",Column(Double1d(tempK)))
        for band in spireBands():
            kColPSrc_BB_table.addColumn(band,Column(kColPSrc_BB[band]))
        return kColPSrc_BB_table

def calcKColESrc(alphaK,src,verbose=False,table=False):
    """
    ========================================================================
    calcKColESrc(alphaK,src,verbose=False,table=False):
        Calculates peak surface brightness colour correction for source with
          power law spectrum and a given (frequency-independent) source profile

    Inputs:
        alphaK:  (scalar/float array) Spectral index(es) to use
        src:     (Source object class) Source object to use
        verbose: (boolean) Set to print more information. Default=True.
        table:   (boolean) Set to output TableDataset. Default=False.

    Outputs:
        If table==False:
            (dict) colour correction factor (one scalar/list per band,
              depending on whether alphaK is scalar/list)
        If table==True:
            (TableDataset) colour correction factor (one column per band,
              one row per alphaK value)
        
    Calculation:
        [See Handbook Eq 5.40-43]
        Calculates KMonE parameter for fully-extended source and pipeline spectrum (alpha=-1)
        Calculates flux conversion to get peak surface brightness for
          partially-resolved source with power law spectrum and profile defined
          by Source object class
        Divides by the KMonE Parameter to get colour correction
        If table is set, returns a TableDataset
        
        KColESrc converts from fully-extended source monochromatic pipeline
          surface brightness to peak surface brightness of partially-resolved source.

    Dependencies:
        sources_dev.Source - Source object class
        sources_dev.SourceProfile - Source profile object class
        SpireHandbookBundle_dev.spireBands() - get list of spire bands
        SpireHandbookBundle_dev.getCal() - retrieve photometer calibration product
        SpireHandbookBundle_dev.getSpireFreq() - retrieve frequency raster
        SpireHandbookBundle_dev.getSpireRefFreq() - retrieve reference frequencies
        calcBeamSrcMonoArea() - calculate area of source multiplied by beam
        herschel.ia.dataset.TableDataset()
    """
    #-----------------------------------------------------------------------
    #print '\nCalculating extended source colour correction parameters over alpha...'
    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #freq=getSpireFreq()
    #spireFilt=getSpireFilt(cal)

    k4E_Tot=sh.calcKMonE()
    print 'k4E_Tot:',k4E_Tot
    try:
        na=len(alphaK)
        aList=True
    except:
        aList=False

    if not aList:
        # alphaK is scalar
        kColESrc = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
        #kBeamK = calcKBeam(alphaK)
        for band in spireBands():
            k4EaTot_x=sh.calcSpireKcorr(sh.getSpireRefFreq()[band], sh.getSpireFreq(),\
             sh.getSpireFilt()[band], BB=False, alpha=alphaK,\
             ext=True, monoArea=calcBeamSrcMonoArea(src,verbose=verbose)[band])[0]/1.e6
            kColESrc[band] = k4EaTot_x / k4E_Tot[band]
        if (verbose): print 'Calculated KColE for alpha=%f and source %s: '%(alphaK,src.key),kColESrc
    else:
        # alphaK is list
        kColESrc = {'PSW': Double1d(na,Double.NaN), 'PMW': Double1d(na,Double.NaN), 'PLW': Double1d(na,Double.NaN)}
        for a in range(na):
            for band in spireBands():
                k4EaTot_x=sh.calcSpireKcorr(sh.getSpireRefFreq()[band], sh.getSpireFreq(),\
                 sh.getSpireFilt()[band], BB=False, alpha=alphaK[a],\
                 ext=True, monoArea=calcBeamSrcMonoArea(src,verbose=verbose)[band])[0]/1.e6
                kColESrc[band][a] = k4EaTot_x / k4E_Tot[band]
            if (verbose): print 'Calculated KColE for alpha=%f and source %s: '%(alphaK[a],src.key),kColESrc["PSW"][a],kColESrc["PMW"][a],kColESrc["PLW"][a]

    if not table:
        #return as is
        return kColESrc
    else:
        #create and returnTableDataset
        kColESrc_table=TableDataset()
        kColESrc_table.setDescription("Partially Extended Surface Brightness Colour Correction (Spectral Index)")
        kColESrc_table.addColumn("alpha",Column(Double1d(alphaK)))
        for band in spireBands():
            kColESrc_table.addColumn(band,Column(kColESrc[band]))
        return kColESrc_table

def calcKColESrc_BB(betaK,tempK,src,verbose=False,table=False):
    """
    ========================================================================
    calcKColESrc(alphaK,src,verbose=False,table=False):
        Calculates peak surface brightness colour correction for source with
          modified black body spectrum and a given (frequency-independent) source profile

    Inputs:
        betaK:   (float) Emmisivity spectral index to use
        tempK:   (scalar/float array) Temperature(s) to use
        src:     (Source object class) Source object to use
        verbose: (boolean) Set to print more information. Default=True.
        table:   (boolean) Set to output TableDataset. Default=False.

    Outputs:
        If table==False:
            (dict) colour correction factor (one scalar/list per band,
              depending on whether tempK is scalar/list)
        If table==True:
            (TableDataset) colour correction factor (one column per band,
              one row per tempK value)
        
    Calculation:
        [See Handbook Eq 5.40-43]
        Calculates KMonE parameter for fully-extended source and pipeline spectrum (alpha=-1)
        Calculates flux conversion to get peak surface brightness for
          partially-resolved source with mod-BB spectrum and profile defined
          by Source object class
        Divides by the KMonE Parameter to get colour correction
        If table is set, returns a TableDataset
        
        KColESrc converts from fully-extended source monochromatic pipeline
          surface brightness to peak surface brightness of partially-resolved source.

    Dependencies:
        sources_dev.Source - Source object class
        sources_dev.SourceProfile - Source profile object class
        SpireHandbookBundle_dev.spireBands() - get list of spire bands
        SpireHandbookBundle_dev.getCal() - retrieve photometer calibration product
        SpireHandbookBundle_dev.getSpireFreq() - retrieve frequency raster
        SpireHandbookBundle_dev.getSpireRefFreq() - retrieve reference frequencies
        calcBeamSrcMonoArea() - calculate area of source multiplied by beam
        herschel.ia.dataset.TableDataset()
    """
    #-----------------------------------------------------------------------
    #print '\nCalculating extended source colour correction parameters over alpha...'
    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #freq=getSpireFreq()
    #spireFilt=getSpireFilt(cal)

    k4E_Tot=sh.calcKMonE()
    print 'k4E_Tot:',k4E_Tot
    try:
        nt=len(tempK)
        tList=True
    except:
        tList=False

    if not tList:
        # alphaK is scalar
        kColESrc_BB = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
        #kBeamK = calcKBeam(alphaK)
        for band in spireBands():
            k4EaTot_x=sh.calcSpireKcorr(sh.getSpireRefFreq()[band], sh.getSpireFreq(),\
             sh.getSpireFilt()[band], BB=True, beta=betaK, temp=tempK,\
             ext=True, monoArea=calcBeamSrcMonoArea(src,verbose=verbose)[band])[0]/1.e6
            kColESrc_BB[band] = k4EaTot_x / k4E_Tot[band]
        if (verbose): print 'Calculated KColE for modBB with T=%f K, beta=%f and source %s: '%(betaK,tempK,src.key),kColESrc_BB
    else:
        # alphaK is list
        kColESrc_BB = {'PSW': Double1d(nt,Double.NaN), 'PMW': Double1d(nt,Double.NaN), 'PLW': Double1d(nt,Double.NaN)}
        for t in range(nt):
            for band in spireBands():
                k4EaTot_x=sh.calcSpireKcorr(sh.getSpireRefFreq()[band], sh.getSpireFreq(),\
                 sh.getSpireFilt()[band], BB=True, beta=betaK, temp=tempK[t],\
                 ext=True, monoArea=calcBeamSrcMonoArea(src,verbose=verbose)[band])[0]/1.e6
                kColESrc_BB[band][t] = k4EaTot_x / k4E_Tot[band]
            if (verbose): print 'Calculated KColE for modBB with T=%f K, beta=%f and source %s: '%(betaK,tempK[t],src.key),kColESrc_BB["PSW"][t],kColESrc_BB["PMW"][t],kColESrc_BB["PLW"][t]

    if not table:
        #return as is
        return kColESrc_BB
    else:
        #create and returnTableDataset
        kColESrc_BB_table=TableDataset()
        kColESrc_BB_table.setDescription("Partially Extended Surface Brightness Colour Correction (Spectral Index)")
        kColESrc_BB_table.meta['beta']=DoubleParameter(betaK,"Emissivity spectral index")
        kColESrc_BB_table.addColumn("Temperature",Column(Double1d(tempK)))
        for band in spireBands():
            kColESrc_BB_table.addColumn(band,Column(kColESrc_BB[band]))
        return kColESrc_BB_table

def semiExtendedTest():
    """
    ================================================================================
    semiExtendedTest():
        Test script to demonstrate use of partially-resolved colour correction
          calculations
        Demonstrates that for very small sources, tends to point source case,
          and for very large sources, tends to fully extended case

    Inputs:
        NONE
        
    Outputs:
        (float array) source widths used
        (list) KColPSrc flux density colour correction for all source widths
                 (one float array per band)
        (list) KColESrc peak surface brightness colour correction for all source
                 widths (one float array per band)
                    
    Calculation:
      - Fixes spectrum at alphaK=2.0 for testing
      - Calculates point-source flux density colour correction for spectrum
      - Calculates fully extended surface brightness colour correction for spectrum
      + For range of source widths, from 0.2arcsec to 10^4 arcsec:
        - Generates Source objects for each width
        - Calculates partially resolved flux density colour correction for
            each source width and the fixed spectrum
        - Calculates partially resolved peak surface brightness colour correction
            for each source width and the the fixed spectrum
      - Puts resulting corrections into lists
    
    Dependencies:
        SpireHandbookBundle_dev.calcKColP() - calculate point source colour correction for power law spectrum
        SpireHandbookBundle_dev.calcKColE() - calculate fully-extended source colour correction for power law spectrum
        calcKColPSrc() - calculate partially resolved flux density correction for power law spectrum
        calcKColESrc() - calculate partially resolved peak surface brightness correction for power law spectrum
    """
    #make Gaussian sources
    alphaK=2.0
    srcWidths=Float1d([0.2,0.5,1.,2.,5.,10.,20.,50.,100.,200.,500.,1000.,2000.,5000.,10000.])
    KColE_full=sh.calcKColE(alphaK,verbose=True)
    KColP_point=sh.calcKColP(alphaK,verbose=True)
    print 'full',KColE_full
    print 'point',KColP_point
    KColE_partialPSW=[]
    KColE_partialPMW=[]
    KColE_partialPLW=[]
    KColP_partialPSW=[]
    KColP_partialPMW=[]
    KColP_partialPLW=[]
    for wid in srcWidths:
        src=srcMod.Source('Gaussian',[wid])
        #calculate flux density colour correction
        KColP_partial=calcKColPSrc(alphaK,src,verbose=True,table=False)
        KColP_partialPSW.append(KColP_partial['PSW'])
        KColP_partialPMW.append(KColP_partial['PMW'])
        KColP_partialPLW.append(KColP_partial['PLW'])
        #calculate surface brightness colour correction
        KColE_partial=calcKColESrc(alphaK,src,verbose=True,table=False)
        KColE_partialPSW.append(KColE_partial['PSW'])
        KColE_partialPMW.append(KColE_partial['PMW'])
        KColE_partialPLW.append(KColE_partial['PLW'])
    KColPSrc=[Float1d(KColP_partialPSW),Float1d(KColP_partialPMW),Float1d(KColP_partialPLW)]
    KColESrc=[Float1d(KColE_partialPSW),Float1d(KColE_partialPMW),Float1d(KColE_partialPLW)]
    return(srcWidths,KColPSrc,KColSrc)
    
def test2():
    srcWidth=1e6
    src=srcMod.Source('Gaussian',srcWidth)

    alphaK=[2.]
    KColE_full=sh.calcKColE(alphaK,verbose=False)
    KColP_point=sh.calcKColP(alphaK,verbose=False)
    print 'point',KColP_point
    print 'full',KColE_full
    KColP_partial=calcKColPSrc(alphaK,src,verbose=False,table=False)
    KColE_partial=calcKColESrc(alphaK,src,verbose=False,table=False)
    print 'partial P',KColP_partial
    print 'partial E',KColE_partial

    print '-----\nBB'
    betaK=1.8
    tempK=[20.]
    srcWidth=5.
    KColE_BB_full=sh.calcKColE_BB(betaK,tempK,verbose=False)
    KColP_BB_point=sh.calcKColP_BB(betaK,tempK,verbose=False)
    print 'point (BB)',KColP_BB_point
    print 'full (BB)',KColE_BB_full
    KColP_BB_partial=calcKColPSrc_BB(betaK,tempK,src,verbose=False,table=False)
    KColE_BB_partial=calcKColESrc_BB(betaK,tempK,src,verbose=False,table=False)
    print 'partial P (BB_',KColP_BB_partial
    print 'partial E (BB)',KColE_BB_partial