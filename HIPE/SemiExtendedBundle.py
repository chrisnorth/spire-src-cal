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
from herschel.share.util import Configuration
from herschel.ia.numeric import Double1d,Float1d,Int1d,String1d
from herschel.ia.dataset import TableDataset,Column
from herschel.ia.toolbox.util import AsciiTableWriterTask, SimpleFitsWriterTask
from herschel.ia.toolbox.util import AsciiTableReaderTask, SimpleFitsReaderTask
asciiTableWriter=AsciiTableWriterTask()
simpleFitsWriter=SimpleFitsWriterTask()
asciiTableReader=AsciiTableReaderTask()
simpleFitsReader=SimpleFitsReaderTask()
#from herschel.spire.ia.cal import SpireCalTask
#spireCal = SpireCalTask()
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


def calcBeamSrcMonoArea(src,verbose=False,forceRecalc=False):
    """
    ========================================================================
    calcBeamSrcMonoArea(src,verbose=False):
        Calculate monochromatic beam areas over frequency raster, convolved with
          source profile
        Stored as entry in global variable beamMonoSrcArea, so only calculated
          once per source

    Inputs:
      src:         (Source class) Source object to multiply beam profile by
                   OR
                   (SourceProfile class) Source profile object to multiply beam profile by
      verbose:     (boolean) print more information. Default=False
      forceRecalc: (boolean) force recalculation of profile, even if already stored

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

    assert type(src)==srcMod.Source or type(src)==srcMod.SourceProfile,\
      'src mus be either Source or SourceProfile object'
    #read in beamProfiles
    beamProfs = sh.getCal().getProduct("RadialCorrBeam")
    beamRad=beamProfs.getCoreCorrectionTable().getColumn('radius').data

    if type(src)==srcMod.Source:
        #create source profile
        
        fromSrc=True
        srcProf=src.calcProfile(beamRad)
        if verbose:print'Creating SourceProfile %s from Source %s'%(srcProf.key,src.key)
        
    elif type(src)==srcMod.SourceProfile:
        fromSrc=False
        srcProf=src.copy()

    if srcProf.key!=None:
        #store in beamMonoSrcArea global variable        
        try:
            beamSrcArea=beamMonoSrcArea[srcProf.key]
            if (verbose):print'Using existing monochromatic areas (key %s)'%srcProf.key
            #exists, only recalculate if forced
            reCalc=forceRecalc
        except:
            #doesn't exist, recalculate
            reCalc=True
            #check if main structure exists
            try:
                beamMonoSrcArea
                #exists, don't create
            except:
                #doesn't exist, must create
                beamMonoSrcArea={}
    else:
        #no key
        reCalc=True
    
    if reCalc:
        if (verbose):print'Calculating monochromatic areas (key %s)'%srcProf.key
        beamSrcArea = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
        
        if 0.5*srcProf.calcFwhm()<beamRad[10]:
            #quite a small source (< 10 steps in radius array)
            beamProfsFine=fineBeam(beamProfs,0.5*srcProf.calcFwhm())
            beamRadFine=beamProfsFine.getCoreCorrectionTable().getColumn('radius').data
            #define srcProfFine to be same length as beamRadFine
            srcProfFine=srcProf.regrid(beamRadFine)
            srcAreaFine=srcProfFine.calcArea()*sh.arcsec2Sr()
            if srcProf.calcFwhm()<beamRad[1] and fromSrc==True:
                #VERY small (scale width < 1 step in radius array)
                if (verbose):print 'Very small source [FWHM=%g arcsec]. Using source area only (ignoring beam)'%srcProf.calcFwhm()
                for band in spireBands():
                    spireFreq=sh.getSpireFreq()
                    beamSrcArea[band]=Double1d(len(spireFreq),srcAreaFine)
            else:
                #quote small source (scale width 1-10 steps in radius array
                if (verbose):print 'Calculating monochromatic beam areas for small [FHWM=%g arcsec] %s SourceProfile'%(srcProf.calcFwhm(),src.key)
                for band in spireBands():
                    #monochromatic beam areas
                    gamma = beamProfsFine.meta['gamma'].double
                    beamSrcArea[band] = beam.spireMonoSrcAreas(sh.getSpireFreq(), beamProfsFine, 
                      sh.getSpireEffFreq()[band], gamma, srcProfFine, band)
        else:
            #calculate area of beam and source
            gamma = beamProfs.meta['gamma'].double
            if (verbose):print 'Calculating monochromatic beam areas for %s SourceProfile'%(srcProf.key)
            if srcProf.checkRadArr(beamRad)==False:
                #regrid to match beam profile
                srcProfRegrid=srcProf.regrid(beamRad)
                if (verbose):print 'Regridding %s SourceProfile to match beam profile'
            else:
                srcProfRegrid=srcProf.copy()
            for band in spireBands():
                #monochromatic beam areas
                beamSrcArea[band] = beam.spireMonoSrcAreas(sh.getSpireFreq(), beamProfs, 
                  sh.getSpireEffFreq()[band], gamma, srcProfRegrid, band)

    if srcProf.key!=None:
        #store in global variable
        beamMonoSrcArea[srcProf.key]=beamSrcArea
    #if (verbose):print 'All beam areas calculated for %s source'%src.key
    return (beamSrcArea)

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

    #make whole array fine-sampled
    maxRadFine=MAX(beamRadIn)
    radFineStep=scaleWidth/10.
    nRadFine=int(maxRadFine/radFineStep)
    print 'fine sampling entire beam profile at %g arcsec for %g arcsec (%d steps)'%(radFineStep,maxRadFine,nRadFine)
    beamRad=Double1d(range(nRadFine)) * radFineStep

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
        
    return(beamProfs)

#-------------------------------------------------------------------------------
# Calculate the effective beam profile and area for SPIRE
def calcSpireSrcBeam(freq0,freq, transm, beamProfs, effFreq, gamma, array, BB=False,temp=20.0,beta=1.8,alpha=-1.0,verbose=False):

    """
    ========================================================================
    Computes an effective beam profile for a given source spectrum
    ***
    N.B. Computing the beam area this way integrates over frequency *then* radius,
      while spireEffArea integrates over radius then frequency.
      The method used here produces areas which area lower by ~0.1%
    ***

    Inputs:
      freq:       (array float) frequency vector [Hz] for which monochromatic
                    beams areas should be calculated
      transm:     (array float) relative spectral response (RSRF) corresponding to freq
                    Note that this should *not* include the aperture efficiency
      beamProfs:  (dataset) PhotRadialCorrBeam object from calibration tree
      effFreq:    (float) effective frequency [Hz] of array
      gamma:      (float) Exponent of powerlaw describing FWHM dependence
                    on frequency
      array:      (string) spire array ('Psw'|'Pmw'|'Plw')
      BB:         (boolean) set to use modified black body spectrum with
                        temperature temp and emissivity index beta
                    OPTIONAL. Default=False
      alpha:      (float) Exponent of power-law sky background model
                        (if BB==False)
                    OPTIONAL. Default=-1
      temp:       (float) Dust/sky temperature (if BB=True)
                    OPTIONAL. Default=20.0
      beta:       (float) Dust/sky spectral index (if BB=True)
                    OPTIONAL. Default=1.8
      verbose:    (boolean) Set to print additional information.
                    OPTIONAL. Default=False.

    Outputs:
      (array float) radialised beam profile

    Calculation:
      Calculates the source spectrum (either modifies black body or power law)
      Loops over beam radius, and calculates scaled radii for full frequency range
      For that radius, gets the values of profile at those scaled radii
      Integrates profile values over frequency, weighted by RSRF

    Dependencies:
      herschel.ia.numeric.toolbox.interp.LinearInterpolator
      herschel.ia.numeric.toolbox.interp.CubicSplineInterpolator
      herschel.ia.numeric.toolbox.integr.TrapezoidalIntegrator

    2014/01/16  C. North  initial version

    """

    h = Constant.H_PLANCK.value
    k = Constant.K_BOLTZMANN.value
    c = Constant.SPEED_OF_LIGHT.value

    #
    # Calculate sky background model
    #
    if BB == 1:
        # Modified black body
        if verbose:
            print 'Using greybody with T=%.2f, beta=%.3f'%(temp,beta)
        fSky  = 2*h * freq**3 / c**2 / (EXP(h*freq/k/temp) - 1.) * freq**beta
        fSky0 = 2*h * freq0**3 / c**2 / (EXP(h*freq/k/temp) - 1.) * freq0**beta
    #
    else:
        # Power-Law spectrum
        if verbose:
            print 'Using power law spectrum with spectral index %.2f'%(alpha)
        fSky  = freq**alpha
        fSky0 = freq0**alpha


    #integrate transm*fSky over frequency for nomalisation
    integrator=TrapezoidalIntegrator(min(freq),max(freq))
    denomInterp=CubicSplineInterpolator(freq,freq0*transm)
    denomInteg=integrator.integrate(denomInterp)

    #get beam radius list from calibration table
    beamRad=beamProfs.getCoreCorrectionTable().getColumn('radius').data
    #get core beam profile from calibration table
    beamCore=beamProfs.getCoreCorrectionTable().getColumn(array).data
    #create interpolation object
    beamCoreInt=CubicSplineInterpolator(beamRad,beamCore)

    #make array for new beam
    nRad=len(beamRad)
    maxRad=max(beamRad)
    effBeam=Float1d(beamRad)
    #loop over radius
    for r in range(nRad):
        #calculate the "scaled" radius for range of frequencies
        radFreq=beamRad[r]*(freq/effFreq)**-gamma
        #ensure it doesn't fo beyong maximum radius
        radFreq[radFreq.where(radFreq > maxRad)]=maxRad
        #compute value beam profile at each scaled radius
        beamCoreFreq=beamCoreInt(radFreq)
        #apply constant beam profile value where appropriate
        beamConstRad=beamProfs.getConstantCorrection(beamRad[r],array)
        isConst=beamCoreFreq.where(beamCoreFreq < beamConstRad)
        beamCoreFreq[isConst]=beamConstRad

        #integrate beamCoreFreq*transm*fSky over frequency
        numInterp=CubicSplineInterpolator(freq,beamCoreFreq*transm*fSky)
        numInteg = integrator.integrate(numInterp)

        #write value into table
        effBeam[r]=numInteg/denomInteg

    #if verbose:print'calculated area:',effBeamArea
    return(effBeam)


def calcEffBeamSrc(alphaK,src,verbose=False,forceRecalc=False):
    """
    ================================================================================
    calcEffBeamSrc(alphaK,src,verbose=False,forceRecalc=False):
        Calculates effective beam for power law spectrum convolved with source profile
        Returns beam profile and stores in global variable

    Inputs:
        alphaK:      (float) power law spectral index to use
        src:         (Source class) Source object to use
                     OR
                     (SourceProfile class) SourceProfile object to use
        verbose:     (boolean) Set to print more information. Default=True.
        forceRecalc: (boolean) Set to force recalculation of profiles
    Outputs:
        (dict SourceProfile) SourceProfile objects containing convolved beam (one per band)
                    
    Calculation:
      + If corresponding key exists in global variable effBeamSrcProfs and forceRecalc=False:
          - Returns previously computed profile
      + If not:
          - Uses SoruceProfile, or computes SoureProfile from Source object
          - Calculates effective beam profile and area for given spectrum (can be slow)
          - Makes maps from source profile and beam map
          - Convolves maps and converts to SourceProfile object
    
    Dependencies:
        spireBands() - list of spire bands
        sources_dev.Source
        sources_dev.SourceProfile
        SpireHandbookBundle_dev.getCal() - retrieve photometer calibration product
        SpireHandbookBundle_dev.spireEffBeam() - calculate effective beam
        sources_dev.convolveProfiles()  -convolve beam profiles
    """
    
    #check type of src
    assert type(src)==srcMod.Source or type(src)==srcMod.SourceProfile,\
      'src must be Source object or SourceProfile object'
      
    #make global variable
    global effBeamSrcProfs

    #get beam profile from cal product
    beamProfs=sh.getCal().getProduct('RadialCorrBeam')
    beamRad=beamProfs.getCoreCorrectionTable().getColumn('radius').data
        
    if type(src)==srcMod.Source:
        #create source profile
        fromSrc=True
        srcProf=src.calcProfile(beamRad)
        if verbose:print'Creating SourceProfile %s from Source %s'%(srcProf.key,src.key)
        
    elif type(src)==srcMod.SourceProfile:
        fromSrc=False
        srcProf=src.copy()
        if verbose:print'Using SourceProfile %s'%(srcProf.key)
    
    
    if srcProf.key!=None:
        effBeamKey='%s_alpha_%g'%(srcProf.key,alphaK)
        try:
            srcConv=effBeamSrcProfs[effBeamKey]
            if verbose:print 'Using existing convolved beam profile %s'%effBeamKey
            reCalc=forceRecalc
            #global variable already defined, so do nothing.
        except:
            if verbose:print 'Calculating convolved beams for %s'%effBeamKey
            try:
                effBeamSrcProfs
                #exists, don't create
            except:
                #doesn't exist, must create
                effBeamSrcProfs={}#not defined, so need to recalculate
            reCalc=True
            
    else:
        if verbose:print 'no key for SourceProfile. Recalculating'
        effBeamKey=None
        reCalc=True
        
    if reCalc:
        gamma=beamProfs.meta['gamma'].double
        maxRad=MAX(beamRad)
        nRad=len(beamRad)
        
        #make blank variable
        srcConv={'PSW':Double.NaN,'PMW':Double.NaN,'PLW':Double.NaN}
        #calculate effective beams (no source)
        effBeams=sh.calcSpireEffBeam(alphaK,verbose=verbose)
        for band in spireBands():
            if verbose:print 'Computing %s effective beam for %s'%(effBeamKey,band)
            effBeam_x=effBeams[band]
            effBeamProfile=srcMod.SourceProfile(beamRad,effBeam_x,key='EffBeam_alpha_%g_%s'%(alphaK,band))
            effBeamNormProfile=effBeamProfile.normArea()
            srcConv[band]=srcMod.convolveProfiles(srcProf,effBeamNormProfile,verbose=verbose,key='%s_%s'%(effBeamKey,band))
        if effBeamKey!=None:
            effBeamSrcProfs[effBeamKey]=srcConv
            
    return(srcConv)

def calcEffBeamSrc_BB(betaK,tempK,src,verbose=False,forceRecalc=False):
    """
    ================================================================================
    calcEffBeamSrc(betaK,tempK,src,verbose=False,forceRecalc=False):
        Calculates effective beam for modified black body convolved with source profile
        Returns beam profile and stores in global variable

    Inputs:
        betaK:       (float) Emmisivity spectral index to use
        tempK:       (scalar/float array) Temperature to use
        src:         (Source class) Source object to use
                     OR
                     (SourceProfile class) SourceProfile object to use
        verbose:     (boolean) Set to print more information. Default=True.
        forceRecalc: (boolean) Set to force recalculation of profiles
    Outputs:
        (dict SourceProfile) SourceProfile objects containing convolved beam (one per band)
                    
    Calculation:
      + If corresponding key exists in global variable effBeamSrcProfs and forceRecalc=False:
          - Returns previously computed profile
      + If not:
          - Uses SoruceProfile, or computes SoureProfile from Source object
          - Calculates effective beam profile and area for given spectrum (can be slow)
          - Makes maps from source profile and beam map
          - Convolves maps and converts to SourceProfile object
    
    Dependencies:
        spireBands() - list of spire bands
        sources_dev.Source
        sources_dev.SourceProfile
        SpireHandbookBundle_dev.getCal() - retrieve photometer calibration product
        SpireHandbookBundle_dev.spireEffBeam() - calculate effective beam
        sources_dev.convolveProfiles()  -convolve beam profiles
    """
    
    #check type of src
    assert type(src)==srcMod.Source or type(src)==srcMod.SourceProfile,\
      'src must be Source object or SourceProfile object'
      
    #make global variable
    global effBeamSrcProfs

    #get beam profile from cal product
    beamProfs=sh.getCal().getProduct('RadialCorrBeam')
    beamRad=beamProfs.getCoreCorrectionTable().getColumn('radius').data
        
    if type(src)==srcMod.Source:
        #create source profile
        fromSrc=True
        srcProf=src.calcProfile(beamRad)
        if verbose:print'Creating SourceProfile %s from Source %s'%(srcProf.key,src.key)
        
    elif type(src)==srcMod.SourceProfile:
        fromSrc=False
        srcProf=src.copy()
        if verbose:print'Using SourceProfile %s'%(srcProf.key)
    
    
    if srcProf.key!=None:
        effBeamKey='%s_beta_%g_temp_%g'%(srcProf.key,betaK,tempK)
        try:
            srcConv=effBeamSrcProfs[effBeamKey]
            if verbose:print 'Using existing convolved beam profile %s'%effBeamKey
            reCalc=forceRecalc
            #global variable already defined, so do nothing.
        except:
            if verbose:print 'Calculating convolved beams for %s'%effBeamKey
            try:
                effBeamSrcProfs
                #exists, don't create
            except:
                #doesn't exist, must create
                effBeamSrcProfs={}#not defined, so need to recalculate
            reCalc=True
            
    else:
        if verbose:print 'no key for SourceProfile. Recalculating'
        effBeamKey=None
        reCalc=True
        
    if reCalc:
        gamma=beamProfs.meta['gamma'].double
        maxRad=MAX(beamRad)
        nRad=len(beamRad)
        
        #make blank variable
        srcConv={'PSW':Double.NaN,'PMW':Double.NaN,'PLW':Double.NaN}
        #calculate effective beams
        effBeams=sh.calcSpireEffBeam_BB(betaK,tempK,verbose=verbose)
        for band in spireBands():
            if verbose:print 'Computing %s effective beam for %s'%(effBeamKey,band)
            effBeam_x=effBeams[band]
            effBeamProfile=srcMod.SourceProfile(beamRad,effBeam_x,key='EffBeam_beta_%g_temp_%g_%s'%(betaK,tempK,band))            
            effBeamNormProfile=effBeamProfile.normArea()
                        
            srcConv[band]=srcMod.convolveProfiles(srcProf,effBeamNormProfile,key='%s_%s'%(effBeamKey,band))
        if effBeamKey!=None:
            effBeamSrcProfs[effBeamKey]=srcConv

    return(srcConv)

def calcKColPSrc(alphaK,src,verbose=False,table=False):
    """
    ========================================================================
    calcKColPSrc(alphaK,src,verbose=False,table=False):
        Calculates flux density colour correction for source with power law
          spectrum and a given (frequency-independent) source profile

    Inputs:
        alphaK:  (scalar/float array) Spectral index(es) to use
        src:     (Source class) Source object to use
                 OR
                 (SourceProfile class) SourceProfile object to use
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

    #check type of src
    assert type(src)==srcMod.Source or type(src)==srcMod.SourceProfile,\
      'src must be Source object or SourceProfile object'
      
    k4P=sh.calcK4P()
    print 'k4P:',k4P
    try:
        na=len(alphaK)
        aList=True
    except:
        aList=False

    #calculate source area analytically (if possible)
    srcArea=src.calcArea()
    print 'source area:',srcArea
    if Double.isNaN(srcArea):
        #calculate SourceProfile from src (if necessary)
        if type(src)==Source:
            #make profile and integrate
            beamRad=sh.getCal().getProduct("RadialCorrBeam").getCoreCorrectionTable().getColumn('radius').data
            srcArea=src.calcProfile(beamRad).calcArea()
    print srcArea
    assert (not Double.isNaN(srcArea)),'Problem calculating source are for %s (%s)'%(src.key,type(src))
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
        src:     (Source class) Source object to use
                 OR
                 (SourceProfile class) SourceProfile object to use
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

    #check type of src
    assert type(src)==srcMod.Source or type(src)==srcMod.SourceProfile,\
      'src must be Source object or SourceProfile object'
      
    k4P=sh.calcK4P()
    print 'k4P:',k4P
    try:
        nt=len(tempK)
        tList=True
    except:
        tList=False

    #calculate source area analytically (if possible)
    srcArea=src.calcArea()
    if Double.isNaN(srcArea):
        #calculate SourceProfile from src (if necessary)
        if type(src)==Source:
            #make profile and integrate
            beamRad=sh.getCal().getProduct("RadialCorrBeam").getCoreCorrectionTable().getColumn('radius').data
            srcArea=src.calcProfile(beamRad).calcArea()

    assert (not Double.isNaN(srcArea)),'Problem calculating source are for %s (%s)'%(src.key,type(src))
    
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
        src:     (Source class) Source object to use
                 OR
                 (SourceProfile class) SourceProfile object to use
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

    #check type of src
    assert type(src)==srcMod.Source or type(src)==srcMod.SourceProfile,\
      'src must be Source object or SourceProfile object'
      
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
        src:     (Source class) Source object to use
                 OR
                 (SourceProfile class) SourceProfile object to use
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

    #check type of src
    assert type(src)==srcMod.Source or type(src)==srcMod.SourceProfile,\
      'src must be Source object or SourceProfile object'
            
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

def calcApCorrSrc(alphaK,src,aperture=[22., 30.,45.],annulus=[60.,90],verbose=False,table=False):
    """
    ================================================================================
    calcApCorr(alphaK,src,aperture=[22., 30.,45.],annulus=[60.,90],verbose=False,table=False):
        Calculates aperture corrections for semi-extended source with power law spectrum
        Returns corrections both with and without background annulus

    Inputs:
        alphaK:   (float) power law spectral index to use
        src:      (Source class) Source object to use
                  OR
                  (SourceProfile class) SourceProfile object to use
        aperture: (list float) source apperture radius (in arcsec) for each band.
                    Default=[22., 30.,45.]
        annulus:  (list float) background annulus radius (in arcsec) for each band
                    Either 2-element list (uses same for each band), or list of 
                    three 2-element lists. Default=[60.,90.]
        verbose:  (boolean) Set to print more information. Default=True.
        table:    (boolean) Set to output TableDataset. Default=False.
    Outputs:
        If table==False:
            (dict) aperture correction INCLUDING background annulus (one scalar/list
              per band, depending on whether alphaK is scalar/list)
            (dict) aperture correction WITHOUT background annulus (one scalar/list
              per band, depending on whether alpha is scalar/list)
        If table==True:
            (TableDataset) aperture correction INCLUDING background annulus
              (one column per band, one row per alphaK value)
            (TableDataset) aperture correction WITHOUT background annulus
              (one column per band, one row per alphaK value)
                    
    Calculation:
        Takes effective beam profile convo
        Integrates profile over aperture and annulus
        Calculates aperture correction factors, including and excluding background annulus
    
    Dependencies:
        spireBands() - list of spire bands
        SpireHandbookBundle_dev.getCal() - retrieve photometer calibration product
        calceffBeamSrc() - calculate effective beam convolved with source
        herschel.ia.numeric.toolbox.interp.CubicSplineInterpolator
        herschel.ia.numeric.toolbox.integr.TrapezoidalIntegrator
        herschel.ia.dataset.TableDataset()
    """
    
    #check type of src
    assert type(src)==srcMod.Source or type(src)==srcMod.SourceProfile,\
      'src must be Source object or SourceProfile object'
    #read aperture from input
    assert type(aperture)==list,'aperture must be 3-element list'
    assert len(aperture)==3,'aperture must be 3-element list'
    apPhotRad={}
    try:
        apPhotRad["PSW"]=aperture[0]
        apPhotRad["PMW"]=aperture[1]
        apPhotRad["PLW"]=aperture[2]
    except:
        'Error reading from aperture: ',aperture
    
    #read annulus from input
    assert type(annulus)==list,'Annulus must be 2 or 3-element list'
    assert len(annulus)==3 or len(annulus)==2,'Annulus must be 2 or 3-element list'
    apPhotBGRad={}
    if len(annulus)==2:
        #use same for each band
        for band in spireBands():
            apPhotBGRad[band]=annulus
    else:
        #use different for each band
        for b in range(3):
            band=spireBands()[b]
            apPhotBGRad[band]=annulus[b]
            assert type(apPhotBGRad[band])==list,'annulus[%d] is not 2-element list: '+str(annulus[b])
            assert len(apPhotBGRad[band])==2,'annulus[%d] is not 2-element list: '+str(annulus[b])

    #check they are valid
    assert len(apPhotRad)==3,'Error reading from aperture: '+str(aperture)
    assert len(apPhotBGRad)==3,'Error reading from annulus: '+str(annulus)
    
    try:
        na=len(alphaK)
        aList=True
    except:
        aList=False
    
    #get beam profile from cal product
    beamProfs=sh.getCal().getProduct('RadialCorrBeam')
    beamRad=beamProfs.getCoreCorrectionTable().getColumn('radius').data
    gamma=beamProfs.meta['gamma'].double
    maxRad=MAX(beamRad)
    #setup integrator over whole beam
    integTot=TrapezoidalIntegrator(0.,maxRad)

    #calculate SourceProfile from src (if necessary)
    if type(src)==srcMod.Source:
        srcProf=src.calcProfile(beamRad)
    else:
        if src.checkRadArr(beamRad):
            srcProf=src.regrid(beamRad)
        else:
            srcProf=src.copy()
        
    #normalised to total area 1.
    srcProfNorm=srcProf.normArea()
    srcArea=srcProf.calcArea(forceNumerical=True)
    
    #create interpolator for computing aperture/annulus size
    sizeInterp = CubicSplineInterpolator(beamRad,2.*PI*beamRad)
    
    if not aList:
        #alphaK is scalar
        apCorrNoBG={'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
        apCorrIncBG={'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
        
        #calculate effective beam
        effBeamsSrc=calcEffBeamSrc(alphaK,src,verbose=True)
        
        for band in spireBands():
            #create integrators
            integAp=TrapezoidalIntegrator(0.,apPhotRad[band])
            integBG=TrapezoidalIntegrator(apPhotBGRad[band][0],apPhotBGRad[band][1])
            integTot=TrapezoidalIntegrator(0.,maxRad)
            #integrate to find size of aperture/annulus
            sizeAp=integAp.integrate(sizeInterp)
            sizeBG=integBG.integrate(sizeInterp)
            
            
            srcConv=effBeamsSrc[band]
            srcConvArea=srcConv.calcArea()
            
            #interpolate convolved profile
            beamInterp = CubicSplineInterpolator(srcConv.radArr,srcConv.profile * 2.*PI*srcConv.radArr)
            
            #if verbose:print '%s aperture/background sizes:'%band,sizeAp,sizeBG
            #integrate profile to find beam area in aperture/annulus
            omegaAp=integAp.integrate(beamInterp)
            omegaBG=integBG.integrate(beamInterp)
            omegaTot=integTot.integrate(beamInterp)	
            #if verbose:print 'aperture/background beam areas:',omegaAp,omegaBG
            #calculate apertur  ecorrections
            apCorrNoBG[band] = \
                srcConvArea/omegaAp
            apCorrIncBG[band] =  \
                srcConvArea/(omegaAp - omegaBG*sizeAp/sizeBG)
        if (verbose):
            print 'Calculated apCorr (noBG) for alpha=%f: '%alphaK,apCorrNoBG["PSW"],apCorrNoBG["PMW"],apCorrNoBG["PLW"]
            print 'Calculated apCorr (incBG) for alpha=%f: '%alphaK,apCorrIncBG["PSW"],apCorrIncBG["PMW"],apCorrIncBG["PLW"]
    else:
        #alphaK is scalar
        apCorrNoBG={'PSW': Double1d(na,Double.NaN), 'PMW': Double1d(na,Double.NaN), 'PLW': Double1d(na,Double.NaN)}
        apCorrIncBG={'PSW': Double1d(na,Double.NaN), 'PMW': Double1d(na,Double.NaN), 'PLW': Double1d(na,Double.NaN)}
        
        for a in range(na):
            #calculate effective beam
            effBeamsSrc=calcEffBeamSrc(alphaK[a],src,verbose=True)
            
            for band in spireBands():
                #create integrators
                integAp=TrapezoidalIntegrator(0.,apPhotRad[band])
                integBG=TrapezoidalIntegrator(apPhotBGRad[band][0],apPhotBGRad[band][1])
                #integrate to find size of aperture/annulus
                sizeAp=integAp.integrate(sizeInterp)
                sizeBG=integBG.integrate(sizeInterp)
                
                srcConv=effBeamsSrc[band]
                srcConvArea=srcConv.calcArea()
                
                #interpolate convolved profile
                beamInterp = CubicSplineInterpolator(srcConv.radArr,srcConv.profile * 2.*PI*srcConv.radArr)
            
                #integrate profile to find beam area in aperture/annulus
                omegaAp=integAp.integrate(beamInterp)
                omegaBG=integBG.integrate(beamInterp)
                omegaTot=integTot.integrate(beamInterp)
                #calculate aperture corrections
                apCorrNoBG[band][a] = \
                    srcConvArea/omegaAp    
                apCorrIncBG[band][a] =  \
                    srcConvArea/(omegaAp - omegaBG*sizeAp/sizeBG)
            if (verbose):
                print 'Calculated apCorr (noBG) for alpha=%f: '%alphaK[a],apCorrNoBG["PSW"][a],apCorrNoBG["PMW"][a],apCorrNoBG["PLW"][a]
                print 'Calculated apCorr (incBG) for alpha=%f: '%alphaK[a],apCorrIncBG["PSW"][a],apCorrIncBG["PMW"][a],apCorrIncBG["PLW"][a]
        
    if not table:
        #return as is
        return(apCorrIncBG,apCorrNoBG)
    else:
        #create and returnTableDataset
        apCorrNoBG_table=TableDataset()
        apCorrIncBG_table=TableDataset()
        apCorrNoBG_table.setDescription("Aperture Correction without background annulus (Spectral Index)")
        apCorrIncBG_table.setDescription("Aperture Correction including background annulus (Spectral Index)")
        apCorrNoBG_table.addColumn("alpha",Column(Double1d(alphaK)))
        apCorrIncBG_table.addColumn("alpha",Column(Double1d(alphaK)))
        for band in spireBands():
            apCorrNoBG_table.addColumn(band,Column(apCorrNoBG[band]))
            apCorrIncBG_table.addColumn(band,Column(apCorrIncBG[band]))
        return (apCorrIncBG_table,apCorrNoBG_table)
    
def calcApCorrSrc_BB(betaK,tempK,src,aperture=[22., 30.,45.],annulus=[60.,90],verbose=False,table=False):
    """
    ================================================================================
    calcApCorr(betaK,tempK,src,aperture=[22., 30.,45.],annulus=[60.,90],verbose=False,table=False):
        Calculates aperture corrections for semi-extended source with mod-BB spectrum
        Returns corrections both with and without background annulus

    Inputs:
        betaK:    (float) Emmisivity spectral index to use
        tempK:    (scalar/float array) Temperature(s) to use
        src:      (Source class) Source object to use
                  OR
                  (SourceProfile class) SourceProfile object to use
        aperture: (list float) source apperture radius (in arcsec) for each band.
                    Default=[22., 30.,45.]
        annulus:  (list float) background annulus radius (in arcsec) for each band
                    Either 2-element list (uses same for each band), or list of 
                    three 2-element lists. Default=[60.,90.]
        verbose:  (boolean) Set to print more information. Default=True.
        table:    (boolean) Set to output TableDataset. Default=False.
    Outputs:
        If table==False:
            (dict) aperture correction INCLUDING background annulus (one scalar/list
              per band, depending on whether alphaK is scalar/list)
            (dict) aperture correction WITHOUT background annulus (one scalar/list
              per band, depending on whether alpha is scalar/list)
        If table==True:
            (TableDataset) aperture correction INCLUDING background annulus
              (one column per band, one row per alphaK value)
            (TableDataset) aperture correction WITHOUT background annulus
              (one column per band, one row per alphaK value)
                    
    Calculation:
        Takes effective beam profile convo
        Integrates profile over aperture and annulus
        Calculates aperture correction factors, including and excluding background annulus
    
    Dependencies:
        spireBands() - list of spire bands
        SpireHandbookBundle_dev.getCal() - retrieve photometer calibration product
        calceffBeamSrc() - calculate effective beam convolved with source
        herschel.ia.numeric.toolbox.interp.CubicSplineInterpolator
        herschel.ia.numeric.toolbox.integr.TrapezoidalIntegrator
        herschel.ia.dataset.TableDataset()
    """
    
    #check type of src
    assert type(src)==srcMod.Source or type(src)==srcMod.SourceProfile,\
      'src must be Source object or SourceProfile object'
    #read aperture from input
    assert type(aperture)==list,'aperture must be 3-element list'
    assert len(aperture)==3,'aperture must be 3-element list'
    apPhotRad={}
    try:
        apPhotRad["PSW"]=aperture[0]
        apPhotRad["PMW"]=aperture[1]
        apPhotRad["PLW"]=aperture[2]
    except:
        'Error reading from aperture: ',aperture
    
    #read annulus from input
    assert type(annulus)==list,'Annulus must be 2 or 3-element list'
    assert len(annulus)==3 or len(annulus)==2,'Annulus must be 2 or 3-element list'
    apPhotBGRad={}
    if len(annulus)==2:
        #use same for each band
        for band in spireBands():
            apPhotBGRad[band]=annulus
    else:
        #use different for each band
        for b in range(3):
            band=spireBands()[b]
            apPhotBGRad[band]=annulus[b]
            assert type(apPhotBGRad[band])==list,'annulus[%d] is not 2-element list: '+str(annulus[b])
            assert len(apPhotBGRad[band])==2,'annulus[%d] is not 2-element list: '+str(annulus[b])

    #check they are valid
    assert len(apPhotRad)==3,'Error reading from aperture: '+str(aperture)
    assert len(apPhotBGRad)==3,'Error reading from annulus: '+str(annulus)
    
    try:
        nt=len(tempK)
        tList=True
    except:
        tList=False
    
    #get beam profile from cal product
    beamProfs=sh.getCal().getProduct('RadialCorrBeam')
    beamRad=beamProfs.getCoreCorrectionTable().getColumn('radius').data
    gamma=beamProfs.meta['gamma'].double
    maxRad=MAX(beamRad)
    #setup integrator over whole beam
    integTot=TrapezoidalIntegrator(0.,maxRad)

    #calculate SourceProfile from src (if necessary)
    if type(src)==srcMod.Source:
        srcProf=src.calcProfile(beamRad)
    else:
        if src.checkRadArr(beamRad):
            srcProf=src.regrid(beamRad)
        else:
            srcProf=src.copy()
        
    #normalised to total area 1.
    srcProfNorm=srcProf.normArea()
    srcArea=srcProf.calcArea(forceNumerical=True)
    
    #create interpolator for computing aperture/annulus size
    sizeInterp = CubicSplineInterpolator(beamRad,2.*PI*beamRad)
    
    if not tList:
        #alphaK is scalar
        apCorrNoBG={'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
        apCorrIncBG={'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
        
        #calculate effective beam
        effBeamsSrc=calcEffBeamSrc_BB(betaK,tempK,src,verbose=True)
        
        for band in spireBands():
            #create integrators
            integAp=TrapezoidalIntegrator(0.,apPhotRad[band])
            integBG=TrapezoidalIntegrator(apPhotBGRad[band][0],apPhotBGRad[band][1])
            integTot=TrapezoidalIntegrator(0.,maxRad)
            #integrate to find size of aperture/annulus
            sizeAp=integAp.integrate(sizeInterp)
            sizeBG=integBG.integrate(sizeInterp)
            
            
            srcConv=effBeamsSrc[band]
            srcConvArea=srcConv.calcArea()
            
            #interpolate convolved profile
            beamInterp = CubicSplineInterpolator(srcConv.radArr,srcConv.profile * 2.*PI*srcConv.radArr)
            
            #if verbose:print '%s aperture/background sizes:'%band,sizeAp,sizeBG
            #integrate profile to find beam area in aperture/annulus
            omegaAp=integAp.integrate(beamInterp)
            omegaBG=integBG.integrate(beamInterp)
            omegaTot=integTot.integrate(beamInterp)	
            #if verbose:print 'aperture/background beam areas:',omegaAp,omegaBG
            #calculate apertur  ecorrections
            apCorrNoBG[band] = \
                srcConvArea/omegaAp
            apCorrIncBG[band] =  \
                srcConvArea/(omegaAp - omegaBG*sizeAp/sizeBG)
        if (verbose):
            print 'Calculated apCorr (noBG) for beta=%f, temp=%f: '%(betaK,tempK),apCorrNoBG["PSW"],apCorrNoBG["PMW"],apCorrNoBG["PLW"]
            print 'Calculated apCorr (incBG) for beta=%f, temp=%f: '%(betaK,tempK),apCorrIncBG["PSW"],apCorrIncBG["PMW"],apCorrIncBG["PLW"]
    else:
        #alphaK is scalar
        apCorrNoBG={'PSW': Double1d(nt,Double.NaN), 'PMW': Double1d(nt,Double.NaN), 'PLW': Double1d(nt,Double.NaN)}
        apCorrIncBG={'PSW': Double1d(nt,Double.NaN), 'PMW': Double1d(nt,Double.NaN), 'PLW': Double1d(nt,Double.NaN)}
        
        for t in range(nt):
            #calculate effective beam
            effBeamsSrc=calcEffBeamSrc_BB(betaK,tempK[t],src,verbose=True)
            
            for band in spireBands():
                #create integrators
                integAp=TrapezoidalIntegrator(0.,apPhotRad[band])
                integBG=TrapezoidalIntegrator(apPhotBGRad[band][0],apPhotBGRad[band][1])
                #integrate to find size of aperture/annulus
                sizeAp=integAp.integrate(sizeInterp)
                sizeBG=integBG.integrate(sizeInterp)
                
                srcConv=effBeamsSrc[band]
                srcConvArea=srcConv.calcArea()
                
                #interpolate convolved profile
                beamInterp = CubicSplineInterpolator(srcConv.radArr,srcConv.profile * 2.*PI*srcConv.radArr)
            
                #integrate profile to find beam area in aperture/annulus
                omegaAp=integAp.integrate(beamInterp)
                omegaBG=integBG.integrate(beamInterp)
                omegaTot=integTot.integrate(beamInterp)
                #calculate aperture corrections
                apCorrNoBG[band][t] = \
                    srcConvArea/omegaAp    
                apCorrIncBG[band][t] =  \
                    srcConvArea/(omegaAp - omegaBG*sizeAp/sizeBG)
            if (verbose):
                print 'Calculated apCorr (noBG) for beta=%f, temp=%f: '%(betaK,tempK[t]),apCorrNoBG["PSW"][t],apCorrNoBG["PMW"][t],apCorrNoBG["PLW"][t]
                print 'Calculated apCorr (incBG) for beta=%f, temp=%f: '%(betaK,tempK[t]),apCorrIncBG["PSW"][t],apCorrIncBG["PMW"][t],apCorrIncBG["PLW"][t]
        
    if not table:
        #return as is
        return(apCorrIncBG,apCorrNoBG)
    else:
        #create and returnTableDataset
        apCorrNoBG_table=TableDataset()
        apCorrIncBG_table=TableDataset()
        apCorrNoBG_table.setDescription("Aperture Correction without background annulus (modified black body, beta=%.2f)"%betaK)
        apCorrIncBG_table.setDescription("Aperture Correction including background annulus (modified black body, beta=%.2f)"%betaK)
        apCorrNoBG_table.addColumn("temperature",Column(Double1d(tempK)))
        apCorrIncBG_table.addColumn("temperature",Column(Double1d(tempK)))
        for band in spireBands():
            apCorrNoBG_table.addColumn(band,Column(apCorrNoBG[band]))
            apCorrIncBG_table.addColumn(band,Column(apCorrIncBG[band]))
        return (apCorrIncBG_table,apCorrNoBG_table)
    
def clear(verbose=False):
    #sets global variables to NoneType
    varList=['beamMonoSrcArea','effBeamSrcProfs']
    delList=[]
    noDelList=[]
    for var in varList:
        exec('global %s'%var)
        try:
            exec('%s=None'%var)
            delList.append(var)
        except:
            noDelList.append(var)
        exec('%s=None'%var)
    if verbose:
        if len(delList)>0:
            print 'Set Global Variables to None:',delList
        if len(noDelList)>0:
            print 'Unable to set Global Variables to None:',noDelList
            
def saveMonoBeams(directory=None,verbose=False):
    #save global variables to file
    global beamMonoSrcArea
    if directory==None:
        directory=Configuration.getProperty('var.hcss.workdir')
    fileMono='SemiExtendedBundle_SAVED_beamMonoSrcArea.csv'
    if beamMonoSrcArea!=None:
        keyList=[]
        fileList=[]
        for key in beamMonoSrcArea:
            keyList.append(key)
            fileKey='%s.fits'%key
            fileList.append(fileKey)
            beamTable=TableDataset()
            beamTable.addColumn('freq',Column(Double1d(sh.getSpireFreq())/1.e9,description='Frequency (GHz)'))
            for band in spireBands():
                beamTable.addColumn(band,Column(Double1d(beamMonoSrcArea[key][band]),description='%s monochromatic area'))
            simpleFitsWriter(product=beamTable,file=os.path.join(directory,fileKey))
        fileTable=TableDataset()
        fileTable.addColumn('key',Column(String1d(keyList)))
        fileTable.addColumn('filename',Column(String1d(fileList)))
        asciiTableWriter(table=fileTable,file=os.path.join(directory,fileMono))
        if verbose:
            print 'Saved beamMonoSrcArea to %s directory. Files listed in %s.'%(directory,fileMono)
    
def loadMonoBeams(directory=None,verbose=False):
    #load beamMonoSrcArea global variables from file
    global beamMonoSrcArea
    try:
        beamMonoSrcArea
    except:
        #make variable
        beamMonoSrcArea={}
    if directory==None:
        directory=Configuration.getProperty('var.hcss.workdir')
    fileMono='SemiExtendedBundle_SAVED_beamMonoSrcArea.csv'
    fileTable=asciiTableReader(file=os.path.join(directory,fileMono))
    keyList=fileTable['key'].data
    fileList=fileTable['filename'].data
    nKey=len(keyList)
    for n in range(nKey):
        key=keyList[n]
        table=simpleFitsReader(file=os.path.join(directory,fileList[n]))
        print table
        print key
        try:
            beamMonoSrcArea[key]
            if verbose:print 'Overwriting %s'%key
            beamMonoSrcArea[key]={'PSW':Double.NaN,'PMW':Double.NaN,'PLW':Double.NaN}
        except:
            if verbose:print 'Reading %s'%key
            beamMonoSrcArea[key]={'PSW':Double.NaN,'PMW':Double.NaN,'PLW':Double.NaN}
        for band in spireBands():
            beamMonoSrcArea[key][band]=table[band].data
    if verbose:print 'BeamMonoSrcArea read from files in %s'%directory
        
def saveEffBeams(directory=None,verbose=False):
    global effBeamSrcProfs
    if directory==None:
        directory=Configuration.getProperty('var.hcss.workdir')
    fileEff='SemiExtendedBundle_SAVED_effBeamSrcProfs.csv'
    if effBeamSrcProfs!=None:
        keyList=[]
        fileList={'PSW':[],'PMW':[],'PLW':[]}
        for key in effBeamSrcProfs:
            keyList.append(key)
            for band in spireBands():
                fileKey='%s_%s.fits'%(key,band)
                fileList[band].append(fileKey)
                effBeamSrcProfs[key][band].saveFits(directory=directory,filename=fileKey,verbose=verbose)
        fileTable=TableDataset()
        fileTable.addColumn('key',Column(String1d(keyList)))
        for band in spireBands():
            fileTable.addColumn(band,Column(String1d(fileList[band]),description='%s files'%band))
        asciiTableWriter(table=fileTable,file=os.path.join(directory,fileEff))
        if verbose:
            print 'Saved effBeamSrcProfs to %s directory. Files listed in %s.'%(directory,fileEff)
        
def loadEffBeams(directory=None,verbose=False):
    #load beamMonoSrcArea global variables from file
    global effBeamSrcProfs
    try:
        effBeamSrcProfs
    except:
        #make variable
        effBeamSrcProfs={}
    if directory==None:
        directory=Configuration.getProperty('var.hcss.workdir')
    fileEff='SemiExtendedBundle_SAVED_effBeamSrcProfs.csv'
    fileTable=asciiTableReader(file=os.path.join(directory,fileEff))
    keyList=fileTable['key'].data
    fileList={}
    for band in spireBands():
        fileList[band]=fileTable[band].data
    nKey=len(keyList)
    for n in range(nKey):
        key=keyList[n]
        print key
        try:
            effBeamSrcProfs[key]
            if verbose:print 'Overwriting %s'%key
            effBeamsSrcProfs[key]={'PSW':Double.NaN,'PMW':Double.NaN,'PLW':Double.NaN}
        except:
            effBeamSrcProfs[key]={'PSW':Double.NaN,'PMW':Double.NaN,'PLW':Double.NaN}
        for band in spireBands():
            table=simpleFitsReader(file=os.path.join(directory,fileList[band][n]))
            effBeamSrcProfs[key][band]=srcMod.SourceProfile().loadFits(directory=directory,filename=fileList[band][n],verbose=verbose)
    if verbose:print 'effBeamSrcProfs read from files in %s'%directory
    
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

    #create arrays for each band, rather than each width
    KColE_partialPSW=[]
    KColE_partialPMW=[]
    KColE_partialPLW=[]
    KColP_partialPSW=[]
    KColP_partialPMW=[]
    KColP_partialPLW=[]
    
    for wid in srcWidths:
        src=srcMod.Source('Gaussian',[wid])
        #calculate colour corrections from source
        KColP_partial=calcKColPSrc(alphaK,src,verbose=True,table=False)
        KColE_partial=calcKColESrc(alphaK,src,verbose=True,table=False)

        #calculate colour corrections from source profile
        srcProf=src.calcProfile(Double1d(range(700)))
        KColP_partial=calcKColPSrc(alphaK,srcProf,verbose=True,table=False)
        KColE_partial=calcKColESrc(alphaK,srcProf,verbose=True,table=False)

        KColP_partialPSW.append(KColP_partial['PSW'])
        KColP_partialPMW.append(KColP_partial['PMW'])
        KColP_partialPLW.append(KColP_partial['PLW'])
        KColE_partialPSW.append(KColE_partial['PSW'])
        KColE_partialPMW.append(KColE_partial['PMW'])
        KColE_partialPLW.append(KColE_partial['PLW'])
        
    #combine into lists for point (surface brightness) and extended (surface brightness)
    KColPSrc=[Float1d(KColP_partialPSW),Float1d(KColP_partialPMW),Float1d(KColP_partialPLW)]
    KColESrc=[Float1d(KColE_partialPSW),Float1d(KColE_partialPMW),Float1d(KColE_partialPLW)]
    return(srcWidths,KColPSrc,KColESrc)

def testApCorr():
    alphaK=[2.0,3.0]
    betaK=1.75
    tempK=[20.,30.,40.]
    srcWidths=Float1d([0.2,0.5,1.,2.,5.,10.,20.,50.,100.,200.,500.,1000.,2000.,5000.,10000.])
    #srcWidths=Float1d([0.2,0.5,1.,10.,100.])
    #srcWidths=Float1d([10,20])
    apCorrNoBG_partialPSW=[]
    apCorrNoBG_partialPMW=[]
    apCorrNoBG_partialPLW=[]
    apCorrIncBG_partialPSW=[]
    apCorrIncBG_partialPMW=[]
    apCorrIncBG_partialPLW=[]
    for wid in srcWidths:
        src=srcMod.Source('Gaussian',[wid])
        
        apCorrBoth=calcApCorrSrc(alphaK,src,verbose=True,table=False)
        apCorrBoth_BB=calcApCorrSrc_BB(betaK,tempK,src,verbose=True,table=False)
        apCorrIncBG_partialPSW.append(apCorrBoth[0]["PSW"])
        apCorrIncBG_partialPMW.append(apCorrBoth[0]["PMW"])
        apCorrIncBG_partialPLW.append(apCorrBoth[0]["PLW"])
        apCorrNoBG_partialPSW.append(apCorrBoth[1]["PSW"])
        apCorrNoBG_partialPMW.append(apCorrBoth[1]["PMW"])
        apCorrNoBG_partialPLW.append(apCorrBoth[1]["PLW"])
        
    apCorrIncBGSrc=[Float1d(apCorrIncBG_partialPSW),Float1d(apCorrIncBG_partialPMW),Float1d(apCorrIncBG_partialPLW)]
    apCorrNoBGSrc=[Float1d(apCorrNoBG_partialPSW),Float1d(apCorrNoBG_partialPMW),Float1d(apCorrNoBG_partialPLW)]
    
    return(srcWidths,apCorrIncBGSrc,apCorrNoBGSrc)

    