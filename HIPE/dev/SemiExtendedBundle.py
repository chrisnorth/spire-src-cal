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
#  Herschel-SPIRE Colour Correction factors 
# 
#  This routine calculates colour correction tables to convert from the standard
#  pipeline flux densities (which are quoted for a nu*F_nu=const. spectrum) into 
#  flux densities for a range of source spectra. These produce monochromatic flux
#  densities at the SPIRE reference wavelengths of 250, 350 and 500 micron. The
#  source spectra include power law spectra, and modified black body spectra with
#  a range of emissivities. Metadata values are produced for the flux conversion
#  parameters applied to the pipeline (the "K4P" and "K4E" parameters).
#
#  Input:
#    RSRF and Aperture efficiency profiles from SPIRE calibration tree
#    Beam profiles and beam model parameters from SPIRE calibration tree
#    Flux conversion parameters from SPIRE calibration tree
#    Name and version of output file
#
#  Output:
#    SCalPhotColorCorrK_point product
#    SCalPhotColorCorrK_extended product
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
#   E. Polehampton   22-10-2013  - First version adapted from Andreas' script - SPCAL-83
#   E. Polehampton   31-10-2013  - update for new input file
#   E. Polehampton   21-01-2014  - add table descriptions and update numbers (SPCAL-93)
#   Chris North      18-02-2014  - updated to use full calculation instead of reading csv files
#   Ivan Valtchanov  15-03-2014  - reformatted to functions and a script bundle to distribute with the handbook
#   Chris North      03-04-2014  - reformatted functions to remove dependence on global variables
#
#===============================================================================
import os
scriptVersionString = "SemiExtendedBundle.py $Revision: 1.0 $"

#-------------------------------------------------------------------------------
# Loading physical and math constants
from herschel.share.unit import *
import sys
print sys.path
execfile('SpireHandbookBundle_CEN.py')
#-------------------------------------------------------------------------------
# Calculate monochromatic beam profile and area at a given frequency (Eq. 5.32) 
def spireMonoBeamSrc(freqx,beamRad,beamProfs,beamConst,effFreq,gamma,srcProf,array):
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
      beamConst: (array float) constant part of beam profile for "array"
                       [passed to prevent repeated calls]
      effFreq:   (float) effective frequency [Hz] of array
      gamma:     (float) Exponent of powerlaw describing FWHM dependence
                   on frequency
      srcProf:   (array float) source density profile, corresponding to
                   radius column in beamProfs
      array:     (string) spire array ('Psw'|'Pmw'|'Plw')

    Outputs:     (list of objects)
                (float) Beam area [arcsec^2] at frequency freqx

    Calculation:
      Scales the core beam profile width as (freqx/effFreq)^gamma.
      Queries the calibration file to generate new core beam profile.
      Uses constant beam profile where it is larger than core profile.
      Multiplies by a source profile
      Integrates over radius to calculate beam area.

    Dependencies:
      herschel.ia.numeric.toolbox.interp.LinearInterpolator
      herschel.ia.numeric.toolbox.integr.TrapezoidalIntegrator
      
    2013/12/19  C. North  initial version

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
    #apply the "constant" beam where appropriate
    #beamConst=beamProfs.getConstantCorrectionTable().getColumn(array).data
    isConst=beamNew.where(beamNew < beamConst)
    beamNew[isConst]=beamConst[isConst]

    #multiply by source Profile
    beamNew = beamNew*srcProf

    ## THIS IS ONLY VALID FOR AREA
    ## FULL PROFILE OF CONVOLUTION REQUIRES PROPER CONVOLUTION

    #integrate to get solid angle (in arcsec^2)    
    beamInterp=LinearInterpolator(beamRad,beamNew * 2. * Math.PI * beamRad)
    integrator=TrapezoidalIntegrator(0,maxRad)
    beamMonoArea=integrator.integrate(beamInterp)

    return (beamMonoArea)

#-------------------------------------------------------------------------------
# Calculate monochromatic beam areas over a range of frequencies
def spireMonoSrcAreas(freq,beamProfs,effFreq,gamma,srcProf,array,freqFact=100):

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
      srcProf:   (array float) source density profile, corresponding to
                   radius column in beamProfs
      freqFact:  (int) Factor by which to reduce size of freq.
                   OPTIONAL. Default=100.

    Outputs:     
                 (array float) Monochromatic Beam area [sr] at frequencies
                    corresponding to freq

    Calculation:
      Generates sparse frequency array of full range
      Uses spireMonoBeam to calculate monochromatic area of beam convolved with 
        a source profile at sparse freqs
      Interpolates to full frequency grid

    Dependencies:
      spireMonoBeam
      herschel.ia.numeric.toolbox.interp.CubicSplineInterpolator
      
    2013/12/19  C. North  initial version

    """

    arcsec2Sr = (Math.PI/(60.*60.*180))**2

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
    #get constant beam profile from calibration table
    beamConst=beamProfs.getConstantCorrectionTable().getColumn(array).data

    # calculate at sparse frequencies
    for fx in range(nNuArea):
        #get corresponding index in full frequency array
        f=iNuArea[fx]
        #populate frequency array
        beamMonoFreqSparse[fx]=freq[f]
        #populate beam area array
        beamMonoAreaSparse[fx]=spireMonoBeamSrc(freq[f],beamRad,beamProfs,beamConst,effFreq,gamma,srcProf,array)

    # interpolate to full frequency array and convert to Sr
    beamInterp=CubicSplineInterpolator(beamMonoFreqSparse,beamMonoAreaSparse)
    beamMonoArea=beamInterp(freq)*arcsec2Sr #in sr
    
    return(beamMonoArea)
    
#-------------------------------------------------------------------------------
#===============================================================================
#=====                  CALCULATE MONOCHROMATIC BEAM AREAS                 =====
#===============================================================================
#-------------------------------------------------------------------------------

def calcBeamSrcMonoAreaGauss(srcWidth,verbose=False):
    #calculate monochromatic beam areas using full or simple beam treatment
    #print '\nCalculating monochromatic beam areas...'

    try:
        arcsec2Sr
    except:
        #define arcsec2Sr
        global arcsec2Sr
        arcsec2Sr = (Math.PI/(60.*60.*180))**2
    
    beamProfs = spireCalPhot.getProduct("RadialCorrBeam")
    beamRad=beamProfs.getCoreCorrectionTable().getColumn('radius').data
    #define srcProf to be same length as beamRad

    srcProf=exp(-beamRad**2/(2.*srcWidth**2))

    gamma = beamProfs.meta['gamma'].double
    beamMonoArea = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
    for band in spireBands:
        if (verbose):print 'Calculating monochromatic beam areas for %.2f arcsec Gaussian source for %s band'%(srcWidth,band)
        #monochromatic beam areas
        beamMonoArea[band] = spireMonoSrcAreas(getSpireFreq(), beamProfs, 
          getSpireEffFreq()[band], gamma, srcProf, band)

    return beamMonoArea


