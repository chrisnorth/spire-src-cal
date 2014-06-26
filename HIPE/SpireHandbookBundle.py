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
#  Herschel-SPIRE Colour Correction factor calculations
# 
#  This routine calculates colour correction tables to convert from the standard
#  pipeline flux densities (which are quoted for a nu*F_nu=const. spectrum) into 
#  flux densities for a range of source spectra. These produce monochromatic flux
#  densities at the SPIRE reference wavelengths of 250, 350 and 500 micron. The
#  source spectra include power law spectra, and modified black body spectra with
#  a range of emissivities. Metadata values are produced for the flux conversion
#  parameters applied to the pipeline (the "K4P" and "K4E" parameters).
#
#  Note that these functions are used to calculate the values in the Handbook.
#
#  Input:
#    SPIRE Calibration Tree, or tree name to import
#    Requires spire_cal_12_2 or later
#
#  Output:
#    Functions which can be used to calculate colour correction parameters
#    Some of them also generate global variables
#
#  Usage:
#    A test script is provided at the end for usage of the functions.
#    See how-to-use_dev.py for examples with more detailed comments.
#
#===============================================================================
# 
#  Herschel-SPIRE Colour Corrections
# 
#  This module provides functions to calculate the colour corrections for a
#  range of source spectra, and for both point sources and extended sources
# 
#  The following global variables are defined and used:
#         spireCalPhot:  [context] Spire Calibration Context for Photometer
#         spireBeamRad:  [float array] radius array used for Spire beams
#         spireFreq:     [float array] frequency raster across all bands
#         spireRefFreq:  [dict] SPIRE reference frequencies (1 scalar per band)
#         spireEffFreq:  [dict] SPIRE effective frequencies (1 scalar per band)
#         spireFiltOnly: [dict] filter profile without aperture efficiency (1 array per band)
#         spireFilt:     [dict] filter profile with aperture efficiency (1 array per band)
#         beamMonoArea:  [dict] SPIRE monochromatic areas (1 array per bands, plus 'beamType')
#         spireEffBeams: [dict] SPIRE effective beam profiles (1 per spectrum)
#
#  The functions are as follows:
#  * spireBands:
#     Defines names for SPIRE bands
#     - Inputs:
#         NONE
#     - Outputs:
#         [string list] SPIRE band names ["PSW","PMW","PLW"]
#
#  * arcsec2Sr:
#     Defines conversion from arcsec^2 so sr
#     - Inputs:
#         NONE
#     - Outputs:
#         [float] conversion factor
#
#  * getCal:
#     Gets calibration Context from pool or file, and/or checks existing cal.
#     If no inputs, reads from default pool
#     - Inputs:
#         cal:     [SpireCal context] Calibration Tree (optional)
#         calTree: [string] Name of calibration tree to read from HSA (optional)
#         calPool: [string] Name of calibration pool (optional)
#         calFile: [string] Name of calibration file (optional)
#         verbose: [boolean] Set to print more information to terminal
#     - Outputs:
#         cal:     [SpireCal context] Spire Calibration Context for Photometer
#     - Global Variables used:
#         spireCalPhot: [context] Spire Calibration Context for Photometer
#
#  * getBeamRad:
#     Gets beam radius array.
#     - Inputs:
#         NONE
#     - Outputs:
#         [float array] radius array used for spire beams
#     - Global variables used:
#         spireBeamRad [float array] radius array used for spire beams
#
#  * getSpireFreq:
#     Gets frequency raster.
#     - Inputs:
#         NONE
#     - Outputs:
#         [float array] frequency raster across all bands
#     - Global variables used:
#         spireFreq [float array] frequency raster across all bands
#
#  * getSpireRefFreq:
#     - Inputs:
#         NONE
#     - Outputs:
#         [float dict] SPIRE reference frequencies (at 250,350,500 um)
#     - Global variables used:
#         spireRefFreq [dict] SPIRE reference frequencies (1 scalar per band)
#
#  * getSpireFilt:
#     Gets SPIRE filter profiles, either RSRF only or RSRF*ApEff
#     - Inputs:
#         rsrfOnly: [boolean] set to only outut RSRF (without Aperture Efficiency). Default=False
#     - Outputs:
#         [float dict] filter profiles with/without aperture efficiency (1 per band)
#     - Global Variables:
#         spireFiltOnly [dict] filter profile without aperture efficiency (1 array per band)
#         spireFilt     [dict] filter profile with aperture efficiency (1 array per band)
#
#  * getSpireEffFreq:
#     Gets SPIRE effective frequencies from calibration tree
#     - Inputs:
#         NONE
#     - Outputs:
#         [dict] SPIRE effective frequencies for three bands (1 scalar per band)
#     - Global Variables used:
#         spireEffFreq: [dict] SPIRE effective frequencies (1 scalar per band)
#
#   * calcBeamMonoArea:
#      Calulates monochromatic beam areas (used for many functions).
#      Calculated for frequencies produced by getFreq
#      - Inputs:
#          beamType: [string] beam type ('Full'|'Simple') to calculate. Default=None
#              If not specified, existing beamType is used
#          verbose:  [boolean] Set to print more information to terminal
#      - Outputs:
#          [dict] SPIRE monochromatic areas (1 array per band + beamType)
#      - Global variables used:
#          beamMonoArea [dict] SPIRE monochromatic areas (1 array per band + beamType)
#
#   * calcSpireEffBeam:
#      Calculates effective beam profile for power law spectrum.
#      - Inputs:
#          alphaK:  [float] power law spectral index (scalar or list)
#          verbose: [boolean] Set to print more information to terminal
#      - Outputs:
#          [dict] SPIRE effective beam profile (1 Double1d per band)
#      - Global variables used:
#          spireEffBeams [dict] SPIRE monochromatic areas (1 per spectrum)
#
#   * calcSpireEffBeam_BB:
#      Calculates effective beam profile for modified black body spectrum.
#      - Inputs:
#          betaK:   [float] modBB emissivity index
#          tempK:   [float] modBB temperature
#          verbose: [boolean] Set to print more information to terminal
#      - Outputs:
#          [dict] SPIRE effective beam profiles (1 Double1d per band)
#      - Global variables used:
#          spireEffBeams [dict] SPIRE monochromatic areas (1 per spectrum)
#
#   * calcOmegaEff:
#      Calculates effective beam area for power law spectrum.
#      - Inputs:
#          alphaK: [float/list] power law spectral index (scalar or list)
#          verbose: [boolean] Set to print more information to terminal
#          table:   [boolean] Set to output TableDataset
#      - Outputs:
#          [dict] SPIRE effective beam areas
#              (if alphaK is list, 1 list per band, otherwise 1 scalar per band)
#              (if 'table' is set, outputs as TableDataset)
#
#   * calcOmegaEff_BB
#      Calculates effective beam area for modified blackbody spectrum.
#      - Inputs:
#          betaK: [float] modBB emissivity index
#          tempK: [float/list] modBB temperature (scalar or list)
#          verbose: [boolean] Set to print more information to terminal
#          table:   [boolean] Set to output TableDataset
#      - Outputs:
#          [dict] SPIRE effective beam areas
#              (if tempK is list, 1 list per band, otherwise 1 scalar per band)
#              (if 'table' is set, outputs as TableDataset)
#
#   * calcKBeam
#      Calculates beam correction factor for power law spectrum
#      - Inputs:
#          alphaK: [float/list] power law spectral index (scalar or list)
#          verbose: [boolean] Set to print more information to terminal
#          table:   [boolean] Set to output TableDataset
#      - Outputs:
#          [dict] beam correction factors
#              (if alphaK is list, 1 list per band, otherwise 1 scalar ber band)
#              (if 'table' is set, outputs as TableDataset)
#
#   * calcKBeam_BB
#      Calculates beam correction factor for modified blackbody spectrum.
#      - Inputs:
#          betaK: [float] modBB emissivity index
#          tempK: [float/list] modBB temperature (scalar or list)
#          verbose: [boolean] Set to print more information to terminal
#          table:   [boolean] Set to output TableDataset
#      - Outputs:
#          [dict] beam correction factors
#              (if tempK is list, 1 list per band, otherwise 1 scalar per band)
#              (if 'table' is set, outputs as TableDataset)
#
#   * calcK4P
#      Calculates K4P calibration parameter (point source flux density, alpha=-1)
#      - Inputs:
#          NONE
#      - Outputs:
#          [dict] K4P parameter (1 per band)
#
#   * calcKMonE
#      Calculates KMonE calibration parameter (extended source surface brightness, alpha=-1)
#      - Inputs:
#          NONE
#      - Outputs:
#          [dict] KMonE parameter for SPIRE bands (1 per band)
#
#   * calcK4E
#      Calculates K4E calibration parameter (extended source flux density, alpha=-1)
#      - Inputs:
#          NONE
#      - Outputs:
#          [dict] K4E parameter for SPIRE bands (1 per band)
#
#   * calcKPtoE
#      Calculates KPtoE calibration parameter (point flux density -> extended surface brightness, alpha=-1)
#      - Inputs:
#          NONE
#      - Outputs:
#          [dict] KPtoE parameter for SPIRE bands (1 per band)
#
#   * calcKColP
#      Calculates KColP colour correction parameter for power law spectrum
#      - Inputs:
#          alphaK:  [float/list] power law spectral index (scalar or list)
#          verbose: [boolean] Set to print more information to terminal
#          table:   [boolean] Set to output TableDataset
#      - Outputs:
#          [dict|TableDataset] KColP colour corrections
#              (if alphaK is list, 1 list per band, otherwise 1 scalar per band)
#              (if 'table' is set, outputs as TableDataset)
#
#   * calcKColP_BB
#      Calculates KColP colour correction parameter for modified blackbody spectrum
#      - Inputs:
#          betaK:   [float] modified blackbody emissivity index
#          tempK:   [float/list] modified black body temperature (scalar or list)
#          verbose: [boolean] Set to print more information to terminal
#          table:   [boolean] Set to output TableDataset
#      - Outputs:
#          [dict|TableDataset] KColP colour correction for SPIRE bands
#              (if tempK is list, 1 list per band, otherwise 1 scalar per band)
#              (if 'table' is set, outputs as TableDataset)
#
#   * calcKColE
#      Calculates KColE colour correction parameter power law spectrum
#      - Inputs:
#          alphaK:  [float/list] power law spectral index (scalar or list)
#          verbose: [boolean] Set to print more information to terminal
#          table:   [boolean] Set to output TableDataset
#      - Outputs:
#          [dict|TableDataset] KColE colour correction for SPIRE bands
#              (if alphaK is list, 1 list per band, otherwise 1 scalar per band)
#              (if 'table' is set, outputs as TableDataset)
#
#   * calcKColE_BB
#      Calculates KColE colour correction parameter for modified blackbody spectrum
#      - Inputs:
#          betaK:   [float] modified blackbody emissivity index
#          tempK:   [float/list] modified black body temperature (scalar or list)
#          verbose: [boolean] Set to print more information to terminal
#          table:   [boolean] Set to output TableDataset
#      - Outputs:
#          [dict|TableDataset] KColE colour correction for SPIRE bands
#              (if tempK is list, 1 list per band, otherwise 1 scalar per band)
#              (if 'table' is set, outputs as TableDataset)
#  
#  Other functions for dealing with the full beam model and colour corrections are also included@
#   * spireEffArea: Calculate effective RSRF-weighted beam area for a given spectrum (power law or modified black-body)
#   * spireMonoBeam: Calculate monochromatic beam profile & area at a given freq
#   * spireMonoAreas: Calculate monochromatic beam areas at range of frequencies
#   * calcSpireKcorr: Calculate K-correction parameters for given spectrum (power law or modified black body)
#
#  Example test functions are provided:
#   * spireColCorrTest
#      Tests the above functions
#      - Inputs:
#          beamType: [string] beamType to use ('Full'|'Simple'). Default=Full.
#      - Outputs:
#          [dict] SPIRE monochromatic areas (1 array per band + beamType)
#          [float array] values of alpha used
#          [float array] calculated values of OmegaEff for alpha
#          [float array] calculated values of KBeam for alpha
#          [float array] calculated values of KColP for alpha
#          [float array] calculated values of KColE for alpha
#   * compFullSimple
#      Tests the calculations for full and simple beams and plots results and comparisons
#      - Inputs:
#          NONE
#      - Outputs:
#          Plots
#===============================================================================
# 
#  Edition History
#   E. Polehampton   22-10-2013  - First version adapted from Andreas' script - SPCAL-83
#   E. Polehampton   31-10-2013  - update for new input file
#   E. Polehampton   21-01-2014  - add table descriptions and update numbers (SPCAL-93)
#   Chris North      18-02-2014  - updated to use full calculation instead of reading csv files
#   Ivan Valtchanov  15-03-2014  - reformatted to functions and a script bundle to distribute with the handbook
#   Chris North      03-04-2014  - reformatted functions to remove dependence on global variables
#   Chris North      08-04-2014  - added global variables back in(!)
#                                  provided support for "Simple" beam model
#                                  added functionality to test scripts
#                                  added import statements so module can be imported
#   Chris North      11-04-2014  - added support for TableDataset output
#   Chris North      24-04-2014  - added aperture correction calculations and more help
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
from herschel.ia.dataset import TableDataset,DoubleParameter,Column
from herschel.ia.numeric.toolbox.interp import LinearInterpolator,CubicSplineInterpolator
from herschel.ia.numeric.toolbox.integr import TrapezoidalIntegrator
from java.lang.Math import PI
from java.lang import Double
from herschel.share.unit import Constant,Temperature,SolidAngle
from herschel.ia.gui.plot import *
import java.awt.Color
from herschel.ia.numeric.toolbox.basic import Floor,Min,Max,Exp
FLOOR=herschel.ia.numeric.toolbox.basic.Floor.PROCEDURE
EXP=herschel.ia.numeric.toolbox.basic.Exp.PROCEDURE
MAX=herschel.ia.numeric.toolbox.basic.Max.FOLDR
MIN=herschel.ia.numeric.toolbox.basic.Min.FOLDR

#-------------------------------------------------------------------------------
# Loading physical and math constants

#scriptVersionString = "SpireHandbookBundle.py $Revision: 1.0 $"

#-------------------------------------------------------------------------------
#===============================================================================
#=====                       DEFINE GLOBAL VARIABLES                       =====
#===============================================================================
#-------------------------------------------------------------------------------

global beamMonoArea,spireEffBeams
## SPIRE monochromatic areas (1 array per bands, plus 'beamType')
#beamMonoArea={}
## SPIRE effective beam peofiles (1 per spectrum)
#spireEffBeams={}

#-------------------------------------------------------------------------------
#===============================================================================
#=====                           SETUP FUNCTIONS                           =====
#===============================================================================
#-------------------------------------------------------------------------------

def spireBands():
    """
    ========================================================================
    spireBands():
        Return list of bands for addessing dicts: ["PSW","PMW",PLW"]
        
    Inputs:
        NONE
        
    Outputs:
        (string list) list of beam names
        
    Calculation:
        ['PSW','PMW','PLW']
    """
    return(['PSW','PMW','PLW'])

def arcsec2Sr():
    """
    ========================================================================
    arcsec2Sr():
        Calculate conversion from square arcsec to steradians (PI/(180*3600))^2
    Inputs:
        NONE
    Outputs:
        (double) conversion from square arcsec to steradians
        
    Calculation:
        conversion is (PI/(180*3600))^2
    """
    return((PI/(60.*60.*180))**2)
    
def getCal(cal=None,calTree=None,calPool=None,calFile=None,verbose=False):
    """
    ========================================================================
    getCal(cal=None,calTree=None,calPool=None,calFile=None,verbose=False):
        Read photometer Calibration tree from appropriate input
        Stored in global variable spireCalPhot so only read once

    Inputs:
        cal:     (SpireCal Context) calibration tree to extract photometer tree
                 from. Default=None
        calTree: (string) name of cal tree to read from HSA. Default=None
        calPool: (string) name of pool to read cal tree from. Default=None
        calFile: (string) filename to read cal tree from. Default=None
        verbose: (boolean) Print extra information. Default=False
    Outputs:
        (PhotCal Context) reference frequencies in Hz (one per band)
                    
    Calculation:
        If spireCalPhot already defined, check it is valid and return
        If not, and no inputs defined, try to read from local pool
        If that fails, try in the following order:
            If cal set, extract PhotCalContext from there
            If calTree set, read from HSA calibration tree
            If calPool set, read from local pool
            If calFile set, read from jar file
        Check spireCalPhot is a valie PhotCalContext
    
    Dependencies:
        herschel.spire.ia.cal.SpireCalTask
    """
    # if cal not defined, read from pool or jarFile
    # if cal is defined, don't read anything new. Just check and return cal
    #print cal,calTree,calPool,calFile,verbose
    
    #set global variables
    global spireCalPhot #Photometer Calibration Context (spireCal.getPhot())

    try:
        spireCalPhot
        #global variable already defined, so do nothing
    except:
        #get spireCalTree
        if cal==None and calTree==None and calPool==None and calFile==None:
            #no import provided. do default action
            if (verbose):print 'Reading from default pool'
            cal=spireCal()
            if (verbose):print 'Calibration read from default pool: %s'%cal.meta["version"].value
        try:
            #try getting from cal
            spireCalPhot=cal.getPhot()
        except:
            if cal==None and calTree!=None:
                #try to read from HSA with calTree
                cal=spireCal(calTree=calTree)
                #try:
                #    cal=spireCal(calTree=calTree)
                #except:
                #    if verbose: print 'unable to read from HSA'
            if cal==None and calPool!=None:
                #try to read from local pool
                cal=spireCal(pool=calPool)
                #try:
                #    cal=spireCal(pool=calPool)
                #except:
                #    if verbose: print 'unable to read from local pool'
            if cal==None and calFile!=None:
                #try to read from jarFile
                cal=spireCal(jarfile=calFile)
                #try:
                #    cal=spireCal(jarfile=calFile)
                #except:
                #    if verbose: print 'unable to read from jar file'
            assert cal!=None,\
                'ERROR: Unable to read calibration context'

            spireCalPhot=cal.getPhot()
        
        assert spireCalPhot.isValid(),'ERROR: Invalid SPIRE Photometer calibration tree'
    
    return(spireCalPhot)

def getBeamRad():
    """
    ========================================================================
    getBeamRad():
        Get beam radius array used for beam profiles
        Stored in global variable spireBeamRad so only calculated once

    Inputs:
      NONE
    Outputs:
        (double array) beam radius array
                    
    Calculation:
        Retrieves radius array from RadialCorrBeam calibration product
    
    Dependencies:
        getCal() - retrieve calibration tree
    """
    global spireBeamRad
    try:
        spireBeamRad
    except:
        beamProfs=getCal().getProduct('RadialCorrBeam')
        spireBeamRad=beamProfs.getCoreCorrectionTable().getColumn('radius').data
    return(spireBeamRad)
            
def getSpireFreq():
    """
    ========================================================================
    getSpireFreq():
        Calculate spire frequencies raster covering all bands
        Stored in global variable spireFreq so only calculated once

    Inputs:
      NONE
    Outputs:
        (double array) reference frequencies in Hz (one per band)
                    
    Calculation:
        Calculates frequency raster in Hz from 300-1800 GHz, with 0.1GHz step
        That range covers all three SPIRE bands, and the 545/857 GHz Planck bands
    
    Dependencies:
        herschel.ia.numeric.toolbox.basic.Floor
    """
    # Frequency raster of common frequency grid spanning all three spire bands

    #define global variable
    global spireFreq #raster of frequencies used throughout functions
    try:
        spireFreq
        #global variable already defined, so do nothing
    except:
        #not defined, so recalculate
        deltaNu = 1.e9        # 0.1 GHz
        nuMin   = 300.e9
        nuMax   = 1800.e9
        nNu     = FLOOR((nuMax-nuMin)/deltaNu)
        spireFreq    = Double1d(range(nNu)) * deltaNu + nuMin
        
    return(spireFreq)
    
def getSpireRefFreq(array=None,verbose=False):
    """
    ========================================================================
    getSpireRefFreq():
        Calculate spire reference frequencies for each band
        Stored in global variable spireRefFreq so only calculated once

    Inputs:
      array:      (string) OPTIONAL. SPIRE array to use ('PSW'|'PMW'|'PLW')
                     Default=None (returns all three bands)
      verbose:    (boolean) print more information. Default=False

    Outputs:
        (double dict) reference frequencies in Hz (one per band)
                    
    Calculation:
        Calculates frequencies corresponding to 250, 350, 500 micron
    
    Dependencies:
        spireBands() - get list of spire bands
        herschel.share.unit.Constant.SPEED_OF_LIGHT
    """

    #define global variable
    global spireRefFreq #spire reference frequencies (at 250, 350, 500 microns)

    if array:
        #check array is valid
        assert (array=='PSW' or array=='PMW' or array=='PLW'),\
          'array must be "PSW"|"PMW"|"PLW" or None'
        bands=[array]
    else:
        #use all bands
        bands=spireBands()
        
    if verbose: print 'Getting SPIRE reference frequencies for ',bands

    try:
        spireRefFreq
        #global variable already defined, so do nothing
    except:
        spireRefFreq = {}
        #not defined, so recalculate
        c = Constant.SPEED_OF_LIGHT.value
        spireRefWl = {"PSW":250.*1e-6, "PMW":350.*1.e-6, "PLW":500.*1.e-6}
        for band in spireBands():
            spireRefFreq[band] = c/spireRefWl[band]
    
    if array:
        return(spireRefFreq[array])
    else:
        return(spireRefFreq)


def getSpireEffFreq(array=None,verbose=False):
    """
    ========================================================================
    getSpireEffFreq():
        Retrieve spire effective frequencies from calibration tree
        Stored in global variable spireEffFreq so only retrieved once

    Inputs:
      array:      (string) OPTIONAL. SPIRE array to use ('PSW'|'PMW'|'PLW')
                     Default=None (returns all three bands)
      verbose:    (boolean) print more information. Default=False

    Outputs:
        if array set:
            (double) effective frequency in Hz for array
        if array not set:
            (double dict) effective frequencies in Hz (one per band)
                    
    Calculation:
        Retrieves metadata freqEffPxw from RadialCorrBeam
    
    Dependencies:
        getCal() - retrieve calibration tree
    """
    #define global variable
    global spireEffFreq #spire Effective frequencies
    
    if array:
        #check array is valid
        assert (array=='PSW' or array=='PMW' or array=='PLW'),\
          'array must be "PSW"|"PMW"|"PLW" or None'
        bands=[array]
    else:
        #use all bands
        bands=spireBands()
        
    if verbose: print 'Getting SPIRE effective frequencies for ',bands

    try:
        spireEffFreq
        #global variable already exists, so do nothing
    except:
        #doesn't exist, so get from calFile
        beamProfs = getCal().getProduct("RadialCorrBeam")
        spireEffFreq = {"PSW":beamProfs.meta['freqEffPsw'].double*1.e9,\
            "PMW":beamProfs.meta['freqEffPmw'].double*1.e9,\
            "PLW":beamProfs.meta['freqEffPlw'].double*1.e9}

    if array:
        return(spireEffFreq[array])
    else:
        return(spireEffFreq)

def getSpireFilt(rsrfOnly=False,array=None,verbose=False):
    """
    ========================================================================
    getSpireFilt(rsrfOnly=False):
        Retrieve spire filter profiles from calibration tree (with/without aperture
        efficiency) and interpolated to frequency raster from getSpireFreq()
        Stored in global spireFilt and spireFiltOnly so only calculated once
        Store filter only in global variable spireFiltOnly
        Store filter*apertureefficiency in global variable spireFilt

    Inputs:
      rsrfOnly:   (boolean) set to retrieve spire filter profiles without
                      aperture efficiency. Default=False (include ap. eff.)
      array:      (string) OPTIONAL. SPIRE array to use ('PSW'|'PMW'|'PLW')
                     Default=None (returns all three bands)
      verbose:    (boolean) print more information. Default=False

    Outputs:
        (float array) Filter profile over frequency raster
           [if rsrfOnly==True: RSRF only, else RSRF*Ap.Eff]
                    
    Calculation:
      + If rsrfOnly=True
        + If global variable spireFiltOnly exists, return that
        + If not:
          - Retrieves photometer calibration product (using getCal())
          - Retrieves filter profile from calibration
          - Retrieves the SPIRE frequency raster (using getSpireFreq())
          - Interpolates filter to frequency raster
      + If rsrfOnly=False
        + If global variable spireFilt exists, return that
        + If not:
          - Retrieves photometer calibration product (using getCal())
          - Retrieves filter profile and aperture efficiency from calibration
          - Retrieves the SPIRE frequency raster (using getSpireFreq())
          - Interpolates filter and ap.eff to frequency raster and multiplies

    Dependencies:
        spireBands() - get list of spire bands
        getCal() - retrieve photometer calibration product
        getSpireFreq() - retrieve frequency raster
        herschel.ia.numeric.toolbox.interp.LinearInterpolator
    """

    #define global variables
    global spireFiltOnly #spire filter profiles (without ap. eff.) over spireFreq
    global spireFilt #spire filter profiles (with ap. eff.) over spireFreq

    if array:
        #check array is valid
        assert (array=='PSW' or array=='PMW' or array=='PLW'),\
          'array must be "PSW"|"PMW"|"PLW" or None'
        bands=[array]
    else:
        #use all bands
        bands=spireBands()
        
    if rsrfOnly:
        if verbose: print 'calculating RSRF for ',bands
    else:
        if verbose: print 'calculating RSRF and Aperture Efficiency for ',bands
    #check spireCalPhot exists
    spireCalPhot=getCal()
    assert spireCalPhot.isValid(), 'Invalid SPIRE Photometer calibration tree. Run getCal()'
    
    for band in bands:

        # spire RSRF only
        try:
            spireFiltOnly[band]
            #global variable spireFiltOnly defined, so don't recalculate
        except:
            try:
                spireFiltOnly
            except:
                spireFiltOnly={}
            #spireFiltOnly not defined, so must calculate
    
           
            #read RSRF and Aperture Efficiency from calibration tree
            rsrf=spireCalPhot.getProduct('Rsrf')
            rsrfVersion=rsrf.getVersion()
    
            #
            spireRsrfFreq   = rsrf.getFrequency()*1e9  # Frequency in Hz
            #indexes of freq in rsrf
            spireFreq=getSpireFreq()
            nNu=len(spireFreq)
            ixR = spireFreq.where((spireFreq>=MIN(spireRsrfFreq)) & (spireFreq<=MAX(spireRsrfFreq)))
            #
            #interpolate to freq array
            #create Rsrf and ApEff interpolation objects
            interpRsrf = LinearInterpolator(spireRsrfFreq, rsrf.getRsrf(band))
            #make arrays for final objects
            spireFiltOnly[band] = Double1d(nNu)
            #interpolate Rsrf to freq array
            spireFiltOnly[band][ixR] = interpRsrf(spireFreq[ixR])

        # spire RSRF * ApEff
        if not rsrfOnly:
            try:
                spireFilt[band]
                #global variable spireFilt defined, so do nothing
            except:
                try:
                    spireFilt
                except:
                    #define spireFilt
                    spireFilt={}
                #spireFilt not defined, so must calculate
                #add in aperture efficiency
                apertureEfficiency = spireCalPhot.getProduct('ApertureEfficiency')
                apertureEfficiencyVersion=apertureEfficiency.getVersion()
                spireApEffFreq = apertureEfficiency.getApertEffTable()["frequency"].data * 1e9 #comes in [GHz]
                #indexes of freq in apEff
                spireFreq=getSpireFreq()
                nNu=len(spireFreq)
                ixA = spireFreq.where((spireFreq>=MIN(spireApEffFreq)) & (spireFreq<=MAX(spireApEffFreq)))
                #interpolate to freq array and apply to RSRF
    
                #create ApEff interpolation objects
                interpAp = LinearInterpolator(spireApEffFreq, apertureEfficiency.getApertEffTable()[band].data)
                #make arrays for final objects
                spireFilt[band] = Double1d(nNu)
                #copy into Rsrf*ApEff array
                spireFilt[band] = spireFiltOnly[band].copy()
                spireFilt[band][ixA] = spireFilt[band][ixA] * interpAp(spireFreq[ixA])

    if rsrfOnly:
        #return spireFiltOnly
        if array:
            return(spireFiltOnly[array])
        else:
            return(spireFiltOnly)
    else:
        #return spireFilt
        if array:
            return(spireFilt[array])
        else:
            return(spireFilt)

def calcBeamMonoArea(beamType=None,array=None,verbose=False):
    """
    ========================================================================
    calcBeamMonoArea(beamType=None,verbose=False):
        Retrieve the monochromatic beam areas across frequency raster, using either
          simple beam treatment or full beam treatment
        Stored in global variable beamMonoArea, so only calculated once, providing
          beamType matches previously calculated version

    Inputs:
      beamType:   (string) type of beam calculate (None|'full'|'simple').
                    Default=None (i.e. use existing, or 'full' if recalculating)
      array:      (string) OPTIONAL. SPIRE array to use ('PSW'|'PMW'|'PLW')
                     Default=None (returns all three bands)
      verbose:    (boolean) print more information. Default=False

    Outputs:     
      (dict) Monochromatic beam areas over frequency raster.
        [one float array per band, plus beamType='full'|'simple']
        
    Calculation:
      [See Handbook Eq 5.32/5.36]
      + If global variable beamMonoArea exists, and beamType matches (or None),
          returns previosly calculated version
      + If not:
          - Retrieves photometer calibration product (using getCal())
          - Extracts radial beam profile
          - Retrieves the SPIRE frequency raster (using getSpireFreq())
          - Retrieves the SPIRE effective frequencies (using getSpireEffFreq())
          + If beamType=simple:
              [See Handbook Eq 5.36]
              - Retrieves Neptune beam areas and scales with frequency
          + If beamType=full of None:
              [See Handbook Eq 5.32]
              - Calculates beam area using full treatment over frequency raster

    Dependencies:
        spireBands() - get list of spire bands
        getCal() - retrieve photometer calibration product
        getSpireFreq() - retrieve frequency raster
        getSpireEffFreq() - retrieve effective frequencies
        spireMonoAreas() - calculate monochromatic beam areas (full treatment)
    """

    #define global variables
    global beamMonoArea #monochromatic beam areas over spireFreq

    if array:
        #check array is valid
        assert (array=='PSW' or array=='PMW' or array=='PLW'),\
          'array must be "PSW"|"PMW"|"PLW" or None'
        bands=[array]
    else:
        #use all bands
        bands=spireBands()
        
    if verbose: print 'Calculating monochromatic beam areas for ',bands

    
    #print docalc
    #Check whether existing beamType matches specified
    try:
        beamMonoArea['beamType']
    except:
        #no beamtype. Reset
        if verbose: print 'Resetting beamMonoArea'
        beamMonoArea={'beamType':None}
        
    if beamMonoArea['beamType']=='Simple' and beamType=='Full':
        if (verbose):print 'Replacing Simple Beams with Full Beams'
        recalc=True
    elif beamMonoArea['beamType']=='Full' and beamType=='Simple':
        if (verbose):print 'Replacing Full Beams with Simple Beams'
        recalc=True
    elif beamMonoArea['beamType']==None:
        if (verbose):print 'Existing beamType undefined, so recalculating'
        recalc=True
    elif beamType==None:
        if (verbose):print 'No beamType specified, so using existing %s beams'%beamMonoArea['beamType']
        recalc=False
        beamType=beamMonoArea['beamType']
    else:
        if (verbose):print 'Using existing %s beams'%beamMonoArea['beamType']
        recalc=False
        beamType=beamMonoArea['beamType']
            
    if recalc:
        beamMonoArea={}
        
    for band in bands:
        try:
            beamMonoArea[band]
            docalc=False
            #global variable already defined, so do nothing.
        except:
            #not defined, so need to recalculate
            docalc=True
            
        if docalc:

            #recalculate beam profiles
            beamProfs = getCal().getProduct("RadialCorrBeam")
            gamma = beamProfs.meta['gamma'].double
    
            if beamType=='Simple':
                beamMonoArea['beamType']='Simple'
                #use simple beam treatment
                neptuneArea = {}
                neptuneArea['PSW']=beamProfs.meta['beamNeptunePswSr'].double
                neptuneArea['PMW']=beamProfs.meta['beamNeptunePmwSr'].double
                neptuneArea['PLW']=beamProfs.meta['beamNeptunePlwSr'].double
                if (verbose):print 'Calculating monochromatic beam areas for %s band (Simple treatment)'%(band)
                beamMonoArea[band]=neptuneArea[band] * (getSpireFreq()/getSpireEffFreq(array=band))**(2.*gamma)
    
        
            else:
                #use full beam treatment
                beamMonoArea['beamType']='Full'
                #monochromatic beam areas
                if (verbose):print 'Calculating monochromatic beam areas for %s band (Full treatment)'%(band)
                beamMonoArea[band] = spireMonoAreas(getSpireFreq(), beamProfs, 
                  getSpireEffFreq(array=band), gamma, band)
        else:
            if verbose:print 'Using existing monochromatic beam areas for %s band (%s treatment)'%(band,beamType)

    if array:
        return beamMonoArea[array]
    else:
        return beamMonoArea

#-------------------------------------------------------------------------------
#===============================================================================
#=====                            BEAM FUNCTIONS                           =====
#===============================================================================
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Calculate the effective beam area for a given spectrum (Eq. 5.34 in SPIRE Handbook v2.5)
def spireEffArea(freq, transm, monoArea, BB=False, temp=20.0, beta=1.8, alpha=-1.0):
    """
    ========================================================================
    spireEffArea(freq, transm, monoArea, BB=False, temp=20.0, beta=1.8, alpha=-1.0):
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
        See Eq. 5.34 in SPIRE Handbook v2.5
      Calculates the source spectrum (either modifies black body or power law)
          Multiplies the monochromatic beam area by RSRF and source spectrum
      Integrates over frequency
      Normalises by integral over frequency of RSRF and source spectrum

    Dependencies:
      herschel.ia.numeric.toolbox.interp.LinearInterpolator
      herschel.ia.numeric.toolbox.integr.TrapezoidalIntegrator

    2013/12/19  C. North  initial version

    """    
    h = Constant.H_PLANCK.value
    k = Constant.K_BOLTZMANN.value
    c = Constant.SPEED_OF_LIGHT.value
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
    numInterp=CubicSplineInterpolator(freq,transm * fSky * monoArea)
    denomInterp=CubicSplineInterpolator(freq,transm * fSky)
    minFreq=MIN(freq)
    maxFreq=MAX(freq)
    integrator=TrapezoidalIntegrator(minFreq,maxFreq)
    numInteg=integrator.integrate(numInterp)
    denomInteg=integrator.integrate(denomInterp)
    effArea = numInteg / denomInteg

    return(effArea)

#-------------------------------------------------------------------------------
# Calculate monochromatic beam profile and area at a given frequency (Eq. 5.32) 
def spireMonoBeam(freqx,beamRad,beamProfs,beamConst,effFreq,gamma,array):
    """
    ========================================================================
    spireMonoBeam(freqx,beamRad,beamProfs,beamConst,effFreq,gamma,array):
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
      array:     (string) spire array ('Psw'|'Pmw'|'Plw')

    Outputs:     (list of objects)
      [0]:       (float) Beam area [arcsec^2] at frequency freqx
      [1]:       (array float) Monochromatic beam profile at frequency freqx

    Calculation:
        [See Handbook Eq 5.32]
          Scales the core beam profile width as (freqx/effFreq)^gamma.
      Queries the calibration file to generate new core beam profile.
      Uses constant beam profile where it is larger than core profile.
      Integrates over radius to calculate beam area.

    Dependencies:
      herschel.ia.numeric.toolbox.interp.LinearInterpolator
      herschel.ia.numeric.toolbox.integr.TrapezoidalIntegrator
      
    2013/12/19  C. North  initial version

    """

    #calculate the "scaled" radius, as nu^-gamma
    radNew=beamRad*(freqx/effFreq)**-gamma
    maxRad=MAX(beamRad)
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

    #integrate to get solid angle (in arcsec^2)
    
    beamInterp=CubicSplineInterpolator(beamRad,beamNew * 2. * PI * beamRad)
    integrator=TrapezoidalIntegrator(0,maxRad)
    beamMonoArea=integrator.integrate(beamInterp)

    return (beamMonoArea,beamNew)

#-------------------------------------------------------------------------------
# Calculate monochromatic beam areas over a range of frequencies
def spireMonoAreas(freq,beamProfs,effFreq,gamma,array,freqFact=100):

    """
    ========================================================================
    spireMonoAreas(freq,beamProfs,effFreq,gamma,array,freqFact=100):
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

    """

    #arcsec2Sr = (PI/(60.*60.*180))**2

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
        beamMonoAreaSparse[fx]=spireMonoBeam(freq[f],beamRad,beamProfs,beamConst,effFreq,gamma,array)[0]

    # interpolate to full frequency array and convert to Sr
    beamInterp=CubicSplineInterpolator(beamMonoFreqSparse,beamMonoAreaSparse)
    beamMonoArea=beamInterp(freq)*arcsec2Sr() #in sr
    
    return(beamMonoArea)


#-------------------------------------------------------------------------------
#===============================================================================
#=====                     COLOUR CORRECTION FUNCTION                      =====
#===============================================================================
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Calculate K-correction parameters for given spectrum & source type
def calcSpireKcorr(freq0, freq, transm, BB=True, temp=20.0, beta=1.8, alpha=-1.0,
  ext=False, monoArea=None):
    """
    ================================================================================
    calcSpireKcorr(freq0, freq, transm, BB=True, temp=20.0, beta=1.8, alpha=-1.0, ext=False, monoArea=None):
        Calculation of the K-correction factor from isophotal flux to a monochromatic 
          flux-density at a given reference frequency (data to be multiplied!)
    
    Inputs:
      freq0:     (float) waveband reference frequency [Hz] for which monochromatic
                   flux-density is given
      freq:      (array float) frequency vector corresponding to RSRF values [Hz]
      transm:    (array float) relative spectral response (RSRF) corresponding to freq
      BB:        (boolean) spectral function to use for source spectrum:
                    'True': a modified black body
                    'False' a power-law with exponent alpha=-1
                    OPTIONAL. Default=False
      temp:      (float) Dust/sky temperature [K] (only for modified black body)
            OPTIONAL. Deafult=20K; 
                    only for modified black body]
      beta:      (float) Dust/sky spectral index (only for modified black body]
            OPTIONAL. Default=1.8
      alpha:     (float) Exponent of power-law sky background model (only for
            power-law spectrum)
            OPTIONAL. Default=-1
      ext:       (boolean) calculating for extended source
                    OPTIONAL. Default=False
      monoArea:  (array float) Monochromatic Beam solid angle [Sr] corresponding
                    to freq.
                    OPTIONAL. Only required if ext=True

    Outputs:
     (list)     [0]: K-correction factor
                [1]: Sky emission at reference fequency (fSky0)
    
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

      N.B. If ext=True, and monoBeam is in units sr, then this procedure outputs
        K-correction factor in [Jy/sr per Jy/beam] and Sky emission in Jy/sr.
        Units will change if a different input unit is used.
    
    Dependencies:
      herschel.ia.numeric.toolbox.interp.LinearInterpolator
      herschel.ia.numeric.toolbox.integr.TrapezoidalIntegrator

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
    2013/12/03  C. North     corrected procedure include area where appropriate
                  NB: if ext=True , note output units
    2014/04/15  C. North     Renamed script to calcSpireKcorr
    
    ================================================================================
    """

    h = Constant.H_PLANCK.value
    k = Constant.K_BOLTZMANN.value
    c = Constant.SPEED_OF_LIGHT.value

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

    if ext == True:
        area = monoArea
    else:
        #don't use area for point sources
        area = Float1d(len(freq))
        area[:] = 1.0

        # monoArea is the monochromatic beam solid angle at effFreq

    # integrate over frequency
    numInterp=CubicSplineInterpolator(freq,transm)
    denomInterp=CubicSplineInterpolator(freq,transm * fSky * area)
    minFreq=MIN(freq)
    maxFreq=MAX(freq)

    integrator=TrapezoidalIntegrator(minFreq,maxFreq)
    numInteg = integrator.integrate(numInterp)
    denomInteg = integrator.integrate(denomInterp)

    kWave = fSky0 * numInteg / denomInteg

    #
    # Return the result as a 2-element array of K-correction and flux at freq0
    return (Double1d([kWave, fSky0]))

#-------------------------------------------------------------------------------
# Calculate the effective beam profile and area for SPIRE
def spireEffBeam(freq, transm, beamProfs, effFreq, gamma, array, BB=False,temp=20.0,beta=1.8,alpha=-1.0,verbose=False):

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
      area:    (float) beam solid angle
      profile: (array float) radialised beam profile

    Calculation:
      Calculates the source spectrum (either modifies black body or power law)
      Loops over beam radius, and calculates scaled radii for full frequency range
      For that radius, gets the values of profile at those scaled radii
      Integrates profile values over frequency, weighted by RSRF and source spectrum

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
    #
    else:
        # Power-Law spectrum
        if verbose:
            print 'Using power law spectrum with spectral index %.2f'%(alpha)
        fSky  = freq**alpha


    #integrate transm*fSky over frequency for nomalisation
    integrator=TrapezoidalIntegrator(min(freq),max(freq))
    denomInterp=CubicSplineInterpolator(freq,transm*fSky)
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

    beamInterp = CubicSplineInterpolator(beamRad,effBeam * 2.*PI*beamRad)
    effBeamAreaInt=TrapezoidalIntegrator(0,maxRad)
    effBeamArea=effBeamAreaInt.integrate(beamInterp)
    effBeamDict = {'area':effBeamArea,'profile':effBeam}
    #if verbose:print'calculated area:',effBeamArea
    return(effBeamDict)
    
#-------------------------------------------------------------------------------
#===============================================================================
#=====                    CALCULATE EFFECTIVE BEAM AREAS                   =====
#===============================================================================
#-------------------------------------------------------------------------------

def calcSpireEffBeam(alphaK,array=None,verbose=False):
    """
    ================================================================================
    calcSpireEffBeam(alphaK,array=None,verbose=False):
        Calculated effective beam area for modified blackbody spectrum.
        Stores result in global variable e

    Inputs:
        alphaK:  (scalar/float array) Spectral index(es) to use
        array:   (string) OPTIONAL. SPIRE array to use ('PSW'|'PMW'|'PLW')
                    Default=None (returns all three bands)
        verbose: (boolean) Set to print more information. Default=True.
    Outputs:
        if array set:
            (double array) effective beam profile for array
        if no array set:
            (dict) effective beam profiles (one per band)
                    
    Calculation:
        [See Handbook Eq 5.32-34]
        Integrates monochromatic beam profile over frequency raster, weighted by
          spire filter profiles and source spectrum
    
    Dependencies:
        spireBands() - list of spire bands
        getSpireFreq() - get frequency raster
        getSpireEffFreq() - get spire effective frequencies
        getCal() - retrieve calibration tree
        getSpireFilt() - get spire Filter profiles
        calcBeamMonoArea() - calculate monochromatic beam areas
        spireEffArea - calculate effective beam area
        herschel.ia.dataset.TableDataset()
    """

    if array:
        #check array is valid
        assert (array=='PSW' or array=='PMW' or array=='PLW'),\
          'array must be "PSW"|"PMW"|"PLW" or None'
        bands=[array]
    else:
        #use all bands
        bands=spireBands()
        
    if verbose: print 'Calculating effective beam profile for power law for ',bands

    key='alpha_%g'%alphaK
    
    global spireEffBeams
    try:
        spireEffBeams[key]
    except:
        try:
            spireEffBeams
        except:
            spireEffBeams={}
        spireEffBeams[key]={}
        
    for band in bands:
        try:
            #see if it currently exists
            effBeam=spireEffBeams[key][band]
            if verbose: print 'Using existing effective %s beam for alpha %.2f'%(band,alphaK)
        except:
            #recalculate effBeam
            if verbose: print 'Calculating effective %s beam for alpha %.2f'%(band,alphaK)
            
            #check global variable exists
            ##now done when module is loaded
            #try:
            #    spireEffBeams
            #except:
            #    #create global variable
            #    spireEffBeams={}
                
            #effBeams={'PSW':Double.NaN,'PMW':Double.NaN,'PLW':Double.NaN}
            beamProfs=getCal().getProduct('RadialCorrBeam')
            beamRad=getBeamRad()
            effBeam=Double1d(beamRad)
            
            freq=getSpireFreq()
            fSky  = freq**alphaK
                    
            effFreq=getSpireEffFreq()
            gamma=beamProfs.meta['gamma'].double
    
            #integrate transm*fSky over frequency for nomalisation
            integrator=TrapezoidalIntegrator(MIN(freq),MAX(freq))

            transm=getSpireFilt(rsrfOnly=True)[band]
            denomInterp=CubicSplineInterpolator(freq,transm*fSky)
            denomInteg=integrator.integrate(denomInterp)
            #get core beam profile from calibration table
            beamCore=beamProfs.getCoreCorrectionTable().getColumn(band).data
            #create interpolation object
            beamCoreInt=CubicSplineInterpolator(beamRad,beamCore)
            
            #make array for new beam
            nRad=len(beamRad)
            maxRad=max(beamRad)
            
            #loop over radius
            for r in range(nRad):
                #calculate the "scaled" radius for range of frequencies
                radFreq=beamRad[r]*(freq/effFreq[band])**-gamma
                #ensure it doesn't fo beyong maximum radius
                radFreq[radFreq.where(radFreq > maxRad)]=maxRad
                #compute value beam profile at each scaled radius
                beamCoreFreq=beamCoreInt(radFreq)
                #apply constant beam profile value where appropriate
                beamConstRad=beamProfs.getConstantCorrection(beamRad[r],band)
                isConst=beamCoreFreq.where(beamCoreFreq < beamConstRad)
                beamCoreFreq[isConst]=beamConstRad
        
                #integrate beamCoreFreq*transm*fSky over frequency
                numInterp=CubicSplineInterpolator(freq,beamCoreFreq*transm*fSky)
                numInteg = integrator.integrate(numInterp)
        
                #write value into table
                effBeam[r]=numInteg/denomInteg

            #store in global variable
            spireEffBeams[key][band]=effBeam

    if array:
        #return single band
        return(spireEffBeams[key][array])
    else:
        #return all bands
        return(spireEffBeams[key])
    
def calcSpireEffBeam_BB(betaK,tempK,array=None,verbose=False):
    """
    ================================================================================
    calcSpireEffBeam_BB(betaK,tempK,array=None,verbose=False):
        Calculated effective beam area for modified blackbody spectrum.
        Stores result in global variable e

    Inputs:
        betaK:   (float) Emmisivity index to use
        tempK:   (scalar/float array) Temperature(s) to use
        array:   (string) OPTIONAL. SPIRE array to use ('PSW'|'PMW'|'PLW')
                    Default=None (returns all three bands)
        verbose: (boolean) Set to print more information. Default=True.
    Outputs:
        if array set:
            (double array) effective beam profile for array
        if no array set:
            (dict) effective beam profiles (one per band)
                    
    Calculation:
        [See Handbook Eq 5.32-34]
        Integrates monochromatic beam profile over frequency raster, weighted by
          spire filter profiles and source spectrum
    
    Dependencies:
        spireBands() - list of spire bands
        getSpireFreq() - get frequency raster
        getSpireEffFreq() - get spire effective frequencies
        getCal() - retrieve calibration tree
        getSpireFilt() - get spire Filter profiles
        calcBeamMonoArea() - calculate monochromatic beam areas
        spireEffArea - calculate effective beam area
        herschel.ia.dataset.TableDataset()
    """

    if array:
        #check array is valid
        assert (array=='PSW' or array=='PMW' or array=='PLW'),\
          'array must be "PSW"|"PMW"|"PLW" or None'
        bands=[array]
    else:
        #use all bands
        bands=spireBands()
        
    if verbose: print 'Calculating effective beam profile for modified black body for ',bands

    key='beta_%g_temp_%g'%(betaK,tempK)
    
    global spireEffBeams
    try:
        spireEffBeams[key]
    except:
        try:
            spireEffBeams
        except:
            spireEffBeams={}
        spireEffBeams[key]={}
    
    for band in bands:
        try:
            effBeam=spireEffBeams[key][band]
            if verbose: print 'Using existing effective %s beam for beta %.2f and temp %.2f'%(band,betaK,tempK)
        except:
            #check global variable exists
            ##now done when module is loaded
            #try:
            #    spireEffBeams
            #except:
            #    #create global variable
            #    spireEffBeams={}
                
            #effBeams={'PSW':Double.NaN,'PMW':Double.NaN,'PLW':Double.NaN}
            beamProfs=getCal().getProduct('RadialCorrBeam')
            beamRad=getBeamRad()
            effBeam=Double1d(beamRad)
            
            h = Constant.H_PLANCK.value
            k = Constant.K_BOLTZMANN.value
            c = Constant.SPEED_OF_LIGHT.value
        
            if verbose: print 'Calculating effective beams for beta %.2f and temp %.2f'%(betaK,tempK)
            freq=getSpireFreq()
            fSky  = 2*h * freq**3 / c**2 / (EXP(h*freq/k/tempK) - 1.) * freq**betaK
                    
            effFreq=getSpireEffFreq()
            gamma=beamProfs.meta['gamma'].double
    
            #integrate transm*fSky over frequency for nomalisation
            integrator=TrapezoidalIntegrator(MIN(freq),MAX(freq))
        
            transm=getSpireFilt(rsrfOnly=True)[band]
            denomInterp=CubicSplineInterpolator(freq,transm*fSky)
            denomInteg=integrator.integrate(denomInterp)

            #get core beam profile from calibration table
            beamCore=beamProfs.getCoreCorrectionTable().getColumn(band).data
            #create interpolation object
            beamCoreInt=CubicSplineInterpolator(beamRad,beamCore)
        
            #make array for new beam
            nRad=len(beamRad)
            maxRad=max(beamRad)
            #loop over radius
            for r in range(nRad):
                #calculate the "scaled" radius for range of frequencies
                radFreq=beamRad[r]*(freq/effFreq[band])**-gamma
                #ensure it doesn't fo beyong maximum radius
                radFreq[radFreq.where(radFreq > maxRad)]=maxRad
                #compute value beam profile at each scaled radius
                beamCoreFreq=beamCoreInt(radFreq)
                #apply constant beam profile value where appropriate
                beamConstRad=beamProfs.getConstantCorrection(beamRad[r],band)
                isConst=beamCoreFreq.where(beamCoreFreq < beamConstRad)
                beamCoreFreq[isConst]=beamConstRad
        
                #integrate beamCoreFreq*transm*fSky over frequency
                numInterp=CubicSplineInterpolator(freq,beamCoreFreq*transm*fSky)
                numInteg = integrator.integrate(numInterp)
        
                #write value into table
                effBeam[r]=numInteg/denomInteg    

            #store in global variable
            spireEffBeams[key][band]=effBeam
        
    if array:
        return(spireEffBeams[key][array])
    else:
        return(spireEffBeams[key])
    
def calcOmegaEff(alphaK,array=None,verbose=False,table=False):
    """
    ================================================================================
    calcOmegaEff(alphaK,array=None,verbose=False,table=False):
        Calculated effective beam area for power law spectrum

    Inputs:
        alphaK:  (scalar/float array) Spectral index(es) to use
        array:   (string) OPTIONAL. SPIRE array to use ('PSW'|'PMW'|'PLW')
                    Default=None (returns all three bands in dict)
        verbose: (boolean) Set to print more information. Default=True.
        table:   (boolean) Set to output TableDataset. Default=False.
    Outputs:
        If table==False:
            if array set:
              (float scalar/list) effective areas in arcsec^2 for array
                (one scalar/list per band, depending on whether alphaK is scalar/list)
            if array not set:
              (dict) effective areas in arcsec^2 (one scalar/list per band,
                depending on whether alphaK is scalar/list)
        If table==True:
            (TableDataset) effective areas in arcsec^2 (one column per band,
              one row per alphaK value)
                    
    Calculation:
        [See Handbook Eq 5.32-34]
        Integrates monochromatic beam areas over frequency raster, weighted by
          spire filter profiles and source spectrum
        If table is set, returns a TableDataset
    
    Dependencies:
        spireBands() - list of spire bands
        getSpireFreq() - get frequency raster
        getCal() - retrieve calibration tree
        getSpireFilt() - get spire Filter profiles
        calcBeamMonoArea() - calculate monochromatic beam areas
        spireEffArea - calculate effective beam area
        herschel.ia.dataset.TableDataset()
    """
    # calculate effective beam area for power law spectrum

    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #freq=getSpireFreq()
    #spireFiltOnly=getSpireFilt(cal,rsrfOnly=True)

    try:
        na=len(alphaK)
        aList=True
    except:
        aList=False

    if array:
        #check array is valid
        assert (array=='PSW' or array=='PMW' or array=='PLW'),\
          'array must be "PSW"|"PMW"|"PLW" or None'
        bands=[array]
    else:
        #use all bands
        bands=spireBands()
        
    if verbose: print 'Calculating effective beam area for power law for ',bands
    
    beamArea = {}
    if not aList:
        # alphaK is a scalar
        for band in bands:
            #pipeline beam areas
            beamArea[band]=spireEffArea(getSpireFreq(), getSpireFilt(rsrfOnly=True,array=band), \
              calcBeamMonoArea(array=band), BB=False, alpha=alphaK)/arcsec2Sr()
        if (verbose): print 'Calculated %s Omega_eff for alpha=%f: '%(band,alphaK),beamArea[band]
    else:
        # tempK is a list
        for band in bands:
            beamArea[band]=Double1d(na,Double.NaN)
            for a in range(na):
                beamArea[band][a]=spireEffArea(getSpireFreq(), getSpireFilt(rsrfOnly=True,array=band),\
                  calcBeamMonoArea(array=band), BB=False, alpha=alphaK[a])/arcsec2Sr()
            if (verbose): print 'Calculated %s Omega_eff for alpha=%f: '%(band,alphaK[a]),beamArea[band][a]

    if not table:
        #return as is
        if array:
            return(beamArea[array])
        else:
            return(beamArea)
    else:
        #create and returnTableDataset
        beamArea_table=TableDataset()
        beamArea_table.setDescription("Beam Solid Angle (Spectral Index)")
        beamArea_table.addColumn("alpha",Column(Double1d(alphaK)))
        for band in bands:
            beamArea_table.addColumn(band,Column(beamArea[band],unit=SolidAngle.STERADIANS,description=''))
        return(beamArea_table)


def calcOmegaEff_BB(betaK,tempK,array=None,verbose=False,table=False):
    """
    ================================================================================
    calcOmegaEff_BB(betaK,tempK,array=None,verbose=False,table=False):
        Calculated effective beam area for modified black-body spectrum

    Inputs:
        betaK:  (float) Emmisivity index to use
        tempK:  (scalar/float array) Temperature(s) to use
        array:   (string) OPTIONAL. SPIRE array to use ('PSW'|'PMW'|'PLW')
                    Default=None (returns all three bands in dict)
        verbose: (boolean) Set to print more information. Default=True.
        table:   (boolean) Set to output TableDataset. Default=False.
    Outputs:
        If table==False:
            if array set:
              (float scalar/list) effective areas in arcsec^2 for array
                (one scalar/list per band, depending on whether tempK is scalar/list)
            if array not set:
              (dict) effective areas in arcsec^2 (one scalar/list per band,
                depending on whether tempK is scalar/list)
        If table==True:
            (TableDataset) effective areas in arcsec^2 (one column per band,
              one row per tempK value)
                    
    Calculation:
        [See Handbook Eq 5.32-34]
        Integrates monochromatic beam areas over frequency raster, weighted by
          spire filter profiles and source spectrum
        If table is set, returns TableDataset (with betaK in metadata)
    
    Dependencies:
        spireBands() - list of spire bands
        getSpireFreq() - get frequency raster
        getCal() - retrieve calibration tree
        getSpireFilt() - get spire Filter profiles
        calcBeamMonoArea() - calculate monochromatic beam areas
        herschel.ia.dataset.TableDataset()
    """    # calculate effective beam area for modified BB spectrum
    # allows multiple temperatures, but only one beta value

    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #freq=getSpireFreq()
    #spireFiltOnly=getSpireFilt(cal,rsrfOnly=True)
    
    try:
        nt=len(tempK)
        tList=True
    except:
        tList=False

    if array:
        #check array is valid
        assert (array=='PSW' or array=='PMW' or array=='PLW'),\
          'array must be "PSW"|"PMW"|"PLW" or None'
        bands=[array]
    else:
        #use all bands
        bands=spireBands()
        
    if verbose: print 'Calculating effective beam area for modified black body for ',bands
        
    beamAreaBB = {}
        
    if not tList:
        # tempK is scalars
        for band in bands:
            #pipeline beam areas
            beamAreaBB[band]=spireEffArea(getSpireFreq(), getSpireFilt(rsrfOnly=True,array=band), \
              calcBeamMonoArea(array=band), BB=True, beta=betaK, temp=tempK)/arcsec2Sr()
        if (verbose): print 'Calculated %s Omega_eff for modBB with beta=%f and T=%f: '%(band,betaK,tempK),beamAreaBB[band]

    else:
        # tempK is a list
        for band in bands:
            beamAreaBB[band]=Double1d(nt,Double.NaN)
            for t in range(nt):
                #pipeline beam areas
                beamAreaBB[band][t]=spireEffArea(getSpireFreq(), getSpireFilt(rsrfOnly=True,array=band), \
                  calcBeamMonoArea(array=band), BB=True, beta=betaK, temp=tempK[t])/arcsec2Sr()
            if (verbose): print 'Calculated %s Omega_eff for modBB with beta=%f and T=%f: '%(band,betaK,tempK[t]),beamAreaBB[band][t]

    if not table:
        #return as is
        if array:
            return(beamAreaBB[array])
        else:
            return beamAreaBB
    else:
        #create and returnTableDataset
        beamAreaBB_table=TableDataset()
        beamAreaBB_table.setDescription("Beam Solid Angle (Modified Black Body, beta=%.2f)"%betaK)
        beamAreaBB_table.meta['beta']=DoubleParameter(betaK,"Emissivity spectral index")
        beamAreaBB_table.addColumn("Temperature",Column(Double1d(tempK),unit=Temperature.KELVIN,description=''))
        for band in bands:
            beamAreaBB_table.addColumn(band,Column(beamAreaBB[band],unit=SolidAngle.STERADIANS,description=''))
        return(beamAreaBB_table)
       
def calcKBeam(alphaK,array=None,verbose=False,table=False):
    """
    ================================================================================
    calcKBeam(alphaK,array=None,verbose=False,table=False):
        Calculated beam correction factor for power law spectrum

    Inputs:
        alphaK:  (scalar/float array) Spectral index(es) to use
        array:   (string) OPTIONAL. SPIRE array to use ('PSW'|'PMW'|'PLW')
                    Default=None (returns all three bands in dict)
        verbose: (boolean) Set to print more information. Default=True.
        table:   (boolean) Set to output TableDataset. Default=False.
    Outputs:
        If table==False:
            if array set:
              (float scalar/list) beam correction factor for array
                (one scalar/list per band, depending on whether alphaK is scalar/list)
            if array not set:
              (dict) beam correction factor (one scalar/list per band,
                depending on whether alphaK is scalar/list)
        If table==True:
            (TableDataset) beam correction factor (one column per band,
              one row per alphaK value)
                    
    Calculation:
        [See Handbook Eq 5.32-39]
        Calculates effective area for pipeline (alpha=-1)
        Calculates effective beam area for input spectrum(s)
        KBeam(spectrum) = Omega_pipeline / Omega(spectrum)
        If table is set, returns a TableDataset
    
    Dependencies:
        spireBands() - list of spire bands
        calcOmegaEff() - calculate effective beam area
        herschel.ia.dataset.TableDataset()
    """
    # Calculate pipeline colour correction parameters
    # cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    try:
        na=len(alphaK)
        aList=True
    except:
        aList=False

    if array:
        #check array is valid
        assert (array=='PSW' or array=='PMW' or array=='PLW'),\
          'array must be "PSW"|"PMW"|"PLW" or None'
        bands=[array]
    else:
        #use all bands
        bands=spireBands()
        
    if verbose: print 'Calculating beam colour corrections for power law for ',bands

    kBeam = {}
    for band in bands:
        beamAreaPip = calcOmegaEff(-1.0,array=band)
        if not aList:
            # alphaK is a scalar
            beamEff = calcOmegaEff(alphaK,array=band)
            #pipeline beam areas
            kBeam[band] = beamAreaPip/beamEff
            if verbose: print 'Calculated %s KBeam for alpha=%f: '%(band,alphaK),kBeam[band]
        else:
            # alphaK is a list
            #kBeam[band] = Double1d(na,Double.NaN)
            beamEff = calcOmegaEff(alphaK,array=band)
            #pipeline beam areas
            kBeam[band] = beamAreaPip/beamEff
            if verbose:
                for a in range(na):
                    print 'Calculated %s KBeam for alpha=%f: '%(band,alphaK[a]),kBeam[band][a]

    if not table:
        #return as is
        if array:
            return(kBeam[array])
        else:
            return kBeam
    else:
        #create and returnTableDataset
        kBeam_table=TableDataset()
        kBeam_table.setDescription("Beam Colour Correction (Spectral Index)")
        kBeam_table.addColumn("alpha",Column(Double1d(alphaK)))
        for band in bands:
            kBeam_table.addColumn(band,Column(kBeam[band]))
        return kBeam_table
#
def calcKBeam_BB(betaK,tempK,array=None,verbose=False,table=False):
    """
    ================================================================================
    calcKBeam_BB(betaK,tempK,array=None,verbose=False,table=False):
        Calculated beam correction factor for modified black body spectrum

    Inputs:
        betaK:   (float) Emmisivity index to use
        tempK:   (scalar/float array) Temperature(s) to use
        array:   (string) OPTIONAL. SPIRE array to use ('PSW'|'PMW'|'PLW')
                    Default=None (returns all three bands in dict)
        verbose: (boolean) Set to print more information. Default=True.
        table:   (boolean) Set to output TableDataset. Default=False.
    Outputs:
        If table==False:
            if array set:
                (scalar/list) beam correction factor for array
                  (scalar/list, depending on whether tempK is scalar/list)
            if array not set:
                (dict) beam correction factor (one scalar/list per band,
                  depending on whether tempK is scalar/list)
        If table==True:
            (TableDataset) beam correction factor (one column per band,
              one row per tempK value)
                    
    Calculation:
        [See Handbook Eq 5.32-39]
        Calculates effective area for pipeline (alpha=-1)
        Calculates effective beam area for input spectrum(s)
        KBeam(spectrum) = Omega_pipeline / Omega(spectrum)
        If table is set, returns a TableDataset (betaK in metadata)
    
    Dependencies:
        spireBands() - list of spire bands
        calcOmegaEff() - calculate effective beam area for power law spectrum
        calcOmegaEff_BB() - calculate effective beam area for mod-BB spectrum
        herschel.ia.dataset.TableDataset()
    """
    # Calculate pipeline colour correction parameters
    # allows multiple temperatures, but only one beta value
    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    try:
        nt=len(tempK)
        tList=True
    except:
        tList=False

    if array:
        #check array is valid
        assert (array=='PSW' or array=='PMW' or array=='PLW'),\
          'array must be "PSW"|"PMW"|"PLW" or None'
        bands=[array]
    else:
        #use all bands
        bands=spireBands()
        
    if verbose: print 'Calculating beam colour corrections for modified black body for ',bands

    kBeamBB = {}
    for band in bands:
        beamAreaPip = calcOmegaEff(-1.0,array=band)
        if not tList:
            # tempK is scalar
            beamEff = calcOmegaEff_BB(betaK,tempK,array=band)
            #pipeline beam areas
            kBeamBB[band] = beamAreaPip/beamEff
            if (verbose): print 'Calculated %s KBeam for modBB with beta=%f and T=%f: '%(band,betaK,tempK),kBeamBB[band]
        else:
            # tempK is a list
            #kBeamBB[band] = Double1d(nt,Double.NaN)
            beamEff = calcOmegaEff_BB(betaK,tempK,array=band)
            kBeamBB[band] = beamAreaPip/beamEff
            if (verbose):
                for t in range(nt): 
                    print 'Calculated %s KBeam for modBB with beta=%f and T=%f: '%(band,betaK,tempK[t]),kBeamBB[band][t]

    if not table:
        #return as is
        if array:
            return(kBeamBB[array])
        else:
            return kBeamBB
    else:
        #create and returnTableDataset
        kBeamBB_table=TableDataset()
        kBeamBB_table.setDescription("Beam Colour Correction (Modified Black Body, beta=%.2f)"%betaK)
        kBeamBB_table.meta['beta']=DoubleParameter(betaK,"Emissivity spectral index")
        kBeamBB_table.addColumn("Temperature",Column(Double1d(tempK),unit=Temperature.KELVIN,description=''))
        for band in bands:
            kBeamBB_table.addColumn(band,Column(kBeamBB[band]))
        return(kBeamBB_table)

#-----------------------------------------------------------------------
#=======================================================================
#=====                 CALCULATE PIPELINE PARAMETERS               =====
#=======================================================================
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Calculate pipeline colour correction parameters
def calcK4P(array=None,verbose=False):
    """
    ================================================================================
    calcK4P(array=None,verbose=False):
        Calculated point souce K4 parameter for pipeline spectrum

    Inputs:
        array: (string) OPTIONAL. SPIRE array to use ('PSW'|'PMW'|'PLW')
                  Default=None (returns all three bands in dict)

    Outputs:
        if array set:
            (double) pipeline point source K4P parameter for array
        if array not set:
            (dict double) pipeline point source K4P parameter (one per band)
                    
    Calculation:
        [See Handbook Eq 5.16]
        Calculates point source K4P parameter for pipeline spectrum (alpha=-1)
    
    Dependencies:
        getCal() - retrieve spire photometer calibration context
        spireBands() - list of spire bands
        getSpireFreq() - get frequency raster
        getSpireRefFreq() - get spire reference frequencies
        getSpireFilt() - get spirefilter profiles
        calcSpireKcorr() - calculate colour correction factor
    """
    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #freq=getSpireFreq()
    #spireFilt=getSpireFilt(cal)

    if array:
        #check array is valid
        assert (array=='PSW' or array=='PMW' or array=='PLW'),\
          'array must be "PSW"|"PMW"|"PLW" or None'
        bands=[array]
    else:
        #use all bands
        bands=spireBands()
        
    if verbose: print 'Calculating K4P for pipeline spectrum for ',bands

    getCal()
    k4P = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
    for band in bands:
        k4P[band] = calcSpireKcorr(getSpireRefFreq(array=band), getSpireFreq(), getSpireFilt(array=band), \
        BB=False, ext=False)[0]
        pass
    if array:
        return k4P[array]
    else:
        return k4P
    
def calcKMonE(array=None,verbose=False):
    """
    ================================================================================
    calcKMonE(array=None,verbose=False):
        Calculated fully- extended source KMonE conversion for pipeline spectrum

    Inputs:
        array: (string) OPTIONAL. SPIRE array to use ('PSW'|'PMW'|'PLW')
                  Default=None (returns all three bands in dict)

    Outputs:
        if array set:
            (double) KMonE conversion for pipeline spectrum for array
        if array not set:
            (dict double) KMonE conversion for pipeline spectrum (one per band)
                    
    Calculation:
        [See Handbook Eq 5.20]
        Calculates fully-extended source KMonE parameter for pipeline spectrum
          (alpha=-1
        KMonE converts broadband RSRF-weighted point source flux density to
          monochromatic fully-extended source surface brightness at reference
          frequencies.
    
    Dependencies:
        getCal() - retrieve spire photometer calibration context
        spireBands() - list of spire bands
        getSpireFreq() - get frequency raster
        getSpireRefFreq() - get spire reference frequencies
        getSpireFilt() - get spirefilter profiles
        calcBeamMonoArea() - calculate monochromatic beam areas
        calcSpireKcorr() - calculate spire flux conversion factor
    """
    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #freq=getSpireFreq()
    #spireFilt=getSpireFilt(cal=cal)

    if array:
        #check array is valid
        assert (array=='PSW' or array=='PMW' or array=='PLW'),\
          'array must be "PSW"|"PMW"|"PLW" or None'
        bands=[array]
    else:
        #use all bands
        bands=spireBands()
        
    if verbose: print 'Calculating KMonE for pipeline spectrum for ',bands

    getCal()
    kMonE = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
    for band in bands:
        kMonE[band] = calcSpireKcorr(getSpireRefFreq(array=band), getSpireFreq(), getSpireFilt(array=band), \
        BB=False, alpha=-1.0, ext=True, monoArea=calcBeamMonoArea(array=band))[0]/1.0e6
        pass

    if array:
        return kMonE[array]
    else:
        return kMonE

def calcK4E(array=None,verbose=False):
    """
    ================================================================================
    calcK4E(array=None,verbose=False):
        Calculated fully- extended source K4E conversion for pipeline spectrum

    Inputs:
        array: (string) OPTIONAL. SPIRE array to use ('PSW'|'PMW'|'PLW')
                  Default=None (returns all three bands in dict)

    Outputs:
        if array set:
            (double) K4E conversion for pipeline spectrum for array
        if array not set:
            (dict double) K4E conversion for pipeline spectrum (one per band)
                    
    Calculation:
        [See Handbook Eq 5.23]
        Calculate KMonE conversion
        Calculate effective beam area for pipeline spectrum
        K4E = KMonE * omega_pipeline
        K4E converts broadband RSRF-weighted point source flux density to
          monochromatic fully-extended source flux density at reference
          frequencies.
    
    Dependencies:
        getCal() - retrieve spire photometer calibration context
        spireBands() - list of spire bands
        calcKMonE() - calculate fully extended
        calcOmegaEff() - calculate effective beam area
    """
    getCal()
    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #spireBands=["PSW","PMW","PLW"]
    
    if array:
        #check array is valid
        assert (array=='PSW' or array=='PMW' or array=='PLW'),\
          'array must be "PSW"|"PMW"|"PLW" or None'
        bands=[array]
    else:
        #use all bands
        bands=spireBands()
        
    if verbose: print 'Calculating K4E for pipeline spectrum for ',bands

    k4E = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
    
    for band in bands:
        kMonE = calcKMonE(array=band)
        omegaEff = calcOmegaEff(-1.0,array=band)
        k4E[band] = kMonE * omegaEff * arcsec2Sr() * 1.0e6
        pass
    if array:
        return k4E[array]
    else:
        return k4E

def calcKPtoE(array=None,verbose=False):
    """
    ================================================================================
    calcKPtoE(array=None,verbose=False):
        Calculated conversion from point source flux density to full-extended source
        surface brightness for pipeline spectrum

    Inputs:
        array: (string) OPTIONAL. SPIRE array to use ('PSW'|'PMW'|'PLW')
                  Default=None (returns all three bands in dict)

    Outputs:
        if array set:
            (double) KPtoE for pipeline spectrum for array
        if array not set:
            (dict double) KPtoE for pipeline spectrum (one per band)
                    
    Calculation:
        [See Handbook Eq 5.22]
        Calculate KMonE conversion
        Calculate K4P conversion
        KPtoE = KMonE / K4P
        KPtoE converts monochromatic point source flux density to
          monochromatic fully-extended source surface brightness at reference
          frequencies.
    
    Dependencies:
        getCal() - retrieve spire photometer calibration context
        spireBands() - list of spire bands
        calcKMonE() - calculate fully extended KMonE parameter
        calcK4P() - calculate point source K4P parameter
    """
    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #spireBands=["PSW","PMW","PLW"]
    
    if array:
        #check array is valid
        assert (array=='PSW' or array=='PMW' or array=='PLW'),\
          'array must be "PSW"|"PMW"|"PLW" or None'
        bands=[array]
    else:
        #use all bands
        bands=spireBands()
        
    if verbose: print 'Calculating KPtoE for pipeline spectrum for ',bands

    kPtoE = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
    
    for band in bands:
        kmonE = calcKMonE(array=band)
        k4p = calcK4P(array=band)
        kPtoE[band] = kmonE/k4p
        
    if array:
        return kPtoE[array]
    else:
        return kPtoE

#-----------------------------------------------------------------------
#=======================================================================
#=====           CALCULATE POINT SOURCE COLOR CORRECTIONS          =====
#=======================================================================
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
def calcKColP(alphaK,array=None,verbose=False,table=False):
    """
    ================================================================================
    calcKColP(alphaK,array=None,verbose=False,table=False):
        Calculated point source colour correction factor for power law spectrum

    Inputs:
        alphaK:  (scalar/float array) Spectral index(es) to use
        array:   (string) OPTIONAL. SPIRE array to use ('PSW'|'PMW'|'PLW')
                    Default=None (returns all three bands in dict)
        verbose: (boolean) Set to print more information. Default=True.
        table:   (boolean) Set to output TableDataset. Default=False.
    Outputs:
        If table==False:
            if array set:
                (scalar/list) colour correction factor for array
                  (scalar/list, depending on whether alphaK is scalar/list)
            if array not set:
                (dict) colour correction factor (one scalar/list per band,
                  depending on whether alphaK is scalar/list)
        If table==True:
            (TableDataset) colour correction factor (one column per band,
              one row per alphaK value)
                    
    Calculation:
        [See Handbook Eq 5.25]
        Calculates K4P parameter pipeline (alpha=-1)
        Calculates point source conversion for input spectrum(s)
        Divides to calculate colour correction factor
        If table is set, returns a TableDataset
        KColP converts from monochromatic point source flux density for
          pipeline spectrum to monochromatic point source flux density for input
          spectrum.
    
    Dependencies:
        spireBands() - list of spire bands
        calcK4P() - calculate K4P parameter
        calcSpireKcorr() - calculate spire flux conversion factor
        herschel.ia.dataset.TableDataset()
    """
    #print '\nCalculating point source colour correction parameters for a given alpha...'

    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #freq=getSpireFreq()
    #spireFilt=getSpireFilt(cal=cal)

    if array:
        #check array is valid
        assert (array=='PSW' or array=='PMW' or array=='PLW'),\
          'array must be "PSW"|"PMW"|"PLW" or None'
        bands=[array]
    else:
        #use all bands
        bands=spireBands()
        
    if verbose: print 'Calculating point source colour corrections for power law for ',bands

    try:
        na=len(alphaK)
        aList=True
    except:
        aList=False
    
    kColP={}
    for band in bands:
        k4P=calcK4P(array=band)
        if not aList:
            # alphaK is scalar
            kConvPsrc=calcSpireKcorr(getSpireRefFreq(array=band),\
               getSpireFreq(), getSpireFilt(array=band), BB=False, alpha=alphaK)[0]
            #point source colour correction for current alpha
            kColP[band] = kConvPsrc/k4P
            if (verbose): print 'Calculated %s KColP for alpha=%f: '%(band,alphaK),kColP[band]
        else:
            # alphaK is list
            kColP[band] = Double1d(na,Double.NaN)

            for a in range(na):
                kConvPsrc=calcSpireKcorr(getSpireRefFreq(array=band),\
                   getSpireFreq(), getSpireFilt(array=band), BB=False, alpha=alphaK[a])[0]
                #point source colour correction for current alpha
                kColP[band][a] = kConvPsrc/k4P
                if (verbose): print 'Calculated %s KColP for alpha=%f: '%(band,alphaK[a]),kColP[band][a]
            
    if not table:
        #return as is
        if array:
            return kColP[array]
        else:
            return kColP
    else:
        #create and returnTableDataset
        kColP_table=TableDataset()
        kColP_table.setDescription("Point Source Colour Correction (Spectral Index)")
        kColP_table.addColumn("alpha",Column(Double1d(alphaK)))
        for band in bands:
            kColP_table.addColumn(band,Column(kColP[band]))
        return kColP_table


def calcKColP_BB(betaK,tempK,array=None,verbose=False,table=False):
    """
    ================================================================================
    calcKColP_BB(betaK,tempK,array=None,verbose=False,table=False):
        Calculated point source colour correction factor for modified Blackbody spectrum

    Inputs:
        betaK:   (float) Emmisivity spectral index to use
        tempK:   (scalar/float array) Temperature(s) to use
        array:   (string) OPTIONAL. SPIRE array to use ('PSW'|'PMW'|'PLW')
                    Default=None (returns all three bands in dict)
        verbose: (boolean) Set to print more information. Default=True.
        table:   (boolean) Set to output TableDataset. Default=False.
    Outputs:
        If table==False:
            if array set:
                (scalar/list) colour correction factor for array
                  (scalar/list, depending on whether tempK is scalar/list)
            if array not set:
                (dict) colour correction factor (one scalar/list per band,
                  depending on whether tempK is scalar/list)
        If table==True:
            (TableDataset) colour correction factor (one column per band,
              one row per tempK value)
                    
    Calculation:
        [See Handbook Eq 5.29]
        Calculates K4P parameter for pipeline (alpha=-1)
        Calculates point source conversion for input spectrum(s)
        Divides to calculate colour correction factor
        If table is set, returns a TableDataset (betaK in metadata)
        KColP converts from monochromatic point source flux density for
          pipeline spectrum to monochromatic point source flux density for input
          spectrum.
    
    Dependencies:print 'Calculated area for alpha=%f: '%alphaK[a],effarea["PSW"][a],effArea["PMW"][a],effArea["PLW"][a]
        spireBands() - list of spire bands
        calcK4P() - calculate K4P parameter
        calcSpireKcorr() - calculate spire flux conversion factor
        herschel.ia.dataset.TableDataset()
    """
    #-----------------------------------------------------------------------
    #print 'Calculating point source colour correction parameters over beta & temp...'
    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #freq=getSpireFreq()
    #spireFilt=getSpireFilt(cal=cal)

    if array:
        #check array is valid
        assert (array=='PSW' or array=='PMW' or array=='PLW'),\
          'array must be "PSW"|"PMW"|"PLW" or None'
        bands=[array]
    else:
        #use all bands
        bands=spireBands()
        
    if verbose: print 'Calculating point source colour corrections for modified black body for ',bands

    try:
        nt=len(tempK)
        tList=True
    except:
        tList=False
        
    kColPBB = {}
    for band in bands:
        k4P=calcK4P(array=band)
        if not tList:
            kConvPsrc=calcSpireKcorr(getSpireRefFreq(array=band),\
                getSpireFreq(), getSpireFilt(array=band), BB=True, beta=betaK, temp=tempK)[0]
            #point source colour correction for current beta,temp
            kColPBB[band] = kConvPsrc/k4P
            if (verbose): print 'Calculated %s KColP for modBB with T=%f K, beta=%f: '%(band,tempK,betaK),kColPBB[band]
        
        else:
            kColPBB[band] = Double1d(nt,Double.NaN)
            for t in range(nt):
                kConvPsrc=calcSpireKcorr(getSpireRefFreq(array=band),\
                    getSpireFreq(), getSpireFilt(array=band), BB=True, beta=betaK, temp=tempK[t])[0]
                #point source colour correction for current beta,temp
                kColPBB[band][t] = kConvPsrc/k4P
                if (verbose): print 'Calculated %s KColP for modBB with T=%f K, beta=%f: '%(band,tempK[t],betaK),kColPBB[band][t]
    
    if not table:
        #return as is
        if array:
            return kColPBB[array]
        else:
            return kColPBB
    else:
        #create and returnTableDataset
        kColPBB_table=TableDataset()
        kColPBB_table.setDescription("Point Source Colour Correction (Modified Black Body, beta=%.2f)"%betaK)
        kColPBB_table.meta['beta']=DoubleParameter(betaK,"Emissivity spectral index")
        kColPBB_table.addColumn("Temperature",Column(Double1d(tempK),unit=Temperature.KELVIN,description=''))
        for band in bands:
            kColPBB_table.addColumn(band,Column(kColPBB[band]))
        return(kColPBB_table)
#-----------------------------------------------------------------------
#=======================================================================
#=====         CALCULATE EXTENDED SOURCE COLOR CORRECTIONS         =====
#=======================================================================
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

def calcKColE(alphaK,array=None,verbose=False,table=False):
    """
    ================================================================================
    calcKColE(alphaK,array=None,verbose=False,table=False):
        Calculated full-extended ource colour correction factor for power law spectrum

    Inputs:
        alphaK:  (scalar/float array) Spectral index(es) to use
        array:   (string) OPTIONAL. SPIRE array to use ('PSW'|'PMW'|'PLW')
                    Default=None (returns all three bands in dict)
        verbose: (boolean) Set to print more information. Default=True.
        table:   (boolean) Set to output TableDataset. Default=False.
    Outputs:
        If table==False:
            if array set:
                (scalar/list) colour correction factor for array
                  (scalar/list, depending on whether alphaK is scalar/list)
            if array not set:
                (dict) colour correction factor (one scalar/list per band,
                  depending on whether alphaK is scalar/list)
        If table==True:
            (TableDataset) colour correction factor (one column per band,
              one row per alphaK value)
                    
    Calculation:
        [See Handbook Eq 5.27]
        Calculates KMonE parameter for pipeline (alpha=-1)
        Calculates full-extended source conversion for input spectrum(s)
        Divides to calculate colour correction factor
        If table is set, returns a TableDataset
        KColE converts from monochromatic fully-extended surface brightness for
          pipeline spectrum to monochromatic fully-extended surface brightness for input
          spectrum.
    
    Dependencies:
        spireBands() - list of spire bands
        calcKMonE() - calculate KMonE parameter
        calcSpireKcorr() - calculate spire flux conversion factor
        herschel.ia.dataset.TableDataset()
    """
    #-----------------------------------------------------------------------
    #print '\nCalculating extended source colour correction parameters over alpha...'
    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #freq=getSpireFreq()
    #spireFilt=getSpireFilt(cal)

    if array:
        #check array is valid
        assert (array=='PSW' or array=='PMW' or array=='PLW'),\
          'array must be "PSW"|"PMW"|"PLW" or None'
        bands=[array]
    else:
        #use all bands
        bands=spireBands()
        
    if verbose: print 'Calculating extended source colour corrections for power law for ',bands

    try:
        na=len(alphaK)
        aList=True
    except:
        aList=False

    kColE={}
    for band in bands:
        k4E_Tot=calcKMonE(array=band)
        if not aList:
            # alphaK is scalar
            k4EaTot_x=calcSpireKcorr(getSpireRefFreq(array=band), getSpireFreq(),\
              getSpireFilt(array=band), BB=False, alpha=alphaK,\
              ext=True, monoArea=calcBeamMonoArea(array=band))[0]/1.e6
            kColE[band] = k4EaTot_x / k4E_Tot
            if (verbose): print 'Calculated %s KColE for alpha=%f: '%(band,alphaK),kColE[band]
        else:
            # alphaK is list
            kColE[band] = Double1d(na,Double.NaN)
            for a in range(na):
                k4EaTot_x=calcSpireKcorr(getSpireRefFreq(array=band), getSpireFreq(),\
                  getSpireFilt(array=band), BB=False, alpha=alphaK[a],\
                  ext=True, monoArea=calcBeamMonoArea(array=band))[0]/1.e6
                kColE[band][a] = k4EaTot_x / k4E_Tot
                if (verbose): print 'Calculated %s KColE for alpha=%f: '%(band,alphaK[a]),kColE[band][a]

    if not table:
        #return as is
        if array:
            return kColE[array]
        else:
            return kColE
    else:
        #create and returnTableDataset
        kColE_table=TableDataset()
        kColE_table.setDescription("Extended Source Colour Correction (Spectral Index)")
        kColE_table.addColumn("alpha",Column(Double1d(alphaK)))
        for band in bands:
            kColE_table.addColumn(band,Column(kColE[band]))
        return kColE_table
#-----------------------------------------------------------------------

def calcKColE_BB(betaK,tempK,array=None,verbose=False,table=False):
    """
    ================================================================================
    calcKColE_BB(betaK,tempK,array=None,verbose=False,table=False):
        Calculated full-extended ource colour correction factor for mod-BB spectrum

    Inputs:
        betaK:   (float) Emmisivity spectral index to use
        tempK:   (scalar/float array) Temperature(s) to use
        array:   (string) OPTIONAL. SPIRE array to use ('PSW'|'PMW'|'PLW')
                    Default=None (returns all three bands in dict)
        verbose: (boolean) Set to print more information. Default=True.
        table:   (boolean) Set to output TableDataset. Default=False.
    Outputs:
        If table==False:
            if array set:
                (scalar/list) colour correction factor for array
                  (scalar/list, depending on whether tempK is scalar/list)
            if array not set:
                (dict) colour correction factor (one scalar/list per band,
                  depending on whether tempK is scalar/list)
        If table==True:
            (TableDataset) colour correction factor (one column per band,
              one row per tempK value)
                    
    Calculation:
        [See Handbook Eq 5.30]
        Calculates KMonE parameter for pipeline (alpha=-1)
        Calculates full-extended source conversion for input spectrum(s)
        Divides to calculate colour correction factor
        If table is set, returns a TableDataset (betaK in metadata)
        KColE converts from monochromatic fully-extended surface brightness for
          pipeline spectrum to monochromatic fully-extended surface brightness for input
          spectrum.
    
    Dependencies:
        spireBands() - list of spire bands
        calcKMonE() - calculate KMonE parameter
        calcSpireKcorr() - calculate spire flux conversion factor
        herschel.ia.dataset.TableDataset()
    """
    #
    #print 'Calculating extended source colour correction parameters over beta & temp...'
    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #freq=getSpireFreq()
    #spireFilt=getSpireFilt(cal=cal)

    if array:
        #check array is valid
        assert (array=='PSW' or array=='PMW' or array=='PLW'),\
          'array must be "PSW"|"PMW"|"PLW" or None'
        bands=[array]
    else:
        #use all bands
        bands=spireBands()
        
    if verbose: print 'Calculating extended source colour corrections for modified black body for ',bands


    try:
        nt=len(tempK)
        tList=True
    except:
        tList=False

    kColEBB = {}
    for band in bands:
        k4E_Tot=calcKMonE(array=band)
        if not tList:
            k4EbTot_x=calcSpireKcorr(getSpireRefFreq(array=band), getSpireFreq(),\
                  getSpireFilt(array=band), BB=True, beta=betaK, temp=tempK,\
                  ext=True, monoArea=calcBeamMonoArea(array=band))[0]/1.e6
            kColEBB[band] = k4EbTot_x / k4E_Tot
            if (verbose): print 'Calculated %s KColE for modBB with T=%f K, beta=%f: '%(band,tempK,betaK),kColEBB[band]
        else:
            kColEBB[band] = Double1d(nt,Double.NaN)
            for t in range(nt):
                k4EbTot_x=calcSpireKcorr(getSpireRefFreq(array=band), getSpireFreq(),\
                      getSpireFilt(array=band), BB=True, beta=betaK, temp=tempK[t],\
                      ext=True, monoArea=calcBeamMonoArea(array=band))[0]/1.e6
                kColEBB[band][t] = k4EbTot_x / k4E_Tot
                if (verbose): print 'Calculated %s KColE for modBB with T=%f K, beta=%f: '%(band,tempK[t],betaK),kColEBB[band][t]

    if not table:
        #return as is
        if array:
            return kColEBB[array]
        else:
            return kColEBB
    else:
        #create and returnTableDataset
        kColEBB_table=TableDataset()
        kColEBB_table.setDescription("Extended Source Colour Correction (Modified Black Body, beta=%.2f)"%betaK)
        kColEBB_table.meta['beta']=DoubleParameter(betaK,"Emissivity spectral index")
        kColEBB_table.addColumn("Temperature",Column(Double1d(tempK),unit=Temperature.KELVIN,description=''))
        for band in spireBands():
            kColEBB_table.addColumn(band,Column(kColEBB[band]))
        return(kColEBB_table)

def calcApCorr(alphaK,aperture=None,annulus=[60.,90],array=None,verbose=False,table=False):
    """
    ================================================================================
    calcApCorr(alphaK,aperture=[22., 30.,45.],annulus=[60.,90],array=None,verbose=False,table=False):
        Calculates aperture corrections for point source with power law spectrum
        Returns corrections both with and without background annulus

    Inputs:
        alphaK:   (float) power law spectral index to use
        aperture: (scalar/list float) source apperture radius (in arcsec) for array
                    or for each band.
                    Default=None (uses default of 22/30/45 for PSW/PMW/PLW)
        annulus:  (list float) background annulus radius (in arcsec) for each band
                    Either 2-element list (uses same for each band), or list of 
                    three 2-element lists. Default=[60.,90.]
        array:    (string) OPTIONAL. SPIRE array to use ('PSW'|'PMW'|'PLW')
                     Default=None (returns all three bands in dict)
        verbose:  (boolean) Set to print more information. Default=True.
        table:    (boolean) Set to output TableDataset. Default=False.
    Outputs:
        If table==False:
            if array set:
                (scalar/list) aperture correction INCLUDING background annulus for array
                  (one scalar/list per band, depending on whether alphaK is scalar/list)
                (scalar/list) aperture correction WITHOUT background annulus for array
                  (scalar/list per band, depending on whether alphaK is scalar/list)
            if array not set:
                (dict) aperture correction INCLUDING background annulus (one scalar/list
                  per band, depending on whether alphaK is scalar/list)
                (dict) aperture correction WITHOUT background annulus (one scalar/list
                  per band, depending on whether alphaK is scalar/list)
        If table==True:
            (TableDataset) aperture correction INCLUDING background annulus
              (one column per band, one row per alphaK value)
            (TableDataset) aperture correction WITHOUT background annulus
              (one column per band, one row per alphaK value)
                    
    Calculation:
        Calculates effective beam profile and area for given spectrum (can be slow)
        Integrates profile over aperture and annulus
        Calculates aperture correction factors, including and excluding background annulus
    
    Dependencies:
        spireBands() - list of spire bands
        getCal() - retrieve photometer calibration product
        herschel.ia.numeric.toolbox.interp.CubicSplineInterpolator
        herschel.ia.numeric.toolbox.integr.TrapezoidalIntegrator
        herschel.ia.dataset.TableDataset()
    """
    
    if array:
        #check array is valid
        assert (array=='PSW' or array=='PMW' or array=='PLW'),\
          'array must be "PSW"|"PMW"|"PLW" or None'
        bands=[array]
    else:
        #use all bands
        bands=spireBands()
        
    if verbose: print 'Calculating aperture corrections for power law for ',bands

    #read aperture from input
    apDef={'PSW':22.,'PMW':30.,'PLW':45.}
    apPhotRad={}
    if array:
        if aperture==None: aperture=apDef[array]
        assert type(aperture)==float,'aperture must be scalar float for single array'
        apPhotRad[array]=aperture
    else:
        if aperture==None: aperture=[apDef['PSW'],apDef['PMW'],apDef['PLW']]
        assert type(aperture)==list,'aperture must be 3-element list'
        assert len(aperture)==3,'aperture must be 3-element list'
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
        for band in bands:
            apPhotBGRad[band]=annulus
    else:
        #use different for each band
        assert array==None,'Cannot set multiple sets of annulus radii for single array'
        for b in range(3):
            band=spireBands()[b]
            apPhotBGRad[band]=annulus[b]
            assert type(apPhotBGRad[band])==list,'annulus[%d] is not 2-element list: '+str(annulus[b])
            assert len(apPhotBGRad[band])==2,'annulus[%d] is not 2-element list: '+str(annulus[b])

    #check they are valid
    assert len(apPhotRad)==len(bands),'Error reading from aperture: '+str(aperture)
    assert len(apPhotBGRad)==len(bands),'Error reading from annulus: '+str(annulus)
    
    try:
        na=len(alphaK)
        aList=True
    except:
        aList=False
    
    #get beam profile from cal product
    beamProfs=getCal().getProduct('RadialCorrBeam')
    beamRad=getBeamRad()
    gamma=beamProfs.meta['gamma'].double
    maxRad=MAX(beamRad)
    #setup integrator over whole beam
    integTot=TrapezoidalIntegrator(0.,maxRad)

    #create interpolator for computing aperture/annulus size
    sizeInterp = CubicSplineInterpolator(beamRad,2.*PI*beamRad)
    apCorrNoBG = {}
    apCorrIncBG = {}
    for band in bands:
        #create integrators
        integAp=TrapezoidalIntegrator(0.,apPhotRad[band])
        integBG=TrapezoidalIntegrator(apPhotBGRad[band][0],apPhotBGRad[band][1])
        integTot=TrapezoidalIntegrator(0.,maxRad)
        #integrate to find size of aperture/annulus
        sizeAp=integAp.integrate(sizeInterp)
        sizeBG=integBG.integrate(sizeInterp)
        if verbose:print '%s aperture/background sizes:'%band,sizeAp,sizeBG
        if not aList:
            #alphaK is scalar
            if verbose:print 'Constructing effective %s beam profile for alpha=%f:'%(band,alphaK)
            effBeam=calcSpireEffBeam(alphaK,array=band,verbose=verbose)
            
            #interpolate effective beam profile
            beamInterp = CubicSplineInterpolator(beamRad,effBeam * 2.*PI*beamRad)
            
            #integrate profile to find beam area in aperture/annulus
            omegaAp=integAp.integrate(beamInterp)
            omegaBG=integBG.integrate(beamInterp)
            omegaTot=integTot.integrate(beamInterp)
            if verbose:print 'aperture/background beam areas:',omegaAp,omegaBG
            #calculate apertur  ecorrections
            apCorrNoBG[band] = \
                omegaTot/omegaAp    
            apCorrIncBG[band] =  \
                omegaTot/(omegaAp - omegaBG*sizeAp/sizeBG)
            if (verbose):
                print 'Calculated %s apCorr (noBG) for alpha=%f: '%(band,alphaK),apCorrNoBG[band]
                print 'Calculated %s apCorr (incBG) for alpha=%f: '%(band,alphaK),apCorrIncBG[band]
        else:
            #alphaK is list
            apCorrNoBG[band] = Double1d(na,Double.NaN)
            apCorrIncBG[band] = Double1d(na,Double.NaN)
        
            for a in range(na):
                if verbose:print 'Constructing effective beam profiles for %f:'%(alphaK[a])
                effBeam=calcSpireEffBeam(alphaK[a],array=band,verbose=verbose)
                

                #interpolate effective beam profile
                beamInterp = CubicSplineInterpolator(beamRad,effBeam * 2.*PI*beamRad)

                #integrate profile to find beam area in aperture/annulus
                omegaAp=integAp.integrate(beamInterp)
                omegaBG=integBG.integrate(beamInterp)
                omegaTot=integTot.integrate(beamInterp)
                #calculate aperture corrections
                apCorrNoBG[band][a] = \
                    omegaTot/omegaAp    
                apCorrIncBG[band][a] =  \
                    omegaTot/(omegaAp - omegaBG*sizeAp/sizeBG)
            if (verbose):
                print 'Calculated %s apCorr (noBG) for alpha=%f: '%(band,alphaK[a]),apCorrNoBG[band][a]
                print 'Calculated %s apCorr (incBG) for alpha=%f: '%(band,alphaK[a]),apCorrIncBG[band][a]
        
    if not table:
        #return as is
        if array:
            return(apCorrIncBG[array],apCorrNoBG[array])
        else:
            return(apCorrIncBG,apCorrNoBG)
    else:
        #create and returnTableDataset
        apCorrNoBG_table=TableDataset()
        apCorrIncBG_table=TableDataset()
        effArea_table=TableDataset()
        apCorrNoBG_table.setDescription("Aperture Correction without background annulus (Spectral Index)")
        apCorrIncBG_table.setDescription("Aperture Correction including background annulus (Spectral Index)")
        effArea_table.setDescription("Effective area from Effective beam (Spectral Index)")
        apCorrNoBG_table.addColumn("alpha",Column(Double1d(alphaK)))
        apCorrIncBG_table.addColumn("alpha",Column(Double1d(alphaK)))
        effArea_table.addColumn("alpha",Column(Double1d(alphaK)))
        for band in bands:
            apCorrNoBG_table.addColumn(band,Column(apCorrNoBG[band]))
            apCorrIncBG_table.addColumn(band,Column(apCorrIncBG[band]))
            effArea_table.addColumn(band,Column(effArea[band]))
        return (apCorrIncBG_table,apCorrNoBG_table)
    
def calcApCorr_BB(betaK,tempK,aperture=None,annulus=[60.,90],array=None,verbose=False,table=False):
    """
    ================================================================================
    calcApCorr(betaK,tempK,aperture=None,annulus=[60.,90],array=None,verbose=False,table=False):
        Calculates aperture corrections for point source with modified black-body spectrum
        Returns corrections both with and without background annulus

    Inputs:
        betaK:   (float) Emmisivity spectral index to use
        tempK:   (scalar/float array) Temperature(s) to use
        aperture: (list float) source apperture radius (in arcsec) for each band.
                    Default=None (uses default of 22/30/45 for PSW/PMW/PLW)
        annulus:  (list float) background annulus radius (in arcsec) for each band
                    Either 2-element list (uses same for each band), or list of 
                    three 2-element lists. Default=[60.,90.]
        array:    (string) OPTIONAL. SPIRE array to use ('PSW'|'PMW'|'PLW')
                     Default=None (returns all three bands in dict)
        verbose:  (boolean) Set to print more information. Default=False.
        table:    (boolean) Set to output TableDataset. Default=False.
    Outputs:
        If table==False:
            if array set:
                (scalar/list) aperture correction INCLUDING background annulus for array
                  (one scalar/list per band, depending on whether tempK is scalar/list)
                (scalar/list) aperture correction WITHOUT background annulus for array
                  (scalar/list per band, depending on whether tempK is scalar/list)
            if array not set:
                (dict) aperture correction INCLUDING background annulus (one scalar/list
                  per band, depending on whether tempK is scalar/list)
                (dict) aperture correction WITHOUT background annulus (one scalar/list
                  per band, depending on whether tempK is scalar/list)
        If table==True:
            (TableDataset) aperture correction INCLUDING background annulus
              (one column per band, one row per tempK value)
            (TableDataset) aperture correction WITHOUT background annulus
              (one column per band, one row per tempK value)
                    
    Calculation:
        Calculates effective beam profile and area for given spectrum (can be slow)
        Integrates profile over aperture and annulus
        Calculates aperture correction factors, including and excluding background annulus
    
    Dependencies:
        spireBands() - list of spire bands
        getCal() - retrieve photometer calibration product
        herschel.ia.numeric.toolbox.interp.CubicSplineInterpolator
        herschel.ia.numeric.toolbox.integr.TrapezoidalIntegrator
        herschel.ia.dataset.TableDataset()
    """
    
    if array:
        #check array is valid
        assert (array=='PSW' or array=='PMW' or array=='PLW'),\
          'array must be "PSW"|"PMW"|"PLW" or None'
        bands=[array]
    else:
        #use all bands
        bands=spireBands()
        
    if verbose: print 'Calculating aperture corrections for modified black body for ',bands

    #read aperture from input
    apDef={'PSW':22.,'PMW':30.,'PLW':45.}
    apPhotRad={}
    if array:
        if aperture==None: aperture=apDef[array]
        assert type(aperture)==float,'aperture must be scalar float for single array'
        apPhotRad[array]=aperture
    else:
        if aperture==None: aperture=[apDef['PSW'],apDef['PMW'],apDef['PLW']]
        assert type(aperture)==list,'aperture must be 3-element list'
        assert len(aperture)==3,'aperture must be 3-element list'
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
        for band in bands:
            apPhotBGRad[band]=annulus
    else:
        #use different for each band
        assert array==None,'Cannot set multiple sets of annulus radii for single array'
        for b in range(3):
            band=spireBands()[b]
            apPhotBGRad[band]=annulus[b]
            assert type(apPhotBGRad[band])==list,'annulus[%d] is not 2-element list: '+str(annulus[b])
            assert len(apPhotBGRad[band])==2,'annulus[%d] is not 2-element list: '+str(annulus[b])

    #check they are valid
    assert len(apPhotRad)==len(bands),'Error reading from aperture: '+str(aperture)
    assert len(apPhotBGRad)==len(bands),'Error reading from annulus: '+str(annulus)
        
    try:
        nt=len(tempK)
        tList=True
    except:
        tList=False
    
    #get beam profile from cal product
    beamProfs=getCal().getProduct('RadialCorrBeam')
    beamRad=getBeamRad()
    gamma=beamProfs.meta['gamma'].double
    maxRad=MAX(beamRad)
    #setup integrator over whole beam
    integTot=TrapezoidalIntegrator(0.,maxRad)

    #create interpolator for computing aperture/annulus size
    sizeInterp = CubicSplineInterpolator(beamRad,2.*PI*beamRad)
    apCorrNoBG_BB = {}
    apCorrIncBG_BB = {}
    for band in bands:
        #create integrators
        integAp=TrapezoidalIntegrator(0.,apPhotRad[band])
        integBG=TrapezoidalIntegrator(apPhotBGRad[band][0],apPhotBGRad[band][1])
        integTot=TrapezoidalIntegrator(0.,maxRad)
        #integrate to find size of aperture/annulus
        sizeAp=integAp.integrate(sizeInterp)
        sizeBG=integBG.integrate(sizeInterp)

        if not tList:
        #tempK is scalar
            if verbose:print 'Constructing effective %s beam profile for beta=%f, Temp=%f:'%(band,betaK,tempK)
            effBeam=calcSpireEffBeam_BB(betaK,tempK,array=band,verbose=verbose)
            
            #interpolate effective beam profile
            beamInterp = CubicSplineInterpolator(beamRad,effBeam * 2.*PI*beamRad)
            
            #integrate profile to find beam area in aperture/annulus
            omegaAp=integAp.integrate(beamInterp)
            omegaBG=integBG.integrate(beamInterp)
            omegaTot=integTot.integrate(beamInterp)
            if verbose:print 'aperture/background beam areas:',omegaAp,omegaBG
            #calculate apertur  ecorrections
            apCorrNoBG_BB[band] = \
                omegaTot/omegaAp    
            apCorrIncBG_BB[band] =  \
                omegaTot/(omegaAp - omegaBG*sizeAp/sizeBG)
            if (verbose):
                print 'Calculated %s apCorr (noBG) for beta %f, %f K: '%(band,betaK,tempK),apCorrNoBG_BB[band]
                print 'Calculated %s apCorr (incBG) for beta %f, %f K: '%(band,betaK,tempK),apCorrIncBG_BB[band]
        else:
            #tempK is list
            apCorrNoBG_BB[band] = Double1d(nt,Double.NaN)
            apCorrIncBG_BB[band] = Double1d(nt,Double.NaN)
        
            for t in range(nt):
                if verbose:print 'Constructing %s effective beam profile beta=%f, Temp=%f:'%(band,betaK,tempK[t])
                effBeam=calcSpireEffBeam_BB(betaK,tempK[t],array=band,verbose=verbose)
                #interpolate effective beam profile
                beamInterp = CubicSplineInterpolator(beamRad,effBeam * 2.*PI*beamRad)
            
                #integrate profile to find beam area in aperture/annulus
                omegaAp=integAp.integrate(beamInterp)
                omegaBG=integBG.integrate(beamInterp)
                omegaTot=integTot.integrate(beamInterp)
                #calculate aperture corrections
                apCorrNoBG_BB[band][t] = \
                    omegaTot/omegaAp    
                apCorrIncBG_BB[band][t] =  \
                    omegaTot/(omegaAp - omegaBG*sizeAp/sizeBG)
                if (verbose):
                    print 'Calculated %s apCorr (noBG) for beta %f, %f K: '%(band,betaK,tempK[t]),apCorrNoBG_BB[band][t]
                    print 'Calculated %s apCorr (incBG) for beta %f, %f K: '%(band,betaK,tempK[t]),apCorrIncBG_BB[band][t]
        
    if not table:
        #return as is
        if array:
            return(apCorrIncBG_BB[band],apCorrNoBG_BB[band])
        else:
            return(apCorrIncBG_BB,apCorrNoBG_BB)
    else:
        #create and returnTableDataset
        apCorrNoBG_BB_table=TableDataset()
        apCorrIncBG_BB_table=TableDataset()
        apCorrNoBG_BB_table.setDescription("Aperture Correction without background annulus (Modified Black Body, beta=%.2f)"%betaK)
        apCorrIncBG_BB_table.setDescription("Aperture Correction including background annulus (Modified Black Body, beta=%.2f)"%betaK)
        apCorrNoBG_BB_table.meta['beta']=DoubleParameter(betaK,"Emissivity spectral index")
        apCorrIncBG_BB_table.meta['beta']=DoubleParameter(betaK,"Emissivity spectral index")
        apCorrNoBG_BB_table.addColumn("Temperature",Column(Double1d(tempK),unit=Temperature.KELVIN,description=''))
        apCorrIncBG_BB_table.addColumn("Temperature",Column(Double1d(tempK),unit=Temperature.KELVIN,description=''))
        for band in bands:
            apCorrNoBG_BB_table.addColumn(band,Column(apCorrNoBG_BB[band]))
            apCorrIncBG_BB_table.addColumn(band,Column(apCorrIncBG_BB[band]))
        return (apCorrIncBG_BB_table,apCorrNoBG_BB_table)

def clear(verbose=False):
    #deletes global variables

    global spireCalPhot,spireBeamRad,spireFreq,spireRefFreq,spireEffFreq
    global spireFilt,spireFiltOnly,beamMonoArea,spireEffBeams
    try:
        del spireCalPhot
        if verbose: print 'Deleted spireCalPhot'
    except:
        if verbose: print 'Cannot delete spireCalPhot'
    try:
        del spireBeamRad
        if verbose: print 'Deleted spireBeamRad'
    except:
        if verbose: print 'Cannot delete spireBeamRad'
    try:
        del spireFreq
        if verbose: print 'Deleted spireFreq'
    except:
        if verbose: print 'Cannot delete spireFreq'
    try:
        del spireRefFreq
        if verbose: print 'Deleted spireRefFreq'
    except:
        if verbose: print 'Cannot delete spireRefFreq'
    try:
        del spireEffFreq
        if verbose: print 'Deleted spireEffFreq'
    except:
        if verbose: print 'Cannot delete spireEffFreq'
    try:
        del spireFilt
        if verbose: print 'Deleted spireFilt'
    except:
        if verbose: print 'Cannot delete spireFilt'
    try:
        del spireFiltOnly
        if verbose: print 'Deleted spireFiltOnly'
    except:
        if verbose: print 'Cannot delete spireFiltOnly'
    try:
        del beamMonoArea
        if verbose: print 'Deleted beamMonoArea'
    except:
        if verbose: print 'Cannot delete beamMonoArea'
    try:
        del spireEffBeams
        if verbose: print 'Deleted spireEffBeams'
    except:
        if verbose: print 'Cannot delete spireEffBeams'


def spireColCorrTest(beamType='Full',array=None,verbose=False):
    """
    ================================================================================
    spireColCorrTest(beamType='Full',verbose=False):
        Test script to demonstrate use of colour correction calculations

    Inputs:
        beamType: (string) beam calculation to use ('Full'|'Simple')
                    Default='Full'
        array:    (string) OPTIONAL. SPIRE array to use ('PSW'|'PMW'|'PLW')
                     Default=None (returns all three bands in dict)
        verbose:  (boolean) Set to print more information. Default=False.
    Outputs:
        (dict double) beamMonoArea for chose beamType (one array per band)
        (double array) alphaK values used in test script
        (dict double) effective areas for alphaK values (one array per band)
        (dict double) beam corrections for alphaK values (one array per band)
        (dict double) KColP for alphaK values (one array per band)
        (dict double) KColE for alphaK values (one array per band)
        (dict double) apCorrIncBG for alphaK values (one array per band)
        (dict double) apCorrNoBG for alphaK values (one array per band)
                    
    Calculation:
        Calculates all the correction parameters
        Returns some for range of alphas.
    
    Dependencies:
        getCal() - get calibration tree
        calcBeamMonoArea() - calculate beam monochromatic area
        calcK4P() - calculate K4P parameter
        calcKMonE() - calculate KMonE parameter
        calcK4E() - calculate K4E parameter
        calcKPtoE() - calculate KPtoE parameter
        calcOmegaEff() - calculate effective beam area for powerl law spectrum
        calcOmegaEff_BB() - calculate effective beam area for mod-BB spectrum
        calcKBeam() - calculate beam correction for power law spectrum
        calcKBeam_BB() - calculate beam correction for power mod-BB spectrum
        calcKColP() - calculate point source colour correction for power law spectrum
        calcKColP_BB() - calculate point source colour correction for mod-BB spectrum
        calcKColE() - calculate fully-extended source colour correction for power law spectrum
        calcKColE_BB() - calculate fully-extended source colour correction for mod-BB spectrum
        calcApCorr() - calculate fully-extended source colour correction for power law spectrum
        calcApCorr_BB() - calculate fully-extended source colour correction for mod-BB spectrum
    """

    if array:
        #check array is valid
        assert (array=='PSW' or array=='PMW' or array=='PLW'),\
          'array must be "PSW"|"PMW"|"PLW" or None'
        bands=[array]
    else:
        #use all bands
        bands=spireBands()
    
    calphot=getCal()
    alphaArr=Float1d(range(-4,5)) #range of alphas to use.
    beamMonoArea=calcBeamMonoArea(beamType=beamType,array=array,verbose=True)

    print 'Testing pipeline parameters'
    print 'K4P=',calcK4P(array=array,verbose=verbose)
    print 'KMonE=',calcKMonE(array=array,verbose=verbose)
    print 'K4E=',calcK4E(array=array,verbose=verbose)
    print 'KPtoE=',calcKPtoE(array=array,verbose=verbose)
    print 'Omega(-1)=',calcOmegaEff(-1.0,array=array,verbose=verbose)
    alphaNep={'PSW':spireCalPhot.getProduct('RadialCorrBeam').meta['alphaNeptunePsw'].double,\
      'PMW':spireCalPhot.getProduct('RadialCorrBeam').meta['alphaNeptunePmw'].double,\
      'PLW':spireCalPhot.getProduct('RadialCorrBeam').meta['alphaNeptunePlw'].double}
    omegaNep={}
    for band in bands:
        omegaNep[band]=calcOmegaEff(alphaNep[band],array=band,verbose=verbose)
    print 'Omega_Nep=',omegaNep
    
    print '\nTesting calcOmegaEff'
    print 'OmegaEff(-2)=',calcOmegaEff(-2.0,array=array,verbose=verbose)
    omegaEff = calcOmegaEff(alphaArr,array=array,verbose=verbose)
    print 'OmegaEff_BB(1.75,20)=',calcOmegaEff_BB(1.75,20.0,array=array,verbose=verbose)
    omegaEff_BB = calcOmegaEff_BB(1.75,[20.0,30.,40.],array=array,verbose=verbose)

    print '\nTesting calcKBeam'
    print 'KBeam(-2)=',calcKBeam(-2.0,array=array,verbose=verbose)
    KBeam = calcKBeam(alphaArr,array=array,verbose=verbose)
    print 'KBeam(1.75,20)=',calcKBeam_BB(1.75,20.0,array=array,verbose=verbose)
    KBeam_BB = calcKBeam_BB(1.75,[20.0,30.,40.],array=array,verbose=verbose)

    print '\nTesting calcKColP'
    print 'KColP(-2)=',calcKColP(-2.0,array=array,verbose=verbose)
    KColP = calcKColP(alphaArr,array=array,verbose=verbose)
    print 'KColP(1.75,20)=',calcKColP_BB(1.75,20.0,array=array,verbose=verbose)
    KColP_BB = calcKColP_BB(1.75,[20.0,30.,40.],array=array,verbose=verbose)

    print '\nTesting calcKColE'
    print 'KColE(-2)=',calcKColE(-2.0,array=array,verbose=verbose)
    KColE = calcKColE(alphaArr,array=array,verbose=verbose)
    print 'KColE(1.75,20)=',calcKColE_BB(1.75,20.0,array=array,verbose=verbose)
    KColE_BB = calcKColE_BB(1.75,[20.0,30.,40.],array=array,verbose=verbose)
    
    print '\nTesting aperture corrections'
    apCorr=calcApCorr(-2.0,array=array,verbose=verbose)
    print 'Ap Corr [incBG] (-2)=',apCorr[0]
    print 'Ap Corr [noBG] (-2)=',apCorr[1]
    apCorr = calcApCorr(alphaArr,array=array,verbose=verbose)
    apCorrIncBG=apCorr[0]
    apCorrNoBG=apCorr[1]
    
    apCorr_BB = calcApCorr_BB(1.75,20.0,array=array,verbose=verbose)
    print 'Ap Corr [incBG] (1.75,20)=',apCorr_BB[0]
    print 'Ap Corr [noBG] (1.75,20)=',apCorr_BB[1]
    apCorr_BB = calcApCorr_BB(1.75,[20.0,30.,40.],array=array,verbose=verbose)
    apCorrIncBG_BB=apCorr_BB[0]
    apCorrNoBG_BB=apCorr_BB[1]
    
    return(beamMonoArea,alphaArr,omegaEff,KBeam,KColP,KColE,apCorrIncBG,apCorrNoBG)

#def compFullSimple():
#    """
#    ================================================================================
#    compFullSimple():
#        Test script to compare effect of full/simple beam treatments
#    Inputs:
#        NONE
#    Outputs:
#        8 x Plots
#                    
#    Calculation:
#        Calculates correction parameters for full beam treatment
#        Calculates correction parameters for simple beam treatment
#        Plots values and full/simple comparisons
#    
#    Dependencies:
#        spireColCorrTest() - test colour correction parameters
#    """
#    #compare full vs single beam treatments:
#    beamMonoArea_F,alphaArr,omegaEff_F,KBeam_F,KColP_F,KColE_F=spireColCorrTest(beamType='Full')
#    beamMonoArea_S,alphaArr,omegaEff_S,KBeam_S,KColP_S,KColE_S=spireColCorrTest(beamType='Simple')
#    pB=PlotXY()
#    pOE=PlotXY()
#    pKB=PlotXY()
#    pKCE=PlotXY()
#
#    pBd=PlotXY()
#    pOEd=PlotXY()
#    pKBd=PlotXY()
#    pKCEd=PlotXY()
#
#    cols={'PSW':java.awt.Color.BLUE,\
#      'PMW':java.awt.Color.GREEN,\
#      'PLW':java.awt.Color.RED}
#    for band in spireBands():
#        pB.addLayer(LayerXY(spireFreq/1.e9,beamMonoArea_F[band],name='Beam Area %s (full)'%band,\
#          color=cols[band]))
#        pB.addLayer(LayerXY(spireFreq/1.e9,beamMonoArea_S[band],name='Beam Area %s (simple)'%band,\
#          color=cols[band],line=Style.DASHED))
#        pBd.addLayer(LayerXY(spireFreq/1.e9,beamMonoArea_S[band]/beamMonoArea_F[band]-1.,name='Beam Area %s (s-f)'%band,\
#          color=cols[band]))
#        
#        pOE.addLayer(LayerXY(alphaArr,omegaEff_F[band],name='OmegaEff %s (full)'%band,\
#          color=cols[band]))
#        pOE.addLayer(LayerXY(alphaArr,omegaEff_S[band],name='OmegaEff %s (simple)'%band,\
#          color=cols[band],line=Style.DASHED))
#        pOEd.addLayer(LayerXY(alphaArr,omegaEff_S[band]/omegaEff_F[band]-1.,name='OmegaEff %s (s-f)'%band,\
#          color=cols[band]))
#    
#        pKB.addLayer(LayerXY(alphaArr,KBeam_F[band],name='KBeam %s (full)'%band,\
#          color=cols[band]))
#        pKB.addLayer(LayerXY(alphaArr,KBeam_S[band],name='KBeam %s (simple)'%band,\
#          color=cols[band],line=Style.DASHED))
#        pKBd.addLayer(LayerXY(alphaArr,KBeam_S[band]/KBeam_F[band]-1.,name='KBeam %s (s-f)'%band,\
#          color=cols[band]))
#    
#        pKCE.addLayer(LayerXY(alphaArr,KColE_F[band],name='KColE %s (full)'%band,\
#          color=cols[band]))
#        pKCE.addLayer(LayerXY(alphaArr,KColE_S[band],name='KColE %s (simple)'%band,\
#          color=cols[band],line=Style.DASHED))
#        pKCEd.addLayer(LayerXY(alphaArr,KColE_S[band]/KColE_F[band]-1.,name='KColE %s (s-g)'%band,\
#          color=cols[band]))
#
#    pB.xaxis.titleText = 'Frequency (GHz)'
#    pB.yaxis.titleText = 'Beam Area (sr)'
#    pBd.xaxis.titleText = 'Frequency (GHz)'
#    pBd.yaxis.titleText = 'Beam Area rel. diff'
#    
#    pOE.xaxis.titleText = 'Spectra Index (alpha)'
#    pOE.yaxis.titleText = 'Effective Beam Area (sr)'
#    pOEd.xaxis.titleText = 'Spectra Index (alpha)'
#    pOEd.yaxis.titleText = 'Effective Beam Area rel. diff'
#    
#    pKB.xaxis.titleText = 'Spectra Index (alpha)'
#    pKB.yaxis.titleText = 'KBeam'
#    pKBd.xaxis.titleText = 'Spectra Index (alpha)'
#    pKBd.yaxis.titleText = 'KBeam rel. diff'
#
#    pKCE.xaxis.titleText = 'Spectra Index (alpha)'
#    pKCE.yaxis.titleText = 'KColE'
#    pKCEd.xaxis.titleText = 'Spectra Index (alpha)'
#    pKCEd.yaxis.titleText = 'KColE rel. diff'

#def testApCorr():
#    #test aperture corrections
#    apCorrIncBGCal=getCal().getProduct('ColorCorrApertureList')[0]['alpha']
#    apCorrNoBGCal=getCal().getProduct('ColorCorrApertureList')[1]['alpha']
#    na=len(apCorrIncBGCal['alpha'].data)
#    alphaK=apCorrIncBGCal['alpha'].data[0:na]
#    result=calcApCorr(alphaK,verbose=True,table=True)
#    apCorrIncBG=result[0]
#    apCorrNoBG=result[1]
#    effArea=result[2]
#    omegaEff=calcOmegaEff(alphaK,verbose=True,table=True)
#    apCorrIncBG_diff=apCorrIncBG.copy()
#    apCorrNoBG_diff=apCorrNoBG.copy()
#    effArea_diff=effArea.copy()
#    for band in spireBands():
#        apCorrIncBG_diff[band].data = apCorrIncBG[band].data/apCorrIncBGCal[band].data[0:na]
#        apCorrNoBG_diff[band].data = apCorrNoBG[band].data/apCorrNoBGCal[band].data[0:na]
#        effArea_diff[band].data = effArea[band].data/omegaEff[band].data
#        
#    return(alphaK,[apCorrIncBGCal,apCorrNoBGCal,omegaEff],[apCorrIncBG,apCorrNoBG,effArea],[apCorrIncBG_diff,apCorrNoBG_diff,effArea_diff])