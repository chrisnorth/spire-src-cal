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
#    A script is provided at the end for usage of the functions
#
#===============================================================================
# 
#  Herschel-SPIRE Colour Corrections
# 
#  This module provides functions to calculate the colour corrections for a
#  range of source spectra, and for both point sources and extended sources
# 
#  The following global variables are defined and used:
#         spireCalPhot: [context] Spire Calibration Context for Photometer
#         spireBands:   [string list] SPIRE band names ["PSW","PMW","PLW"]
#         spireFreq:    [float array] frequency raster across all bands
#         spireRefFreq: [dict] SPIRE reference frequencies (1 scalar per band)
#         spireEffFreq: [dict] SPIRE effective frequencies (1 scalar per band)
#         spireFiltOnly:[dict] filter profile without aperture efficiency (1 array per band)
#         spireFilt:    [dict] filter profile with aperture efficiency (1 array per band)
#         beamMonoArea: [dict] SPIRE monochromatic areas (1 array per bands, plus 'beamType')
#         arcsec2Sr:    [float] conversion fron square arcseconds to steradians
#
#  The functions are as follows:
#  * getCal:
#     Gets calibration Context from pool or file, and/or checks existing cal
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
#         spireBands:   [string list] SPIRE band names ["PSW","PMW","PLW"]
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
#          arcsec2Sr: [float] conversion fron square arcseconds to steradians
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
#          tempK: [float/list)] modBB temperature (scalar or list)
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
from herschel.ia.dataset import TableDataset,DoubleParameter,Column
spireCal = SpireCalTask()
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
#=====                           INPUT PARAMETERS                          =====
#===============================================================================
#-------------------------------------------------------------------------------
#
## need the calibration tree 
##
#cal=spireCal(pool='spire_cal_12_2')
### alternatively, read from jarFile
##cal=spireCal(jarFile='spire_cal_12_2.jar')
#
##set verbose to print more information during processing
#verbose = True

#-------------------------------------------------------------------------------
#===============================================================================
#=====                       DEFINE GLOBAL VARIABLES                       =====
#===============================================================================
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#===============================================================================
#=====                           SETUP FUNCTIONS                           =====
#===============================================================================
#-------------------------------------------------------------------------------

def getCal(cal=None,calTree=None,calPool=None,calFile=None,verbose=False):
    # if cal not defined, read from pool or jarFile
    # if cal is defined, don't read anything new. Just check and return cal
    #print cal,calTree,calPool,calFile,verbose
    
    #set global variables
    global spireCalPhot #Photometer Calibration Context (spireCal.getPhot())
    global spireBands # list of SPIRE band names

    try:
        spireCalPhot
        #global variable already defined, so do nothing
    except:
        #get spireCalTree
        if cal==None and calTree==None and calPool==None and calFile==None:
            #no import provided. do default action
            defaultPool='spire_cal_12_2'
            print 'Reading from default pool: %s'%defaultPool
            cal=spireCal(pool=defaultPool)
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

    #define spireBands
    try:
        spireBands
    except:
        spireBands=['PSW','PMW','PLW']
    
    return(spireCalPhot)
    
def getSpireFreq():
    # Frequency raster of common frequency grid spanning all three spire bands

    #define spireBands
    try:
        spireBands
    except:
        spireBands=['PSW','PMW','PLW']
        
    #define global variable
    global spireFreq #raster of frequencies used throughout functions
    try:
        spireFreq
        #global variable already defined, so do nothing
    except:
        #not defined, so recalculate
        deltaNu = 0.1e9        # 0.1 GHz
        nuMin   = 300.e9
        nuMax   = 1800.e9
        nNu     = FLOOR((nuMax-nuMin)/deltaNu)
        spireFreq    = Double1d(range(nNu)) * deltaNu + nuMin
        
    return(spireFreq)
    
def getSpireRefFreq():

    #define global variable
    global spireRefFreq #spire reference frequencies (at 250, 350, 500 microns)

    try:
        spireRefFreq
        #global variable already defined, so do nothing
    except:
        #not defined, so recalculate
        spireRefFreq = {}
        c = Constant.SPEED_OF_LIGHT.value
        spireRefWl = {"PSW":250.*1e-6, "PMW":350.*1.e-6, "PLW":500.*1.e-6}
        spireBands=["PSW","PMW","PLW"]
        for band in spireBands:
            spireRefFreq[band] = c/spireRefWl[band]
     
    return(spireRefFreq)


def getSpireEffFreq():
    #get SPIRE Effective Frequencies from metadata

    #define global variable
    global spireEffFreq #spire Effective frequencies

    try:
        spireEffFreq
        #global variable already exists, so do nothing
    except:
        #doesn't exist, so get from calFile
        beamProfs = spireCalPhot.getProduct("RadialCorrBeam")
        spireEffFreq = {"PSW":beamProfs.meta['freqEffPsw'].double*1.e9,\
            "PMW":beamProfs.meta['freqEffPmw'].double*1.e9,\
            "PLW":beamProfs.meta['freqEffPlw'].double*1.e9}
    
    return(spireEffFreq)

def getSpireFilt(rsrfOnly=False):
    
    #define global variables
    global spireFiltOnly #spire filter profiles (without ap. eff.) over spireFreq
    global spireFilt #spire filter profiles (with ap. eff.) over spireFreq

    try:
        spireFiltOnly
        #global variable spireFiltOnly defined, so don't recalculate
    except:
        #spireFiltOnly not defined, so must calculate

        #check spireCalPhot exists
        assert spireCalPhot.isValid(), 'Invalid SPIRE Photometer calibration tree. Run getCalPhot()'

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
        # spire RSRF only
        spireFiltOnly={}
        #interpolate to freq array
    
        for band in spireBands:
            #create Rsrf and ApEff interpolation objects
            interpRsrf = LinearInterpolator(spireRsrfFreq, rsrf.getRsrf(band))
            #make arrays for final objects
            spireFiltOnly[band] = Double1d(nNu)
            #interpolate Rsrf to freq array
            spireFiltOnly[band][ixR] = interpRsrf(spireFreq[ixR])

    try:
        spireFilt
        #global variable spireFilt defined, so do nothing
    except:
        #spireFilt not defined, so must calculate
        #add in aperture efficiency
        apertureEfficiency = spireCalPhot.getProduct('ApertureEfficiency')
        apertureEfficiencyVersion=apertureEfficiency.getVersion()
        spireApEffFreq = apertureEfficiency.getApertEffTable()["frequency"].data * 1e9 #comes in [GHz]
        #indexes of freq in apEff
        spireFreq=getSpireFreq()
        nNu=len(spireFreq)
        ixA = spireFreq.where((spireFreq>=MAX(spireApEffFreq)) & (spireFreq<=MAX(spireApEffFreq)))
        # spire RSRF * ApEff
        spireFilt={}
        #interpolate to freq array and apply to RSRF
        for band in spireBands:
            #create ApEff interpolation objects
            interpAp = LinearInterpolator(spireApEffFreq, apertureEfficiency.getApertEffTable()[band].data)
            #make arrays for final objects
            spireFilt[band] = Double1d(nNu)
            #copy into Rsrf*ApEff array
            spireFilt[band] = spireFiltOnly[band].copy()
            spireFilt[band][ixA] = spireFilt[band][ixA] * interpAp(spireFreq[ixA])

    if rsrfOnly:
        #return spireFiltOnly
        return(spireFiltOnly)
    else:
        #return spireFilt
        return(spireFilt)

def calcBeamMonoArea(beamType=None,verbose=False):
    #calculate monochromatic beam areas using full or simple beam treatment
    #print '\nCalculating monochromatic beam areas...'

    #define global variables
    global arcsec2Sr #conversion from acrsec^2 to steradians
    global beamMonoArea #monochromatic beam areas over spireFreq

    try:
        arcsec2Sr
    except:
        #define global variable arcsec2Sr
        arcsec2Sr = (PI/(60.*60.*180))**2
    
    try:
        beamMonoArea['PSW'][0]
        docalc=False
        #global variable already defined, so do nothing.
    except:
        #not defined, so need to recalculate
        docalc=True
        recalc=False
    
    #print docalc
    #Check whether existing beamType matches specified
    if not docalc:
        if beamMonoArea['beamType']=='Simple' and beamType=='Full':
            if (verbose):print 'Replacing Simple Beams with Full Beams'
            recalc=True
        elif beamMonoArea['beamType']=='Full' and beamType=='Simple':
            if (verbose):print 'Replacing Full Beams with Simple Beams'
            recalc=True
        elif beamType==None:
            if (verbose):print 'No beamType specified, so using existing %s beams'%beamMonoArea['beamType']
            recalc=False
        else:
            if (verbose):print 'Using existing %s beams'%beamMonoArea['beamType']
            recalc=False
    if recalc or docalc:

        #recalculate beam profiles
        beamProfs = spireCalPhot.getProduct("RadialCorrBeam")
        gamma = beamProfs.meta['gamma'].double
        beamMonoArea = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}

        if beamType=='Simple':
            beamMonoArea['beamType']='Simple'
            #use simple beam treatment
            neptuneArea = {}
            neptuneArea['PSW']=beamProfs.meta['beamNeptunePswSr'].double
            neptuneArea['PMW']=beamProfs.meta['beamNeptunePmwSr'].double
            neptuneArea['PLW']=beamProfs.meta['beamNeptunePlwSr'].double
            for band in spireBands:
                if (verbose):print 'Calculating monochromatic beam areas for %s band (Simple treatment)'%band
                beamMonoArea[band]=neptuneArea[band] * (getSpireFreq()/getSpireEffFreq()[band])**(2.*gamma)

    
        else:
            #use full beam treatment
            beamMonoArea['beamType']='Full'
            for band in spireBands:
                if (verbose):print 'Calculating monochromatic beam areas for %s band (Full treatment)'%band
                #monochromatic beam areas
                beamMonoArea[band] = spireMonoAreas(getSpireFreq(), beamProfs, 
                  getSpireEffFreq()[band], gamma, band)

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
    numInterp=LinearInterpolator(freq,transm * fSky * monoArea)
    denomInterp=LinearInterpolator(freq,transm * fSky)
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
    
    beamInterp=LinearInterpolator(beamRad,beamNew * 2. * PI * beamRad)
    integrator=TrapezoidalIntegrator(0,maxRad)
    beamMonoArea=integrator.integrate(beamInterp)

    return (beamMonoArea,beamNew)

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

    """

    arcsec2Sr = (PI/(60.*60.*180))**2

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
    beamMonoArea=beamInterp(freq)*arcsec2Sr #in sr
    
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
    Calculation of the K-correction factor from isophotal flux to a monochromatic 
    flux-density at a given reference frequency (data to be multiplied!)
    This routine is needed by hpXcalColorCorr.py
    
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
    numInterp=LinearInterpolator(freq,transm)
    denomInterp=LinearInterpolator(freq,transm * fSky * area)
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
#===============================================================================
#=====                    CALCULATE EFFECTIVE BEAM AREAS                   =====
#===============================================================================
#-------------------------------------------------------------------------------

def calcOmegaEff(alphaK,verbose=False,table=False):
    # calculate effective beam area for power law spectrum

    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #freq=getSpireFreq()
    #spireFiltOnly=getSpireFilt(cal,rsrfOnly=True)

    try:
        na=len(alphaK)
        aList=True
    except:
        aList=False

    if not aList:
        # alphaK is a scalar
        beamArea = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
        for band in spireBands:
            #pipeline beam areas
            beamArea[band]=spireEffArea(getSpireFreq(), getSpireFilt(rsrfOnly=True)[band], \
              calcBeamMonoArea()[band], BB=False, alpha=alphaK)/arcsec2Sr
        if (verbose): print 'Calculated Omega_eff for alpha=%f: '%alphaK,beamArea
    else:
        # tempK is a list
        beamArea = {'PSW': Double1d(na,Double.NaN),'PMW': Double1d(na,Double.NaN), 'PLW': Double1d(na,Double.NaN)}
        for a in range(na):
            for band in spireBands:
                beamArea[band][a]=spireEffArea(getSpireFreq(), getSpireFilt(rsrfOnly=True)[band],\
                  calcBeamMonoArea()[band], BB=False, alpha=alphaK[a])/arcsec2Sr
            if (verbose): print 'Calculated Omega_eff for alpha=%f: '%alphaK[a],beamArea["PSW"][a],beamArea["PMW"][a],beamArea["PLW"][a]

    if not table:
        #return as is
        return beamArea
    else:
        #create and returnTableDataset
        beamArea_table=TableDataset()
        beamArea_table.setDescription("Beam Solid Angle (Spectral Index)")
        beamArea_table.addColumn("alpha",Column(Double1d(alphaK)))
        for band in spireBands:
            beamArea_table.addColumn(band,Column(beamArea[band],unit=SolidAngle.STERADIANS,description=''))
        return(beamArea_table)


def calcOmegaEff_BB(betaK,tempK,verbose=False,table=False):
    # calculate effective beam area for modified BB spectrum
    # allows multiple temperatures, but only one beta value

    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #freq=getSpireFreq()
    #spireFiltOnly=getSpireFilt(cal,rsrfOnly=True)
    
    try:
        nt=len(tempK)
        tList=True
    except:
        tList=False

    if not tList:
        # tempK is scalars
        beamAreaBB = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
        for band in spireBands:
            #pipeline beam areas
            beamAreaBB[band]=spireEffArea(getSpireFreq(), getSpireFilt(rsrfOnly=True)[band], \
              calcBeamMonoArea()[band], BB=True, beta=betaK, temp=tempK)/arcsec2Sr
        if (verbose): print 'Calculated Omega_eff for modBB with beta=%f and T=%f: '%(betaK,tempK),beamAreaBB

    else:
        # tempK is a list
        beamAreaBB = {'PSW': Double1d(nt,Double.NaN), 'PMW': Double1d(nt,Double.NaN), 'PLW': Double1d(nt,Double.NaN)}
        for t in range(nt):
            for band in spireBands:
                #pipeline beam areas
                beamAreaBB[band][t]=spireEffArea(getSpireFreq(), getSpireFilt(rsrfOnly=True)[band], \
                  calcBeamMonoArea()[band], BB=True, beta=betaK, temp=tempK[t])/arcsec2Sr
            if (verbose): print 'Calculated Omega_eff for modBB with beta=%f and T=%f: '%(betaK,tempK[t]),beamAreaBB["PSW"][t],beamAreaBB["PMW"][t],beamAreaBB["PLW"][t]

    if not table:
        #return as is
        return beamAreaBB
    else:
        #create and returnTableDataset
        beamAreaBB_table=TableDataset()
        beamAreaBB_table.setDescription("Beam Solid Angle (Modified Black Body, beta=%.2f)"%betaK)
        beamAreaBB_table.meta['beta']=DoubleParameter(betaK,"Emissivity spectral index")
        beamAreaBB_table.addColumn("Temperature",Column(Double1d(tempK),unit=Temperature.KELVIN,description=''))
        for band in spireBands:
            beamAreaBB_table.addColumn(band,Column(beamAreaBB[band],unit=SolidAngle.STERADIANS,description=''))
        return(beamAreaBB_table)
       
def calcKBeam(alphaK,verbose=False,table=False):
    # Calculate pipeline colour correction parameters
    # cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    try:
        na=len(alphaK)
        aList=True
    except:
        aList=False

    beamAreaPip = calcOmegaEff(-1.0)

    if not aList:
        # alphaK is a scalar
        kBeam = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
        beamEff = calcOmegaEff(alphaK)
        for band in spireBands:
            #pipeline beam areas
            kBeam[band] = beamAreaPip[band]/beamEff[band]
        if verbose: print 'Calculated KBeam for alpha=%f: '%alphaK,kBeam
    else:
        # alphaK is a list
        kBeam = {'PSW': Double1d(na,Double.NaN), 'PMW': Double1d(na,Double.NaN), 'PLW': Double1d(na,Double.NaN)}
        beamEff = calcOmegaEff(alphaK)
        for a in range(na):
            for band in spireBands:
                #pipeline beam areas
                kBeam[band][a] = beamAreaPip[band]/beamEff[band][a]
            if verbose: print 'Calculated KBeam for alpha=%f: '%alphaK[a],kBeam["PSW"][a],kBeam["PMW"][a],kBeam["PLW"][a]

    if not table:
        #return as is
        return kBeam
    else:
        #create and returnTableDataset
        kBeam_table=TableDataset()
        kBeam_table.setDescription("Beam Colour Correction (Spectral Index)")
        kBeam_table.addColumn("alpha",Column(Double1d(alphaK)))
        for band in spireBands:
            kBeam_table.addColumn(band,Column(kBeam[band]))
        return kBeam_table
#
def calcKBeam_BB(betaK,tempK,verbose=False,table=False):
    # Calculate pipeline colour correction parameters
    # allows multiple temperatures, but only one beta value
    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    try:
        nt=len(tempK)
        tList=True
    except:
        tList=False

    beamAreaPip = calcOmegaEff(-1.0)
    if not tList:
        # tempK is scalar
        kBeamBB = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
        beamEff = calcOmegaEff_BB(betaK,tempK)
        for band in spireBands:
            #pipeline beam areas
            kBeamBB[band] = beamAreaPip[band]/beamEff[band]
        if (verbose): print 'Calculated KBeam for modBB with beta=%f and T=%f: '%(betaK,tempK),kBeamBB
    else:
        # tempK is a list
        kBeamBB = {'PSW': Double1d(nt,Double.NaN), 'PMW': Double1d(nt,Double.NaN), 'PLW': Double1d(nt,Double.NaN)}
        beamEff = calcOmegaEff_BB(betaK,tempK)
        for t in range(nt):
            for band in spireBands:
                #pipeline beam areas
                kBeamBB[band][t] = beamAreaPip[band]/beamEff[band][t]
            if (verbose): print 'Calculated KBeam for modBB with beta=%f and T=%f: '%(betaK,tempK[t]),kBeamBB["PSW"][t],kBeamBB["PMW"][t],kBeamBB["PLW"][t]

    if not table:
        #return as is
        return kBeamBB
    else:
        #create and returnTableDataset
        kBeamBB_table=TableDataset()
        kBeamBB_table.setDescription("Beam Colour Correction (Modified Black Body, beta=%.2f)"%betaK)
        kBeamBB_table.meta['beta']=DoubleParameter(betaK,"Emissivity spectral index")
        kBeamBB_table.addColumn("Temperature",Column(Double1d(tempK),unit=Temperature.KELVIN,description=''))
        for band in spireBands:
            kBeamBB_table.addColumn(band,Column(kBeamBB[band]))
        return(kBeamBB_table)

#-----------------------------------------------------------------------
#=======================================================================
#=====                 CALCULATE PIPELINE PARAMETERS               =====
#=======================================================================
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Calculate pipeline colour correction parameters
def calcK4P():
    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #freq=getSpireFreq()
    #spireFilt=getSpireFilt(cal)

    k4P = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
    for band in spireBands:
        k4P[band] = calcSpireKcorr(getSpireRefFreq()[band], getSpireFreq(), getSpireFilt()[band], \
        BB=False, ext=False)[0]
        pass
    return k4P
    
def calcKMonE():
    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #freq=getSpireFreq()
    #spireFilt=getSpireFilt(cal=cal)

    kMonE = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
    for band in spireBands:
        kMonE[band] = calcSpireKcorr(getSpireRefFreq()[band], getSpireFreq(), getSpireFilt()[band], \
        BB=False, alpha=-1.0, ext=True, monoArea=calcBeamMonoArea()[band])[0]/1.0e6
        pass
    return kMonE

def calcK4E():

    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #spireBands=["PSW","PMW","PLW"]
    
    k4E = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
    kMonE = calcKMonE()
    omegaEff = calcOmegaEff(-1.0)
    for band in spireBands:
        k4E[band] = kMonE[band] * omegaEff[band] * arcsec2Sr * 1.0e6
        pass
    return k4E

def calcKPtoE():
    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #spireBands=["PSW","PMW","PLW"]
    
    kPtoE = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
    kmonE = calcKMonE()
    k4p = calcK4P()
    for band in spireBands:
        kPtoE[band] = kmonE[band]/k4p[band]
    return kPtoE

#-----------------------------------------------------------------------
#=======================================================================
#=====           CALCULATE POINT SOURCE COLOR CORRECTIONS          =====
#=======================================================================
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
def calcKColP(alphaK,verbose=False,table=False):
    #print '\nCalculating point source colour correction parameters for a given alpha...'

    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #freq=getSpireFreq()
    #spireFilt=getSpireFilt(cal=cal)

    k4P=calcK4P()

    try:
        na=len(alphaK)
        aList=True
    except:
        aList=False
    
    if not aList:
        # alphaK is scalar
        kColP = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
        for band in spireBands:
            kConvPsrc=calcSpireKcorr(getSpireRefFreq()[band],\
               getSpireFreq(), getSpireFilt()[band], BB=False, alpha=alphaK)[0]
            #point source colour correction for current alpha
            kColP[band] = kConvPsrc/k4P[band]
        if (verbose): print 'Calculated KColP for alpha=%f: '%alphaK,kColP
    else:
        # alphaK is list
        kColP = {'PSW': Double1d(na,Double.NaN), 'PMW': Double1d(na,Double.NaN), 'PLW': Double1d(na,Double.NaN)}

        for a in range(na):
            for band in spireBands:
                kConvPsrc=calcSpireKcorr(getSpireRefFreq()[band],\
                   getSpireFreq(), getSpireFilt()[band], BB=False, alpha=alphaK[a])[0]
                #point source colour correction for current alpha
                kColP[band][a] = kConvPsrc/k4P[band]
            if (verbose): print 'Calculated KColP for alpha=%f: '%alphaK[a],kColP["PSW"][a],kColP["PMW"][a],kColP["PLW"][a]
            
    if not table:
        #return as is
        return kColP
    else:
        #create and returnTableDataset
        kColP_table=TableDataset()
        kColP_table.setDescription("Point Source Colour Correction (Spectral Index)")
        kColP_table.addColumn("alpha",Column(Double1d(alphaK)))
        for band in spireBands:
            kColP_table.addColumn(band,Column(kColP[band]))
        return kColP_table


def calcKColP_BB(betaK,tempK,verbose=False,table=False):
    #-----------------------------------------------------------------------
    #print 'Calculating point source colour correction parameters over beta & temp...'
    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #freq=getSpireFreq()
    #spireFilt=getSpireFilt(cal=cal)

    k4P=calcK4P()

    try:
        nt=len(tempK)
        tList=True
    except:
        tList=False
    
    if not tList:
        kColPBB = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
        for band in spireBands:
            kConvPsrc=calcSpireKcorr(getSpireRefFreq()[band],\
                getSpireFreq(), getSpireFilt()[band], BB=True, beta=betaK, temp=tempK)[0]
            #point source colour correction for current beta,temp
            kColPBB[band] = kConvPsrc/k4P[band]
        if (verbose): print 'Calculated KColP for modBB with T=%f K, beta=%f: '%(tempK,betaK),kColPBB
    else:
        kColPBB = {'PSW': Double1d(nt,Double.NaN), 'PMW': Double1d(nt,Double.NaN), 'PLW': Double1d(nt,Double.NaN)}
        for t in range(nt):
            for band in spireBands:
                kConvPsrc=calcSpireKcorr(getSpireRefFreq()[band],\
                    getSpireFreq(), getSpireFilt()[band], BB=True, beta=betaK, temp=tempK[t])[0]
                #point source colour correction for current beta,temp
                kColPBB[band][t] = kConvPsrc/k4P[band]
            if (verbose): print 'Calculated KColP for modBB with T=%f K, beta=%f: '%(tempK[t],betaK),kColPBB["PSW"][t],kColPBB["PMW"][t],kColPBB["PLW"][t]
    
    if not table:
        #return as is
        return kColPBB
    else:
        #create and returnTableDataset
        kColPBB_table=TableDataset()
        kColPBB_table.setDescription("Point Source Colour Correction (Modified Black Body, beta=%.2f)"%betaK)
        kColPBB_table.meta['beta']=DoubleParameter(betaK,"Emissivity spectral index")
        kColPBB_table.addColumn("Temperature",Column(Double1d(tempK),unit=Temperature.KELVIN,description=''))
        for band in spireBands:
            kColPBB_table.addColumn(band,Column(kColPBB[band]))
        return(kColPBB_table)
#-----------------------------------------------------------------------
#=======================================================================
#=====         CALCULATE EXTENDED SOURCE COLOR CORRECTIONS         =====
#=======================================================================
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

def calcKColE(alphaK,verbose=False,table=False):
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
        kColE = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}
        #kBeamK = calcKBeam(alphaK)
        for band in spireBands:
            k4EaTot_x=calcSpireKcorr(getSpireRefFreq()[band], getSpireFreq(),\
             getSpireFilt()[band], BB=False, alpha=alphaK,\
             ext=True, monoArea=calcBeamMonoArea()[band])[0]/1.e6
            kColE[band] = k4EaTot_x / k4E_Tot[band]
        if (verbose): print 'Calculated KColE for alpha=%f: '%alphaK,kColE
    else:
        # alphaK is list
        kColE = {'PSW': Double1d(na,Double.NaN), 'PMW': Double1d(na,Double.NaN), 'PLW': Double1d(na,Double.NaN)}
        for a in range(na):
            for band in spireBands:
                k4EaTot_x=calcSpireKcorr(getSpireRefFreq()[band], getSpireFreq(),\
                 getSpireFilt()[band], BB=False, alpha=alphaK[a],\
                 ext=True, monoArea=calcBeamMonoArea()[band])[0]/1.e6
                kColE[band][a] = k4EaTot_x / k4E_Tot[band]
            if (verbose): print 'Calculated KColE for alpha=%f: '%alphaK[a],kColE["PSW"][a],kColE["PMW"][a],kColE["PLW"][a]

    if not table:
        #return as is
        return kColE
    else:
        #create and returnTableDataset
        kColE_table=TableDataset()
        kColE_table.setDescription("Extended Source Colour Correction (Spectral Index)")
        kColE_table.addColumn("alpha",Column(Double1d(alphaK)))
        for band in spireBands:
            kColE_table.addColumn(band,Column(kColE[band]))
        return kColE_table
#-----------------------------------------------------------------------

def calcKColE_BB(betaK,tempK,verbose=False,table=False):
    #
    #print 'Calculating extended source colour correction parameters over beta & temp...'
    #cal=getCal(cal=cal,calPool=calPool,calFile=calFile)
    #freq=getSpireFreq()
    #spireFilt=getSpireFilt(cal=cal)

    k4E_Tot=calcKMonE()

    try:
        nt=len(tempK)
        tList=True
    except:
        tList=False

    if not tList:
        kColEBB = {'PSW': Double.NaN, 'PMW': Double.NaN, 'PLW': Double.NaN}

        for band in spireBands:
            k4EbTot_x=calcSpireKcorr(getSpireRefFreq()[band], getSpireFreq(),\
                  getSpireFilt()[band], BB=True, beta=betaK, temp=tempK,\
                  ext=True, monoArea=calcBeamMonoArea()[band])[0]/1.e6
            kColEBB[band] = k4EbTot_x / k4E_Tot[band]
        if (verbose): print 'Calculated KColE for modBB with T=%f K, beta=%f: '%(tempK,betaK),kColEBB
    else:
        kColEBB = {'PSW': Double1d(nt,Double.NaN), 'PMW': Double1d(nt,Double.NaN), 'PLW': Double1d(nt,Double.NaN)}
        for t in range(nt):
            for band in spireBands:
                k4EbTot_x=calcSpireKcorr(getSpireRefFreq()[band], getSpireFreq(),\
                      getSpireFilt()[band], BB=True, beta=betaK, temp=tempK[t],\
                      ext=True, monoArea=calcBeamMonoArea()[band])[0]/1.e6
                kColEBB[band][t] = k4EbTot_x / k4E_Tot[band]
            if (verbose): print 'Calculated KColE for modBB with T=%f K, beta=%f: '%(tempK[t],betaK),kColEBB["PSW"][t],kColEBB["PMW"][t],kColEBB["PLW"][t]

    if not table:
        #return as is
        return kColEBB
    else:
        #create and returnTableDataset
        kColEBB_table=TableDataset()
        kColEBB_table.setDescription("Extended Source Colour Correction (Modified Black Body, beta=%.2f)"%betaK)
        kColEBB_table.meta['beta']=DoubleParameter(betaK,"Emissivity spectral index")
        kColEBB_table.addColumn("Temperature",Column(Double1d(tempK),unit=Temperature.KELVIN,description=''))
        for band in spireBands:
            kColEBB_table.addColumn(band,Column(kColEBB[band]))
        return(kColEBB_table)
#
# test for some parameters
#
def spireColCorrTest(beamType='Full'):
    
    calphot=getCal(calPool='spire_cal_12_2')
    alphaArr=Float1d(range(-4,5)) #range of alphas to use.
    print beamType
    beamMonoArea=calcBeamMonoArea(beamType=beamType,verbose=True)

    print 'Testing pipeline parameters'
    print 'K4P=',calcK4P()
    print 'KMonE=',calcKMonE()
    print 'K4E=',calcK4E()
    print 'KPtoE=',calcKPtoE()
    print 'Omega(-1)=',calcOmegaEff(-1.0)
    omegaNepPsw=calcOmegaEff(spireCalPhot.getProduct('RadialCorrBeam').meta['alphaNeptunePsw'].double)['PSW']
    omegaNepPmw=calcOmegaEff(spireCalPhot.getProduct('RadialCorrBeam').meta['alphaNeptunePmw'].double)['PMW']
    omegaNepPlw=calcOmegaEff(spireCalPhot.getProduct('RadialCorrBeam').meta['alphaNeptunePlw'].double)['PLW']
    print 'Omega_Nep=',[omegaNepPsw,omegaNepPmw,omegaNepPlw]
    print '\nTesting calcOmegaEff'
    print 'OmegaEff(-2)=',calcOmegaEff(-2.0)
    omegaEff = calcOmegaEff(alphaArr)
    print 'OmegaEff_BB(1.75,20)=',calcOmegaEff_BB(1.75,20.0)
    omegaEff_BB = calcOmegaEff_BB(1.75,[20.0,30.,40.])

    print '\nTesting calcKBeam'
    print 'KBeam(-2)=',calcKBeam(-2.0)
    KBeam = calcKBeam(alphaArr)
    print 'KBeam(1.75,20)=',calcKBeam_BB(1.75,20.0)
    KBeam_BB = calcKBeam_BB(1.75,[20.0,30.,40.])

    print '\nTesting calcKColP'
    print 'KColP(-2)=',calcKColP(-2.0)
    KColP = calcKColP(alphaArr)
    print 'KColP(1.75,20)=',calcKColP_BB(1.75,20.0)
    KColP_BB = calcKColP_BB(1.75,[20.0,30.,40.])

    print '\nTesting calcKColE'
    print 'KColE(-2)=',calcKColE(-2.0)
    KColE = calcKColE(alphaArr)
    print 'KColE(1.75,20)=',calcKColE_BB(1.75,20.0)
    KColE_BB = calcKColE_BB(1.75,[20.0,30.,40.])
    
    return(beamMonoArea,alphaArr,omegaEff,KBeam,KColP,KColE)

def compFullSimple():
    #compare full vs single beam treatments:
    beamMonoArea_F,alphaArr,omegaEff_F,KBeam_F,KColP_F,KColE_F=spireColCorrTest(beamType='Full')
    beamMonoArea_S,alphaArr,omegaEff_S,KBeam_S,KColP_S,KColE_S=spireColCorrTest(beamType='Simple')
    pB=PlotXY()
    pOE=PlotXY()
    pKB=PlotXY()
    pKCE=PlotXY()

    pBd=PlotXY()
    pOEd=PlotXY()
    pKBd=PlotXY()
    pKCEd=PlotXY()

    cols={'PSW':java.awt.Color.BLUE,\
      'PMW':java.awt.Color.GREEN,\
      'PLW':java.awt.Color.RED}
    for band in spireBands:
        pB.addLayer(LayerXY(spireFreq/1.e9,beamMonoArea_F[band],name='Beam Area %s (full)'%band,\
          color=cols[band]))
        pB.addLayer(LayerXY(spireFreq/1.e9,beamMonoArea_S[band],name='Beam Area %s (simple)'%band,\
          color=cols[band],line=Style.DASHED))
        pBd.addLayer(LayerXY(spireFreq/1.e9,beamMonoArea_S[band]/beamMonoArea_F[band]-1.,name='Beam Area %s (s-f)'%band,\
          color=cols[band]))
        
        pOE.addLayer(LayerXY(alphaArr,omegaEff_F[band],name='OmegaEff %s (full)'%band,\
          color=cols[band]))
        pOE.addLayer(LayerXY(alphaArr,omegaEff_S[band],name='OmegaEff %s (simple)'%band,\
          color=cols[band],line=Style.DASHED))
        pOEd.addLayer(LayerXY(alphaArr,omegaEff_S[band]/omegaEff_F[band]-1.,name='OmegaEff %s (s-f)'%band,\
          color=cols[band]))
    
        pKB.addLayer(LayerXY(alphaArr,KBeam_F[band],name='KBeam %s (full)'%band,\
          color=cols[band]))
        pKB.addLayer(LayerXY(alphaArr,KBeam_S[band],name='KBeam %s (simple)'%band,\
          color=cols[band],line=Style.DASHED))
        pKBd.addLayer(LayerXY(alphaArr,KBeam_S[band]/KBeam_F[band]-1.,name='KBeam %s (s-f)'%band,\
          color=cols[band]))
    
        pKCE.addLayer(LayerXY(alphaArr,KColE_F[band],name='KColE %s (full)'%band,\
          color=cols[band]))
        pKCE.addLayer(LayerXY(alphaArr,KColE_S[band],name='KColE %s (simple)'%band,\
          color=cols[band],line=Style.DASHED))
        pKCEd.addLayer(LayerXY(alphaArr,KColE_S[band]/KColE_F[band]-1.,name='KColE %s (s-g)'%band,\
          color=cols[band]))

    pB.xaxis.titleText = 'Frequency (GHz)'
    pB.yaxis.titleText = 'Beam Area (sr)'
    pBd.xaxis.titleText = 'Frequency (GHz)'
    pBd.yaxis.titleText = 'Beam Area rel. diff'
    
    pOE.xaxis.titleText = 'Spectra Index (alpha)'
    pOE.yaxis.titleText = 'Effective Beam Area (sr)'
    pOEd.xaxis.titleText = 'Spectra Index (alpha)'
    pOEd.yaxis.titleText = 'Effective Beam Area rel. diff'
    
    pKB.xaxis.titleText = 'Spectra Index (alpha)'
    pKB.yaxis.titleText = 'KBeam'
    pKBd.xaxis.titleText = 'Spectra Index (alpha)'
    pKBd.yaxis.titleText = 'KBeam rel. diff'

    pKCE.xaxis.titleText = 'Spectra Index (alpha)'
    pKCE.yaxis.titleText = 'KColE'
    pKCEd.xaxis.titleText = 'Spectra Index (alpha)'
    pKCEd.yaxis.titleText = 'KColE rel. diff'
