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
# Classes:
#   Source(type,params): a source type and parameters
#   SourceProfile(radArr,profile): a radial source profile
#   RftProfile(kArr,rftProfile): radial Fourier Transform of a source profile
#
# Other functions:
#   calcJ0(z): calculate 0th order Bessel function J0(z)
#   j0Arg(theta): calculate the argument used in integral of the J0 Bessel function
#   rftArg(r): calculate the argument of the RFT integral
#
# Source Class:
#   An object containing a parameterised source description
#   Initialisation:
#     Source(type,params): define source type based on source type and parameters
#   Contains:
#     type: source type (string): GAUSSIAN|POWERLAW|LINEAR|CONSTANT|LINCONST
#     params: the parameters of the source type
#     name: the name of the source type
#     key: unique ID based on type and params
#     availableTypes: dictionary of available types and default parameters
#   Methods:
#     check(): check type and params are correct, and use default params if necessary
#     setKey(): set unique key
#     calcProfile(radArr): produce SourceProfile object corresponding to radArr
#     initAvailableTypes(): initialise the available source types
#
# SourceProfile Class:
#   An object containing a source profile
#   Initialisation:
#     SourceProfile(): create empty object
#     SourceProfile(radArr,profile): create object containing radArr, profile
#   Contains:
#     radArr: radius vector
#     nRad: length of radius vector
#     maxRad: max value of radArr
#     profile: radial profile corresponding to radArr
#   Methods:
#     setRadArr(radArr): set radArr,nRad,maxRad
#     generate(Source,radArr): calculate profile from Source corresponding to radArr
#     rad2k: calculate k-space vector from radArr
#     calcRft: produce RftProfile object from radArr, profile
#
# RftProfile Class:
#   An object containing a radial Fourier transform of a source profile
#   Initialisation:
#     RftProfile(): create empty object
#     RftProfile(kArr,rftProfile): create object containing kArr, rftProfile
#   Contains:
#     kArr: radius vector
#     nK: length of kArr vector
#     maxK: max value of kArr
#     rftProfile: radial Fourier transform profile corresponding to kArr
#   Methods:
#     setKArr(kArr): set radArr,nRad,maxRad
#     generate(Source,kArr): calculate radial Fourier transform profile from Source and kArr
#     k2rad: calculate radius array from kArr
#     calcInvRft: produce SourceProfile object from kArr, rftProfile
#
#-------------------------------------------------------------------------------
#===============================================================================
#=====                      IMPORT HIPE & JAVA MODULES                     =====
#===============================================================================
#-------------------------------------------------------------------------------

import os
import herschel
from herschel.share.util import Configuration
from herschel.ia.numeric import Double1d,Float1d,Int1d,Double2d
from herschel.ia.dataset import TableDataset,Column
from herschel.ia.dataset import DoubleParameter,LongParameter,StringParameter,BooleanParameter
from herschel.ia.dataset.image import SimpleImage
from herschel.ia.dataset.image.wcs import Wcs
from herschel.ia.toolbox.image import TransposeTask,CropTask,CircleHistogramTask,ImageConvolutionTask
from herschel.ia.toolbox.util import AsciiTableWriterTask, SimpleFitsWriterTask
from herschel.ia.toolbox.util import AsciiTableReaderTask, SimpleFitsReaderTask
crop=CropTask()
transpose=TransposeTask()
circleHistogram=CircleHistogramTask()
imageConvolution=ImageConvolutionTask()
asciiTableWriter=AsciiTableWriterTask()
simpleFitsWriter=SimpleFitsWriterTask()
asciiTableReader=AsciiTableReaderTask()
simpleFitsReader=SimpleFitsReaderTask()
from herschel.ia.numeric.toolbox.basic import Floor,Min,Max,Exp,Cos,Abs,Sqrt,Log
FLOOR=herschel.ia.numeric.toolbox.basic.Floor.PROCEDURE
EXP=herschel.ia.numeric.toolbox.basic.Exp.PROCEDURE
COS=herschel.ia.numeric.toolbox.basic.Cos.PROCEDURE
MAX=herschel.ia.numeric.toolbox.basic.Max.FOLDR
MIN=herschel.ia.numeric.toolbox.basic.Min.FOLDR
ABS=herschel.ia.numeric.toolbox.basic.Abs.FUNCTION
SQRT=herschel.ia.numeric.toolbox.basic.Sqrt.PROCEDURE
LOG=herschel.ia.numeric.toolbox.basic.Log.PROCEDURE
from java.lang import Double,Float,String,Integer
from java.lang.Math import PI
from herschel.ia.numeric.toolbox import RealFunction
from herschel.ia.numeric.toolbox.integr import TrapezoidalIntegrator
from herschel.ia.numeric.toolbox.interp import CubicSplineInterpolator,LinearInterpolator

#-------------------------------------------------------------------------------
#===============================================================================
#=====                          DEFINE Source CLASS                        =====
#===============================================================================
#-------------------------------------------------------------------------------
class Source(object):
    
    def __init__(self,srcTypeIn,paramsIn=None,verbose=False):
        #set available types
        self.initAvailableTypes()
        #set variables
        self.type=srcTypeIn
        self.params=paramsIn or []
        #check against available tmplates
        self.check(verbose=verbose)
        self.setKey()

    def setKey(self):
        #set unique key for source based on parameters
        key=self.type
        for par in self.params:
            key='%s_%g'%(key,par)
        self.key=key
        return(key)
    
    def check(self,verbose=False):
        #capitalise type and truncate to 8 letters
        typeNew=self.type[:8].upper()
        if typeNew != self.type:
            if (verbose):print 'Converting to upper case and truncating: %s -> %s'%(self.type,typeNew)
            self.type=typeNew


        #make sure parameters are in Double1d
        try:
            #see if len() works (doesn't work for scalars)
            len(self.params)
            #turn into Double1d
            self.params=Double1d(self.params)
        except:
            #convert scalar to 1-element Float1d
            self.params=Double1d([self.params])

        #set the typeTemplate (checks type and nParams are valid)
        self.setTemplate(verbose=verbose)
        
        self.nParams=len(self.params)
        
    def setTemplate(self,verbose=False):
        #set typeTemplate based on available SourceTypes
        if not self.availableTypes:
            #initialise available types
            self.initAvailableTypes()
        srcFound=False
        #loop over available types
        for typename in self.availableTypes:
            typex=self.availableTypes[typename]
            if self.type==typex.type:
                srcFound=True
                self.typeTemplate=self.availableTypes[typename]
        #error if source not found
        assert srcFound,'Unknown source type: %s'%self.type
        
        #set typeDesc
        self.typeDesc=self.typeTemplate.desc
        #check number of parameters is valid
        self.checkNParams()

        #apply default parameters from Template
        self.applyTemplateDefaults(verbose=verbose)

    def checkNParams(self):
        #check number of parameters is valud
        nParam=len(self.params)
        validNParam= (nParam >= self.typeTemplate.minPar and nParam <= self.typeTemplate.maxPar)

        assert (nParam >= self.typeTemplate.minPar),\
          '%s source requires at least %d parameters'%(self.typeTemplate.name,self.typeTemplate.minPar)
        assert (nParam <= self.typeTemplate.maxPar),\
          '%s source has maximum %d parameters'%(self.typeTemplate.name,self.typeTemplate.maxPar)
        return(validNParam)

    def applyTemplateDefaults(self,verbose=False):
        #apply default parameters from typeTemplate
        
        #get current number of parameters
        nParamIn=len(self.params)
        
        #get total number of parameters
        nParamOut=self.typeTemplate.maxPar
        paramsOut=Double1d(nParamOut)
        #loop over parameters
        for p in range(nParamOut):
            #if parameter is given
            if p < nParamIn:
                #parameter given
                if Double.isNaN(self.params[p]):
                    #if NaN use default value
                    paramsOut[p] = self.typeTemplate.defaults[p]
                    if (verbose):print 'Using default value for "%s": %f'%\
                        (self.typeTemplate.paramNames[p],paramsOut[p])
                else:
                    #use provided value
                    paramsOut[p] = self.params[p]
            else:
                #use default value
                paramsOut[p] = self.typeTemplate.defaults[p]
                if (verbose):print 'Using default value for "%s": %f'%\
                    (self.typeTemplate.paramNames[p],paramsOut[p])
        #set new params
        self.params=paramsOut
        #set paramNames
        self.paramNames=self.typeTemplate.paramNames
        #update nParams
        self.nParams=len(self.params)
        #set zeroValue
        self.zeroVal=self.typeTemplate.calcZeroVal(self.params)
        #set scaleWidth
        self.scaleWidth=self.typeTemplate.calcScaleWidth(self.params)

    def calcProfile(self,radArr):
        #run calcProfile associated with typeTemplate
        srcProf=self.typeTemplate.calcProfile(radArr,self.params)
        newProf=SourceProfile(radArr,srcProf,None,self,key='%s_Profile'%self.key)
        return(newProf)

    def calcArea(self):
        #run calcArea associated with typeTemplate
        srcArea=self.typeTemplate.calcArea(self.params)
        return(srcArea)
        
    def calcFwhm(self):
        #run calcArea associated with typeTemplate
        srcFwhm=self.typeTemplate.calcFwhm(self.params)
        return(srcFwhm)
    
    def initAvailableTypes(self):
        self.availableTypes={}
        self.availableTypes["GAUSSIAN"]=GaussianSourceType()
        self.availableTypes["EXPONENT"]=ExponentialSourceType()
        self.availableTypes['LINEAR']=LinearSourceType()
        self.availableTypes['POWERLAW']=PowerLawSourceType()
        self.availableTypes['CONSTANT']=ConstantSourceType()
        self.availableTypes['TOPHAT']=TopHatSourceType()
        self.availableTypes['LINCONST']=LinConstSourceType()

    def listSrcTypes(self):
        #list all source types
        print '====\nSource Types\n===='
        for typeName in self.availableTypes:
            print 'SourceType %s:'%(typeName)
            print str(self.availableTypes[typeName])
            print '-----'

    def toString(self):
        return(self.__str__())
        
    def __str__(self):
        #print helpful string
        result='%s'%self.__class__
        result=result+'\nKey=%s'%(self.key)
        result=result+"\nSourceType: %s (%s)"%(self.type,self.typeDesc)
        result=result+"\nnParams: %d"%self.nParams
        for p in range(self.nParams):
            result=result+"\n  params[%d] (%s): %g"%(p,self.paramNames[p],self.params[p])
        result=result+'\nscaleWidth=%g'%self.scaleWidth
        result=result+'\nzeroVal=%g'%self.zeroVal
        return(result)
        
#-------------------------------------------------------------------------------
#===============================================================================
#=====                        DEFINE SourceType CLASSES                    =====
#===============================================================================
#-------------------------------------------------------------------------------
class SourceType(object):
    def __init__(self,name=type,minPar=0,paramNames=None,defaults=None):
        self.desc=type
        self.minPar=0
        self.paramNames=[]
        self.defaults=[]
        self.maxPar=0
        if paramNames:
            for p in range(len(paramNames)):
                #first minPar are required, rest are not
                if p<minPar:
                    required=True
                else:
                    required=False
                self.addParameter(paramNames[p],defaults[p],required=required)
            
    def addParameter(self,paramName,default,required=False):
        print paramName,default,required
        if required:
            if self.minPar<self.maxPar:
                print 'cannot make required parameter, as there are already non-required ones'
            else:
                self.minPar+=1
        self.paramNames.append(paramName)
        self.defaults.append(default)
        self.maxPar=len(self.paramNames)

    def calcArea(self,params=None):
        print 'Warning: no analytical value for area of %s source'%self.type
        return(Double.NaN)
        
    def calcFwhm(self,params=None):
        print 'Warning: no analytical value for FHWM of %s source'%self.type
        return(Double.NaN)
        
    def __str__(self):
        result=self.__class__
        result="\nType name: %s"%self.desc
        result=result+"\nParameters: %d (%d required)"%(self.maxPar,self.minPar)
        for p in range(self.maxPar):
            result=result+"\n  %d: %s (default %g)"%(p,self.paramNames[p],self.defaults[p])
        return(result)
        
    
class NoneSourceType(SourceType):
    def __init__(self):
        self.type='NONE'
        self.desc='NONE source type'
        self.minPar=0
        self.paramNames=[]
        self.defaults=[]
        self.maxPar=len(self.paramNames)
        
    def calcProfile(self,radArr,params):
        return(Double1d(len(radArr),Double.NaN))
        
    def calcArea(self,params=None):
        return(Double.NaN)
        
    def calcZeroVal(self,params):
        return(Double.NaN)
        
    def calcFwhm(self,params):
        return(Double.NaN)
    
    def calcScaleWidth(self,params):
        return(Double.NaN)
    
class GaussianSourceType(SourceType):
    def __init__(self):
        self.type='GAUSSIAN'
        self.desc='Gaussian profile'
        self.minPar=1
        self.paramNames=["Source width","Peak value"]
        self.defaults=[100.,1.0]
        self.maxPar=len(self.paramNames)
    
    def calcProfile(self,radArr,params):
        srcWidth=params[0]
        srcPeak=params[1]
        srcProf=srcPeak * EXP(-radArr**2/(2.*srcWidth**2))
        return(srcProf)
        
    def calcArea(self,params=None):
        #integrate Gaussian
        assert params!=None,'must provide parameters'
        srcWidth=params[0]
        srcPeak=params[1]
        srcArea=srcPeak * 2*PI*srcWidth**2
        return(srcArea)
        
    def calcZeroVal(self,params):
        zeroVal=params[1]
        return(zeroVal)
        
    def calcFwhm(self,params):
        srcWidth=params[0]
        srcFwhm=2. * srcWidth * SQRT(2.*LOG(2.))
        return(srcFwhm)
        
    def calcScaleWidth(self,params):
        scaleWidth=params[0]
        return(scaleWidth)

class ExponentialSourceType(SourceType):
    def __init__(self):
        self.type='EXPONENT'
        self.desc='Exponential profile'
        self.minPar=1
        self.paramNames=["Scale width","Peak value"]
        self.defaults=[100.,1.0]
        self.maxPar=len(self.paramNames)
    
    def calcProfile(self,radArr,params):
        srcWidth=params[0]
        srcPeak=params[1]
        srcProf=srcPeak * EXP(-radArr/srcWidth)
        return(srcProf)
        
    def calcArea(self,params=None):
        #integrate exponential
        assert params!=None,'must provide parameters'
        srcWidth=params[0]
        srcPeak=params[1]
        if srcWidth > 0:
            srcArea=srcPeak * 2.*PI*1.*srcWidth**2
        else:
            print 'Warning: no analytical value for area of EXPONENT source with negative scale width'
            srcArea=Double.NaN
        return(srcArea)
        
    def calcZeroVal(self,params):
        zeroVal=params[1]
        return(zeroVal)
        
    def calcFwhm(self,params):
        srcWidth=params[0]
        srcFwhm=2.*srcWidth * LOG(2.)
        return(srcFwhm)
        
    def calcScaleWidth(self,params):
        scaleWidth=params[0]
        return(scaleWidth)
        
class LinearSourceType(SourceType):
    def __init__(self):
        self.type='LINEAR'
        self.desc='Linear profile'
        self.minPar=1
        self.paramNames=["Gradiant","Value at r=0"]
        self.defaults=[-0.01,1.0]
        self.maxPar=len(self.paramNames)
        
    def calcProfile(self,radArr,params):
        print params
        srcGrad=params[0]
        srcZero=params[1]
        srcProf=srcZero + srcGrad * radArr
        return(srcProf)
        
    def calcZeroVal(self,params):
        zeroVal=params[1]
        return(zeroVal)
        
    def calcScaleWidth(self,params):
        #calculate half-width-half-max
        scaleWidth=ABS(params[1]/params[0])
        return(scaleWidth)
    
class PowerLawSourceType(SourceType):
    def __init__(self):
        self.type='POWERLAW'
        self.desc='Power law profile'
        self.minPar=2
        self.paramNames=["Spectral index","Scale Radius","min radius to extend to","value below min radius"]
        self.defaults=[-1.,100.,1.e-3,Double.NaN]
        self.maxPar=len(self.paramNames)
        
    def calcProfile(self,radArr,params):
        srcIdx=params[0]
        srcScalRad=params[1]
        srcMinRad=params[2]
        srcValMinRad=params[3]
        srcProf=Double1d(len(radArr))
        srcProf[radArr.where(radArr>=srcMinRad)]=(radArr[radArr.where(radArr>=srcMinRad)]/srcScalRad)**srcIdx
        srcProf[radArr.where(radArr<srcMinRad)]=srcValMinRad
        return(srcProf)
        
    def calcZeroVal(self,params):
        zeroVal=params[3]
        return(zeroVal)
        
    def calcFwhm(self,params):
        srcIdx=params[0]
        srcScalRad=params[1]
        srcMinRad=params[2]
        srcValMinRad=params[3]
        if srcValMinRad==Double.NaN:
            #calculate from minRad
            srcFwhm = 2.*srcMinRad / 2.**(1./srcIdx)
        else:
            srcFwhm = 2. * (srcValMinRad/2)**(1./srcIdx) * srcScalRad
        return(srcFwhm)
        
    def calcScaleWidth(self,params):
        #return scale radius
        scaleWidth=params[1]
        return(scaleWidth)	
        
class ConstantSourceType(SourceType):
    def __init__(self):
        self.type='CONSTANT'
        self.desc='Constant profile'
        self.minPar=0
        self.paramNames=["Value"]
        self.defaults=[1.]
        self.maxPar=len(self.paramNames)
        
    def calcProfile(self,radArr,params):
        srcVal=params[0]
        srcProf=Double1d(len(radArr),srcVal)
        return(srcProf)
        
    def calcZeroVal(self,params):
        zeroVal=params[0]
        return(zeroVal)
        
    def calcFwhm(self,params):
        return(Double.POSITIVE_INFINITY)
        
    def calcScaleWidth(self,params):
        #set as infinity
        scaleWidth=Double.POSITIVE_INFINITY
        return(scaleWidth)

class TopHatSourceType(SourceType):
    def __init__(self):
        self.type='TOPHAT'
        self.desc='Top Hat profile'
        self.minPar=1
        self.paramNames=["Width","Value"]
        self.defaults=[100.,1.]
        self.maxPar=len(self.paramNames)
        
    def calcProfile(self,radArr,params):
        srcWidth=params[0]
        srcVal=params[1]
        srcProf=Double1d(len(radArr))
        srcProf[radArr.where(radArr<=srcWidth)]=srcVal
        return(srcProf)
        
    def calcArea(self,params):
        srcWidth=params[0]
        srcVal=params[1]
        area = srcVal * PI * srcWidth**2
        return(area)
        
    def calcZeroVal(self,params):
        zeroVal=params[1]
        return(zeroVal)
        
    def calcFwhm(self,params):
        srcFwhm=2.*params[0]
        return(srcFwhm)
        
    def calcScaleWidth(self,params):
        #set as infinity
        scaleWidth=params[0]
        return(scaleWidth)

class LinConstSourceType(SourceType):
    def __init__(self):
        self.type='LINCONST'
        self.desc='Linear profile with limits'
        self.minPar=1
        self.paramNames=["Gradiant","Value at r=0","Min Limit","Max Limit"]
        self.defaults=[-0.01,1.0, 0.,Double.POSITIVE_INFINITY]
        self.maxPar=len(self.paramNames)
        
    def calcProfile(self,radArr,params):
        srcGrad=params[0]
        srcZero=params[1]
        srcMin=params[2]
        srcMax=params[3]
        srcProf=srcZero + srcGrad * radArr
        if MAX(srcProf) > srcMax:
            srcProf[srcProf.where(srcProf > srcMax)]=srcMax
        if MIN(srcProf) < srcMin:
            srcProf[srcProf.where(srcProf < srcMin)]=srcMin
        return(srcProf)
        
    def calcArea(self,params=None):
        assert params!=None,'must provide parameters'
            
        srcGrad=params[0]
        srcZero=params[1]
        srcMin=params[2]
        srcMax=params[3]
        if srcMin==0 and srcGrad<0:
            #profile -> 0 at high r
            srcArea=(PI/3)*srcZero**3/srcGrad**2
            if srcMax < srcZero:
                #account for central region being uniform
                srcArea=srcArea - (PI/3)*(srcZero-srcMax)**3/srcGrad**2
        else:
            #no analytical solution
            if srcMin!=0:
                print 'Warning: no analytical value for area of LINCONST source with min < 0'
            if srcGrad>=0:
                print 'Warning: no analytical value for area of LINCONST source with gradient >= 0'
            srcArea=Double.NaN
            
        return(srcArea)
        
    def calcFwhm(self,params):
        srcGrad=params[0]
        srcZero=params[1]
        srcMin=params[2]
        srcMax=params[3]
        if srcGrad<0:
            srcFwhm=-(srcZero-srcMin)/srcGrad
        else:
            #no analytical solution
            print 'Warning: no analytical value for FWHM of LINCONST source with gradient >= 0'
            srcFwhm=Double.NaN
        return(srcFwhm)
        
    def calcZeroVal(self,params):
        zeroVal=params[1]
        if zeroVal < params[2]: zeroVal=params[2]
        if zeroVal > params[3]: zeroVal=params[3]
        return(zeroVal)
        
    def calcScaleWidth(self,params):
        #calculate half-width-half-max
        if params[0] < 0:
            scaleWidth=-(params[1]-params[2])/params[0]
        else:
            scaleWidth=(params[3]-params[1])/params[0]
        return(scaleWidth)
    
#-------------------------------------------------------------------------------
#===============================================================================
#=====                       DEFINE SourceProfile CLASS                    =====
#===============================================================================
#-------------------------------------------------------------------------------
class SourceProfile(object):
    def __init__(self,radArr=None,profile=None,error=None,originator=None,key=None):
        #set radArr
        self.setRadArr(radArr)
        #set profile
        self.setProfile(profile)
        self.setError(error)
        #check radArr and profile are compatible
        self.check()
        #set original Source object
        self.setOrig(originator)
        self.setKey(key)

    def setRadArr(self,radArr=None):
        #set radArr, nRad and maxRad
        self.radArr=radArr or None
        if self.radArr:
            self.nRad=len(radArr)
            self.maxRad=MAX(radArr)
        else:
            self.nRad=0
            self.maxRad=Double.NaN
            
        return(self)

    def setProfile(self,profile=None):
        #set error
        if self.radArr==None:
            #no radArr
            assert profile==None,\
              'cannot set profile without radArr'
            self.profile=None
        else:
            #radArr exists
            if profile==None:
                #make blank array
                self.profile=Double1d(self.nRad)
            else:
                assert len(profile)==self.nRad,\
                  'profile array must be same length as rad array'
                self.profile=profile
        return(self)
        
    def setError(self,error=None):
        #set error
        if self.radArr==None:
            #no radArr
            assert error==None,\
              'cannot set error without radArr'
            self.error=None
        else:
            #radArr exists
            if error==None:
                #make black array
                self.error=Double1d(self.nRad)
            else:
                assert len(error)==self.nRad,\
                  'error array must be same length as rad array'
                self.error=error
        return(self)

    def setKey(self,key,fallbackKey=None):
        if key==None:
            newKey=fallbackKey
            #print 'setting fallback key %s'%fallbackKey
        else:
            newKey=key
            #print 'setting key %s'%key
        self.key=newKey
        return(self)

    def checkRadArr(self,radArr2):
        #check whether radius array matches another
        if self.nRad==len(radArr2) and self.maxRad==MAX(radArr2):
            #arrays have same max and length
            #take difference of arrays
            radDiff=self.radArr - radArr2
            if MIN(radDiff)==0. and MAX(radDiff)==0.:
                #rad arrays match
                radMatch=True
            else:
                #rad arrays don't exactly match
                radMatch=False
        else:
            #arrays have different length and/or max value
            radMatch=False
        return(radMatch)

    def check(self):
        #check radArr and profile are both/neither None
        #check both same length
        if self.radArr==None:
            assert (self.profile==None),\
              'must set both radArr and profile, or neither'
        else:
            assert (self.profile!=None),\
              'must set both radArr and profile, or neither'
            assert len(self.radArr)==len(self.profile),\
              'radArr and profile must be of same lengths'%(len(radArr),len(profile))
        return(self)

    def copy(self):
        return(SourceProfile(radArr=self.radArr,profile=self.profile,error=self.error,originator=self.origSrc,key=self.key))

    def makeTable(self):
        table=TableDataset()
        table.setDescription('Radial Profile Table')
        table.addColumn('radius',Column(self.radArr,description='radius array'))
        table.addColumn('profile',Column(self.profile,description='radial profile'))
        table.addColumn('error',Column(self.error,description='error on radial profile'))
        table.meta['nRad']=LongParameter(self.nRad,description='length of radius array')
        table.meta['maxRad']=DoubleParameter(self.maxRad,description='maximum radius')
        table.meta['key']=StringParameter(self.key,description='unique key')
        table.meta['fwhm']=DoubleParameter(self.calcFwhm(),description='source FWHM')
        table.meta['area']=DoubleParameter(self.calcArea(),description='source area')
        table.meta['fromSrc']=BooleanParameter(self.fromSrc,description='is profile based on Source')
        if self.fromSrc:
            table.meta['origSrc']=StringParameter(self.origSrc.key,description='original Source key')
            table.meta['srcType']=StringParameter(self.origSrc.type,description='source type')
            for n in range(self.origSrc.nParams):
                table.meta['srcPar_%d'%n]=DoubleParameter(Double(self.origSrc.params[n]),\
                  description=self.origSrc.paramNames[n])
        self.table=table
        return(table)
        
    def saveFits(self,directory=None,filename=None,verbose=False):
        #save table in FITS file
        
        #make table object
        table=self.makeTable()
        #set directory and filename
        if directory==None:
            directory=Configuration.getProperty('var.hcss.workdir')
        if self.key==None:
            assert filename!=None,'No key in SourceProfile. Filename must be provided'
        if filename==None:
            filename='%s.fits'%self.key
        simpleFitsWriter(table,os.path.join(directory,filename))
        if verbose:print 'Written to %s'%os.path.join(directory,filename)
        
    def saveAscii(self,directory=None,filename=None,verbose=False):
        #save table in ASCII CSV file
        
        #make table object
        table=self.makeTable()
        #set directory and filename
        if directory==None:
            directory=Configuration.getProperty('var.hcss.workdir')
        if self.key==None:
            assert filename!=None,'No key in SourceProfile. Filename must be provided'
        if filename==None:
            filename='%s.csv'%self.key
        asciiTableWriter(table,os.path.join(directory,filename))
        if verbose:print 'Written to %s'%os.path.join(directory,filename)

    def loadFits(self,directory=None,filename=None,verbose=False):
        #read from a FITS file
        if directory==None:
            directory=Configuration.getProperty('var.hcss.workdir')
        assert filename!=None,'Must provide filename'
        table=simpleFitsReader(os.path.join(directory,filename))
        self.setRadArr(table['radius'].data)
        self.setProfile(table['profile'].data)
        self.setError(table['error'].data)
        self.setKey(table.meta['key'].string)
        self.clearOrig()
        if verbose:
            print 'Reading profile from %s'%(os.path.join(directory,filename))
        return self
        
    def loadAscii(self,directory=None,filename=None,verbose=False):
        #read from a FITS file
        if directory==None:
            directory=Configuration.getProperty('var.hcss.workdir')
        assert filename!=None,'Must provide filename'
        table=asciiTableReader(os.path.join(directory,filename))
        self.setRadArr(table['radius'].data)
        self.setProfile(table['profile'].data)
        self.setError(table['error'].data)
        set.setKey(table.meta['key'].string)
        if verbose:
            print 'Reading profile from %s'%(os.path.join(directory,filename))
        return self
        self.clearOrig()

    def generate(self,src,radArr,key=None):
        #generate profile from Source object
        self.setRadArr(radArr)
        self.setProfile(src.calcProfile(radArr))
        #make blank error array
        self.setError(None)
        self.setOrig(src)
        self.setKey(key,'%s_Profile'%src.key)
        return(self)

    def setOrig(self,src):
        #set whether profile comes from Source object
        if src:
            self.fromSrc=True
            self.origSrc=src
        else:
            self.fromSrc=False
            self.origSrc=NoneSourceType()
            
    def clearOrig(self):
        #clear Source object originator
        self.fromSrc=False
        self.origSrc=NoneSourceType()

    def regrid(self,radNew,key=None):
        #interpolate profile to a new radial grid
        #create interpolation objects
        profInterp=CubicSplineInterpolator(Double1d(self.radArr),Double1d(self.profile))
        errInterp=CubicSplineInterpolator(Double1d(self.radArr),Double1d(self.error))
        
        #get valid range of radii
        minRad=min(self.radArr)
        maxRad=self.maxRad
        iR=radNew.where(radNew>=minRad and radNew <= maxRad)
        nRadNew=len(radNew)
        #make new profile and error arrays
        profNew=Double1d(nRadNew)
        errNew=Double1d(nRadNew)
        profNew[iR]=profInterp(radNew[iR])
        errNew[iR]=errInterp(radNew[iR])
        newProf=self.copy()
        newProf.setRadArr(radNew)
        newProf.setProfile(profNew)
        newProf.setError(errNew)
        if self.key==None:
            newProf.setKey(key)
        else:
            newProf.setKey(key,'%s_regrid'%self.key)
        
        return(newProf)
        
    def normProf(self,key=None):
        #normalise such that profile=1 at r=0
        prof0=self.profile[0]
        newProf=self.copy()
        newProf.profile=self.profile/prof0
        newProf.error=self.error/prof0
        newProf.clearOrig()
        if self.key==None:
            newProf.setKey(key)
        else:
            newProf.setKey(key,'%s_normProf'%(self.key))
        return(newProf)
        
    def normArea(self,key=None):
        #normalise such that area=1
        area0=self.calcArea(forceNumerical=True)
        #print 'area0=',area0
        newProf=self.copy()
        newProf.profile=self.profile/area0
        newProf.error=self.error/area0
        newProf.clearOrig()
        #print 'new area=',newProf.calcArea()
        if self.key==None:
            newProf.setKey(key)
        else:
            newProf.setKey(key,'%s_normArea'%(self.key))
        return(newProf)

    def calcArea(self,forceNumerical=False,forceAnalytical=False):
        #print 'calc area'
        if forceNumerical:
            #use the numerial integration method
            doNum=True
        elif forceAnalytical:
            #calculate analytically from the Source originator
            doNum=False
            if self.fromSrc:
                srcArea=self.origSrc.calcArea()
            else:
                srcArea=Double.NaN
        else:
            #calculate from source if possible, else use numerical integration
            if self.fromSrc:
                #try to get analytical area from origSrc
                srcArea=self.origSrc.calcArea()
                if Double.isNaN(srcArea):
                    #can't calculate analytically
                    doNum=True
                else:
                    doNum=False
            else:
                 doNum=True

        if doNum:
            #calculate area numerically
            #print 'numerical integration'
            profInterp=CubicSplineInterpolator(self.radArr,2.*PI*self.radArr*self.profile)
            integrator=TrapezoidalIntegrator(0,self.maxRad)
            srcArea=integrator.integrate(profInterp)
            
        return(srcArea)
        
    def calcFwhm(self,forceNumerical=False,forceAnalytical=False):
        #print 'calc area'
        if forceNumerical:
            doNum=True
        elif forceAnalytical:
            doNum=False
            if self.fromSrc:
                srcFwhm=self.origSrc.calcFwhm()
            else:
                srcFwhm=Double.NaN
        else:
            if self.fromSrc:
                #try to get analytical area from origSrc
                srcFwhm=self.origSrc.calcFwhm()
                if Double.isNaN(srcFwhm):
                    doNum=True
                else:
                    doNum=False
            else:
                 doNum=True

        if doNum:
            #calculate area numerically
            try:
                #try interpolation (only works if profile is monatonic
                profInterp=CubicSplineInterpolator(self.profile,self.radArr)
                profHalf=self.profile[0]/2.
                srcFwhm=2.*profInterp(profHalf)
            except:
                #step through manually
                halfFound=False
                profHalf=self.profile[0]/2.
                r=0
                profThis=self.profile[0]
                radThis=self.radArr[0]
                while not halfFound and r<self.nRad-1:
                    r=r+1
                    profPrev=profThis
                    radPrev=radThis
                    radThis=self.radArr[r]
                    profThis=self.profile[2]
                    if profThis<=profHalf:
                        halfFound=True
                if not halfFound:
                    srcFwhm=Double.NaN
                else:
                    srcFwhm = 2. * (radPrev + (profHalf-profPrev)*(radThis-radPrev)/(profThis-profPrev))

        return(srcFwhm)
        
    def add(self,prof2,err2=None,key=None):
        #add profile to another
        newProf=self.copy()
        if type(prof2)==SourceProfile:
            #interpolate new prof to same as existing profile
            prof2Regrid=prof2.regrid(self.radArr)
            #add profiles
            newProf.profile = self.profile + prof2Regrid.profile
            #add errors in quadrature
            newProf.error = SQRT(self.error**2 + prof2Regrid.error**2)
            #set key
            if self.key==None or prof2.key==None:
                newProf.setKey(key)
            else:
                newProf.setKey(key,'%s_ADD_%s'%(self.key,prof2.key))
        else:
            #check whether prof2 is scalar of array
            try:
                #len doesn't work on scalar
                len(prof2)
                isScal=False
            except:
                isScal=True
            
            if not isScal:
                #check vector lengths
                assert len(prof2)==len(self.profile),\
                  'Added profile must scalar of of same length as SourceProfile'
                  
                assert len(err2)==len(self.profile),\
                  'Added profile must scalar of of same length as SourceProfile'
            
            #add profiles
            newProf.profile = self.profile + prof2
            #set key
            newProf.setKey(key)
            
            if error!=None:
            #check whether err2 is scalar of array
                try:
                    #len doesn't work on scalar
                    len(prof2)
                    isScal=False
                except:
                    isScal=True
                
                if not isScal:
                    #check vector lengths
                    assert len(prof2)==len(self.profile),\
                      'Added profile must scalar of of same length as SourceProfile'
                      
                    assert len(err2)==len(self.profile),\
                      'Added profile must scalar of of same length as SourceProfile'
        #no longer based on SourceType
        newProf.clearOrig()
        
        return(newProf)
        
    def mult(self,prof2,key=None):
        #multiply profile by another
        newProf=self.copy()
        
        if type(prof2)==SourceProfile:
            #interpolate new prof to same as existing profile
            prof2Regrid=prof2.regrid(self.radArr)
            #add profiles
            newProf.profile = self.profile * prof2Regrid.profile
            #add relative errors in quadrature
            newProf.error = newProf.profile * SQRT((self.error/self.profile)**2 + (prof2Regrid.error/prof2Regrid.profile)**2)
            #set key
            if self.key==None or prof2.key==None:
                newProf.setKey(key)
            else:
                newProf.setKey(key,'%s_MULT_%s'%(self.key,prof2.key))
        else:
            #check whether prof2 is scalar of array
            try:
                #len doesn't work on scalar
                len(prof2)
                isScal=False
            except:
                isScal=True
            
            if not isScal:
                #check vector lengths
                assert len(prof2)==len(self.profile),\
                  'Multiplied profile must scalar of of same length as SourceProfile'
            
            #multiply profiles and errors
            newProf.profile = self.profile * prof2
            #set key
            newProf.setKey(key)
            #set error
            newProf.error = self.error * prof2
        #no longer based on SourceType
        newProf.clearOrig()
        
        return(newProf)
        
    def makeImage(self):
        #make image of source
        #calculate image size
        nXQuad=self.nRad
        nYQuad=nXQuad
        quadArr=Double2d(nXQuad,nYQuad)
        quadErr=Double2d(nXQuad,nYQuad)
        #print quadArr.dimensions
        #create profile interpolation
        profInterp=CubicSplineInterpolator(Double1d(self.radArr),Double1d(self.profile))
        errInterp=CubicSplineInterpolator(Double1d(self.radArr),Double1d(self.error))
        #fill first quadrant
        yList=range(nYQuad)
        #print MIN(yList),MAX(yList)
        for x in range(nXQuad):
            radList=SQRT(x**2. + Double1d(yList)**2.)
            #print x
            inBeam=radList.where(radList <= self.maxRad)
            #print inBeam
            quadArr[x,inBeam]=profInterp(radList[inBeam])
            quadErr[x,inBeam]=errInterp(radList[inBeam])
            
        #fill other quadrants
        nXIm=2*self.nRad - 1
        nYIm=nXIm
        cXIm=self.nRad
        cYIm=self.nRad
        #imageArr=Double2d(nXIm,nYIm)
        #imageErr=Double2d(nXIm,nYIm)
        image=SimpleImage()
        image.setImage(Double2d(nXIm,nYIm))
        image.setError(Double2d(nXIm,nYIm))
        #print 'image=',image
        #top-right quarter
        image['image'].data[cXIm-1:nXIm,cYIm-1:nYIm]=quadArr
        image['error'].data[cXIm-1:nXIm,cYIm-1:nYIm]=quadErr
        #print 'top-right done'
        imageQuad=SimpleImage()
        imageQuad.setImage(quadArr)
        imageQuad.setError(quadErr)
        #top-left-quarter
        transImage=transpose(imageQuad,TransposeTask.FLIP_VERTICAL)
        image['image'].data[0:cXIm,cYIm-1:nYIm]=transImage.getImage()
        image['error'].data[0:cXIm,cYIm-1:nYIm]=transImage.getError()
        #print 'bottom-right done'
        #bottom-right quarter
        transImage=transpose(imageQuad,TransposeTask.FLIP_HORIZONTAL)
        image['image'].data[cXIm-1:nXIm,0:cYIm]=transImage.getImage()
        image['error'].data[cXIm-1:nXIm,0:cYIm]=transImage.getError()
        #print 'top-left done'
        #bottom-left quarter
        transImage=transpose(imageQuad,TransposeTask.FLIP_ANTIDIAGONAL)
        image['image'].data[0:cXIm,0:cYIm]=transImage.getImage()
        image['error'].data[0:cXIm,0:cYIm]=transImage.getImage()
        #print 'bottom-left done'
        #print image['image'].data.dimensions
        #image=SimpleImage()
        #image.setImage(imageArr)
        
        return(image)
        
    
    def __str__(self):
        result='%s'%self.__class__
        result=result+'\nKey=%s'%self.key
        result=result+'\nnRad=%d'%self.nRad
        if self.nRad > 0:
            result=result+'\nmaxRad=%g'%self.maxRad
            if self.nRad<=5:
                result=result+'\nradArr: '+self.radArr
                result=result+'\nprofile: '+self.profile
                result=result+'\nerror: '+self.error
            else:
                result=result+'\nradArr: [%g, %g, ... %g, %g]'%(self.radArr[0],self.radArr[1],self.radArr[-2],self.radArr[-1])
                result=result+'\nprofile: [%g, %g, ... %g, %g]'%(self.profile[0],self.profile[1],self.profile[-2],self.profile[-1])
                result=result+'\error: [%g, %g, ... %g, %g]'%(self.error[0],self.error[1],self.error[-2],self.error[-1])
            result=result+'\nfromSrc: %r'%self.fromSrc
            if self.fromSrc:
                result=result+'\norigSrc: %s'%self.origSrc.key
        return(result)

def convolveProfiles(profile1,profile2,key=None,verbose=False):
    assert type(profile1)==SourceProfile,'profile1 must be SourceProfile object'
    assert type(profile2)==SourceProfile,'profile2 must be SourceProfile object'
    
    #calculate profile resolutions
    res1=profile1.maxRad/(profile1.nRad-1)
    res2=profile2.maxRad/(profile2.nRad-1)
    if res1<res2:
        if(verbose):print 'regridding profile2 from %g to %g arcsec steps'%(res2,res1)
        #different resolutions, so regrid profile2
        radArr2New=Double1d(int(profile2.maxRad/res1))*res1
        profile2Regrid=profile2.regrid(radArr2New)
        profile1Regrid=profile1.copy()
    elif res2<res1:
        #different resolutions, so regrid profile2
        if(verbose):print 'regridding profile1 from %g to %g arcsec steps'%(res1,res2)
        radArr1New=Double1d(int(profile1.maxRad/res2))*res2
        profile1Regrid=profile1.regrid(radArr1New)
        profile2Regrid=profile2.copy()
    else:
        #same resolution
        if(verbose):print 'no regridding necessary'
        profile1Regrid=profile1.copy()
        profile2Regrid=profile2.copy()

    profile2Regrid=profile2Regrid.normArea()
    #make images of profiles
    image1=profile1Regrid.makeImage()
    #print image1
    image2=profile2Regrid.makeImage()
    #print image2
    #convolve images
    convolvedImage=imageConvolution(image1,image2)
    #compute profile of convolved image)
    convolvedProfile=image2ProfCirc(convolvedImage)
    #set key
    if profile1.key==None or profile2.key==None:
        convolvedProfile.setKey(key)
    else:
        convolvedProfile.setKey(key,'%s_CONV_%s'%(profile1.key,profile2.key))
    return(convolvedProfile)

def image2ProfCirc(imageIn):
    #assumes an image is square, centred and circularly symmetric
    nxIn=int(imageIn['image'].meta['naxis1'].long)
    nyIn=int(imageIn['image'].meta['naxis2'].long)
    cxIn=nxIn/2
    cyIn=nyIn/2
    radArr=Double1d(range(nxIn-cxIn))
    profile=imageIn.getImage()[cxIn:nyIn,cyIn]
    if imageIn.hasError():
        error=imageIn.getError()[cxIn:nyIn,cyIn]
    else:
        error=None
    return(SourceProfile(radArr,profile,error))

def map2Prof(mapIn,raCtr=None,decCtr=None,xCtr=None,yCtr=None,maxRadArcsec=700.,dRadArcsec=None,key=None,verbose=False):
    #general script to take a position in a map and create a radial source profile
    assert type(mapIn)==SimpleImage,'mapIn must be a SimpleImage'
    if raCtr!=None and decCtr!=None:
        useWcs=True
        ctrProvided=True
        if xCtr!=None and yCtr!=None:
            #both types of coords set
            print 'Both (raCtr,decCtr) and (xCtr,yCtr) set. Using (raCtr,decCtr). (xCtr,yCtr) might by overwritten.'
    elif xCtr!=None and yCtr!=None:
        useWcs=False
        ctrProvided=True
    else:
        if(verbose):print 'no (raCtr,decCtr) or (xCtr,yCtr) provided. Using centre pixel'
        ctrProvided=False

    #extract wcs info
    try:
        wcs=mapIn.getWcs()
    except:
        wcs=None
    assert type(wcs)==Wcs,'Cannot extract WCS from mapIn:'
    naxis1=wcs.getNaxis1()
    naxis2=wcs.getNaxis2()
    #assume map has square pixels
    mapPix=abs(wcs.getCdelt1())
    mapPixArcsec=mapPix*3600.
    mapXRadArcsec=naxis1*mapPix
    mapYRadArcsec=naxis2*mapPix
    
    #set up radius array
    if dRadArcsec==None:
        #use map pixel size as default
        dRadArcsec=mapPixArcsec
    nRad=int(maxRadArcsec/dRadArcsec)
    if verbose:
        print 'computing radial profile to %g arcsec, with step of %g arcsec'%(maxRadArcsec,dRadArcsec)
    #make radius array
    radArr=Double1d(range(nRad))*dRadArcsec
    maxRadDeg=maxRadArcsec/3600.
    #make profile and error arrays
    radProf=Double1d(nRad)
    errProf=Double1d(nRad)

    #crop map as appropriate
    if useWcs:
        [xCtr,yCtr]=wcs.getPixelCoordinates(raCtr,decCtr)
        assert xCtr<naxis1 and xCtr>0 and yCtr<naxis2 and yCtr>0,'centre pixel not in centre of map: [%g,%g]'%(xCtr,yCtr)
    else:
        if not ctrProvided:
            #user central pixel
            xCtr=int(naxis1/2)
            yCtr=int(naxis2/2)
            if verbose:print 'using centre pixel: (x,y)=(%g,%g) (Ra,Dec)=(%g,%g)'%(xCtr,yCtr,raCtr,decCtr)
        #get RA,Dec of centre
        [raCtr,decCtr]=wcs.getWorldCoordinates(xCtr,yCtr)
    validCrop=False
    redCrop=False
    cropRadPix=maxRadDeg/mapPix
    while validCrop==False:
        #set rows
        row1=int(xCtr - cropRadPix)
        row2=int(xCtr + cropRadPix)
        column1=int(yCtr - cropRadPix)
        column2=int(yCtr + cropRadPix)
        if row1 <0:
            #row1 is negative
            #decrease crop radius and start again
            validCrop=False
            redCrop=True
            cropRadPix=cropRadPix + row1
            if verbose:print 'crop row1<0. Decreasing crop radius to %g'%cropRadPix
        elif column1<0:
            #column1 is negative
            #decrease crop radius and start again
            validCrop=False
            redCrop=True
            cropRadPix=cropRadPix + column1
            if verbose:print 'crop column1<0. Decreasing crop radius to %g'%cropRadPix
        elif row2>=naxis1:
            #row2 > naxis1
            #decrease crop radius and start again
            validCrop=False
            redCrop=True
            cropRadPix=cropRadPix - (row2-naxis1)
            if verbose:print 'crop row2>naxis1. Decreasing crop radius to %g'%cropRadPix
        elif column2>=naxis1:
            #column2 > naxis2
            #decrease crop radius and start again
            validCrop=False
            redCrop=True
            cropRadPix=cropRadPix - (column2-naxis2)
            if verbose:print 'crop column2>naxis2. Decreasing crop radius to %g'%cropRadPix
        else:
            validCrop=True
            
    cropRadDeg=cropRadPix*mapPix
    cropRadArcsec=cropRadDeg*3600.
    #compute radius of cropped region (since cropRadPix may have changed)
    if redCrop:
        if verbose:print 'max radius reduced to %g arcsec due to edge of map'%cropRadPix
        
    nCropRad=int(cropRadArcsec/dRadArcsec)
    
    #crop map
    mapCrop=crop(image=mapIn,row1=row1,column1=column1,row2=row2, column2=column2)
    [xCtrCrop,yCtrCrop]=mapCrop.getWcs().getPixelCoordinates(raCtr,decCtr)
    
    #get first histogram
    mapMax=max(mapCrop.getImage())
    hist0=circleHistogram(image=mapCrop,lowCut=0.,highCut=mapMax,\
          bins=1000, centerX=xCtrCrop,centerY=yCtrCrop,\
          radiusArcsec=radArr[0]+dRadArcsec/2.)
    histVals=hist0.getValues()
    histThis=hist0.getFrequencies()
    #number of pixels used
    nPixRad=sum(histThis)
    if nPixRad>0:
        #calculate mean
        meanRad=sum(histThis*histVals)/nPixRad
        #calcuate standard deviation (simple)
        sdRad=SQRT(sum(histThis * histVals**2)/nPixRad - meanRad**2)
        radProf[0]=meanRad
        errProf[0]=sdRad
    else:
        #get value and error from central pixel
        radProf[0]=mapCrop.getIntensity(xCtrCrop,yCtrCrop)
        errProf[0]=mapCrop.getError()[int(xCtrCrop),int(yCtrCrop)]
        
    #loop over radii in cropped map
    for r in range(1,nCropRad):
        #compute new max radius
        rad=radArr[r]+dRadArcsec/2.
        #copy previous histogram
        histPrev=histThis.copy()
        #calculate new histogram
        histThis=circleHistogram(image=mapCrop,lowCut=0.,highCut=mapMax,\
          bins=1000, centerX=xCtrCrop,centerY=yCtrCrop,radiusArcsec=rad).getFrequencies()
        #take difference to get histogram for annulus
        histDiff=histThis-histPrev
        nPixRad=sum(histDiff)
        if nPixRad>0:
            #calculate mean
            meanRad=sum(histDiff*histVals)/nPixRad
            #calculate standard deviation (simple way)
            sdRad=SQRT(sum(histDiff * histVals**2)/nPixRad - meanRad**2)
            radProf[r]=meanRad
            errProf[r]=sdRad
        #print r,rad,sum(histThis),sum(histDiff),radProf[r]
        
    #make profile
    profMap=SourceProfile(radArr,radProf,errProf)
    #set key
    profMap.setKey(key)
    return([profMap,mapCrop])

    