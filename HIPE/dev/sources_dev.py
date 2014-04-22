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

import herschel
from herschel.ia.numeric import Double1d,Float1d,Int1d
from herschel.ia.numeric.toolbox.basic import Floor,Min,Max,Exp,Cos,Abs
FLOOR=herschel.ia.numeric.toolbox.basic.Floor.PROCEDURE
EXP=herschel.ia.numeric.toolbox.basic.Exp.PROCEDURE
COS=herschel.ia.numeric.toolbox.basic.Cos.PROCEDURE
MAX=herschel.ia.numeric.toolbox.basic.Max.FOLDR
MIN=herschel.ia.numeric.toolbox.basic.Min.FOLDR
ABS=herschel.ia.numeric.toolbox.basic.Abs.FUNCTION
from java.lang import Double,Float,String
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
        return(SourceProfile(radArr,srcProf,self))

    def calcArea(self):
        #run calcArea associated with typeTemplate
        srcArea=self.typeTemplate.calcArea(self.params)
        return(srcArea)
    
    def initAvailableTypes(self):
        self.availableTypes={}
        self.availableTypes["GAUSSIAN"]=GaussianSourceType()
        self.availableTypes['LINEAR']=LinearSourceType()
        self.availableTypes['POWERLAW']=PowerLawSourceType()
        self.availableTypes['CONSTANT']=ConstantSourceType()
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
    
    def __str__(self):
        result=self.__class__
        result="\nType name: %s"%self.desc
        result=result+"\nParameters: %d (%d required)"%(self.maxPar,self.minPar)
        for p in range(self.maxPar):
            result=result+"\n  %d: %s (default %g)"%(p,self.paramNames[p],self.defaults[p])
        return(result)
        
    
class NoneSourceType:
    def __init__(self):
        self.type='NONE'
        self.desc='NONE source type'
        self.minPar=0
        self.paramNames=[]
        self.defaults=[]
        self.maxPar=len(self.paramNames)
        
    def calcProfile(self,radArr,params):
        return(Double1d(len(radArr),Double.NaN))
        
    def calcArea(self,params):
        return(Double.NaN)
        
    def calcZeroVal(self,params):
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
        
    def calcArea(self,params):
        #integrate Gaussian
        srcWidth=params[0]
        srcPeak=params[1]
        srcArea=srcPeak * 2*PI*srcWidth**2
        return(srcArea)
        
    def calcZeroVal(self,params):
        zeroVal=params[1]
        return(zeroVal)
        
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
        
    def calcScaleWidth(self,params):
        #set as infinity
        scaleWidth=Double.POSITIVE_INFINITY
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
        
    def calcArea(self,params):
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
            print 'Warning: no analytical value for area of LINCONST source'
            srcArea=Double.NaN
            
        return(srcArea)
            
    def calcZeroVal(self,params):
        zeroVal=params[1]
        if zeroVal < params[2]: zeroVal=params[2]
        if zeroVal > params[3]: zeroVal=params[3]
        
    def calcScaleWidth(self,params):
        #calculate half-width-half-max
        if params[0] < 0:
            scaleWidth=(params[1]-params[2])/params[0]
        else:
            scaleWidth=(params[3]-params[1])/params[0]
        return(scaleWidth)
    
#-------------------------------------------------------------------------------
#===============================================================================
#=====                       DEFINE SourceProfile CLASS                    =====
#===============================================================================
#-------------------------------------------------------------------------------
class SourceProfile(object):
    def __init__(self,radArr=None,profile=None,originator=None):
        #set radArr
        self.setRadArr(radArr)
        #set profile
        self.profile=profile or None
        #check radArr and profile are compatible
        self.check()
        #set original Source object
        self.setOrig(originator)

    def setRadArr(self,radArr=None):
        #set radArr, nRad and maxRad
        self.radArr=radArr or None
        if self.radArr:
            self.nRad=len(radArr)
            self.maxRad=MAX(radArr)
        else:
            self.nRad=0
            self.maxRad=Double.NaN

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

    def copy(self):
        new=SourceProfile(radArr=self.radArr,profile=self.profile,originator=self.origSrc)
        return(new)
        
    def generate(self,src,radArr):
        #generate from Source object
        self.setradArr(radArr)
        self.profile=src.calcProfile(radArr)
        self.setOrig(src)

    def setOrig(self,src):
        if src:
            self.fromSrc=True
            self.origSrc=src
        else:
            self.fromSrc=False
            self.origSrc=NoneSourceType()
                        
    def calcArea(self,forceNumerical=False,forceAnalytical=False):
        #print 'calc area'
        if forceNumerical:
            doNum=True
        elif forceAnalytical:
            doNum=False
            if self.fromSrc:
                srcArea=self.origSrc.calcArea()
            else:
                srcArea=Double.NaN()
        else:
            if self.fromSrc:
                #try to get analytical area from origSrc
                srcArea=self.origSrc.calcArea()
                if Double.isNaN(srcArea):
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
        
    def add(self,prof2):
        #add profile to another
        
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
        
        #multiply profiles
        self.profile = self.profile + prof2
        #no longer based on SourceType
        self.origSrc=NoneSourceType()
        self.fromSrc=False
        
        return(self)
        
    def mult(self,prof2):
        #multiply profile by another
        
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
        
        #multiply profiles
        self.profile = self.profile * prof2
        #no longer based on SourceType
        self.origSrc=NoneSourceType()
        self.fromSrc=False
        
        return(self)

    def rad2k(self):
        #make k-array
        kArr=Double1d(self.nRad)
        nK=len(kArr)
        for k in range(len(kArr)-1):
            #print k,radArr[(nRad-1)-k]
            kArr[k+1]=2*PI/self.radArr[(self.nRad-1)-k]
        return(kArr)

    def calcRft(self):
        #compute Radial-Fourier Transform (RFT) of source profile
        #make k-array
        kArr=self.rad2k()
        nK=len(kArr)
        rftProfile=Double1d(nK)
        integrator=TrapezoidalIntegrator(0.,self.maxRad)
        profInterp=CubicSplineInterpolator(self.radArr,self.profile)
        for k in range(nK):
            rftArg_k=rftArg(kArr[k],profInterp)
            print 'rftArg calculated:',k,kArr[k]
            rftProfile[k]=integrator.integrate(rftArg_k)
            print 'integral calculated:',rftProfile[k]
            #print k,kArr[k],rftProfile[k]
        return(RftProfile(kArr,rftProfile))
        
    def __str__(self):
        result='%s'%self.__class__
        result=result+'\nnRad=%d'%self.nRad
        if self.nRad > 0:
            result=result+'\nmaxRad=%g'%self.maxRad
            if self.nRad<=5:
                result=result+'\nradArr: '+self.radArr
                result=result+'\nprofile: '+self.profile
            else:
                result=result+'\nradArr: [%g, %g, ... %g, %g]'%(self.radArr[0],self.radArr[1],self.radArr[-2],self.radArr[-1])
                result=result+'\nprofile: [%g, %g, ... %g, %g]'%(self.profile[0],self.profile[1],self.profile[-2],self.profile[-1])
            result=result+'\nfromSrc: %r'%self.fromSrc
            if self.fromSrc:
                result=result+'\norigSrc: %s'%self.origSrc.key
        return(result)

#-------------------------------------------------------------------------------
#===============================================================================
#=====                       DEFINE RftProfile CLASS                     =====
#===============================================================================
#-------------------------------------------------------------------------------
class RftProfile(object):
    def __init__(self,kArr=None,rtfProfile=None):
        self.setKArr(kArr)
        self.rft=rftProfile or None
        
    def setKArr(kArr=None):
        self.kArr=kArr or None
        if self.kArr:
            self.nK=len(kArr)
            self.maxK=MAX(self.kArr)

    def generate(self,src,kArr):
        self.kArr=kArr
        self.nK=len(self.kArr)
        
        radArr=self.k2rad()
        realProfile=Profile().generate(src,radArr)
        self.rft=realProfile.calcRft()
    
    def k2rad(self):
        #make k-array
        radArr=Double1d(self.nK)
        nR=len(radArr)
        for r in range(nR-1):
            #print k,radArr[(nRad-1)-k]
            radArr[r+1]=2*PI/self.kArr[(self.nK-1)-r]
        return(radArr)
        
    def calcInvRft(self):
        #compute Radial-Fourier Transform (RFT) of source profile
        #make k-array
        radArr=self.k2rad()
        realProfile=Double1d(self.nK)
        integrator=TrapezoidalIntegrator(0.,self.maxK)
        rftInterp=CubicSplineInterpolator(self.kArr,self.rft)
        for r in range(self.nR):
            profArg_x=rftArg(self.radArr[r],rftInterp)
            realProfile[r]=integrator.integrate(profArg_x)
            print r,radArr[r],realProfile[r]
        return(radArr,realProfile)

#def testJ0():
#    #make test arrays of 0th order Bessel functions
#    xarr=Float1d(range(0,2000))/100.
#    nX=len(xarr)
#    J0arr=Double1d(nX)
#    for x in range(nX):
#        J0arr[x]=calcJ0(xarr[x])
#    return(xarr,J0arr)
    
class j0arg(RealFunction):
    #returns COS(z * COS(theta)) for use in Bessel functions
    def __init__(self,z):
        self.z=z
    def calc(self,theta):
        return COS(self.z*COS(theta))        

def calcJ0(z):
    #calculate J0 Bessel function at value z (z=kR)
    #J= (1/pi) * int_0^pi{cos(z cos(x)) dx}
    integrator=TrapezoidalIntegrator(0.,PI)
    j0Arg_z=j0arg(z)
    j0=integrator.integrate(j0Arg_z) / PI
    return(j0)
    
class rftArg(RealFunction):
    #returns J0(k*r) * radInt(r)
    def __init__(self,k,radInt):
        self.k=k
        self.radInt=radInt
    def calc(self,r):
        return (calcJ0(self.k*r)*self.radInt(r))
    