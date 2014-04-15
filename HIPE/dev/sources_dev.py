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
from herschel.ia.numeric.toolbox.basic import Floor,Min,Max,Exp,Cos
FLOOR=herschel.ia.numeric.toolbox.basic.Floor.PROCEDURE
EXP=herschel.ia.numeric.toolbox.basic.Exp.PROCEDURE
COS=herschel.ia.numeric.toolbox.basic.Cos.PROCEDURE
MAX=herschel.ia.numeric.toolbox.basic.Max.FOLDR
MIN=herschel.ia.numeric.toolbox.basic.Min.FOLDR
from java.lang import Double,Float,String
from java.lang.Math import PI
from herschel.ia.numeric.toolbox import RealFunction
from herschel.ia.numeric.toolbox.integr import TrapezoidalIntegrator
from herschel.ia.numeric.toolbox.interp import CubicSplineInterpolator


class Source:
    
    def __init__(self,srcTypeIn,paramsIn):
        #set available types
        self.initAvailableTypes()

        #set variables
        self.type=srcTypeIn
        self.params=paramsIn
        self.check(verbose=True)
        self.setKey()

    def setKey(self):
        #set key of source
        key=self.type
        for par in self.params:
            key='%s_%g'%(key,par)
        self.key=key
        return(key)
    
    #def save(self):
    #    Sources[self.key]=self
        
    def check(self,verbose=False):
        typeNew=self.type[:8].upper()
        if typeNew != self.type:
            print 'Converting to upper case and truncating: %s -> %s'%(self.type,typeNew)
        print self.params
        try:
            #see if len() works (doesn't work for scalars)
            len(self.params)
            paramsIn=Double1d(self.params)
        except:
            #convert scalar to 1-element Float1d
            paramsIn=Double1d([self.params])
        #find number of parameters
        nParamIn=len(paramsIn)
        print paramsIn
        #check number of parameters against available Types
        srcFound=False
        for typename in self.availableTypes:
            typex=self.availableTypes[typename]
            if typeNew==typename[:8]:
                srcFound=True
                maxParReq=typex['maxPar']
                minParReq=typex['minPar']
                assert nParamIn>=minParReq and nParamIn<=maxParReq,\
                  '%s source must have %d--%d parameters. %d provided'%(typename,minParReq,maxParReq,nParamIn)
                if (verbose):
                    if nParamIn < maxParReq:
                        print '%d parameters provided. Using %d defaults'%(nParamIn,maxParReq-nParamIn)
                    else:
                        print '%d parameters provided.'%(nParamIn)
                        
                self.type=typename
                self.nParam=self.availableTypes[self.type]['maxPar']
                self.defaults=self.availableTypes[self.type]['defaults']
                self.paramNames=self.availableTypes[self.type]['paramNames']
                self.typeDesc=self.availableTypes[self.type]['name']

        #error if source not found
        assert srcFound,'Unknown source type: %s'%self.type
        
        #fill in default parameters
        paramsOut=Double1d(self.nParam)
        for p in range(maxParReq):
            if p <= nParamIn-1:
                if Double.isNaN(paramsIn[p]):
                    #use default value
                    paramsOut[p] = self.defaults[p]
                    if (verbose):print 'Using default value for "%s": %f'%\
                        (self.paramNames[p],self.defaults[p])
                else:
                    #use provided value
                    paramsOut[p] = paramsIn[p]
            else:
                #use default value
                paramsOut[p] = self.defaults[p]
                if (verbose):print 'Using default value for "%s": %f'%\
                    (self.paramNames[p],self.defaults[p])
        self.params=paramsOut

    def calcProfile(self,radArr):
        if self.type=='GAUSSIAN':
            #generate profile for Gaussian source
            srcWidth=self.params[0]
            srcPeak=self.params[1]
        
            srcProf=srcPeak * EXP(-radArr**2/(2.*srcWidth**2))
            
        elif self.type=='LINEAR':
            #generate profile for Linear source
            srcGrad=self.params[0]
            srcZero=self.params[1]
        
            srcProf=srcZero + srcGrad * radArr
            
        elif self.type=='POWERLAW':
            #generate profile for Power Law source
            srcIdx=self.params[0]
            srcScalRad=self.params[1]
            srcMinRad=self.params[2]
            srcValMinRad=self.params[3]

            srcProf=Double1d(len(radArr))
            srcProf[radArr.where(radArr>=srcMinRad)]=(radArr[radArr.where(radArr>=srcMinRad)]/srcScalRad)**srcIdx
            srcProf[radArr.where(radArr<srcMinRad)]=srcValMinRad
            
        elif self.type=='CONSTANT':
            #generate profile for Constant source
            
            srcVal=self.params[0]
            srcProf=Double1d(len(radArr),srcVal)
            
        elif self.type=='LINCONST':
            #generate profile for Linear-Constant (linear with max and/or min limits)
            srcGrad=self.params[0]
            srcZero=self.params[1]
            srcMin=self.params[2]
            srcMax=self.params[3]
        
            srcProf=srcZero + srcGrad * radArr
            if MAX(srcProf) > srcMax:
                srcProf[srcProf.where(srcProf > srcMax)]=srcMax
            if MIN(srcProf) < srcMin:
                srcProf[srcProf.where(srcProf < srcMin)]=srcMin
            
        return(SourceProfile(radArr,srcProf))

    def initAvailableTypes(self):
        self.availableTypes={}
        self.availableTypes["GAUSSIAN"]={
          "name":"Gaussian",
          'minPar':1,
          'maxPar':2,
          'paramNames':["Source width","Peak value"],
          'defaults':[100.,1.0]}

        self.availableTypes['LINEAR']={
          "name":"Linear",
          'minPar':1,
          'maxPar':2,
          'paramNames':["Gradiant","Value at r=0"],
          'defaults':[-0.01,1.0]}
        
        self.availableTypes['POWERLAW']={
          "name":"Power law",
          'minPar':2,
          'maxPar':4,
          'paramNames':["Spectral index","Scale Radius","min radius to extend to","value below min radius"],
          'defaults':[-1.,100.,1.e-3,Double.NaN]}
            
        self.availableTypes['CONSTANT']={
          "name":"Constant",
          'minPar':0,
          'maxPar':1,
          'paramNames':["Value"],
          'defaults':[1.]}
            
        self.availableTypes['LINCONST']={
          "name":"LinConst",
          'minPar':1,
          'maxPar':4,
          'paramNames':["Gradiant","Value at r=0","Min Limit","Max Limit"],
          'defaults':[-0.01,1.0, 0.,Double.POSITIVE_INFINITY]}

    def listSrcTypes(self):
        #list all source types
        print '====\nSource Types\n===='
        for typeName in self.availableTypes:
            src=self.availableTypes[typeName]
            print '%s (%s):'%(typeName,src['name'])
            print '  %d--%d parameters:'%(src['minPar'],src['maxPar'])
            for p in range(src['maxPar']):
                print '    %d: %s (default=%g)'%(p,src['paramNames'][p],src['defaults'][p])
            print '-----'

class SourceProfile:
    def __init__(self,radArr=None,profile=None):

        self.setRadArr(radArr)
        self.profile=profile or None

    def setRadArr(self,radArr=None):
        self.radArr=radArr or None
        if self.radArr:
            self.nRad=len(radArr)
            self.maxRad=MAX(radArr)

    def generate(self,src,radArr):
        #generate from Source object
        self.radArr=radArr
        self.profile = src.calcProfile(radArr)
    
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

class RftProfile:
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
    