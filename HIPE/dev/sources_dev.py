import herschel
from herschel.ia.numeric import Double1d,Float1d,Int1d
from herschel.ia.numeric.toolbox.basic import Floor,Min,Max,Exp
FLOOR=herschel.ia.numeric.toolbox.basic.Floor.PROCEDURE
EXP=herschel.ia.numeric.toolbox.basic.Exp.PROCEDURE
MAX=herschel.ia.numeric.toolbox.basic.Max.FOLDR
MIN=herschel.ia.numeric.toolbox.basic.Min.FOLDR
from java.lang import Double,Float,String

global p
p=Float1d(3)

#possible Source Types and number of parameters required    
global srcTypes,sources
srcGaussSpec={
  "name":"Gaussian",
  'minPar':1,
  'maxPar':2,
  'paramNames':["Source width","Peak value"],
  'defaults':[100.,1.0],
  'profFunction':'src2ProfGauss'}
srcLinearSpec={
  "name":"Linear",
  'minPar':1,
  'maxPar':2,
  'paramNames':["Gradiant","Value at r=0"],
  'defaults':[-0.01,1.0],
  'profFunction':'src2ProfLinear'}
srcPowerLawSpec={
  "name":"Power law",
  'minPar':2,
  'maxPar':4,
  'paramNames':["Spectral index","Scale Radius","min radius to extend to","value below min radius"],
  'defaults':[-1.,100.,1.e-3,Double.NaN],
  'profFunction':'src2ProfPowerLaw'}
srcConstantSpec={
  "name":"Constant",
  'minPar':0,
  'maxPar':1,
  'paramNames':["Value"],
  'defaults':[1.],
  'profFunction':'src2ProfConst'}
srcLinConstSpec={
  "name":"LinConst",
  'minPar':1,
  'maxPar':4,
  'paramNames':["Gradiant","Value at r=0","Min Limit","Max Limit"],
  'defaults':[-0.01,1.0, 0.,Double.POSITIVE_INFINITY],
  'profFunction':'src2ProfLinConst'}

#make source types
srcTypes={
  'GAUSSIAN':srcGaussSpec,
  'LINEAR':srcLinearSpec,
  'POWERLAW':srcPowerLawSpec,
  'CONSTANT':srcConstantSpec,
  'LINCONST':srcLinConstSpec}

#define empty dict for sources
sources={}

def listSrcTypes():
    #list all source types
    print '====\nSource Types\n===='
    for typeName in srcTypes:
        src=srcTypes[typeName]
        print '%s (%s):'%(typeName,src['name'])
        print '  %d--%d parameters:'%(src['minPar'],src['maxPar'])
        for p in range(src['maxPar']):
            print '    %d: %s (default=%g)'%(p,src['paramNames'][p],src['defaults'][p])
        print '-----'

def checkSrc(srcTypeIn,paramsIn=None,verbose=False):
    #check source is valid, and fill in default values

    if paramsIn==None:
        #only one parameter, so assume it's a source object
        srcTypeOut=srcTypeIn['type']
        params=srcTypeIn['params']
        if (verbose): print 'Using %s source'%(srcTypeOut)
    else:
        #assume srcType and params passed separately
        #take first 8 characters, and make uppercase
        srcTypeOut=srcTypeIn[:8].upper()
        if srcTypeOut != srcTypeIn:
            print 'Converting to upper case and truncating: %s -> %s'%(srcTypeIn,srcTypeOut)
        
        #check if parameter(s) are in a Double1d, and if not, put them in a Double1d
        try:
            #see if len() works (doesn't work for scalars)
            len(paramsIn)
            params=Double1d(paramsIn)
        except:
            #convert scalar to 1-element Float1d
            params=Double1d([paramsIn])

    #find number of parameters
    nParam=len(params)

    #check number of parameters against list
    srcFound=False
    for typename in srcTypes:
        typex=srcTypes[typename]
        if srcTypeOut[:8]==typename[:8]:
            srcFound=True
            maxParReq=typex['maxPar']
            minParReq=typex['minPar']
            assert nParam>=minParReq and nParam<=maxParReq,\
              '%s source must have %d--%d parameters. %d provided'%(srcTypeOut,minParReq,maxParReq,nParam)
            if (verbose):
                if nParam < maxParReq:
                    print '%d parameters provided. Using %d defaults'%(nParam,maxParReq-nParam)
                else:
                    print '%d parameters provided.'%(nParam)
                    
    #error if source not found
    assert srcFound,'Unknown source type: %s'%srcTypeIn

    #fill in defaults
    defaults=srcTypes[srcTypeOut]['defaults']
    paramsOut=Double1d(maxParReq)
    for p in range(maxParReq):
        if p <= nParam-1:
            if Double.isNaN(params[p]):
                #use default value
                paramsOut[p] = defaults[p]
                if (verbose):print 'Using default value for "%s": %f'%\
                    (srcTypes[srcTypeOut]['paramNames'][p],defaults[p])
            else:
                #use provided value
                paramsOut[p] = params[p]
        else:
            #use default value
            paramsOut[p] = defaults[p]
            if (verbose):print 'Using default value for "%s": %f'%\
                (srcTypes[srcTypeOut]['paramNames'][p],srcTypes[srcTypeOut]['defaults'][p])

    return(srcTypeOut,paramsOut)

def makeSrc(srcTypeIn,paramsIn):
    #generate source object specific source type & params

    #check nParams valid and fill in defaults
    (srcType,params) = checkSrc(srcTypeIn,paramsIn,verbose=True)

    source={'type':srcType,'params':params,'paramNames':srcTypes[srcType]['paramNames']}
    
    return(source)

def setSrc(srcTypeIn,paramsIn=None,keyName=None):
    #set either source object or sourceType & Params as object in global variable

    #check nParams valid and fill in defaults
    (srcType,params) = checkSrc(srcTypeIn,paramsIn,verbose=True)

    #determine key name
    if String.isInstance(keyName):
        #use provided key name
        key=keyName
    else:
        #initialise key
        key=srcType
        for par in params:
            key='%s_%g'%(key,par)

    #add to sources global variable
    sources[key]={'type':srcType,'params':params,'paramNames':srcTypes[srcType]['paramNames']}

    return(key)

def getSrc(key):
    #get source type from global variable
    try:
        #make empty dict
        source=sources[key]
    except:
        print 'ERROR: Source key not found: %s'%key
        #make empty dict
        source={}

    return(source)

def src2Prof(key,radArr,verbose=False):
    #generate profile for given source type

    #get source from sources global variable
    source=getSrc(key)
    srcFound=False
    for typeName in srcTypes:
        if source['type']==typeName:
            exec('srcProf=%s(source,radArr,verbose=verbose)'%srcTypes[typeName]['profFunction'])
            srcFound=True
    assert srcFound, 'source key not found: %s'%key
    
    return(srcProf)

def src2ProfGauss(source,radArr,verbose=False):
    srcWidth=source['params'][0]
    srcPeak=source['params'][1]

    srcProf=srcPeak * EXP(-radArr**2/(2.*srcWidth**2))
    if (verbose):print 'Made profile for %s (MIN=%g; MAX=%g)'%(source['type'],MIN(srcProf),MAX(srcProf))
    return(srcProf)
    
def src2ProfLinear(source,radArr):
    srcGrad=source['params'][0]
    srcZero=source['params'][1]

    srcProf=srcZero + srcGrad * radArr
    
    return(srcProf)
    
def src2ProfPowerLaw(source,radArr):
    srcIdx=source['params'][0]
    srcScalRad=source['params'][1]
    srcMinRad=source['params'][2]
    srcValMinRad=source['params'][3]
    print srcIdx,srcScalRad,srcMinRad,srcValMinRad
    srcProf=Double1d(len(radArr))
    srcProf[radArr.where(radArr>=srcMinRad)]=(radArr[radArr.where(radArr>=srcMinRad)]/srcScalRad)**srcIdx
    srcProf[radArr.where(radArr<srcMinRad)]=srcValMinRad
    
    return(srcProf)
    
def src2ProfConst(source,radArr):
    srcVal=source['params'][0]

    srcProf=srcVal
    
    return(srcProf)

    
def src2ProfLinConst(source,radArr):
    srcGrad=source['params'][0]
    srcZero=source['params'][1]
    srcMin=source['params'][2]
    srcMax=source['params'][3]


    srcProf=srcZero + srcGrad * radArr
    if MAX(srcProf) > srcMax:
        srcProf[srcProf.where(srcProf > srcMax)]=srcMax
    if MIN(srcProf) < srcMin:
        srcProf[srcProf.where(srcProf < srcMin)]=srcMin
    return(srcProf)