library(XDE)
data(expressionSetList)

###Initializing slots common to all Models
commonParams <- XDE:::Parameters(object=expressionSetList, clinicalCovariate="adenoVsquamous")
commonObject <- XDE:::startingValues(commonParams, expressionSetList, "adenoVsquamous")



###Initializing Model C

objectC=as(commonObject, "ParametersC")
objectC@betag=rep(1.0, 1500) 
objectC@nNeighbour=rep(1L, 1500)
objectC@neighbour=rep(1L, 1500)
objectC@alpha=rep(1.0, 1500)
objectC@beta=rep(1.0, 1500)

setMethod("initialize", signature(.Object="ParametersC"), function(.Object, ...){
  object=callNextMethod(.Object, ...)
  object@betag -...
  object@nNeighbour -...
  object@neighbour -...
  object@alpha -...
  object@beta -...
    return(object)
})

#paramsC <- new("ParametersC")

#dms <- dims(paramsC)
#objectC <- XDE:::startingValues(paramsC,
#                                expressionSetList,
#                                "adenoVsquamous")

##
#orig <- objectC
#dms2 <- dims(objectC)
#one.delta <- FALSE
#all.equal(dms, dms2)
##


##Updates for all models
objectC <- modelA_B_MIupdates(objectC)

##C-specific updates
objectC <- rupdateAlpha(object=objectC,
                        nTry=5L,
                        epsilon=0.1,
                        dryrun=FALSE)

objectC <- rupdateBeta(object=objectC,
                       nTry=5L,
                       epsilon=0.1,
                       dryrun=FALSE)

objectC <- rupdateBetag(object=objectC,
                        nTry=5L,
                        epsilon=0.1,
                        dryrun=FALSE)

##Model C update common between C and D
objectC <- rupdateDeltaDDelta_MRF2(object=objectC,
                                   nTry=5L,
                                   dryrun=FALSE,
                                   one.delta=FALSE)

##Model C update common between C and D
objectC <- rupdateDelta_MRF2(object=objectC,
                             nTry=5L,
                             dryrun=FALSE,
                             one.delta=FALSE)




###Initializing Model D


##paramsD <- XDE:::ParametersD(object=expressionSetList,
##                             clinicalCovariate="adenoVsquamous")
##paramsD <- new("ParametersD")
##objectD=as(commonObject, "ParametersC")
##tmp=XDE::Parameters(expressionSetList, clinicalCovariate="adenoVsquamous")
##objD=as(tmp, "ParametersD")


objectD=as(commonObject, "ParametersD")
objectD@nNeighbour=rep(1L, 1500)
objectD@neighbour=rep(1L, 1500)
objectD@alpha=rep(1.0, 1500)
objectD@beta=rep(1.0, 1500)


setMethod("initialize", signature(.Object="ParametersD"), function(.Object, ...){
object=callNextMethod(.Object, ...)
object@nNeighbour -...
object@neighbour -...
object@alpha -...
object@beta -...
return(object)
})

#dms <- dims()
#objectD <- XDE:::startingValues(paramsD,expressionSetList,"adenoVsquamous")
#orig <- objectD
#dms2 <- dims(objectD)
#one.delta <- FALSE
#all.equal(dms, dms2)
                    
##common updates
objectD <- modelA_B_MIupdates(objectD) 

##D-specifc updates
objectD <- rupdateAlpha_onedelta(object=objectD,
                        nTry=5L,
                        epsilon=0.1,
                        dryrun=FALSE)

objectD <- rupdateBeta_onedelta(object=objectD,
                       nTry=5L,
                       epsilon=0.1,
                       dryrun=FALSE)

##Model C update common between C and D
objectD <- rupdateDeltaDDelta_MRF2(object=objectD,
                                   nTry=5L,
                                   dryrun=FALSE,
                                   one.delta=TRUE) 

##Model D update common between C and D
objectD <- rupdateDelta_MRF2(object=objectD,
                             nTry=5L,
                             dryrun=FALSE,
                             one.delta=TRUE) 



