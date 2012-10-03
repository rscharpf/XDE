library(XDE)
data(expressionSetList)

### Initializing slots common to all models
commonParams <- XDE:::Parameters(object=expressionSetList, clinicalCovariate="adenoVsquamous")
## need to specify: sumS
## not specified by HT: c2, gamma2, a, b, l, t, lambda, theta, ...
## initialize clique
##
commonObject <- XDE:::startingValues(commonParams, expressionSetList, "adenoVsquamous")





### Initializing MI

objectMI <- commonObject ##no specific updates for MI
#paramsMI <- XDE:::Parameters(object=expressionSetList,
#                           clinicalCovariate="adenoVsquamous")
#dms <- dims(paramsMI)
#objectMI <- XDE:::startingValues(paramsMI,
#                                expressionSetList,
#                                "adenoVsquamous")
#orig <- objectMI
#dms2 <- dims(objectMI)
#one.delta <- FALSE
#all.equal(dms, dms2)
##common updates (no specific updates for MI)
objectMI <- modelA_B_MIupdates(objectMI)



### Initializing MII

##paramsMII <- XDE:::ParametersMII(object=expressionSetList,
##                              clinicalCovariate="adenoVsquamous")
##dms <- dims(paramsMII)
##objectMII <- XDE:::startingValues(paramsMII,
##                                 expressionSetList,
##                                 "adenoVsquamous")
##orig <- objectMII
##dms2 <- dims(objectMII)
##one.delta <- FALSE
##all.equal(dms, dms2)

objectMII=as(commonObject, "ParametersMII")
objectMII@Omega=rep(1.0, 1500)
objectMII@zeta=rep(1.0, 1500)
objectMII@D=rep(1.0, 1500)
objectMII@nClique=rep(1L, 1500)
objectMII@oldClique=rep(1L, 1500)
objectMII@nOldComponents=rep(1L, 1500)
objectMII@nNewComponents=rep(1L, 1500)
objectMII@oldComponents=rep(1L, 1500)


setMethod("initialize", signature(.Object="ParametersMII"), function(.Object, ...){
  object=callNextMethod(.Object, ...)
  object@Omega -...
  object@zeta -...
  object@D -...
  object@nClique -...
  object@oldClique -...
  object@nOldComponents -...
  object@nNewComponents -...
  object@oldComponents -...
  return(object)
})

##paramsMII <- new("ParametersMII")
##tmp=XDE:::Parameters(expressionSetList, clinicalCovariate="adenoVsquamous")
##objMII=as(tmp, "ParametersMII")

#dms <- dims(commonParams)
#objectMII <- XDE:::startingValues(paramsMII,expressionSetList,"adenoVsquamous")
orig <- objectMII
#dms2 <- dims(objectMII)
one.delta <- FALSE
#all.equal(dms, dms2)

##common updates
objectMII <- modelA_B_MIupdates(objectMII)


##Model II-specific updates

objectMII <- rupdateSigma2_HyperInverseWishart(object=objectMII,
                                             nTry=5L,
                                             epsilon=0.1,
                                             dryrun=FALSE)

objectMII <- rupdateLSigma2_HyperInverseWishart(object=objectMII,
                                              nTry=5L,
                                              epsilon=0.1,
                                              dryrun=FALSE)

objectMII <- rupdateTSigma2_HyperInverseWishart(object=objectMII,
                                              nTry=5L,
                                              epsilon=0.1,
                                              dryrun=FALSE)

objectMII <- rupdateOmega_HyperInverseWishart(object=objectMII,
                                            dryrun=FALSE)

#objectMII <- rupdateDDeltaStar_HyperInverseWishart(object=objectMII,
#                                                dryrun=FALSE)

#objectMII <- rupdateTau2RDDeltaStar_HyperInverseWishart(object=objectMII,
#                                                      nTry=5L,
#                                                      epsilon=0.1,
#                                                     dryrun=FALSE)

#objectMII <- rupdateBDDeltaStar_HyperInverseWishart(object=objectMII,
#                                                  nTry=5L,
#                                                 epsilon=0.1,
#                                                  dryrun=FALSE)

#objectMII <- rupdateRDDeltaStar_HyperInverseWishart(object=objectMII,
#                                                  nTry=5L,
#                                                  epsilon=0.1,
#                                                  dryrun=FALSE)

