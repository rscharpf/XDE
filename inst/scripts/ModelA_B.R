library(XDE)
data(expressionSetList)

## initializing A
paramsA <- XDE:::Parameters(object=expressionSetList,
                           clinicalCovariate="adenoVsquamous")
##No additional slots for A
dms <- dims(paramsA)
objectA <- XDE:::startingValues(paramsA,
                                expressionSetList,
                                "adenoVsquamous")
orig <- objectA
dms2 <- dims(objectA)
one.delta <- FALSE
all.equal(dms, dms2)

##objectA updates (no specific updates for MA)
trace(modelA_B_MIupdates, browser)
objectA <- modelA_B_MIupdates(objectA)


## initializing B
paramsB <- XDE:::Parameters(object=expressionSetList,
                             clinicalCovariate="adenoVsquamous")
##No additional slots for B
dms <- dims(paramsB)
objectB <- XDE:::startingValues(paramsB,
                                expressionSetList,
                                "adenoVsquamous")
orig <- objectB
dms2 <- dims(objectB)
one.delta <- FALSE
all.equal(dms, dms2)

##objectB updates (no specific updates for MA)
objectB <- modelA_B_MIupdates(objectB)

