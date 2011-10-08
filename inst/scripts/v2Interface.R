library(XDE)
data(expressionSetList)
## some initial values
params <- XDE:::Parameters(object=expressionSetList,
			   clinicalCovariate="adenoVsquamous")

dms <- dims(params)
params2 <- XDE:::startingValues(params,
			  expressionSetList,
			  "adenoVsquamous")
dms2 <- dims(params2)
all.equal(dms, dms2)
a <- b <- matrix(NA, 100, 3)
for(i in 1:50){
	## could remove the seed argument -- haakon updates this automatically
	##trace(XDE:::rupdateANu, browser)
	params2 <- XDE:::rupdateANu(object=params2,
				    nTry=5L,
				    epsilon=0.1,
				    dryrun=FALSE)
	params2 <- XDE:::rupdateBDDelta(object=params2,
					nTry=5L,
					epsilon=0.1,
					dryrun=FALSE)
	a[i,] <- params.v2[["a"]]
	b[i,] <- params.v2[["b"]]
}

















