library(XDE)
data(expressionSetList)
## some initial values
params <- XDE:::Parameters(object=expressionSetList,
			   clinicalCovariate="adenoVsquamous")

dms <- dims(params)
params2 <- startingValues(params,
			  expressionSetList,
			  "adenoVsquamous")
dms2 <- dims(params2)
all.equal(dms, dms2)
a <- b <- matrix(NA, 100, 3)
for(i in 1:50){
	## could remove the seed argument -- haakon updates this automatically
	trace(XDE:::rupdateANu, browser)
	params.v2 <- XDE:::rupdateANu(object=params2,
				      nTry=5L,
				      epsilon=0.1,
				      dryrun=FALSE)
	params.v2 <- XDE:::rupdateBDDelta(object=expressionSetList,
					  hyper.params=hyper.params,
					  params=params.v2,
					  nTry=5L,
					  seed=as.integer(ceiling(runif(1,0,100))),
					  epsilon=0.1,
					  dryrun=FALSE)
	a[i,] <- params.v2[["a"]]
	b[i,] <- params.v2[["b"]]
}

















