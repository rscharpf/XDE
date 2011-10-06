library(XDE)
data(expressionSetList)
hyper.params <- XDE:::HyperParams(object=expressionSetList,
				  clinicalCovariate="adenoVsquamous")
## params.v2 <- XDE:::Params(hyper.params)  ## not working -- usually invalid results
params.xde <- new("XdeParameter", esetList=expressionSetList,
		  phenotypeLabel="adenoVsquamous", one.delta=FALSE)
params.v2 <- as(params.xde, "Params")
a <- b <- matrix(NA, 100, 3)
for(i in 1:50){
	## could remove the seed argument -- haakon updates this automatically
	params.v2 <- XDE:::rupdateANu(object=expressionSetList,
				      hyper.params=hyper.params,
				      params=params.v2,
				      nTry=5L,
				      seed=as.integer(ceiling(runif(1,0,100))),
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

















