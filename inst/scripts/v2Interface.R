library(XDE)
data(expressionSetList)
## some initial values
## one.delta is 'model B'
initial.paramsA <- XDE:::Parameters(object=expressionSetList,
				    clinicalCovariate="adenoVsquamous",
				    one.delta=TRUE)
##
## there should probably be a model specification here (model A or model B)
##
## There should be a C entry point to get initial values for all
## parameters under model A or model B... could probably use the 'v1'
## interface for this
##
paramsA <- XDE:::startingValues(initial.paramsA, expressionSetList,
				"adenoVsquamous", one.delta=TRUE)

## ideally, this should all be one step
initial.paramsB <- XDE:::Parameters(object=expressionSetList, clinicalCovariate="adenoVsquamous",
				    one.delta=FALSE)
paramsB <- XDE:::startingValues(initial.paramsB, expressionSetList, "adenoVsquamous", one.delta=FALSE)
orig <- params2
dms2 <- dims(params2)
one.delta <- FALSE
all.equal(dms, dms2)
for(i in 1:50){
	##
	## An entry point for each of the C routines in Rintervace_v2.cpp
	##
	## Updates that can be used by all model variants
	##
	a=paramsA@a
	paramsA <- rupdateANu(object=paramsA, nTry=5L,
			      epsilon=0.1, dryrun=FALSE)
	paramsA <- rupdateTau2RhoNu(object=paramsA,
				    nTry=5L, epsilon=0.1, dryrun=FALSE)
	paramsA <- rupdateNu(object=paramsA, dryrun=FALSE)
	paramsA <- rupdateGamma2(object=paramsA, dryrun=FALSE)
	paramsA <- rupdateRhoNu(object=paramsA,
				nTry=5L,
				epsilon=0.1,
				dryrun=FALSE)
	paramsA <- rupdateRhoGamma2(object=paramsA,
				    nTry=5L,
				    epsilon=0.1,
				    dryrun=FALSE)
	paramsA <- rupdatePhi(object=paramsA,
				 nTry=5L,
				 epsilon=0.1,
				 dryrun=FALSE)
	paramsA <- rupdateTheta(object=paramsA,
				 nTry=5L,
				 epsilon=0.1,
				 dryrun=FALSE)
	paramsA <- rupdateLambda(object=paramsA,
				 nTry=5L,
				 epsilon=0.1,
				 dryrun=FALSE)
	paramsA <- rupdateT(object=paramsA,
				 nTry=5L,
				 epsilon=0.1,
				 dryrun=FALSE)
	paramsA <- rupdateL(object=paramsA,
				 nTry=5L,
				 epsilon=0.1,
				 dryrun=FALSE)
	paramsA <- rupdateLambdaPhi(object=paramsA,
				      nTry=5L,
				      epsilon=0.1,
				      dryrun=FALSE)
	paramsA <- rupdateThetaPhi(object=paramsA,
				      nTry=5L,
				      epsilon=0.1,
				      dryrun=FALSE)
	##
	## finished: updates that can be used by all model variants
	##
	##
	## updates only for model A
	##
	paramsA <- rupdateXi(object=paramsA, dryrun=FALSE, one.delta=TRUE)
	##
	## finished: updates that can be used only for model A
	##
	## updates only for model B
	paramsB <- rupdateXi(object=paramsB, dryrun=FALSE, one.delta=FALSE)

	params2 <- rupdateBDDelta(object=params2,
				  nTry=5L,
				  epsilon=0.1,
				  dryrun=FALSE)
	## stop

	params2 <- rupdateTau2RDDelta(object=params2,
				      nTry=5L,
				      epsilon=0.1,
				      dryrun=FALSE)

	## NaN's appear in nu after the update...?
	params2 <- rupdateDDelta(object=params2,
				 dryrun=FALSE)

	params2 <- rupdateC2(object=params2,
			     nTry=5L,
			     dryrun=FALSE)

	params2 <- rupdateC2DDelta(object=params2,
				   nTry=5L,
				   dryrun=FALSE)



	params2 <- rupdateGamma2Nu(object=params2,
				   nTry=5L,
				   epsilon=0.1,
				   dryrun=FALSE)

	params2 <- rupdateRDDelta(object=params2,
				 nTry=5L,
				 epsilon=0.1,
				 dryrun=FALSE)

	params2 <- rupdateRC2(object=params2,
				 nTry=5L,
				 epsilon=0.1,
				 dryrun=FALSE)





	params2 <- rupdateSigma2(object=params2,
				    nTry=5L,
				    epsilon=0.1,
				    dryrun=FALSE)










	if(!one.delta){
		params2 <- rupdateXi(object=params2,
				     dryrun=FALSE,
				     one.delta=FALSE)
		params2 <- rupdateDeltaDDelta(object=params2,
					      nTry=5L,
					      dryrun=FALSE,
					      one.delta=FALSE)
	}
	if(one.delta){
		## must make sure that the deltas
		## for a gene are the same

		params2 <- rupdateDeltaDDelta(object=params2,
					      nTry=5L,
					      dryrun=FALSE,
					      one.delta=TRUE)
	}
	params2 <- rupdateLSigma2(object=params2,
				  nTry=5L,
				  epsilon=0.1,
				  dryrun=FALSE)

	params2 <- rupdateTSigma2(object=params2,
				      nTry=5L,
				      epsilon=0.1,
				      dryrun=FALSE)





	if(FALSE){
		## add appropriate slots to Parameters object
		params2 <- rupdateDeltaDDelta_MRF(object=params2,
						  nTry=5L,
						  epsilon=0.1,
						  dryrun=FALSE)

		params2 <- rupdateAlpha_MRF(object=params2,
					    nTry=5L,
					    epsilon=0.1,
					    dryrun=FALSE)
		params2 <- rupdateBeta_MRF(object=params2,
					   nTry=5L,
					   epsilon=0.1,
					   dryrun=FALSE)
		params2 <- rupdateBetag_MRF(object=params2,
					    nTry=5L,
					    epsilon=0.1,
					    dryrun=FALSE)
	}
}

















	rupdate <- function(object,
			    type=c("ANu",
			           "BDDelta",
			           "Tau2Rho",
			           "Tau2DDelta"),
			    nTry=10L,
			    epsilon=0.1,
			    dryrun=FALSE){
		update.type <- paste("rupdate", type, sep="")
		.C(update.type, ...)
