setClassUnion("NULLorNumeric", c("NULL", "numeric"))
setClassUnion("NULLorMatrix", c("NULL", "matrix"))
validityMethod <- function(object){
	validExpressionSet <- sapply(object, function(eset) validObject(eset))
	if(any(!validExpressionSet)) {
		return("Each element in ExpressionSetList must be a valid ExpressionSet")
	}
	##featureNames should be identical for all elements of the list.
	fnList <- sapply(object, featureNames)
	if(class(fnList) != "matrix") {
		return("all objects in ExpressionSetList must have the same number of features")
	}

	f <- function(features){
		length(unique(features)) == 1
	}

	identicalFeatureNames <- apply(fnList, 1, f)
	if(any(!identicalFeatureNames)) {
		return("featureNames must be in the same order for each element in the list")
	}
	return(TRUE)
}
setClass("ExpressionSetList", representation("list"),
         validity=validityMethod)




##output of MCMC -- only the posterior averages
setClass("XdeMcmc", ##contains="ExpressionSetList",
         representation(studyNames="character",
                        featureNames="character",
                        iterations="numeric",
                        directory="character",
                        seed="integer",
                        output="numeric",
                        lastMcmc="environment",
                        posteriorAvg="NULLorMatrix",
                        bayesianEffectSize="NULLorMatrix"))

##Parameters needed to fit XDE, including starting values
setClass("XdeParameter",
         representation(updates = "numeric",
                        tuning = "numeric",
                        hyperparameters = "numeric",
                        output = "numeric",
                        iterations = "numeric",
                        burnin="logical",
                        seed = "integer",
                        notes = "character",
                        firstMcmc="environment",
                        showIterations="logical",
                        specifiedInitialValues="logical",
                        directory="character",
                        phenotypeLabel="character",
                        verbose="logical",
                        studyNames="character",
			one.delta="logical"))


##setClass("HyperParams", contains="list")
##setClass("Params", contains="list")
##setClass("foo", contains="list",
##	 representation(A="numeric",
##			B="numeric"))
##setValidity("foo", function(object){
##	c("A", "B") %in% names(object)
##}


setClass("Parameters",
	 representation(seed="integer",
			data="numeric",##"environment",
			phenodata="integer",
			G="integer",
			Q="integer",
			S="integer",
			alphaA="numeric",
			alphaB="numeric",
			betaA="numeric",
			betaB="numeric",
			pA0="numeric",
			pA1="numeric",
			pB0="numeric",
			pB1="numeric",
			nuR="numeric",
			nuRho="numeric",
			alphaXi="numeric",
			betaXi="numeric",
			c2Max="numeric",
			alphaEta="numeric",
			betaEta="numeric",
			pOmega0="numeric",
			lambdaOmega="numeric",
			lambdaKappa="numeric",
			gamma2="numeric",
			c2="numeric",
			tau2Rho="numeric",
			tau2R="numeric",
			a="numeric",
			b="numeric",
			l="numeric",
			t="numeric",
			lambda="numeric",
			theta="numeric",
			phi="numeric",
			sigma2="numeric",
			r="numeric",
			rho="numeric",
			nu="numeric",
			delta="numeric",
			Delta="numeric",
			xi="numeric"))
setGeneric("pnames", function(object) standardGeneric("pnames"))

setValidity("Parameters", function(object){
	if(length(exprs(object)) != object@G*sum(object@S))
		return("data should have length G*sum(S)")
})##nipun: slots are common to Models A,B and MI 


##nipun: ParametersD
setClass("ParametersD",contains= "Parameters",
         representation(nNeighbour="integer",
                        neighbour="integer",
                        alpha="numeric",
                        beta="numeric"))
##setGeneric("pnamesII", function(object) standardGeneric("pnames")) ##:?

setValidity("ParametersD", function(object){
  if(length(exprs(object)) != object@G*sum(object@S))
    return("data should have length G*sum(S)")
})



##nipun: Parameters specific to model C
setClass("ParametersC",contains= "ParametersD",
         representation(betag="numeric"))
##setGeneric("pnamesII", function(object) standardGeneric("pnames")) ##??

setValidity("ParametersC", function(object){
  if(length(exprs(object)) != object@G*sum(object@S))
    return("data should have length G*sum(S)")
})


##nipun: Parameters specific to model II
setClass("ParametersMII",contains= "Parameters",
         representation(Omega="numeric",
                        zeta = "numeric",
                        D = "numeric",
                        nClique="integer",
                        oldClique="integer",
                        nOldComponents="integer",
                        nNewComponents="integer",
                        oldComponents="integer"))
##setGeneric("pnamesII", function(object) standardGeneric("pnames")) ##??

setValidity("ParametersMII", function(object){
  if(length(exprs(object)) != object@G*sum(object@S))
    return("data should have length G*sum(S)")
})