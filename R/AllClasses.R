setClassUnion("NULLorNumeric", c("NULL", "numeric"))
setClassUnion("NULLorMatrix", c("NULL", "matrix"))

validityMethod <- function(object){
  validExpressionSet <- sapply(object, function(eset) validObject(eset))
  if(any(!validExpressionSet)) {
    warning("Each element in ExpressionSetList must be a valid ExpressionSet")
    return(FALSE)
  }
  ##featureNames should be identical for all elements of the list.
  fnList <- sapply(object, featureNames)
  if(class(fnList) != "matrix") {
    warning("all objects in ExpressionSetList must have the same number of features")
    return(FALSE)
  } 

  f <- function(features){
    length(unique(features)) == 1
  }
  
  identicalFeatureNames <- apply(fnList, 1, f)
  if(any(!identicalFeatureNames)) {
    warning("featureNames must be in the same order for each element in the list")
    return(FALSE)
  } else return(TRUE)
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





