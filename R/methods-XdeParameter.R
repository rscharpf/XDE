##The chain in the initialization object always has length 1 for every parameter.
##The chain lengths of the output object has length iterations*output (if writeToFile is TRUE)
setMethod("initialize", "XdeParameter",
          function(.Object,
                   esetList,
                   updates=.updatesDefault(),
                   tuning=.tuningDefault(),
                   hyperparameters=.hyperparametersDefault(P=length(esetList)),
                   output=.outputDefault(),
                   iterations=as.integer(1000),
                   burnin=FALSE,
                   seed=runif(1, 0, 1e6), ##Default is NULL
                   genes=integer(0),
                   studies=integer(),
                   firstMcmc=new.env(),
                   specifiedInitialValues=FALSE,
                   directory=".",
                   phenotypeLabel="",
                   showIterations=TRUE,
                   verbose=FALSE,
		   one.delta=TRUE,
                   studyNames=names(esetList)){
		  if(missing(esetList)) stop("Must provide an ExpressionSetList in order to set up default values for XDE")
		  if(class(esetList) != "ExpressionSetList"){
			  esetList <- as(esetList, "ExpressionSetList")
		  }
            ##check that phenotypeLabel is in each ExpressionSet
            if(!(all(sapply(esetList, function(x, label){ label %in% varLabels(x)}, label=phenotypeLabel)))){
              stop("supplied phenotypeLabel must be present in all ExpressionSets")
            } else{
              .Object@phenotypeLabel <- phenotypeLabel
            }
            if(length(studyNames) != length(esetList)){
              studyNames <- paste("study", 1:length(esetList), sep="")
            }
            .Object@studyNames <- studyNames
            .Object@updates <- updates            
            ##Directory to write files to must exist
            if(!burnin){
		    if(!file.exists(directory)){
			    print(paste(directory, "does not exist. Creating new one"))
			    dir.create(directory, recursive=TRUE)
		    }
            } else{
		    output <- c(1, rep(0, 21))
		    names(output) <- .parameterNames()              
            }
            .Object@tuning <- tuning
            .Object@hyperparameters <- .hyperparametersDefault(length(esetList))
            .Object@output <- output
            .Object@seed <- as.integer(seed)
            .Object@directory <- paste(directory, "/", sep="")
            .Object@verbose <- verbose
	    .Object@one.delta <- one.delta
	    
            ##if firstIteration is not supplied
            if(length(firstMcmc) == 0){
		    .Object@specifiedInitialValues <- FALSE
		    .Object@iterations <- 1
		    .Object@showIterations <- FALSE
		    chain.length <- rep(1, length(2:length(output)))
		    ##only need to store the results from one iteration
		    ##Nothing will be written to file if burnin is TRUE
		    .Object@burnin <- TRUE
		    firstMcmc(.Object) <- .chainInitialize(object=esetList,
							   chain.length=chain.length,
							   verbose=verbose)
		    ##get starting values by simulating from the prior
		    firstMcmc(.Object) <- lastMcmc(xde(paramsMcmc=.Object, esetList=esetList))
            } else{
		    firstMcmc(.Object) <- firstMcmc
            }
            .Object@specifiedInitialValues <- TRUE            
            .Object@burnin <- burnin            
            .Object@iterations <- iterations
            .Object@showIterations <- showIterations
            .Object
          })


setMethod("burnin", "XdeParameter", function(object) object@burnin)
setReplaceMethod("burnin", c("XdeParameter", "logical"), function(object, value){
	object@burnin <- value
	if(value) {
		object@output[2:22] <- rep(0, 21)
	} else {
		object@output[2:22] <- rep(1, 21)
	}
	names(object@output) <- .parameterNames()    
	object
})


setMethod("firstMcmc", "XdeParameter",  function(object) object@firstMcmc$vars)
setReplaceMethod("firstMcmc", c("XdeParameter", "environment"),
                 function(object, value){
                   object@firstMcmc <- value
                   object
                 })
setReplaceMethod("firstMcmc", c("XdeParameter", "list"),
                 function(object, value){
                   eenv <- new.env()
                   assign("vars", value, eenv)
                   names(eenv$vars) <- names(value)
                   object@firstMcmc <- eenv
                   object
                 })

setMethod("directory", "XdeParameter", function(object) object@directory)
setReplaceMethod("directory", "XdeParameter",
                 function(object, value){
			 fname <- paste(value, "/", sep="")
			 if(!file.exists(fname)){
				 print("Directory does not exist. Creating new directory")
				 dir.create(path=fname, recursive=TRUE)
			 } 
			 object@directory <- fname
			 object
                 })

##setMethod("featureData", "XdeParameter", function(object) object@featureData)
##setReplaceMethod("featureData", "XdeParameter", function(object, value){
##  object@featureData <- value
##  object
##})

setMethod("hyperparameters", "XdeParameter", function(object) object@hyperparameters)
setReplaceMethod("hyperparameters", "XdeParameter", function(object, value) {
  object@hyperparameters <- value
  object
})

setMethod("iterations", "XdeParameter", function(object) object@iterations)
setReplaceMethod("iterations", c("XdeParameter", "numeric"),
                 function(object, value){
                   object@iterations <- as.integer(value)
                   object
                 })
setReplaceMethod("iterations", c("XdeParameter", "integer"),
                 function(object, value){
                   object@iterations <- value
                   object
                 })

setMethod("output", "XdeParameter", function(object) object@output)
setReplaceMethod("output", "XdeParameter", function(object, value){
  object@output <- value
  names(object@output) <- .parameterNames()
  ##if output is nonzero, then values are written to file
  if(any(value[2:length(value)] > 0)) object@burnin <- FALSE
  object
})

setMethod("phenotypeLabel", "XdeParameter", function(object) object@phenotypeLabel)
setReplaceMethod("phenotypeLabel", c("XdeParameter", "character"), 
  function(object, value){
    object@phenotypeLabel <- value
    object
    })

setMethod("savedIterations", "XdeParameter", function(object) object@iterations/output(object)[1])

setMethod("show", "XdeParameter",
          function(object){
            cat("Instance of", class(object), "\n\n")
            cat("hyperparameters:\n")
            print(hyperparameters(object)[1:8])
            cat("...\n\n")
            cat("updates (frequency of updates per MCMC iteration):\n")
            print(updates(object)[1:5])
            cat("... \n\n")
            cat("tuning (the epsilon for Metropolis-Hastings proposals): \n")
            print(tuning(object)[1:5])
            cat("...\n\n")
            cat("output (parameters to save (0 = not saved, 1 = saved to log file): \n")
            print(output(object)[2:6])
            cat("...\n\n")
            cat("iterations: ", iterations(object), "\n")
            cat("thin: ", thin(object), "\n")                        
            cat("seed: ", seed(object), "\n")
            cat("notes: ", object@notes, "\n")
            cat("firstMcmc: \n")
            str(firstMcmc(object)[1:5])
            cat("...\n\n")
            cat("showIterations: ", showIterations(object), "\n")
            cat("specifiedInitialValues: ", object@specifiedInitialValues, "\n")
            cat("directory (where to save the MCMC chains): ", directory(object), "\n")
            cat("phenotypeLabel: ", phenotypeLabel(object), "\n")
            cat("studyNames: ", studyNames(object), "\n")
	    cat("one.delta: ", object@one.delta, "\n")
          })
            
setMethod("showIterations", "XdeParameter", function(object) object@showIterations)
setReplaceMethod("showIterations", "XdeParameter", function(object, value){ object@showIterations <- value; return(object)})
setMethod("seed", "XdeParameter", function(object) object@seed)
setReplaceMethod("seed", c("XdeParameter", "integer"),
                 function(object, value){
                   object@seed <- value
                   object
                 })
setReplaceMethod("seed", c("XdeParameter", "numeric"),
                 function(object, value){
                   object@seed <- as.integer(value)
                   object
                 })

setMethod("studyNames", "XdeParameter", function(object) object@studyNames)
setReplaceMethod("studyNames", "XdeParameter",
                 function(object, value){
                   object@studyNames <- value
                   object
                 })

setMethod("thin", "XdeParameter", function(x, ...) output(x)[1])
setReplaceMethod("thin", c("XdeParameter", "numeric"),
                 function(x, value){
                   output(x)[1] <- value
                   x
                 })

setMethod("tuning", "XdeParameter", function(object) object@tuning)
setReplaceMethod("tuning", "XdeParameter",
                 function(object, value){
                   object@tuning <- value
                   names(object@tuning) <- .parameterNames()[4:20]
                   object
                 })

setMethod("updates", "XdeParameter", function(object)  object@updates)
setReplaceMethod("updates", "XdeParameter",
                 function(object, value){
                   object@updates <- value
                   names(object@updates) <- .parameterNames()[4:20]
                   object
                 })






