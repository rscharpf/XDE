setMethod("initialize", "XdeMcmc",
          function(.Object,
                   studyNames="",
                   featureNames=character(),
                   iterations=numeric(),
                   seed=integer(),
                   output=numeric(),
                   directory="./",
                   lastMcmc=new.env(),
                   posteriorAvg=NULL,
                   bayesianEffectSize=NULL){
            .Object@studyNames <- studyNames
            .Object@featureNames <- featureNames
            .Object@iterations <- iterations
            .Object@seed <- seed
            .Object@output <- output
            .Object@directory <- directory
            .Object@lastMcmc <- lastMcmc
            .Object@posteriorAvg <- posteriorAvg
            .Object@bayesianEffectSize <- bayesianEffectSize
            .Object
    })

setMethod("$", "XdeMcmc", function(x, name) {
	if(!(name %in% .parameterNames()[output(x) != 0])){
		stop(paste("Parameter", name, "not saved.  See output for XdeParameter object"))
	}
	mcmc <- scan(paste(directory(x), "/", name, ".log", sep = ""))

	##---------------------------------------------------------------------------
	##Parameters indexed by gene and study
	##
	## If ITER iterations were executed, there are ITER lines written to the file. 
	## Each line is a big vector in the following order:
	## gene1 study 1, gene 1 study2, ... , gene1 study P, gene 2 study 1, ...
	## Or,
	##iter1:     nu_11, nu_12, ..., nu_1P, nu_21, ..., nu_2P, nu_G1, ..., nu_GP
	##iter2:     nu_11, nu_12, ..., nu_1P, nu_21, ..., nu_2P, nu_G1, ..., nu_GP
	##
	## for the nu parameter.
	##
	## For reading in arrays, the 'left subscript moves the fastest'
	##
	## scan() reads in the data by rows.
	## Therefore, the fastest moving index is study, followed by genes and then the mcmc iteration. 
	##
	##  mcmc <- array(mcmc, dim=c(# studies, # genes, # iterations)) is equiv to
	## For iteration 1:
	## mcmc[1, 1, 1] <- nu_11
	## mcmc[2, 1, 1] <- nu_12
	## mcmc[P, 1, 1] <- nu_1P
	## mcmc[1, 2, 1] <- nu_21
	## mcmc[2, 2, 1] <- nu_22
	## ...
	## mcmc[P, 2, 1] <- nu_2P
	## ...
	## Example:
	##	tmp <- outer(1:10, 1:3, FUN="*")
	##	tmp2 <- as.numeric(t(tmp))
	##	tmp2 <- rbind(tmp2, tmp2)
	##      ##rows are iterations
	##	write.table(tmp2, file="tmp2.txt", sep=" ", quote=FALSE, row.names=FALSE,col.names=FALSE)
	##	x=scan("tmp2.txt")
	##	y=array(x, dim=c(3, 10,  2))
	##	y <- aperm(y)
	##---------------------------------------------------------------------------	
	if(name %in% c("DDelta", "nu", "phi", "sigma2", "delta", "probDelta")){
		I <- iterations(x)
		mcmc <- array(mcmc,
			      dim=c(length(studyNames(x)),
			      nrow(x),
			      iterations(x)),
			      dimnames=list(studyNames(x),
			      featureNames(x),
			      as.character(1:I)))
		##Now transform this to an array with dimensions:
		## iteration x features x study
		mcmc <- aperm(mcmc)
		return(mcmc)
	}
	##--------------------------------------------------

	##Parameters indexed by study
	if(name %in% c("a", "b", "lambda", "l", "tau2", "t", "theta", "xi", "tau2R", "tau2Rho")){
		nr <- iterations(x)
		nc <- length(studyNames(x))
		if(length(mcmc) < nr*nc) byrow <- FALSE else byrow <- TRUE
		mcmc <- matrix(mcmc, nrow=iterations(x), ncol=length(studyNames(x)), byrow=byrow)
		colnames(mcmc) <- studyNames(x)
		return(mcmc)    
	}

  ##################################################
	##Parameters not indexed by gene or study
	if(name %in% c("c2", "gamma2")){
		mcmc <- matrix(mcmc, ncol=1, byrow=TRUE)
		if(name == "c2"){
			##check for xi-inflation
			if(median(mcmc, na.rm=TRUE) < 0.01){
				warning("c2 values are small.  The estimated proportion of differentially expressed genes in a study (given by the parameter xi) is likely inflated.  See the discussion in the XDE manuscript for additional details")
			}
		}
		return(mcmc)    
	}

  ##################################################  
	##Correlation parameters
	if(name %in% c("r", "rho")){
		S <- length(studyNames(x))
		mcmc <- matrix(mcmc, ncol=S*(S-1)/2, byrow=TRUE)
		return(mcmc)    
	}

	##---------------------------------------------------------------------------
	##Diagnostics
	if(name == "acceptance"){
		mcmc <- matrix(mcmc, ncol = 18, nrow = iterations(x), byrow = TRUE)
		##Drop parameters updated by gibbs
		##mcmc <- mcmc[, -c(1, 2, 5, 6, 10)]
		colnames(mcmc) <- c("nu", "Delta", "a", "b", "c2", "gamma2",
				    "r", "rho", "delta", "xi", "sigma2",
				    "t", "l", "phi", "theta", "lambda", "tauR", "tau2Rho")
		return(mcmc)
	}
	if(name == "potential"){
		mcmc <- matrix(mcmc, ncol=20, byrow=TRUE)
		colnames(mcmc) <- c("posterior", "likelihood_potential",
				    .parameterNames()[4:21])
		return(mcmc)
	}

	if(name == "diffExpressed"){
		mcmc <- matrix(mcmc, ncol=3, byrow=FALSE)
		colnames(mcmc) <- c("diffExpressed", "concordant", "discordant")
		return(mcmc)
	}
  
	stop(paste("Supplied name not found.  Should be one of", .parameterNames()[2:22]))
	##  eval(substitute(pData(x)$NAME_ARG, list(NAME_ARG=name)))
})

##setMethod("[", "XdeMcmc", function(x, i, j, ..., drop = FALSE){
##	if (missing(drop)) drop <- FALSE
##	if (missing(i) && missing(j))
##	{
##		if (length(list(...))!=0)
##			stop("specify genes or samples to subset; use '",
##			     substitute(x), "$", names(list(...))[[1]],
##			     "' to access phenoData variables")
##		return(x)
##	}
##	if (!missing(j)){
##		f <- function(x, j){
##			x <- x[, j]
##		}
##		x <- lapply(x, f, j)
##	}
##	if(!missing(i)){
##		f <- function(x, i){
##			x <- x[i, ]
##		}
##		x <- lapply(x, f, i)
##	}
##	as(x, "ExpressionSetList")
##})

setMethod("bayesianEffectSize", "XdeMcmc", function(object) object@bayesianEffectSize)

setReplaceMethod("bayesianEffectSize", c("XdeMcmc", "matrix"),
                 function(object, value){
			 object@bayesianEffectSize <- value
			 object
		 })

setMethod("calculateBayesianEffectSize", "XdeMcmc",
          function(object){
		  D <- .standardizedDelta(object)
		  averageDelta <- apply(D, c(2, 3), "mean")
		  averageDelta
          })

setMethod("directory", "XdeMcmc", function(object) object@directory)
setMethod("featureNames", "XdeMcmc", function(object) object@featureNames)
setMethod("iterations", "XdeMcmc", function(object) object@iterations)
setMethod("lastMcmc", "XdeMcmc", function(object) object@lastMcmc)
setMethod("nrow", "XdeMcmc", function(x) length(featureNames(x)))
setMethod("output", "XdeMcmc", function(object) object@output)

##setMethod("calculatePosteriorAvg", "XdeMcmc", function(object, NCONC=2, NDIFF=2){

setMethod("posteriorAvg", "XdeMcmc", function(object) object@posteriorAvg)
setReplaceMethod("posteriorAvg", c("XdeMcmc", "matrix"),
		 function(object, value){
			 object@posteriorAvg <- value
			 object
		 })

setMethod("seed", "XdeMcmc", function(object) object@seed)
setMethod("studyNames", "XdeMcmc", function(object) object@studyNames)
	  
setMethod(".standardizedDelta", "XdeMcmc",
          function(object, deltaProduct=TRUE){
		  delta <- object$delta
		  sqrt.c2 <- sqrt(object$c2)
		  sigma <- sqrt(object$sigma2)
		  if("tau2R" %in% names(output(object))){
			  separateTau <- TRUE
		  } else separateTau <- FALSE
		  if(separateTau){
			  tau <- sqrt(object$tau2R)
			  ##tauRho <- sqrt(object$tau2Rho)
		  } else{
			  tau <- sqrt(object$tau2)
		  }
		  b <- object$b
		  dDelta <- object$DDelta
		  s <- array(NA, dim=dim(sigma))
		  ##dDelta <- Delta
		  P <- dim(sigma)[3]
		  ##For each platform...
		  for(p in 1:P){
			  s[, , p] <- (sigma[, , p]^b[,p])*tau[, p]*as.numeric(sqrt.c2)
		  }
		  dDelta <- dDelta * delta
		  D <- dDelta/s
		  D
          })

setMethod("show", "XdeMcmc",
          function(object){
		  cat("Instance of", class(object), "\n")
		  cat("Study names:", studyNames(object), "\n")
		  cat("Number of features:" , length(featureNames(object)), "\n\n")
		  ##cat("phenotypeLabel", phenotypeLabel(object), "\n\n")
		  cat("seed (a new seed):", seed(object), "\n\n")            
		  cat("iterations (saved):", iterations(object), "\n\n")
		  out <- as.logical(output(object)[2:length(output(object))])
		  names(out) <- .parameterNames()[2:length(output(object))]
		  cat("output (parameters saved):\n")
		  print(out)
		  cat("\n")
		  cat("directory (of saved log files):", directory(object), "\n\n")
		  cat("lastMcmc:", class(lastMcmc(object)), "\n\n")
		  if(!is.null(posteriorAvg(object))){
			  cat("posteriorAvg:\n")
			  str(posteriorAvg(object))
			  cat("\n")
		  } else cat("posteriorAvg: NULL\n")
		  if(!is.null(bayesianEffectSize(object))){
			  cat("bayesianEffectSize:\n")
			  str(bayesianEffectSize(object))
			  cat("\n\n")
		  } else cat("bayesianEffectSize: NULL\n")
		  return()
          })

setMethod("plot", "XdeMcmc",
          function(x, y, op, verbose=FALSE, ...){
		  if(missing(y)){
			  getLogs <- function(object){
				  params <- output(object)[output(object) == 1]
				  params <- params[!(names(params) %in% c("nu", "phi", "DDelta", "delta", "sigma2", "diffExpressed"))]
				  names(params)
			  }
			  y <- getLogs(x)
		  }
		  params <- lapply(lapply(as.list(y), function(name, object) eval(substitute(object$NAME_ARG, list(NAME_ARG=name))), object=x), as.ts)
		  names(params) <- y
		  tracefxn <- function(x, name)  plot(x, plot.type="single", col=1:ncol(x), ylab=name)
		  mapply(tracefxn, params, name=names(params))			  
          })
