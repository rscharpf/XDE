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

  ##################################################
	##Parameters indexed by gene and study
	if(name %in% c("DDelta", "nu", "phi", "sigma2", "delta", "probDelta")){
		I <- iterations(x)
		mcmc <- array(mcmc,
			      dim=c(length(studyNames(x)),
			      nrow(x),
			      iterations(x)),
			      dimnames=list(studyNames(x),
			      featureNames(x),
			      as.character(1:I)))
		return(mcmc)
	}

  
	##--------------------------------------------------
	##Parameters indexed by gene 
##	if(name %in% "probDelta"){
##		mcmc <- matrix(mcmc, ncol=iterations(x))
##		rownames(mcmc) <- featureNames(x)
##		return(mcmc)
##	}

  ##################################################
	##Parameters indexed by study
	if(name %in% c("a", "b", "lambda", "l", "tau2", "t", "theta", "xi")){
		mcmc <- matrix(mcmc, nc=length(studyNames(x)), byrow=TRUE)
		colnames(mcmc) <- studyNames(x)
		return(mcmc)    
	}

  ##################################################
	##Parameters not indexed by gene or study
	if(name %in% c("c2", "gamma2")){
		mcmc <- matrix(mcmc, nc=1, byrow=TRUE)
		return(mcmc)    
	}

  ##################################################  
	##Correlation parameters
	if(name %in% c("r", "rho")){
		S <- length(studyNames(x))
		mcmc <- matrix(mcmc, nc=S*(S-1)/2, byrow=TRUE)
		return(mcmc)    
	}

  ##################################################    
  ##Diagnostics
  if(name == "acceptance"){
    mcmc <- matrix(mcmc, nc = 17, nr = iterations(x), byrow = TRUE)
    ##Drop parameters updated by gibbs
    mcmc <- mcmc[, -c(1, 2, 5, 6, 10)]
    colnames(mcmc) <- c("nu", "Delta", "a", "b", "c2", "gamma2",
                        "r", "rho", "delta", "xi", "sigma2",
                        "t", "l", "phi", "theta", "lambda", "tau")
    return(mcmc)
  }
  if(name == "potential"){
    mcmc <- matrix(mcmc, ncol=19, byrow=TRUE)
    colnames(mcmc) <- .parameterNames()[4:22]
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

setMethod("calculatePosteriorAvg", "XdeMcmc", function(object, NCONC=2, NDIFF=2){
##	postAvg <- object@posteriorAvg
##	if(nrow(postAvg) == 0){
	D <- object$DDelta
	d <- object$delta
	dD <- sign(d*D)
	f <- function(x){
		##define an indicator for concordance
		##indicator for discordance
		tmp <- t(rbind(colSums(x > 0), colSums(x < 0)))
		colnames(tmp) <- c("#up", "#down")
		discordant <- (tmp[, 1] * tmp[, 2]) != 0
		concordant <- (tmp[, 1] * tmp[, 2]) == 0  ##first require no discordant signs
		##Finally, let the number of concordant studies be user-specified
		concordant <- concordant * (rowSums(tmp) >= NCONC)
		diffExpr <- rowSums(abs(tmp)) >= NDIFF
		indicators <- matrix(c(concordant,
				       discordant,
				       diffExpr), ncol=3, byrow=FALSE)
	}
	I <- array(NA, c(dim(dD)[2], 3, dim(dD)[3]))
	dimnames(I) <- list(featureNames(object),
			    c("concordant", "discordant", "diffExpressed"),
			    paste("iterations", 1:dim(dD)[3], sep="_"))
	for(i in 1:dim(dD)[3]){
		I[, , i] <- f(dD[, , i])
	}
	concordant.avg <- rowMeans(I[, 1, ])
	discordant.avg <- rowMeans(I[, 2, ])
	diffExpressed.avg <- rowMeans(I[, 3, ])
	X <- cbind(concordant.avg, discordant.avg, diffExpressed.avg)
	colnames(X) <- c("concordant", "discordant", "diffExpressed")
	rownames(X) <- featureNames(object)
	X
})
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
		  ##C <- sqrt(c2Log(object))
		  sigma <- sqrt(object$sigma2)
		  tau <- sqrt(object$tau2)
		  b <- object$b
		  dDelta <- object$DDelta
		  b <- t(b)
		  tau <- t(tau)
		  s <- array(NA, dim=dim(sigma))
		  ##dDelta <- Delta
		  P <- dim(sigma)[1]
		  ##For each platform...
		  for(i in 1:P){
			  s[i, , ] <- t(apply(as.matrix(sigma[i, , ]), 1, "^", b[i, ]))
			  s[i, , ] <- t(apply(as.matrix(s[i, , ]), 1, "*", tau[i, ]))
			  s[i, , ] <- t(apply(as.matrix(s[i, , ]), 1, "*", sqrt.c2))
##			  if(deltaProduct) dDelta[i, ,] <- as.matrix(dDelta[i, , ]) * delta
			  if(deltaProduct) dDelta <- dDelta * delta
		  }
		  D <- dDelta/s
		  D <- aperm(D)
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

##See help on plot.mcmc
##... specificies arguments to par
setMethod("plot", "XdeMcmc",
          function(x, y, op, verbose=FALSE, ...){
            require(coda) || stop("Method plot for class XdeMcmc requires coda package")
            def.par <- par(no.readonly=TRUE)
            on.exit(def.par)
            
            if(missing(op)){
              par(mfrow=c(2, 6), las=1, mar=c(0.5, 3, 2, 1), oma=c(3, 1, 1, 2))
            } else {
              par(mfrow=op$mfrow, las=1, mar=op$mar, oma=op$oma)
            }
            plotFxn <- function(x, xaxt="n", ylim=c(0, 1), ...){
              x <- as.mcmc(x)
              for(i in 1:ncol(x)){
                if(i == 1){
                  coda:::plot.mcmc(x[, i], density=FALSE, col=i,
                                   smooth=FALSE, ylim=ylim,
                                   auto.layout=FALSE, xaxt=xaxt, ...)
                } else {
                  coda:::plot.mcmc(x[, i], density=FALSE, col=i, add=TRUE,
                                   main=NA, smooth=FALSE, auto.layout=FALSE, ylim=ylim, ...)
                }
              }
            }

            if(output(x)["potential"] == 1){
              pot <- x$potential[, 2, drop=FALSE]
              ylim <- range(pot)
              ##            ylim[2] <- quantile(pot, probs=0.995)
              plotFxn(pot, main="potential", ylim=ylim, ...)
            } else {
              if(verbose) print("potential not saved")
            }

            if(output(x)["a"] == 1){
              a <- x$a
              plotFxn(a, main="a", ...)
            } else{
              if(verbose) print("a not saved")
            }

            if(output(x)["b"] == 1){
              b <- x$b
              plotFxn(b, main="b", ...)
            }

            if(output(x)["l"] == 1){
              l <- x$l
              plotFxn(l, main="l", ...)
            }

            if(output(x)["t"] == 1){
              t <- x$t
              plotFxn(t, main="t", ...)
            }

            if(output(x)["gamma2"] == 1){
              gamma2 <- as.mcmc(x$gamma2)
              coda:::plot.mcmc(gamma2, density=FALSE, col=1, main=expression(gamma^2), smooth=FALSE, auto.layout=FALSE, xaxt="s", ...)              
            }

            if(output(x)["c2"] == 1){            
              c2 <- as.mcmc(x$c2)
              coda:::plot.mcmc(c2, density=FALSE, col=1, main=expression(c^2), smooth=FALSE, auto.layout=FALSE, xaxt="s", ...)
            }

            if(output(x)["tau2"] == 1){                        
              tau2 <- x$tau2
              plotFxn(tau2, main=expression(tau^2), xaxt="s", ylim=range(tau2), ...)
            }

            if(output(x)["xi"] == 1){                        
		    xi <- as.mcmc(x$xi)
		    plotFxn(xi, main=expression(xi), xaxt="s", ylim=c(0, 1), ...)
##		    coda:::plot.mcmc(xi, density=FALSE, col=1, main=expression(xi), smooth=FALSE, auto.layout=FALSE, ...)
            }

            if(output(x)["rho"] == 1){                        
              rho <- x$rho
              plotFxn(rho, main=expression(rho), xaxt="s", ...)
            }

            if(output(x)["r"] == 1){                        
              r <- x$r
              plotFxn(r, main=expression(r), xaxt="s", ...)
            }
          })
