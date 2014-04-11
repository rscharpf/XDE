##setValidity("Params", function(object){
##	nms <- names(object)
##	p.nms <- parameterNames()
##	if(!all(p.nms %in% nms)){
##		message("Missing one or more parameters required for the class:")
##		return(setdiff(p.nms, nms))
##	}
##	TRUE
##})
## initialize an instance of the class Params
##setMethod("Params", signature(object="HyperParams"),
Parameters <- function(object,
		       G=nrow(object),
		       QQ=length(object),
		       S=sapply(object, ncol),
		       clinicalCovariate,
		       psi=phenotype(object, clinicalCovariate),
		       one.delta,
		       MRF=FALSE,
		   alphaA=1.0,
		   betaA=1.0,
		   alphaB=1.0,
		   betaB=1.0,
		   pA0=0.1,
		   pA1=0.1,
		   pB0=0.1,
		   pB1=0.1,
		   nuR=1.0+QQ,
		   nuRho=1.0+QQ,
		   alphaXi=1.0,
		   betaXi=1.0,
		   c2Max=10.0,
		   alphaEta=1.0,
		   betaEta=1.0,
		   pOmega0=0.1,
		   lambdaOmega=1.0,
		   lambdaKappa=1.0,
		   seed=123L,
		   ...){
	if(missing(one.delta)) stop("one.delta should be logical, but is missing")
		  ##require(mvtnorm)
##		  if(!missing(object)){
##			  stopifnot(is(object, "ExpressionSetList"))
##			  Gs <- sapply(object, nrow)
##			  stopifnot(length(unique(Gs)) == 1)
##			  G <- Gs[[1]]
##			  QQ <- length(object)
##			  S <- sapply(object, ncol)
##			  psi <- phenotype(object, clinicalCovariate)
##			  stopifnot(length(vectorize(psi)) == sum(S))
##		  } else {
##		  if(missing(object)){
##			  anymissing <- missing(G) | missing(QQ) | missing(S)
##			  stopifnot(!anymissing)
##			  psi <- NULL
##		  }
	stopifnot(!missing(object))
	psi <- vectorize(psi)
	xi <- rep(0.5, QQ)
	gamma2 <- 0.5*0.5
	c2 <- 0.5*0.5
	tau2Rho <- rep(1.0, QQ)
	tau2R <- rep(1.0,QQ)
	a <- b <- rep(NA, QQ)
	l <- rep(1.0, QQ)
	t <- rep(0.5^2, QQ)
	sigma2 <- matrix(NA, G, QQ)
	## simulate variances from gamma
	param2 <- l/t
	param1 <- l*param2
	for(q in seq(length=QQ)){
		u <- runif(1)
		if(u < pA0){
			a[q] <- 0
		} else {
			if(u < pA0+pA1){
				a[q] <- 1
			} else {
				a[q] <-runif(1)
			}
		}
		u <- runif(1)
		if(u < pB0){
			b[q] <- 0
		} else {
			if(u < pB0+pB1){
				b[q] <- 1
			} else {
				b[q] <-runif(1)
			}
		}
		## in R, mean is shape * scale and variance is shape*scale^2
		## in Haakon's formulation, mean = l and variance = t
		## => shape=l^2/t and scale =t/l
		##  (this means that param2 is the rate and param1 is the shape
		sigma2[, q] <- rgamma(G, shape=param1[q], rate=param2[q])
	}
	## the platform index changes the fastest, then gene
	##sigma2 <- as.numeric(t(sigma2))
	##
	## initialize phi
	lambda <- rep(1.0, QQ)
	theta <- rep(0.1*0.1, QQ)
	param2 <- lambda/theta
	param1 <- lambda*theta
	phi <- matrix(NA, G, QQ)
	for(q in seq(length=QQ)){
		phi[, q] <- rgamma(G, shape=param1[q], rate=param2[q])
	}
	##r <- rho <- rep(0.5, QQ*(QQ-1)/2)
	r <- rho <- matrix(NA, QQ, QQ)
	diag(r) <- diag(rho) <- 1
	rho[upper.tri(rho)] <- r[upper.tri(r)] <- 0.5
	rho[lower.tri(rho)] <- r[lower.tri(r)] <- 0.5
	## nu
	nu <- matrix(NA, G, QQ)
	for(g in seq(length=G)){
		Sig <- makeSigma(Q=QQ, gamma2=gamma2, tau2=tau2Rho, a=a,
				 sigma2=sigma2[g, ], r=rho)
		## simulate from multivariate normal
		nu[g, ] <- rmvnorm(1, sigma=Sig)
	}
	Delta <- matrix(NA, G, QQ)
	for(g in seq(length=G)){
		R <- makeSigma(Q=QQ, gamma2=c2, tau2=tau2Rho, a=b,
			       sigma2=sigma2[g, ], r=r)
		## simulate from multivariate normal
		Delta[g, ] <- rmvnorm(1, sigma=Sig)
	}
	##---------------------------------------------------------------------------
	##
	## Simulating delta
	##
	## There are 4 prior models for delta:
	##    model A: one delta
	##    model B: study-specific delta
	##    model C: MRF study-specific delta
	##    model D: MRF one delta
	##
	##---------------------------------------------------------------------------
	## Delta in the Rinterface.cpp file there are different
	## functions relating to each of the now four prior models for
	##
	## delta
	##
	## model A:
	##         updateXi_onedelta
	##         updateDeltaDDelta_onedelta
	##
	## model B:
	##         updateXi
	##         updateDeltaDDelta
	##
	## model C:
	##         updateAlpha_MRF
	##         updateBeta_MRF
	##         updateBetag_MRF
	##         updateDeltaDDelta_MRF
	##
	##
	## model D:
	##         updateAlpha_MRF_onedelta
	##         updateBeta_MRF_onedelta
	##         updateDeltaDDelta_MRF_onedelta
	##
	alpha.mrf <- -0.2  ## in (-Inf, +Inf)
	beta.mrf <- 3.0    ## > 0
	if(one.delta){
		delta <- rbinom(G, 1, 0.2)
		## This results in Cholesky problems
		##Delta[delta==0, ] <- 0
	} else{
		delta <- matrix(rbinom(G*QQ, 1, 0.2), G, QQ)
	}
	##	if(!MRF){
	##		if(one.delta){
	##			## model A
	##			##delta <- rep(0L, G)
	##			delta <- rbinom(G, 1, 0.2)
	##		} else {
	##			## model B (default)
	##			delta <- matrix(0L, G, QQ)
	##		}
	##	} else {
	##		## 1. specify alpha.mrf, beta.mrf
	##		## 2. specify |N_g|
	##		## 3. sample |N_g| indices from {1, ..., G} to form N_g for each g
	##		## 4. simulate delta_g or delta_gp from Bernoulli(prob), with
	##		##    prob = exp (...)/(1+exp(...))
	##		##    ... =
	##		##
	##		if(one.delta){
	##			delta <- rbinom(G, 1, 0.2)
	##			## model D
	##			## a subset of genes are differentially expressed
	##			##    the probability of a gene's differential expression
	##			##    depends on other genes in the network and is assumed
	##			##    to be the same across platforms.
	##			## Let's assume 20% of the genes are differentially expressed
	##
	##			##
	##			##
	##		} else {
	##			## model C
	##			delta <- matrix(NA, G, Q)
	##		}
	##	}
	##		  hyper.params <- data.frame(
	rho <- rho[upper.tri(rho)]
	r <- r[upper.tri(r)]
	new("Parameters",
	    seed=seed,
	    data=expressionVector(object),
	    phenodata=psi,
	    Q=as.integer(QQ),
	    G=as.integer(G),
	    S=as.integer(S),
	    alphaA=alphaA,
	    alphaB=alphaB,
	    betaA=betaA,
	    betaB=betaB,
	    pA0=pA0,
	    pA1=pA1,
	    pB0=pB0,
	    pB1=pB1,
	    nuR=nuR,
	    nuRho=nuRho,
	    alphaXi=alphaXi,
	    betaXi=betaXi,
	    c2Max=c2Max,
	    alphaEta=alphaEta,
	    betaEta=betaEta,
	    pOmega0=pOmega0,
	    lambdaOmega=lambdaOmega,
	    lambdaKappa=lambdaKappa,
	    gamma2=gamma2,
	    c2=c2,
	    tau2Rho=tau2Rho,
	    tau2R=tau2R,
	    a=a,
	    b=b,
	    l=l,
	    t=t,
	    lambda=lambda,
	    theta=theta,
	    phi=vectorize(phi),
	    sigma2=vectorize(sigma2),
	    r=r,
	    rho=rho,
	    nu=vectorize(nu),
	    delta=vectorize(delta),
	    Delta=vectorize(Delta),
	    xi=vectorize(xi))
#		  res <- as(res, "Params")
##		  return(res)
	  }

setMethod("show", signature(object="Parameters"),
	  function(object){
		  cat("class Parameters: elements of class include\n")
		  nms <- slotNames(object)
		  L <- rep(NA, length(nms))
		  headers <- vector("list", length(nms))
		  for(i in seq_len(length(L))){
			  obj <- eval(substitute(object@NAME_ARG, list(NAME_ARG=nms[i])))
			  L[i] <- length(obj)
			  ##headers <- head(obj)
		  }
		  cls <- sapply(nms, function(x, object){
			  class(eval(substitute(object@NAME_ARG, list(NAME_ARG=x))))
		  }, object=object)
		  size <- paste("[1:", L, "]", sep="")
		  res <- data.frame(nms=nms, cls=substr(cls, 1,3), size=size)
		  row.names(res) <- NULL
		  colnames(res) <- c("slots", "type", "size")
		  print(head(res))
		  cat("\n...\n")
		  colnames(res)<-
		  print(tail(res))
	  })

setMethod("exprs", signature(object="Parameters"),
	  function(object) object@data)

setMethod("vectorize", signature(object="integer"), function(object) return(object))
setMethod("vectorize", signature(object="numeric"), function(object) return(object))
##setMethod("vectorize", signature(object="Params"), function(object){
##	## gene index changes fastest, then study index
##	object$nu <- vectorize(object$nu)
##	object$delta <- vectorize(object$delta)
##	object$Delta <- vectorize(object$Delta)
##	object$phi <- vectorize(object$phi)
##	if(is(object$rho, "matrix")){
##		rho <- object$rho
##		rho <- rho[upper.tri(rho)]
##		object$rho <- rho
##	}
##	object$sigma2 <- vectorize(object$sigma2)
##	return(object)
##})
##setMethod("show", signature(object="Params"), function(object) print(ls.str(object)))
setAs("XdeParameter", "Parameters", function(from){
	x <- firstMcmc(from)
	S <- length(x[["Xi"]])
	G <- length(x[["Nu"]])/S
	res <- list(gamma2=x[["Gamma2"]],
		    c2=x[["C2"]],
		    tau2Rho=x[["Tau2Rho"]],
		    tau2R=x[["Tau2R"]],
		    a=x[["A"]],
		    b=x[["B"]],
		    l=x[["L"]],
		    t=x[["T"]],
		    lambda=x[["Lambda"]],
		    theta=x[["Theta"]],
		    phi=matrix(x[["Phi"]], G, S),
		    sigma2=matrix(x[["Sigma2"]],G,S),
		    r=x[["R"]],
		    rho=x[["Rho"]],
		    delta=matrix(x[["Delta"]], G, S),
		    nu=matrix(x[["Nu"]], G, S),
		    Delta=matrix(x[["DDelta"]], G, S),
		    xi=x[["Xi"]])
	res <- as(res, "Params")
	return(res)
})
setMethod("$", "Parameters", function(x, name){
	slot(x, name)
	##eval(substitute(x@NAME_ARG, list(NAME_ARG=name)))
})

setReplaceMethod("$", "Parameters", function(x, name, value){
	slot(x, name) <- value
	##eval(substitute(x@NAME_ARG, list(NAME_ARG=name)))
})

setMethod("pnames", "Parameters", function(object){
	slotNames(object)
})

setMethod("[[", "Parameters", function(x,i,j,...){
	if(!is.character(i)) return(x)
	if(length(i) > 1) stop("invalid subsetting")
	slot(x, i)
})

setReplaceMethod("[[", "Parameters", function(x,i,j,...,value){
	if(!is.character(i)) return(x)
	if(length(i) > 1) stop("invalid subsetting")
	slot(x, i) <- value
	x
})

setMethod("dims", "Parameters", function(object){
	nms <- slotNames(object)
	L <- length(nms)
	res <- rep(NA, L)
	for(i in seq_len(L)) res[i] <- length(slot(object, nms[i]))
	names(res) <- nms
	return(res)
})

startingValues <- function(params, expressionSetList, phenotypeLabel, one.delta){
	if(missing(one.delta)) stop("one.delta should be logical but is missing")
	params.xde <- firstMcmc(new("XdeParameter", esetList=expressionSetList,
				    phenotypeLabel=phenotypeLabel, one.delta=one.delta))
	nms <- names(params.xde)
	firstletter <- tolower(substr(nms, 1, 1))
	nms2 <- paste(firstletter, substr(nms,2, nchar(nms)),sep="")
	nms2[nms2=="dDelta"] <- "Delta"
	all(nms2%in% slotNames(params))
	names(params.xde) <- nms2
	for(i in seq_along(nms2)){
		nm <- nms2[i]
		params[[nm]] <- params.xde[[nm]]
	}
	return(params)
}



          ##nipun: for ParametersD

          ParametersD <- function(object,...){
            tmp <- Parameters(object,...)
            object <- as(tmp, "ParametersD")
            object@alpha <- 1.0
            object@beta <- 1.0

          }


          ##nipun: for ParametersC

          ParametersC <- function(object,...){
            tmp <- ParametersD(object,...)
            object <- as(tmp, "ParametersC")
            object@betag <- 1.0
          }



          ##nipun: for ParametersMII

          ParametersMII <- function(object,...){
            tmp <- Parameters(object,...)
            object <- as(tmp, "ParametersMII")
            object@Omega <- 1.0
          }



