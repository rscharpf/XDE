HyperParams <- function(object, G, Q, S, seed, clinicalCovariate, ...){
	if(!missing(object)){
		stopifnot(is(object, "ExpressionSetList"))
		Gs <- sapply(object, nrow)
		stopifnot(length(unique(Gs)) == 1)
		G <- Gs[[1]]
		Q <- length(object)
		S <- sapply(object, ncol)
		psi <- phenotype(object, clinicalCovariate)
		stopifnot(length(vectorize(psi)) == sum(S))
	} else {
		anymissing <- missing(G) | missing(Q) | missing(S)
		stopifnot(!anymissing)
		psi <- NULL
	}
	if(missing(seed)) seed <- 123L

	betaA <- alphaA <- 1.0
	pA0 <- pA1 <- 0.1

	##
	betaB <- alphaB <- 1.0
	pB0 <- pB1 <- 0.1

	## degrees of freedom for correlation of Delta
	## p.8 Scharpf et al
	nuR <- 1.0+Q
	## degrees of freedom for correlation of Nu
	nuRho <- 1.0+Q

	betaXi <- alphaXi <- 1.0

	c2Max <- 10.0

	betaEta <- alphaEta <- 1.0

	pOmega0 <- 0.1
	lambdaOmega <- 1.0
	lambdaKappa <- 1.0

	res <- list(
		   G=G,
		    Q=Q,
		    S=S,
		    alphaA=alphaA,
		    betaA=betaA,
		    alphaB=alphaB,
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
		    seed=seed,
		    psi=psi)
	res <- as(res, "HyperParams")
	return(res)
}

##setMethod("show", signature(object="HyperParams"),
##	  function(object){
##		  cat("class HyperParams: elements of class include\n")
##		  print(ls.str(object))
##	  })



setMethod("vectorize", signature(object="HyperParams"),
	  function(object){
		  object$S <- as.integer(object$S)
		  object$psi <- as.integer(object$psi)
		  return(object)
	  })
