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

setMethod("show", signature(object="HyperParams"),
	  function(object){
		  cat("class HyperParams: elements of class include\n")
		  print(ls.str(object))
	  })

setMethod("Params", signature(object="HyperParams"),
	  function(object, one.delta=FALSE, MRF=FALSE){
		  ##require(mvtnorm)
		  G <- object$G
		  Q <- object$Q
		  pA0 <- object$pA0
		  pA1 <- object$pA1
		  pB0 <- object$pB0
		  pB1 <- object$pB1
		  alphaXi <- object$alphaXi
		  betaXi <- object$betaXi

		  gamma2 <- 0.5*0.5
		  c2 <- 0.5*0.5
		  tau2Rho <- rep(1.0, Q)
		  tau2R <- rep(1.0,Q)
		  a <- b <- rep(NA, Q)
		  l <- rep(1.0, Q)
		  t <- rep(0.5^2, Q)
		  sigma2 <- matrix(NA, G, Q)
		  ## simulate variances from gamma
		  param2 <- l/t
		  param1 <- l*param2
		  for(q in seq(length=Q)){
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
		  lambda <- rep(1.0, Q)
		  theta <- rep(0.1*0.1, Q)
		  param2 <- lambda/theta
		  param1 <- lambda*theta
		  phi <- matrix(NA, G, Q)
		  for(q in seq(length=Q)){
			  phi[, q] <- rgamma(G, shape=param1[q], rate=param2[q])
		  }
		  ##r <- rho <- rep(0.5, Q*(Q-1)/2)
		  r <- rho <- matrix(NA, Q, Q)
		  diag(r) <- diag(rho) <- 1
		  rho[upper.tri(rho)] <- r[upper.tri(r)] <- 0.5
		  rho[lower.tri(rho)] <- r[lower.tri(r)] <- 0.5
		  ## nu
		  nu <- matrix(NA, G, Q)
		  for(g in seq(length=G)){
			  Sig <- makeSigma(Q=Q, gamma2=gamma2, tau2=tau2Rho, a=a,
					   sigma2=sigma2[g, ], r=rho)
			  ## simulate from multivariate normal
			  nu[g, ] <- rmvnorm(1, sigma=Sig)
		  }
		  Delta <- matrix(NA, G, Q)
		  for(g in seq(length=G)){
			  R <- makeSigma(Q=Q, gamma2=c2, tau2=tau2Rho, a=b,
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
			  Delta[delta==0, ] <- 0
		  } else{
			  delta <- matrix(rbinom(G*Q, 1, 0.2), G, Q)
			  Delta <- Delta*delta
		  }
		  ##	if(!MRF){
		  ##		if(one.delta){
		  ##			## model A
		  ##			##delta <- rep(0L, G)
		  ##			delta <- rbinom(G, 1, 0.2)
		  ##		} else {
		  ##			## model B (default)
		  ##			delta <- matrix(0L, G, Q)
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
		  res <- list(
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
			      phi=phi,
			      sigma2=sigma2,
			      r=r,
			      rho=rho,
			      delta=delta,
			      nu=nu,
			      Delta=Delta)
		  res <- as(res, "Params")
		  return(res)
	  })

setMethod("vectorize", signature(object="HyperParams"),
	  function(object){
		  object$S <- as.integer(object$S)
		  object$psi <- vectorize(object$psi)
		  return(object)
	  })
