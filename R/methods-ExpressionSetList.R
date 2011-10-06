setMethod("dim", "ExpressionSetList", function(x) sapply(x, dim))

setMethod("rowttests", "ExpressionSetList", function(x, fac, tstatOnly=FALSE){
	sapply(x, function(x) rowttests(x, fac, tstatOnly)$statistic)
})

setMethod("lapply", "ExpressionSetList", function(X, FUN, ...){
	  X <- lapply(as(X, "list"), FUN, ...)
	  if(all(sapply(X, class) == "ExpressionSet")){
		  X <- as(X, "ExpressionSetList")
	  }
	  return(X)
	  })

setMethod("[", "ExpressionSetList", function(x, i, j, ..., drop = FALSE){
            if (missing(drop)) drop <- FALSE
            if (missing(i) && missing(j))
              {
                 if (length(list(...))!=0)
                   stop("specify genes or samples to subset; use '",
                        substitute(x), "$", names(list(...))[[1]],
                        "' to access phenoData variables")
                 return(x)
               }
            if (!missing(j)){
              f1 <- function(x, j){
                x <- x[, j]
              }
              x <- lapply(x, f1, j)
            }
            if(!missing(i)){
              f2 <- function(x, i){
                x <- x[i, ]
              }
              x <- lapply(x, f2, i)
            }
            as(x, "ExpressionSetList")
          })

setAs("list", "ExpressionSetList",
      function(from, to){
        new("ExpressionSetList", from)
      })

setMethod("featureNames", "ExpressionSetList", function(object) featureNames(object[[1]]))
setMethod("nrow", "ExpressionSetList", function(x) nrow(x[[1]]))
setMethod("nSamples", "ExpressionSetList", function(x)  unlist(lapply(x, function(x) ncol(x))))

setMethod(".pca", "ExpressionSetList",
          function(object, x, ...){
		  P <- length(object)
		  N <- nSamples(object)
		  G <- nrow(object)
		  principal.components <- function(x, P, N, G){
			  ##1. obtain the effect sizes in the P studies
			  ##   - these are stored in 'x'

			  ##2. center the effect sizes by their means
			  meancen <- colMeans(x)
			  xcen <- sweep(x, 2, meancen)

			  ##3. Fit principal components to w1, w2, and w3 using
			  ##   covariance, saving a1, a2 and a3, the loadings of the first
			  ##   principal component.
			  fit <- princomp(xcen, cor = FALSE)
			  a <- fit$loadings[1:P, 1]

			  ##4. Estimate the combined effect: (see p.9 of Garrett-Mayer
			  ##   paper)
			  weight <- sum(a * sqrt(N))
			  combined.effect <- matrix(NA, ncol=P, nrow=G)
			  for(i in 1:P){
				  combined.effect[,i] <- a[i] * sqrt(N[i]) * xcen[, i]
			  }
			  score <- apply(combined.effect, 1, sum)
			  score <- score/weight
			  return(score)
		  }
		  score <- principal.components(x=x, P=P, N=N, G=G)
		  return(score)
          })

setMethod("pData", "ExpressionSetList",
          function(object){
            f <- function(x) pData(x)
            lapply(object, f)
          })

setMethod("standardizeSamples", "ExpressionSetList",
          function(object, ...){
            standardizeColumns <- function(object){
              colSds <- function(x) rowSds(t(x))
              sds <- colSds(exprs(object))
              exprs(object) <- t(apply(exprs(object), 1, "/", sds))
              mns <- colMeans(exprs(object))
              exprs(object) <- sweep(exprs(object), 2, mns)
              object
            }
            object <- lapply(object, standardizeColumns)
            object
          })

setMethod("geneCenter", "ExpressionSetList",
          function(object){
		  mean.center <- function(x){
			  ##xx.c <- rowMeans(exprs(x))
			  exprs(x) <- exprs(x) - rowMeans(exprs(x))
			  ##exprs(x) <- sweep(exprs(x), 1, xx.c)
			  return(x)
		  }
		  object <- lapply(object, mean.center)
		  object
	  })

setMethod("studyCenter", "ExpressionSetList",
          function(object){
            mean.center <- function(x){
              overallMean <- mean(exprs(x))
              exprs(x) <- exprs(x) - overallMean
              return(x)
            }
            object <- lapply(object, mean.center)
            as(object, "ExpressionSetList")
          })

setMethod("zeroNu", "ExpressionSetList",
          function(object, phenotypeLabel, one.delta=FALSE, ...){
		  if(missing(phenotypeLabel)) stop("must specify phenotypeLabel (character string)")
		  object <- geneCenter(object)

		  ##Note that by zeroing the nu's we must specify initial
		  ##values and can not draw random samples from the prior
		  ##(could result in illegal values)
		  ##            params <- new("XdeParameter", phenotypeLabel=phenotypeLabel, expressionSetList=object)
		  params <- new("XdeParameter", esetList=object, phenotypeLabel=phenotypeLabel,
				one.delta=one.delta)
		  firstMcmc <- empiricalStart(object, zeroNu=TRUE, phenotypeLabel=phenotypeLabel, one.delta=one.delta)
		  firstMcmc$A <- rep(0, length(object))
		  firstMcmc$Rho <- rep(0, choose(length(object), 2))
		  firstMcmc$Tau2Rho <- rep(1, length(object))
		  up <- updates(params)
		  up["nu"] <- 0
		  up["a"] <- 0
		  up["gamma2"] <- 0
		  up["rhoAndGamma2"] <- 0
		  up["tau2Rho"] <- 0
		  updates(params) <- up
		  firstMcmc(params) <- firstMcmc
		  params
          })

setMethod("goodnessOfFit", c("ExpressionSetList", "XdeMcmc"),
	  function(object,
		   fit,
		   phenotype,
		   firstIteration=1,
		   lastIteration=iterations(fit),
		   by=2){
		  stopifnot(phenotype %in% varLabels(object[[1]]))
		  psi <- lapply(object, function(x) pData(x)[, match(phenotype, varLabels(object))])
		  nus <- fit$nu
		  Delta <- fit$DDelta
		  delta <- fit$delta
		  sigma2 <- fit$sigma2
		  phi <- fit$phi
		  computeGOF(object,
			     nus=nus,
			     Delta=Delta,
			     delta=delta,
			     sigma2=sigma2,
			     phi=phi,
			     firstIteration=firstIteration,
			     lastIteration=lastIteration,
			     by=by)
	  })


##setMethod("getHyperparameters", signature(object="ExpressionSetList"),


##setMethod("getHyperparameters", signature(object="NULL"),
##	  function(object, G, Q, S, ...){
##		  anymissing <- missing(G) | missing(Q) | missing(S)
##		  stopifnot(!anymissing)
##
##		  betaA <- alphaA <- 1.0
##		  pA0 <- pA1 <- 0.1
##
##		  ##
##		  betaB <- alphaB <- 1.0
##		  pA0 <- pB0 <- 0.1
##
##		  nuR <- 1.0+Q
##		  nuRho <- 1.0+Q
##
##		  betaXi <- alphaXi <- 1.0
##
##		  c2Max <- 10.0
##
##		  betaEta <- alphaEta <- 1.0
##
##		  pOmega0 <- 0.1
##		  lambdaOmega <- 1.0
##		  lambdaKappa <- 1.0
##
##		  res <- list(G=G,
##			      Q=Q,
##			      S=S,
##			      alphaA=alphaA,
##			      betaA=betaA,
##			      pA0=pA0,
##			      pA1=pA1,
##			      nuR=nuR,
##			      nuRho=nuRho,
##			      alphaXi=alphaXi,
##			      betaXi=betaXi,
##			      c2Max=c2Max,
##			      alphaEta=alphaEta,
##			      betaEta=betaEta,
##			      pOmega0=pOmega0,
##			      lambdaOmega=lambdaOmega,
##			      lambdaKappa=lambdaKappa)
##		  return(res)
##	  })
