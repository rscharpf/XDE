##setMethod("chain.initialize", "XdeSet",
.chainInitialize <- function(object,
                             chain.length,
                             verbose){
              names(chain.length) <- c("potential",
                                        "acceptance",
                                        "nu",
                                        "deltaDelta",
                                        "a",
                                        "b",
                                        "c2",
                                        "gamma2",
                                        "r",
                                        "rho",
                                        "delta",
                                        "xi",
                                        "sigma2",
                                        "t",
                                        "l",
                                        "phi",
                                        "theta",
                                        "lambda",
                                        "tau2R",
				       "tau2Rho",
                                        "probDelta",
                                        "diffExpressed")

            ##These values are placeholders to create a chain of the
            ##right dimension.  If called, it is assumed that default
            ##starting values for the chain will be used.
            G <- nrow(object)
            QQ <- length(object)
            potential <- rep(0, 20*chain.length[["potential"]])
	      ##acceptance <- rep(0, 17*chain.length[["acceptance"]])
	      acceptance <- rep(0, 18*chain.length[["acceptance"]])
            Nu <- rep(0, G*QQ*chain.length[["nu"]])
            DDelta <- rep(0, G*QQ*chain.length[["deltaDelta"]])
            A <- rep(0, QQ*chain.length[["a"]])
            B <- rep(0, QQ*chain.length[["b"]])
            C2 <- rep(0, chain.length[["c2"]])
            Gamma2 <- rep(0, chain.length[["gamma2"]])
            Rho <- rep(0, QQ*(QQ-1)/2*chain.length[["rho"]])
            R <- rep(0, QQ*(QQ-1)/2*chain.length[["r"]])
            Delta <- rep(0, G*QQ*chain.length[["delta"]])
            Xi <- rep(0, QQ*chain.length[["xi"]])
            Sigma2 <- rep(0, QQ * G*chain.length[["sigma2"]])
            T <- rep(0, QQ*chain.length[["t"]])
            L <- rep(0, QQ*chain.length[["l"]])
            Phi <- rep(0, QQ * G*chain.length[["phi"]])
            Theta <- rep(0, QQ*chain.length[["theta"]])
            Lambda <- rep(0, QQ*chain.length[["lambda"]])
	      Tau2R <- Tau2Rho <- rep(0, QQ*chain.length[["tau2R"]])
            ProbDelta <- rep(0, G*QQ*chain.length[["probDelta"]])
            Diffexpressed <- rep(0, 3*G)
            eenv <- new.env()
            x <- list(Potential=potential,
                      Acceptance=acceptance,
                      Nu=Nu,
                      DDelta=DDelta,
                      A=A,
                      B=B,
                      C2=C2,
                      Gamma2=Gamma2,
                      R=R,
                      Rho=Rho,
                      Delta=Delta,
                      Xi=Xi,
                      Sigma2=Sigma2,
                      T=T,
                      L=L,
                      Phi=Phi,
                      Theta=Theta,
                      Lambda=Lambda,
                      Tau2R=Tau2R,
		      Tau2Rho=Tau2Rho,
                      ProbDelta=ProbDelta,
                      Diffexpressed=Diffexpressed)
              assign("vars", x, eenv)
              eenv
      }

calculatePosteriorAvg <- function(object, NCONC=2, NDIFF=1, burnin=0){
	if(class(object) != "XdeMcmc") stop("object must be of class XdeMcmc")
	if(iterations(object) < 10) stop("Too few iterations to compute posterior means")
	ii <- which((1:iterations(object)) > burnin)
	if(length(ii) < 1) stop("Not enough iterations after burnin to compute the posterior means")
	D <- object$DDelta
	d <- object$delta
	dD <- sign(d*D)
	##For each mcmc iteration, assess concordance, discordance
	myf <- function(x){
		##define an indicator for concordance
		##indicator for discordance
		tmp <- cbind(rowSums(x>0), rowSums(x < 0))
		##tmp <- t(rbind(colSums(x > 0), colSums(x < 0)))
		colnames(tmp) <- c("#up", "#down")
		discordant <- (tmp[, 1] * tmp[, 2]) != 0

		##1) no discordant signs
		concordant <- (tmp[, 1] * tmp[, 2]) == 0
		##2) let the number of concordant studies be user-specified
		concordant <- concordant * (rowSums(tmp) >= NCONC)
		diffExpr <- rowSums(abs(tmp)) >= NDIFF
		indicators <- matrix(c(concordant,
				       discordant,
				       diffExpr), ncol=3, byrow=FALSE)
	}
##	I <- array(NA, dim=dim(dD)
	I <- array(NA, c(dim(dD)[2], 3, dim(dD)[1]),
		   dimnames=list(featureNames(object),
		   c("concordant", "discordant", "diffExpressed"),
		   paste("iter", 1:dim(dD)[1], sep="_")))
	IT <- dim(dD)[1]
	for(i in 1:IT){
		I[, , i] <- myf(dD[i, , ])
	}
	concordant.avg <- rowMeans(I[, "concordant", ])
	discordant.avg <- rowMeans(I[, "discordant", ])
	diffExpressed.avg <- rowMeans(I[, "diffExpressed", ])
	X <- cbind(concordant.avg, discordant.avg, diffExpressed.avg)
	colnames(X) <- c("concordant", "discordant", "diffExpressed")
	rownames(X) <- featureNames(object)
	X
}


.parameterNames <- function(){
	c("thin",
	  "potential",
	  "acceptance",
	  "nu",
	  "DDelta",
	  "a",
	  "b",
	  "c2",
	  "gamma2",
	  "r",
	  "rho",
	  "delta",
	  "xi",
	  "sigma2",
	  "t",
	  "l",
	  "phi",
	  "theta",
	  "lambda",
	  "tau2R",
	  "tau2Rho",
	  "probDelta",
	  "diffExpressed")
}

.permutePhenotype <- function(object, seed, ...){
	if(missing(seed)){
		seed <- sample(0:100000, 1)
		set.seed(seed)
		print(paste("Seed", seed))
	}
	print(paste("permuting phenotypeLabel", phenotypeLabel(object)))
	object@data <- lapply(object, .permuteClassLabels, phenotype=phenotypeLabel(object))
	object
}

.permuteClassLabels <- function(object, phenotype, ...){
	pData(object)[, phenotype] <- permute(pData(object)[, phenotype])
	object
}

ssStatistic <- function(statistic=c("t", "sam", "z")[1],
                        phenotypeLabel,
                        esetList, ...){
	if(!(statistic == "t" | statistic == "sam" | statistic == "z")){
		stop("only t and sam study-specific statistics provided")
	}

	##-----------------------------------------------------------
	##Wrapper for Welch t-statistic
##	multtestWrapper <- function(eset, classlabel, ...){
##		require(multtest) || stop("Package multtest not available")
##		column <- grep(classlabel, colnames(pData(eset)))
##		if(length(column) == 0) { warning("Invalid classlabel.  Using the first column in phenoData"); classlabel <- pData(eset)[, 1]}
##		if(length(column) > 1)  {warning("More than 1 phenotype has the class label.  Using the first"); classlabel <- pData(eset)[, column[1]]}
##		if(length(column) == 1) classlabel <- pData(eset)[, column]
##		X <- exprs(eset)
##		stat <- mt.teststat(X=X, classlabel=classlabel, test="t",
##				    nonpara="n")
##		stat
##	}
	tt <- rowttests(esetList, phenotypeLabel)

	sam.wrapper <- function(eset, classlabel, method, returnQ, gene.names){
		require(siggenes) || stop("Package siggenes not available")
		column <- grep(classlabel, colnames(pData(eset)))
		if(length(column) == 0) { warning("Invalid classlabel.  Using the first column in phenoData"); classlabel <- pData(eset)[, 1]}
		if(length(column) > 1)  {warning("More than 1 phenotype has the class label.  Using the first"); classlabel <- pData(eset)[, column[1]]}
		if(length(column) == 1) classlabel <- pData(eset)[, column]
		X <- exprs(eset)
		stat <- sam(data=X, cl=classlabel, method=method,
			    gene.names=gene.names)
		if(returnQ) stat <- stat@q.value else stat <- stat@d
		stat
	}

	z.wrapper <- function(object, classlabel, useREM=TRUE){
		##object is a list of ExpressionSets
		require(GeneMeta) || stop("Package GeneMeta is not available")
		pDataList <- function(eset, covariateName){
			1-(pData(eset)[, grep(covariateName, colnames(pData(eset)))])
		}
		classes <- lapply(object, pDataList, classlabel)
		z <- zScores(object, classes=classes, useREM=useREM)
		##    z <- data.frame(z[, grep("zSco", colnames(z))])
		z <- z[, grep("zSco", colnames(z))]
		return(z)
	}

	##---------------------------------------------------------------------------
	##Wrapper for SAM-statistic
	stat <- switch(statistic,
##		       t=sapply(esetList, multtestWrapper, classlabel=phenotypeLabel),
		       t=tt,
		       sam=sapply(esetList, sam.wrapper, classlabel=phenotypeLabel,
		       returnQ=FALSE, method="d.stat", gene.names=featureNames(esetList), ...),
		       z=z.wrapper(object=esetList, classlabel=phenotypeLabel, useREM=TRUE))
	rownames(stat) <- featureNames(esetList)
	return(stat)
}

symbolsInteresting <- function(rankingStatistic, percentile=0.9,
                               colors=c("grey50", "royalblue"),
                               symbols=c(".", "o"),
                               size=c(3, 1),
                               background=c("white", "grey70")){
  thr <- quantile(rankingStatistic, probs=percentile)
##  avgC <- postAvg[, "concordant"]
  i <- rankingStatistic > thr
  N <- length(rankingStatistic)
  pch <- rep(symbols[1], N)
  j <- i[order(i)]
  pch[j] <- symbols[2]
  col <- rep(colors[1], N)
  col[j] <- colors[2]
  cex <- rep(size[1], N)
  cex[j] <- size[2]
  bg <- rep(background[1], N)
  bg[j] <- background[2]
  list(order=order(i), pch=pch, col=col, bg=bg, cex=cex)
}


##statistic should be the standardized Delta
xsScores <- function(statistic, N){
	if(class(statistic) != "matrix") stop("study-specific statistics must be provided as a matrix")
	if(missing(N)) stop("Must specify the sample size of each study")

	scores <- matrix(NA, nrow(statistic), ncol=3)
	if(!is.null(rownames(statistic))) rownames(scores) <- rownames(statistic)
	colnames(scores) <- c("diffExpressed", "concordant", "discordant")

	##------------------------------------------------
	##Differential expression
	scores[, 1] <- .pca(abs(statistic), N=N)

	##------------------------------------------------
	##Concordance
	scores[, 2] <- abs(.pca(statistic, N=N))

	##------------------------------------------------
	##Discordance
	tmp <- cbind(rowSums(statistic > 0), rowSums(statistic < 0))
	I <- which((tmp[, 1] * tmp[, 2]) == 0)  ##concordant signs
##	I <- abs(rowSums(sign(statistic))) < ncol(statistic)
##	I[!I] <- -1
	scores[, 3] <- .pca(abs(statistic), N)
	##penalize concordance
	scores[I, 3] <- scores[I, 3] * -1
	return(scores)
}


.outputDefault <- function(thin=1, burnin=FALSE){
##  if(!burnin) out <- rep(1, 22) else out <- rep(0, 22)
	out <- rep(1, 23)
	out[1] <- thin
	out[23] <- 1  ##by default, keep a running average of the posterior statistics attached to the object (do not write to file)
	names(out) <- .parameterNames()
	return(out)
}
##Number of platforms (P)
.hyperparametersDefault <- function(P){
  ##Setting p_a^0, p_a^1, p_b^0, and p_b^1 = 0.1 (instead of 0) allows
  ##a and b to be identical to zero and 1 (10/23/2006 HT)
  hyperparameters <- c(1.0, 1.0, 0.1, 0.1,
                       1.0, 1.0, 0.1, 0.1,
                       P+1,
                       P+1,
                       1.0, 1.0,
                       50.0)
  names(hyperparameters) <- c("alpha.a", "beta.a", "p0.a", "p1.a",
                              "alpha.b", "beta.b", "p0.b", "p1.b",
                              "nu.r",
                              "nu.rho",
                              "alpha.xi", "beta.xi",
                              "c2max")
  hyperparameters
}



##.getPar <- function(object, ...){
##  as.list(...)
##}




###########################################################################
##Number of updates per iterations set as per Haakon's 10/23/2006 e-mail
###########################################################################
.updatesDefault <- function(){
  updates <- c(1, # nu Gibbs
                          1, # Delta Gibbs
                          3, # a
                          3, # b
                          1, # c2 Gibbs
                          1, # gamma2 Gibbs
                          3, # r and c2
                          3, # rho and gamma2
                          1, # delta
                          1, # xi Gibbs
                          1, # sigma2
                          1, # t and l
                          1, # l
                          1, # phi
                          1, # theta and lambda
                          1, # lambda
                          1, # tau2R
	       1) #tau2Rho
  names(updates) <- c("nu",
                      "Delta",
                      "a",
                      "b",
                      "c2",
                      "gamma2",
                      "rAndC2",
                      "rhoAndGamma2",
                      "delta",
                      "xi",
                      "sigma2",
                      "tAndL",
                      "l",
                      "phi",
                      "thetaAndLambda",
                      "lambda",
                      "tau2R",
		      "tau2Rho")
  updates
}

###########################################################################
##Tuning parameters set as per Haakon's 10/23/2006 e-mail
###########################################################################
.tuningDefault <- function(){
	##10/23/2006: HT
	tuning <- c(0.01,  #1. nu Gibbs
		    0.01,  #2. Delta Gibbs
		    0.04,  #3. a
		    0.04,  #4. b
		    0.01,  #5. c2 Gibbs
		    0.01,  #6. gamma2 Gibbs
		    0.01,  #7. r and c2
		    0.01,  #8. rho and gamma2
		    0.01,  #9. delta
		    0.01,  #10. xi Gibbs
		    0.50,  #11. sigma2
		    0.10,  #12. t and l
		    0.04,  #13. l
		    0.40,  #14. phi
		    0.10,  #15. theta and lambda
		    0.02,  #16. lambda
		    0.04,  #17. tau2R
		    0.04) #18. tau2Rho
  names(tuning) <- c("nu", "Delta",  "a", "b", "c2", "gamma2",
                     "rAndC2", "rhoAndGamma2", "delta", "xi", "sigma2",
                     "tAndL", "l", "phi", "thetaAndLambda", "lambda", "tau2R", "tau2Rho")
  return(tuning)
}



empiricalStart <- function(object, zeroNu=FALSE,
			   phenotypeLabel, one.delta=FALSE, T_THRESH=4){
	if(length(grep(phenotypeLabel, varLabels(object[[1]]))) < 1) stop("phenotypeLabel must be in varLabels")
	potential <- rep(0, 19)
	acceptance <- rep(0, 17)

	##order should be platforms, genes for nu, phi, sigma2, delta, and DDelta
	if(!zeroNu){
		meanExpression <- function(eset) rowMeans(exprs(eset))
		Nu <- sapply(object, meanExpression)
		Rho <- cor(Nu)
		##order should be 1_2, 1_3, ..., 1_P, 2_3, ..., 2_P, ...
		##R <- R[upper.tri(R)]  ##Not in the right order!
		cors <- list()
		if(nrow(Rho) > 1){
			for(i in 1:(nrow(Rho)-1)) cors[[i]] <- Rho[i, -(1:i)]
		} else cors[[1]] <- Rho
		Rho <- unlist(cors)

		NuVars <- apply(Nu, 2, "var")
		Gamma2 <- mean(NuVars)
		S <- length(NuVars)
		tau2Rho <- NuVars/(prod(NuVars))^(1/S)
		Nu <- as.vector(t(Nu))
		##the prior for gamma2 is an improper uniform distribution on (0, Inf)
	} else{
		Nu <- as.vector(matrix(0, nrow(object), length(object)))
		Rho <- rep(0, choose(length(object), 2))
		tau2Rho <- rep(0, length(object))
		Gamma2 <- 0
	}
	averageDifference <- function(eset, classLabel){
		## limma would be an easy way to fit a linear model to each gene
		phenoColumn <- grep(classLabel, names(pData(eset)))
		pheno <- pData(eset)[, phenoColumn]
		mn1 <- rowMeans(exprs(eset)[, pheno == 1])
		mn0 <- rowMeans(exprs(eset)[, pheno == 0])
		avdif <- mn1-mn0
		avdif
	}
	DDelta <- sapply(object, averageDifference, classLabel=phenotypeLabel)
	require(genefilter)
	tt <- rowttests(object, phenotypeLabel)
	delta <- abs(tt) > T_THRESH
	if(one.delta){
		tmp <- ifelse(rowMeans(delta) > 0.5, 1, 0)
		delta <- matrix(tmp, nrow=length(tmp), ncol=length(object), byrow=FALSE)
	}
	if(any(colSums(delta) == 0)){
		##if none of the genes are differentially expressed in
		##some of the studies, set R at 0
		R <- rep(0, choose(S,2))
	} else 	{
		R <- cor(delta*DDelta)
		##order should be 1_2, 1_3, ..., 1_P, 2_3, ..., 2_P, ...
		##R <- R[upper.tri(R)]  ##Not in the right order!
		cors <- list()
		if(nrow(R) > 1){
			for(i in 1:(nrow(R)-1)) cors[[i]] <- R[i, -(1:i)]
		} else cors[[1]] <- R
		R <- unlist(cors)
	}
	DeltaVars <- apply(DDelta*delta, 2, "var")
	C2 <- mean(DeltaVars)
	S <- length(DeltaVars)
	if(any(colSums(delta)==0)){
		tau2R <- rep(1, S)
	} else {
		tau2R <- DeltaVars/(prod(DeltaVars))^(1/S)
	}
	##0.1 seems reasonable,  or the proportion of t-statistics > X
	Xi <- apply(delta, 2, "mean")
	Sigma2 <- sapply(lapply(object, exprs), rowVars)
	DDelta <- as.numeric(t(DDelta))
	delta <- as.numeric(t(delta))
	##Should we start at zero, 1, or somewhere in the middle?
	A <- rep(1, length(object))
	B <- rep(1, length(object))
	##Gamma with mean 1
	##L and T are mean and variance, respectively
	L <- as.numeric(colMeans(Sigma2))
	T <- as.numeric(apply(Sigma2, 2, var))
	Sigma2 <- as.numeric(t(Sigma2))
	##sigma2*phi is variance of expression values for patients with binary phenotype = 0
	##sigma2/phi is variance of expression values for patients with binary phenotype = 1
	Phi <- rep(1, (length(object) * nrow(object)))
	##Lambda and theta are the mean and variance for phi (study-specific)
	Lambda <- rep(1, length(object))
	Theta <- rep(0.05, length(object))
	##Tau2 is the relative scale of sigma across platforms
	##- the product is constrained to be one
	##- tau2 is shared by both nu and Delta
	##Tau2 <- rep(1, length(object))
	##Tau2 <- tau2.Delta * tau2.Nu
	linInitialValues <- list(Nu=Nu, DDelta=DDelta,
				 A=A, B=B, C2=C2, Gamma2=Gamma2,
				 R=R, Rho=Rho, Delta=delta,
				 Xi=Xi, Sigma2=Sigma2, T=T, L=L,
				 Phi=Phi, Theta=Theta, Lambda=Lambda,
				 Tau2R=tau2R,
				 Tau2Rho=tau2Rho)
	linInitialValues
}





computeGOF <- function(object,
		       nus,
		       Delta,
		       delta,
		       sigma2,
		       phi,
		       firstIteration,
		       lastIteration,
		       by=2){
	results <- list()
	avgM.k <- list()
	m.k <- 0
	##already thinned by 10
	I <- seq(firstIteration, lastIteration, by=by)
	NN <- length(I)
	for(p in seq(along=object)){
		cat("Study ", p, "\n")
		for(i in I){
			if(i %% 10 == 0) cat(i, " ")
			nus.1 <- nus[i, , ]
			Delta.1 <- Delta[i, , ]
			delta.1 <- delta[i, , ]
			sigma2.1 <- sigma2[i, , ]
			phi.1 <- phi[i, , ]
			sigma.1 <- sqrt(sigma2.1/phi.1)
			sigma.0 <- sqrt(sigma2.1*phi.1)
			##	trace(goodnessOfFit, browser)
			by <- 1/round(ncol(object[[p]])^0.4, 1)
			qtls <- qnorm(seq(0, 1, by=by), mean=0, sd=1)
			results[[p]] <- .goodnessOfFit(x=exprs(object[[p]]),
						      nus=nus.1[,p],
						      delta=delta.1[,p],
						      Delta=Delta.1[,p],
						      sigma.0=sigma.0[, p],
						      sigma.1=sigma.1[, p],
						      psi=psi[[p]],
						      qtls=qtls)
			m.k <- results[[p]]$m.k + m.k
		}
		m.k <- m.k/NN
		avgM.k[[p]] <- m.k
		m.k <- 0
		cat("\n")
	}
	pvals <- R.b <- list()
	for(p in 1:4){
		n.p.k <- ncol(object[[p]])*1/ncol(avgM.k[[p]])
		tmp <- ((avgM.k[[p]] - n.p.k)/sqrt(n.p.k))^2
		R.b[[p]] <- rowSums(tmp)
		pvals[[p]] <- -log10(1-pchisq(R.b[[p]], df=ncol(tmp)-1))
	}
	goodnessOfFit_results <- list(R.b, pvals)
	goodnessOfFit_results
}

.goodnessOfFit <- function(x, nus, delta, Delta, sigma.0, sigma.1, psi, qtls){
	Psi <- matrix(psi, nrow=nrow(x), ncol=ncol(x), byrow=TRUE)
	mus <- nus + delta * (2*Psi - 1)*Delta
	sigma <- sigma.0*(1-Psi) + sigma.1 * Psi
	m.k <- matrix(NA, nrow(x), length(qtls[-1]))

	##x <- exprs(object[[p]]) - mean(as.numeric(exprs(object[[p]])))
	x <- x-mean(as.numeric(x))
	z <- (x-mus)/sigma
	##determine bins by the inverse distribution function evaluated at theta.
	##
	##because we have a different mean and variance for psi=0,
	##psi=1, I standardized the x-values by theta and then
	##calculated the inverse quantiles
	R.b <- rep(NA, nrow(z))
	for(i in 1:nrow(z)){
		m.k[i, ] <- as.numeric(table(cut(as.numeric(z[i, ]), breaks=qtls, labels=seq(along=qtls[-1]))))
		n.p.k <- ncol(z)*1/length(qtls[-1])
		R.b[i] <- sum(((m.k[i, ] - n.p.k)/sqrt(n.p.k))^2)
	}
	pvals <- -log10(1-pchisq(R.b, df=ncol(m.k)-1))
	pvals[is.infinite(pvals)] <- 16
	require(genefilter)
	sds0 <- rowSds(z[, psi==0])
	sds1 <- rowSds(z[, psi==1])
	n0 <- sum(psi==0)
	n1 <- sum(psi==1)
	sds <- (sds0*(n0-1) + sds1*(n1-1))/(n0+n1-2)
	return(list(pvals=pvals, R.b=R.b, sds=sds, z=z, m.k=m.k))
}

makeSigma <- function(G, Q, gamma2, tau2, a, sigma2, r){
	sigma <- sqrt(gamma2*tau2*exp(a*log(sigma2)))
	Sigma <- sigma %*% t(sigma)
	Sigma <- Sigma*r
	## check diagonals
	if(FALSE){
		correct.answer <- gamma2*tau2*exp(a*log(sigma2))
		answer2 <- diag(Sigma)
		all.equal(correct.answer, answer2)
	}
	return(Sigma)
}

##getParameters <- function(hyper.parameters, one.delta=FALSE, MRF=FALSE){
##Params <- function(hyper.parameters, one.delta=FALSE, MRF=FALSE){






