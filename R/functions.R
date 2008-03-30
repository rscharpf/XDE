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
                                        "tau2",
                                        "probDelta",
                                        "diffExpressed")
              
            ##These values are placeholders to create a chain of the
            ##right dimension.  If called, it is assumed that default
            ##starting values for the chain will be used.
            G <- nrow(object)
            QQ <- length(object)
            potential <- rep(0, 19*chain.length[["potential"]])
            acceptance <- rep(0, 17*chain.length[["acceptance"]])
            Nu <- rep(0, G*QQ*chain.length[["nu"]])
            DDelta <- rep(0, G*QQ*chain.length[["deltaDelta"]])
            A <- rep(0, QQ*chain.length[["a"]])
            B <- rep(0, QQ*chain.length[["b"]])
            C2 <- rep(0, chain.length[["c2"]])
            Gamma2 <- rep(0, chain.length[["gamma2"]])
            Rho <- rep(0, QQ*(QQ-1)/2*chain.length[["rho"]])
            R <- rep(0, QQ*(QQ-1)/2*chain.length[["r"]])              
##            Delta <- rep(0, G*chain.length[["delta"]])
            Delta <- rep(0, G*QQ*chain.length[["delta"]])
            Xi <- rep(0, QQ*chain.length[["xi"]])
            Sigma2 <- rep(0, QQ * G*chain.length[["sigma2"]])
            T <- rep(0, QQ*chain.length[["t"]])
            L <- rep(0, QQ*chain.length[["l"]])
            Phi <- rep(0, QQ * G*chain.length[["phi"]])
            Theta <- rep(0, QQ*chain.length[["theta"]])
            Lambda <- rep(0, QQ*chain.length[["lambda"]])
            Tau2 <- rep(0, QQ*chain.length[["tau2"]])
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
                      Tau2=Tau2,
                      ProbDelta=ProbDelta,
                      Diffexpressed=Diffexpressed)
              assign("vars", x, eenv)
              eenv
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
    "tau2",
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
            ##object is an ExpressionSet (or BreastSet)
  require(gtools) || stopifnot("gtools package not available")
  pData(object)[, phenotype] <- permute(pData(object)[, phenotype])
  object
}

ssStatistic <- function(statistic=c("t", "sam", "z")[1],
                        phenotypeLabel,
                        esetList, ...){
  if(!(statistic == "t" | statistic == "sam" | statistic == "z")){
    stop("only t and sam study-specific statistics provided")
  }

  ###########################################################################
  ##Wrapper for Welch t-statistic
  multtestWrapper <- function(eset, classlabel, ...){
    require(multtest) || stop("Package multtest not available")    
    column <- grep(classlabel, colnames(pData(eset)))
    if(length(column) == 0) { warning("Invalid classlabel.  Using the first column in phenoData"); classlabel <- pData(eset)[, 1]}
    if(length(column) > 1)  {warning("More than 1 phenotype has the class label.  Using the first"); classlabel <- pData(eset)[, column[1]]}
    if(length(column) == 1) classlabel <- pData(eset)[, column]
    X <- exprs(eset)
    stat <- mt.teststat(X=X, classlabel=classlabel, test="t",
                        nonpara="n")
    stat
  }

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
    options(warn=-1)
    z <- zScores(object, classes=classes, useREM=useREM)
##    z <- data.frame(z[, grep("zSco", colnames(z))])
    z <- z[, grep("zSco", colnames(z))]
    return(z)
  }

  ###########################################################################
  ##Wrapper for SAM-statistic
  stat <- switch(statistic,
                 t=sapply(esetList, multtestWrapper, classlabel=phenotypeLabel),
                 sam=sapply(esetList, sam.wrapper, classlabel=phenotypeLabel,
                   returnQ=FALSE, method="d.stat", gene.names=featureNames(esetList), ...),
                 z=z.wrapper(object=esetList, classlabel=phenotypeLabel, useREM=TRUE))
  rownames(stat) <- featureNames(esetList)
  return(stat)
}

.integ.norm<-function(x){
  xm<-mean(x,na.rm = TRUE)
  ss<-sum((x-xm)*(x-xm),na.rm = TRUE)
  return((x-xm)/sqrt(ss))
}

##.integ.cal <- .integ.cal ##From MergeMaid.R
.integ.cal<-function(mat1,mat2,method=NULL){
  rn<-rownames(mat1)
  integ<-rep(0,nrow(mat1))

  m1<- t(apply(mat1,1,.integ.norm))
  m2<- t(apply(mat2,1,.integ.norm))

  sziicor=function(index){
    x<-m1[index,]
    y<-m2[index,]
  
    A<-m1[-index,]
    B<-m2[-index,]
    nammx=function(indx){
      aa=A[indx,]
      sum(aa*x,na.rm = TRUE)
    }
    rx=sapply(1:(nrow(m1)-1),nammx)
  
    nammy=function(indy){
      bb=B[indy,]
      sum(bb*y,na.rm = TRUE)
    }
    ry=sapply(1:(nrow(m2)-1),nammy)
  
    integ[index]<-cor(rx,ry,use="pairwise.complete.obs")
  }
  ans=sapply(1:(nrow(mat1)),sziicor)
  names(ans)=rn
  return(ans)
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
  
  ##################################################
  ##Differential expression
  ##################################################  
  scores[, 1] <- .pca(abs(statistic), N=N)

  ##################################################
  ##Concordance
  ##################################################  
  scores[, 2] <- abs(.pca(statistic, N=N))


  ##################################################
  ##Discordance
  ##################################################    
  I <- abs(rowSums(sign(statistic))) < ncol(statistic)
  I[!I] <- -1
  scores[, 3] <- .pca(abs(statistic), N)

  return(scores)
}


.outputDefault <- function(thin=1, burnin=FALSE){
##  if(!burnin) out <- rep(1, 22) else out <- rep(0, 22)
	out <- rep(1, 22)
	out[1] <- thin
	out[22] <- 1  ##by default, keep a running average of the posterior statistics attached to the object (do not write to file)
	names(out) <- c("thin",
			"potential",
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
			"tau2",
			"p.delta",
			"diffExpressed")
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
                          1) # tau2    
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
                      "tau2")
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
              0.04)  #17. tau2                
  names(tuning) <- c("nu", "Delta",  "a", "b", "c2", "gamma2",
                     "rAndC2", "rhoAndGamma2", "delta", "xi", "sigma2",
                     "tAndL", "l", "phi", "thetaAndLambda", "lambda", "tau2")
  return(tuning)
}


empiricalStart <- function(object, zeroNu=FALSE){
	if(length(grep(phenotypeLabel(object), varLabels(object[[1]]))) < 1) stop("phenotypeLabel must be in varLabels")
	potential <- rep(0, 19)
	acceptance <- rep(0, 17)
  

	if(!zeroNu){
		meanExpression <- function(eset) rowMeans(exprs(eset))  
		Nu <- sapply(object, meanExpression)
		Rho <- cor(Nu)[upper.tri(cor(Nu))]
		Nu <- as.vector(Nu)
		Gamma2 <- var(Nu)    
	} else{
		Nu <- as.vector(matrix(0, nrow(object), length(object)))
		Rho <- rep(0, choose(length(object), 2))
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
	DDelta <- sapply(object, averageDifference, classLabel=phenotypeLabel(object))
	R <- c(cor(DDelta[, 1], DDelta[, 2]), cor(DDelta[, 1], DDelta[,3]), cor(DDelta[,2], DDelta[,3]))        

	quant <- quantile(rowMeans(DDelta), probs=c(0.1, 0.9))
	Delta <- ifelse(rowMeans(DDelta) < quant[1] | rowMeans(DDelta) > quant[2], 1, 0)
	Delta <- as.vector(Delta)
	DDelta <- as.vector(DDelta)
	C2 <- var(DDelta)

	##Should we start at zero or somewhere in the middle
	A <- rep(0.5, length(object))
	B <- rep(0.5, length(object))

	##0.2 seems reasonable,  or the proportion of t-statistics > X
	Xi <- 0.2
  
	##Start at the empirical standard deviation
	varianceExpression <- function(eset){
		apply(exprs(eset), 1, var)
	}
	Sigma2 <- sapply(object, varianceExpression)

	##Gamma with mean 1
	##L and T are mean and variance, respectively
	##
	L <- as.numeric(colMeans(Sigma2))
	T <- as.numeric(apply(Sigma2, 2, var))
	Sigma2 <- as.vector(Sigma2)

	##Lambda and theta are the mean and variance
	##respectively (as per Haakon's 10/23/2006 e-mail)
	Lambda <- rep(1, length(object))
	Theta <- rep(0.001, length(object))

  
	Phi <- rep(1, (length(object) * nrow(object)))
	Tau2 <- rep(1, length(object))
	ProbDelta <- Delta
	eenv <- new.env()
	linInitialValues <- list(#potential=potential, 
					#acceptance=acceptance, 
				 Nu=Nu, DDelta=DDelta, 
				 A=A, B=B, C2=C2, Gamma2=Gamma2, 
				 R=R, Rho=Rho, Delta=Delta,
				 Xi=Xi, Sigma2=Sigma2, T=T, L=L, 
				 Phi=Phi, Theta=Theta, Lambda=Lambda, 
				 Tau2=Tau2, ProbDelta=ProbDelta)
	linInitialValues
}












