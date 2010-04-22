xde <- function(paramsMcmc, esetList, outputMcmc, batchSize=NULL, NCONC=2, center=TRUE){
	if(is.null(batchSize)){
		fit <- xdeFit(paramsMcmc=paramsMcmc,
			      esetList=esetList,
			      outputMcmc=outputMcmc,
			      center=center)
	} else{
		## Tabulate running average of posterior means
		## Remove log files after each batch
		print("Log files automatically removed after each batch")
		I <- iterations(paramsMcmc)/batchSize
		iterations(paramsMcmc) <- batchSize
		for(i in 1:I){
			cat("Batch ", i, "\n")
			if(i > 1){
				firstMcmc(paramsMcmc) <- lastMcmc(fit)					
				fit <- xdeFit(paramsMcmc=paramsMcmc,
					      esetList=esetList,
					      center=center)			
				pa.last <- calculatePosteriorAvg(fit, NCONC=NCONC)
				pa.cum <- (pa.last*batchSize + pa.cum*(i-1)*batchSize)/(i*batchSize)
				
				bes.last <- calculateBayesianEffectSize(fit)
				bes.cum <- (bes.last*batchSize + bes.cum*(i-1)*batchSize)/(i*batchSize)
			} else{
				fit <- xdeFit(paramsMcmc=paramsMcmc,
					      esetList=esetList,
					      center=center)						
				pa.cum <- calculatePosteriorAvg(fit, NCONC=NCONC)
				bes.cum <- calculateBayesianEffectSize(fit)
			}
			unlink(directory(paramsMcmc), recursive=TRUE)			
		}
		bayesianEffectSize(fit) <- bes.cum
		posteriorAvg(fit) <- pa.cum
	}
	return(fit)
}

xdeFit <- function(paramsMcmc, esetList, outputMcmc, center=TRUE){
	if(!file.exists(directory(paramsMcmc))){
		dir.create(directory(paramsMcmc), recursive=TRUE)
	}
	if(class(esetList) != "ExpressionSetList") {
		esetList <- as(esetList, "ExpressionSetList")
	}
	##stop("argument esetList must have class ExpresssionSetList")
	if(class(paramsMcmc) != "XdeParameter") stop("paramsMcmc must be an object of class XdeParameter")
	##if outputMcmc is present, use the seed stored in outputMcmc for
	##beginning the MCMC (prevents simulating from the same seed as the
	##previous run).  Furthermore, start the chain using the value of
	##chain stored in outputMcmc (by default, this would be the last
	##iteration from the previous run).
	if(!missing(outputMcmc)){
		if(class(outputMcmc) == "XdeMcmc"){
			seed(paramsMcmc) <- outputMcmc@seed
			firstMcmc(paramsMcmc) <- lastMcmc(outputMcmc)
		} else {
			stop("Argument 'outputMcmc' must be an object of class XdeMcmc")
		}
	}
	if(center){
		studyMean <- function(x) abs(mean(as.numeric(exprs(x)))) > 0.1
		overall.means <- sapply(esetList, studyMean)
		if(any(overall.means)){
			print("Replacing expression matrices in esetList with mean-centered matrices...")
			esetList <- studyCenter(esetList)
		}
	}
	##Convert expression to vector
	f <- function(x)  x <- as.vector(t(exprs(x)))
	expression.vector <- as.numeric(unlist(lapply(esetList, f)))
	if(any(is.na(expression.vector))){
		stop("Missing values in the expression data")
	}
	##Convert phenoData to a vector 
	convertToVector <- function(x, phenotypeLabel){
		x <- pData(x)
		if(dim(x)[2] > 1 & is.null(phenotypeLabel)){
			warning("More than one phenoData covariate -- using the first")
			x <- x[, 1]
		}
		if(dim(x)[2] > 1 & !is.null(phenotypeLabel)){
			var.i <- match(phenotypeLabel, colnames(x))
			x <- x[, var.i]
		}
		return(x)
	}
	psi <- as.numeric(unlist(lapply(esetList, convertToVector, phenotypeLabel=paramsMcmc@phenotypeLabel)))
	if(paramsMcmc@verbose) cat("Covariate to use in phenoData: ", phenotypeLabel(paramsMcmc), "\n")
	if(paramsMcmc@verbose) print(paste("writing chain files to", directory(paramsMcmc)))
	firstMcmc <- firstMcmc(paramsMcmc)
	ZZ <- .xdeparameterlist(paramsMcmc=paramsMcmc,
				esetList=esetList,
				firstMcmc=firstMcmc,
				expression.vector=expression.vector,
				psi=psi)
	TRY <- tryCatch(results <- .C("xdeLIN_main",
				      seed = ZZ[["seed"]],
				      showIterations = ZZ[["showIterations"]],
				      nIt = ZZ[["nIt"]],
				      oneDelta=ZZ[["one.delta"]],  ##if zero, prior model A (multiple delta) is used				
				      Q = ZZ[["Q"]],
				      G = ZZ[["G"]],
				      S = ZZ[["S"]],
				      x = ZZ[["x"]],
				      Psi = ZZ[["Psi"]],
				      specifiedInitialValues = ZZ[["specifiedInitialValues"]],
				      Nu = ZZ[["Nu"]],
				      DDelta = ZZ[["DDelta"]], 
				      A = ZZ[["A"]],
				      B = ZZ[["B"]],
				      C2 = ZZ[["C2"]],
				      Gamma2 = ZZ[["Gamma2"]],
				      R = ZZ[["R"]],
				      Rho = ZZ[["Rho"]],
				      Delta = ZZ[["Delta"]],
				      Xi = ZZ[["Xi"]],
				      Sigma2 = ZZ[["Sigma2"]],
				      T = ZZ[["T"]],
				      L = ZZ[["L"]],
				      Phi = ZZ[["Phi"]],
				      Theta = ZZ[["Theta"]],
				      Lambda = ZZ[["Lambda"]],
				      Tau2R = ZZ[["Tau2R"]],
				      Tau2Rho=ZZ[["Tau2Rho"]],
				      param = ZZ[["param"]],
				      nUpdate = ZZ[["nUpdate"]],
				      epsilon = ZZ[["epsilon"]],
				      output = ZZ[["output"]],
				      writeToFile = ZZ[["writeToFile"]],
				      directory = ZZ[["directory"]],
				      valuePotential = ZZ[["valuePotential"]],
				      valueAcceptance = ZZ[["valueAcceptance"]],
				      valueNu = ZZ[["valueNu"]],
				      valueDDelta = ZZ[["valueDDelta"]],
				      valueA = ZZ[["valueA"]],
				      valueB = ZZ[["valueB"]],
				      valueC2 = ZZ[["valueC2"]],
				      valueGamma2 = ZZ[["valueGamma2"]],
				      valueR = ZZ[["valueR"]],
				      valueRho = ZZ[["valueRho"]],
				      valueDelta = ZZ[["valueDelta"]],
				      valueXi = ZZ[["valueXi"]],
				      valueSigma2 = ZZ[["valueSigma2"]],
				      valueT = ZZ[["valueT"]],
				      valueL = ZZ[["valueL"]],
				      valuePhi = ZZ[["valuePhi"]],
				      valueTheta = ZZ[["valueTheta"]],
				      valueLambda = ZZ[["valueLambda"]],
				      valueTau2R = ZZ[["valueTau2R"]],
				      valueTau2Rho=ZZ[["valueTau2Rho"]],
				      valueProbDelta = ZZ[["valueProbDelta"]],
				      valueDiffexpressed=ZZ[["valueDiffexpressed"]],
				      writeDiffexpressedTofile=ZZ[["writeDiffexpressedTofile"]]),                                          
                  error=function(e) NULL)
	if(is.null(TRY)){
		stop("Problem in C wrapper to xdeLIN_main. Check arguments in .fit")
	}
	##The data can be stored more efficiently in environments?
	##If writeToFile or burnin is TRUE, we only collect the last
	##iteration of the chain
	index <- match(c("Nu",
			 "DDelta",
			 "A",
			 "B",
			 "C2",
			 "Gamma2",
			 "R",
			 "Rho",
			 "Delta",
			 "Xi",
			 "Sigma2",
			 "T",
			 "L",
			 "Phi",
			 "Theta",
			 "Lambda",
			 "Tau2R",
			 "Tau2Rho"), names(results))
	lastMcmc <- results[index]
	eenv <- new.env()
	eenv[["vars"]] <- lastMcmc
	xmcmc <- new("XdeMcmc",
		     studyNames=studyNames(paramsMcmc),
		     featureNames=featureNames(esetList),
		     iterations=savedIterations(paramsMcmc)[[1]],
		     directory=directory(paramsMcmc),
		     seed=results[[1]],
		     output=output(paramsMcmc),
		     lastMcmc=eenv,
		     posteriorAvg=NULL,
		     bayesianEffectSize=NULL)
##	if(!burnin(paramsMcmc)){
##		##paramsMcmc has directory where chains are stored, number of iterations, etc
##		##i <- match(c("c2", "sigma2", "tau2", "b"), names(output(paramsMcmc)))
##		##xmcmc@posteriorAvg <- xmcmc$diffExpressed
##		##if(nrow(xmcmc@posteriorAvg) == ZZ$G){
##		##	rownames(xmcmc@posteriorAvg) <- featureNames(esetList)
##		##}
##	} 
	xmcmc
}


.xdeparameterlist <- function(paramsMcmc,
			      esetList,
			      firstMcmc,
                              expression.vector,
			      psi){
	ZZ <-list(seed = as.integer(seed(paramsMcmc)),#1
		  showIterations = as.integer(paramsMcmc@showIterations),
		  nIt = as.integer(iterations(paramsMcmc)),
		  one.delta=as.integer(paramsMcmc@one.delta),## if TRUE, "prior model B" (multiple delta) is used
		  Q = as.integer(length(esetList)),
		  G = as.integer(nrow(esetList)),
		  S = as.integer(nSamples(esetList)),
		  x = as.double(expression.vector),
		  Psi = as.integer(psi),
		  specifiedInitialValues = as.integer(paramsMcmc@specifiedInitialValues),
		  ##Nu: a vector of GP values.  Must be of the right length even if specifiedInitialValues is FALSE
		  Nu = as.double(firstMcmc$Nu),#11
		  DDelta = as.double(firstMcmc$DDelta),
		  A = as.double(firstMcmc$A),
		  B = as.double(firstMcmc$B),
		  C2 = as.double(firstMcmc$C2),
		  Gamma2 = as.double(firstMcmc$Gamma2),#16
		  R = as.double(firstMcmc$R),
		  Rho = as.double(firstMcmc$Rho),
		  Delta = as.integer(firstMcmc$Delta),
		  Xi = as.double(firstMcmc$Xi),#20
		  Sigma2 = as.double(firstMcmc$Sigma2), #21
		  T = as.double(firstMcmc$T),
		  L = as.double(firstMcmc$L),
		  Phi = as.double(firstMcmc$Phi),
		  Theta = as.double(firstMcmc$Theta), #25
		  Lambda = as.double(firstMcmc$Lambda),#26
##		  Tau2 = as.double(firstMcmc$Tau2),#27  ##decoupling of taus
		  Tau2R=as.double(firstMcmc$Tau2R),
		  Tau2Rho=as.double(firstMcmc$Tau2Rho),
		  param = as.double(hyperparameters(paramsMcmc)),#28
		  nUpdate = as.integer(paramsMcmc@updates),#29
		  epsilon = as.double(paramsMcmc@tuning),#30
		  output = as.integer(paramsMcmc@output),#31 (should have length 22)
		  ##do not write to file if running burnin
		  writeToFile = as.integer(!paramsMcmc@burnin), #32
		  directory = as.character(paramsMcmc@directory), #,#33
		  valuePotential = as.double(firstMcmc$Potential),#34
		  valueAcceptance = as.double(firstMcmc$Acceptance),#35
		  valueNu = as.double(firstMcmc$Nu),#36
		  valueDDelta = as.double(firstMcmc$DDelta),#37
		  valueA = as.double(firstMcmc$A),#38
		  valueB = as.double(firstMcmc$B),#39
		  valueC2 = as.double(firstMcmc$C2),#40
		  valueGamma2 = as.double(firstMcmc$Gamma2),#41
		  valueR = as.double(firstMcmc$R),#42
		  valueRho = as.double(firstMcmc$Rho),#43
		  valueDelta = as.integer(firstMcmc$Delta),#44
		  valueXi = as.double(firstMcmc$Xi),#45
		  valueSigma2 = as.double(firstMcmc$Sigma2),#46
		  valueT = as.double(firstMcmc$T),#47
		  valueL = as.double(firstMcmc$L),#48
		  valuePhi = as.double(firstMcmc$Phi),#49
		  valueTheta = as.double(firstMcmc$Theta),#50
		  valueLambda = as.double(firstMcmc$Lambda),#51
		  valueTau2R = as.double(firstMcmc$Tau2R),#52
		  valueTau2Rho = as.double(firstMcmc$Tau2Rho),#52		  
		  valueProbDelta = as.double(firstMcmc$ProbDelta),#53
		  ##posterior mean of differential expression,  concordant differential expression, and discordant differential expression            
		  valueDiffexpressed = as.double(firstMcmc$Diffexpressed), ##54
		  writeDiffexpressedTofile = as.integer(FALSE)) ##55 -- 
	ZZ
}
