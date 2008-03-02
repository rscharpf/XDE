xde <- function(paramsMcmc, esetList, outputMcmc){
  if(class(esetList) != "ExpressionSetList") stop("argument esetList must have class ExpresssionSetList")
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
  
  studyMean <- function(x) abs(mean(as.vector(exprs(x)))) > 0.1
  overall.means <- sapply(esetList, studyMean)
  if(any(overall.means)){
    if(paramsMcmc@verbose) print("Replacing expression matrices in esetList with mean-centered matrices...")
    esetList <- studyCenter(esetList)
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
                                seed = ZZ[[1]],
                                showIterations = ZZ[[2]],
                                nIt = ZZ[[3]],
                                Q = ZZ[[4]],
                                G = ZZ[[5]],
                                S = ZZ[[6]],
                                x = ZZ[[7]],
                                Psi = ZZ[[8]],
                                specifiedInitialValues = ZZ[[9]],
                                Nu = ZZ[[10]],
                                DDelta = ZZ[[11]], 
                                A = ZZ[[12]],
                                B = ZZ[[13]],
                                C2 = ZZ[[14]],
                                Gamma2 = ZZ[[15]],
                                R = ZZ[[16]],
                                Rho = ZZ[[17]],
                                Delta = ZZ[[18]],
                                Xi = ZZ[[19]],
                                Sigma2 = ZZ[[20]],
                                T = ZZ[[21]],
                                L = ZZ[[22]],
                                Phi = ZZ[[23]],
                                Theta = ZZ[[24]],
                                Lambda = ZZ[[25]],
                                Tau2 = ZZ[[26]],
                                param = ZZ[[27]],
                                nUpdate = ZZ[[28]],
                                epsilon = ZZ[[29]],
                                output = ZZ[[30]],
                                writeToFile = ZZ[[31]],
                                directory = ZZ[[32]],
                                valuePotential = ZZ[[33]],
                                valueAcceptance = ZZ[[34]],
                                valueNu = ZZ[[35]],
                                valueDDelta = ZZ[[36]],
                                valueA = ZZ[[37]],
                                valueB = ZZ[[38]],
                                valueC2 = ZZ[[39]],
                                valueGamma2 = ZZ[[40]],
                                valueR = ZZ[[41]],
                                valueRho = ZZ[[42]],
                                valueDelta = ZZ[[43]],
                                valueXi = ZZ[[44]],
                                valueSigma2 = ZZ[[45]],
                                valueT = ZZ[[46]],
                                valueL = ZZ[[47]],
                                valuePhi = ZZ[[48]],
                                valueTheta = ZZ[[49]],
                                valueLambda = ZZ[[50]],
                                valueTau2 = ZZ[[51]],
                                valueProbDelta = ZZ[[52]],
                                valueDiffexpressed=ZZ[[53]],
                                writeDiffexpressedTofile=ZZ[[54]]),                                          
                  error=function(e) NULL)
  if(is.null(TRY)){
    stop("Problem in C wrapper to xdeLIN_main. Check arguments in .fit")
  }
  ##The data can be stored more efficiently in environments?

  ##If writeToFile or burnin is TRUE, we only collect the last
  ##iteration of the chain
  lastMcmc <- list(results[[10]],
                        results[[11]],
                        results[[12]],
                        results[[13]],
                        results[[14]],
                        results[[15]],
                        results[[16]],
                        results[[17]],
                        results[[18]],
                        results[[19]],
                        results[[20]],
                        results[[21]],
                        results[[22]],
                        results[[23]],
                        results[[24]],
                        results[[25]],
                        results[[26]])
  eenv <- new.env()
  assign("vars", lastMcmc, eenv)
  names(eenv$vars) <- names(results)[10:26]
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
  if(!burnin(paramsMcmc)){
    ##paramsMcmc has directory where chains are stored, number of iterations, etc
    i <- match(c("c2", "sigma2", "tau2", "b"), names(output(paramsMcmc)))
    xmcmc@posteriorAvg <- xmcmc$diffExpressed
    rownames(xmcmc@posteriorAvg) <- featureNames(esetList)
  } 
  xmcmc
}


.xdeparameterlist <- function(paramsMcmc, esetList, firstMcmc,
                              expression.vector, psi){
  ZZ <-list(seed = as.integer(seed(paramsMcmc)),#1
            showIterations = as.integer(paramsMcmc@showIterations),
            nIt = as.integer(iterations(paramsMcmc)),
            Q = as.integer(length(esetList)),
            G = as.integer(nrow(esetList)),
            S = as.integer(nSamples(esetList)),
            x = as.double(expression.vector),
            Psi = as.integer(psi),
            specifiedInitialValues = as.integer(paramsMcmc@specifiedInitialValues),
            ##Nu: a vector of GP values.  Must be of the right length even if specifiedInitialValues is FALSE
            Nu = as.double(firstMcmc$Nu),#10
            DDelta = as.double(firstMcmc$DDelta),
            A = as.double(firstMcmc$A),
            B = as.double(firstMcmc$B),
            C2 = as.double(firstMcmc$C2),
            Gamma2 = as.double(firstMcmc$Gamma2),#15
            R = as.double(firstMcmc$R),
            Rho = as.double(firstMcmc$Rho),
            Delta = as.integer(firstMcmc$Delta),
            Xi = as.double(firstMcmc$Xi),#19
            Sigma2 = as.double(firstMcmc$Sigma2), #20
            T = as.double(firstMcmc$T),
            L = as.double(firstMcmc$L),
            Phi = as.double(firstMcmc$Phi),
            Theta = as.double(firstMcmc$Theta), #24
            Lambda = as.double(firstMcmc$Lambda),#25
            Tau2 = as.double(firstMcmc$Tau2),#26
            param = as.double(hyperparameters(paramsMcmc)),#27
            nUpdate = as.integer(paramsMcmc@updates),#28
            epsilon = as.double(paramsMcmc@tuning),#29
            output = as.integer(paramsMcmc@output),#30 (should have length 22)
            ##do not write to file if running burnin
            writeToFile = as.integer(!paramsMcmc@burnin), #31
            directory = as.character(paramsMcmc@directory), #,#32
            valuePotential = as.double(firstMcmc$Potential),#33
            valueAcceptance = as.double(firstMcmc$Acceptance),#34
            valueNu = as.double(firstMcmc$Nu),#35
            valueDDelta = as.double(firstMcmc$DDelta),#36
            valueA = as.double(firstMcmc$A),#37
            valueB = as.double(firstMcmc$B),#38
            valueC2 = as.double(firstMcmc$C2),#39
            valueGamma2 = as.double(firstMcmc$Gamma2),#40
            valueR = as.double(firstMcmc$R),#41
            valueRho = as.double(firstMcmc$Rho),#42
            valueDelta = as.integer(firstMcmc$Delta),#43
            valueXi = as.double(firstMcmc$Xi),#44
            valueSigma2 = as.double(firstMcmc$Sigma2),#45
            valueT = as.double(firstMcmc$T),#46
            valueL = as.double(firstMcmc$L),#47
            valuePhi = as.double(firstMcmc$Phi),#48
            valueTheta = as.double(firstMcmc$Theta),#49
            valueLambda = as.double(firstMcmc$Lambda),#50
            valueTau2 = as.double(firstMcmc$Tau2),#51
            valueProbDelta = as.double(firstMcmc$ProbDelta),#52
            ##posterior mean of differential expression,  concordant differential expression, and discordant differential expression            
            valueDiffexpressed = as.double(firstMcmc$Diffexpressed), ##53
            writeDiffexpressedTofile = as.integer(FALSE)) ##54 -- indicator of whether to write the posterior mean to file (1=yes, 0 = no)
  ZZ
}
