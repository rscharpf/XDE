setMethod("xapply", "ExpressionSetList", function(X, FUN, ...){
	class(X) <- "list"
	X <- as(lapply(X, FUN, ...), "ExpressionSetList")
	X <- as(X, "ExpressionSetList")
})

##object is a LinearXdeSet or loglinearXdeSet
setMethod(".integrativeCorrelationFilter", "ExpressionSetList",
          function(object, fdrCut, ...){
            require(gtools) || stop("gtools package not available")
            require(MergeMaid) || stop("MergeMaid package not available")
            if(missing(fdrCut)) stop("must specify fdrCut (FDR threshold)")
            ##contains a fix so that intCor will work

	    ##check that this works-- can comment out
##	    gS <- matrix(1, nrow=nrow(object), ncol=length(object))
##	    obj <- new("mergeExpressionSet",
##		       data=object,
##		       geneStudy=gS,
##		       notes="")
	    
            ij <- combinations(length(object), 2)
            exprSetList <- lapply(object, exprs)

            print("calculating integrative correlation coefficients")      
            iC <- function(ij, objList){
              i <- ij[1]
              j <- ij[2]
              ic <- XDE:::.integ.cal(objList[[i]], objList[[j]])
            }
            emp <- apply(ij, 1, iC, exprSetList)
  
            ##This is the approximation to the null (does fine if you have
            ##thousands of genes)
            print("approximating the null distribution for the integrative correlation...")
            ##This function will not work if there are no rownames on geneStudy ...should check for this
##            setAs("ExpressionSetList", "mergeExpressionSet",
##                  function(from, to){
##                    gS <- matrix(1, nrow=nrow(from[[1]]), ncol=length(from))
##                    rownames(gS) <- featureNames(from)
##                    new("mergeExpressionSet",
##                        data=from,
##                        geneStudy=gS,
##                        notes="")
##                  })
	    gS <- matrix(1, nrow=nrow(object), ncol=length(object))
	    rownames(gS) <- featureNames(object)
	    obj <- new("mergeExpressionSet",
		       data=object,
		       geneStudy=gS,
		       notes="")
##            obj <- as(object, "mergeExpressionSet")
            null <- intcorDens(obj)
            null <- do.call("cbind", null)            

            ##Find the cutoff, c,  at which
            ##(NULL_c)/(NULL_c + EMP_c) = 0.3
            ##c is a quantile of the integrative correlation coefficients
            ##Higher c indicates higher reproducibility
            ##Lowering the fdrCut increases the value of c
            findCuts <- function(x1, x2, probs=seq(0.6, 0.99, by=0.005), fdrCut){
              k <- 1
              quantiles <- quantile(x2, probs)
              fdrhat <- 1
              while(fdrhat > fdrCut & k <= length(quantiles)){
                p.null <- sum(x1 > quantiles[k])/length(x1)
                p.emp <- sum(x2 > quantiles[k])/length(x2)
                fdrhat <- p.null/(p.null + p.emp)
                print(fdrhat)
                if(is.na(fdrhat)) { warning("0 in fdr denominator"); break()}
                if(k+1 < length(quantiles)) k <- k+1 else break()
              }
              quantiles[k]
            }
            N <- ncol(emp)
            cut <- rep(NA, N)
            for(i in 1:N)   cut[i] <- findCuts(null[, i], emp[, i], fdrCut=fdrCut)
            print(paste("integrative correlation coefficient threshold:", cut))
  
            ##Use a conservative threshold that keeps more genes
            cuts <- matrix(NA, nrow=nrow(emp), ncol=ncol(emp))
            for(i in 1:N){
              cuts[, i] <- emp[, i] > cut[i]
            }
            
            ##Take union of gene passing cut
            print("keep the union of genes passing the threshold")            
            union <- apply(cuts, 1, "any")
            nG <- sum(union)
            print(paste(nG, "genes remaining after first IC"))
            object <- object[union, ]
            object
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
              f <- function(x, j){
                x <- x[, j]
              }
              x <- lapply(x, f, j)
            }
            if(!missing(i)){
              f <- function(x, i){
                x <- x[i, ]
              }
              x <- lapply(x, f, i)
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

##setMethod("qqplot", signature(x = "ExpressionSetList", y = "character"),
##          function(x, y, plot.it = TRUE, xlab = deparse(substitute(x)), 
##                   ylab = deparse(substitute(y)), filename = NULL,
##                   mar = c(3, 3, 2, 2),
##                   oma = rep(2, 4), ...){
##            object <- x
##            statistic <- y
##            xlab <- statistic[1]
##            ylab <- statistic[2]
##
##            x <- ssStatistic(statistic=xlab,
##                             xdeParams=xdeParams,
##                             esetList=esetList)
##            y <- ssStatistic(statistic=ylab,
##                             xdeParams=xdeParams,
##                             esetList=esetList)
####            if(length(statistic) < 2) stop("must specify two statistics")
####            getStats <- function(object, method){
####                switch(method,
####                       t = pca(object, "t.test"),
####                       sam = pca(object, "sam"),
####                       d = dPosteriorMean(object))
####              }
##            x <- getStats(object, statistic[1])
##            y <- getStats(object, statistic[2])
##            if(!is.null(filename)) postscript(filename)
##            par(las = 1, mar = mar, oma = oma)
##            qqplot(x, y, xlab = xlab, ylab = ylab)
##            abline(0, 1)
##            if(!is.null(filename)) dev.off()
##          })

##pca <- function(object, x, ...){
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
            require(genefilter) || stop("package genefilter not available")
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
              xx.c <- rowMeans(exprs(x))
              exprs(x) <- sweep(exprs(x), 1, xx.c)
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
          function(object, ...){
            object <- geneCenter(object)
            
            ##Note that by zeroing the nu's we must specify initial
            ##values and can not draw random samples from the prior
            ##(could result in illegal values)
            
            params <- new("XdeParameter", phenotypeLabel=phenotypeLabel, expressionSetList=object)
            firstMcmc <- empiricalStart(object, zeroNu=TRUE)
            firstMcmc$A <- rep(0, length(object))
            firstMcmc$Rho <- rep(0, choose(length(object), 2))
            up <- updates(params)
            up["nu"] <- 0
            up["a"] <- 0
            up["gamma2"] <- 0
            up["rhoAndGamma2"] <- 0
            updates(params) <- up
            firstMcmc(params) <- firstMcmc
            params
          })

