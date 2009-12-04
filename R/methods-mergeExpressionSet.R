##object is a LinearXdeSet or loglinearXdeSet
setMethod("integrativeCorrelationFilter", "mergeExpressionSet",
          function(object, fdrCut, ...){
		  ##require(gtools) || stop("gtools package not available")
		  ##require(MergeMaid) || stop("MergeMaid package not available")
		  if(missing(fdrCut)) stop("must specify fdrCut (FDR threshold)")
		  emp <- pairwise.cors(intCor(object))
		  ##This is the approximation to the null (does fine if you have
		  ##thousands of genes)
		  message("approximating the null distribution for the integrative correlation...")
		  null <- intcorDens(object)
		  null <- do.call(cbind, null)            
		  ##Find the cutoff, c,  at which
		  ##(NULL_c)/(NULL_c + EMP_c) = fdrcut
		  ##Lowering the fdrCut applies a more stringent integrative correlation filter
		  findCuts <- function(x1, x2, probs=seq(0.6, 0.99, by=0.005), fdrCut){
			  k <- 1
			  quantiles <- quantile(x2, probs)
			  fdrhat <- 1
			  while(fdrhat > fdrCut & k <= length(quantiles)){
				  p.null <- sum(x1 > quantiles[k])/length(x1)
				  p.emp <- sum(x2 > quantiles[k])/length(x2)
				  fdrhat <- p.null/(p.null + p.emp)
				  ##print(fdrhat)
				  if(is.na(fdrhat)) { warning("0 in fdr denominator"); break()}
				  if(k+1 < length(quantiles)) k <- k+1 else break()
			  }
			  quantiles[k]
		  }
		  N <- ncol(emp)
		  cut <- rep(NA, N)  ##a threshold for each pairwise combination of studies
		  for(i in 1:N)   cut[i] <- findCuts(null[, i], emp[, i], fdrCut=fdrCut)
		  print(paste("integrative correlation coefficient threshold:", round(cut, 2)))
		  
		  ##Define an indicator matrix for each gene x study-pair
		  cuts <- matrix(NA, nrow=nrow(emp), ncol=ncol(emp))
		  for(i in 1:N){
			  cuts[, i] <- emp[, i] > cut[i]
		  }
		  return(cuts)
          })
