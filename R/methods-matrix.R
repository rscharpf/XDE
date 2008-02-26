setMethod("pairs", "matrix",
          function(x, ...){
            op <- par(no.readonly=TRUE)
            on.exit(par(op))
            require(geneplotter) || stop("Package geneplotter is not available")
            require(RColorBrewer) || stop("Package RColorBrewer is not available")
            par(mfrow=c(1,1), las=1)
            panel <- function(...){
              par(new="TRUE")
              colramp <- colorRampPalette(c("white", brewer.pal(9, "Greys")))              
              smoothScatter(..., colramp=colramp, nrpoints=0, xaxt="n", yaxt="n")
              x <- list(...)
              r <- round(cor(x[[1]], x[[2]], method="spearman"), 2)
              legend("topleft", legend=as.character(r), bty="n", cex=1.4)              
            }
            graphics:::pairs(x, panel=panel, upper.panel=NULL, ...)            
          })

##N is a vector of sample sizes
setMethod(".pca", "matrix",
          function(object, N=NULL, ...){
            if(length(N) != ncol(object)) stop("The sample size in each study must be specified")
            ##1. obtain the effect sizes in the P studies
            ##   - these are stored in 'object'
            ##2. center the effect sizes by their means
            meancen <- colMeans(object)
            xcen <- sweep(object, 2, meancen)
  
            ##3. Fit principal components to w1, w2, and w3 using
            ##   covariance, saving a1, a2 and a3, the loadings of the first
            ##   principal component.
            fit <- princomp(xcen, cor = FALSE)
            a <- fit$loadings[1:ncol(object), 1]
  
            ##4. Estimate the combined effect: (see p.9 of Garrett-Mayer
            ##   paper)
            weight <- sum(a * sqrt(N))
            combined.effect <- matrix(NA, ncol=ncol(object), nrow=nrow(object))
            for(i in 1:ncol(object)){
              combined.effect[,i] <- a[i] * sqrt(N[i]) * xcen[, i]
            }
            score <- apply(combined.effect, 1, sum)
            score <- score/weight
            return(score)
          })
