library(XDE)
data(expressionSetList)
##hyper.params <- XDE:::getHyperparameters(G=3592, Q=3, S=c(100,50,75))
hyper.params <- HyperParams(object=expressionSetList, clinicalCovariate="adenoVsquamous")
params <- Params(hyper.params)

## note that some values have changed, but nAccept is 'zero'
## Also, the seed has not changed
trace(rupdateANu, browser)
res <- rupdateANu(expressionSetList,
		  hyper.params,
		  params,
		  nTry=100L,
		  epsilon=0.01)
if(FALSE){
nus <- sapply(expressionSetList, function(object) rowMeans(exprs(object)))
rho <- cor(nus)
vars <- apply(nus,2, var)
##nus <- XDE:::vectorize(nus)
##rho <- rho[upper.tri(rho)]
gamma2 <- mean(vars)
psi <- phenotype(expressionSetList, "adenoVsquamous")
d <- matrix(NA, nrow(expressionSetList), length(expressionSetList))
for(i in seq_along(expressionSetList)){
	p0 <- which(psi[,i]==0)
	mn0 <- rowMeans(exprs(expressionSetList[[i]])[, p0])
	mn1 <- rowMeans(exprs(expressionSetList[[i]])[, -p0])
	d[,i] <- mn1-mn0
}
##Delta <- unlist(d)
Delta <- d
library(genefilter)
tt <- rowttests(expressionSetList, "adenoVsquamous")
delta <- abs(tt) > 2.4
delta <- matrix(as.integer(delta), nrow(delta), ncol(delta))
##delta <- as.integer(delta)
params$delta <- delta
params$Delta <- Delta
params$gamma2 <- gamma2
params$nu <- nus
params$r <- params$rho <- rho
}
res <- rupdateBDDelta(expressionSetList,
		     hyper.params,
		     params,
		     nTry=100L,
		     epsilon=.2)
















