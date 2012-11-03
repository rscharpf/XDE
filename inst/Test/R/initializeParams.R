initializeParams <- function(G=3592L, nIt=1L, Q=3L, S=c(100L, 50L, 75L), seed=123L, x, ...){
	simulateExpression <- if(missing(x)) TRUE else FALSE
	## #########################################
	## Not sure this is the right position
	clique <- rep(0,G)
	rU <- runif(G)
	for (g in 2:G){
		if (rU[g]<0.95){
			clique[g] <- clique[g-1]}
		else clique[g] <- clique[g-1]+1
	}
	if(length(unique(clique)) < 2) stop("must have at least 2 cliques")
	nNewInClique <- tabulate(clique+1)
	oldComponents <- vector("list", length(nNewInClique))
	oldClique <- nNewInClique
	oldClique[1] <- -1
	nTotalInClique <- nNewInClique
	rU <- runif(length(nNewInClique)-1)
	rU2 <- runif(length(nNewInClique)-1)
	for (c in 2:length(nNewInClique)){
		oldClique[c] =  round(rU[c-1]*(c-1))
		nOld = round(rU2[c-1]*rU2[c-1]*(nTotalInClique[oldClique[c]+1]))
		nTotalInClique[c] = nTotalInClique[c] + nOld
		oldComponents[[c]] = rep(0,nOld)
		nSampled = 0
		while(nSampled < nOld) {
			oldComponents[[c]][nSampled+1] = round(runif(1)*nTotalInClique[oldClique[c]+1])
			existAlready = 0
			if (nSampled>1){
				for (i in 1:(nSampled)){
					if (oldComponents[[c]][i]==oldComponents[[c]][nSampled+1]){
						existAlready = 1
					}
				}
			}
			if (existAlready==0) {nSampled=nSampled+1}
		}
	}
	##
	nClique <- max(clique) + 1
	nOldClique <- rep(0,nClique)
	for (i in 2:nClique){
		nOldClique[i] = length(oldComponents[[i]])
	}
	gamma2 <- 0.5^2
	tau2Rho <- rep(1.0, Q)
	psi <- as.integer(runif(sum(S)) > 0.5)
	if(!simulateExpression) stop("C function not set up to handle input expression data yet")
	tmp <- .C("initializeParams",
		  nIt=nIt,
		  seedR=as.integer(seed),
		  G=as.integer(G),
		  Q=as.integer(Q),
		  S=as.integer(S),
		  psi=as.integer(psi),
		  alphaA=1.0,
		  betaA=1.0,
		  pA0=0.1,
		  pA1=0.1,
		  alphaB=1.0,
		  betaB=1.0,
		  pB0=0.1,
		  pB1=0.1,
		  nuR=1.0,
		  nuRho=1.0,
		  alphaXi=1.0,
		  betaXi=1.0,
		  c2Max=10.0,
		  sigma2=rep(0.0, G*Q), ## G * Q
		  tau2Rho=as.numeric(tau2Rho),
		  gamma2=as.numeric(gamma2),
		  tau2R=as.numeric(tau2Rho),
		  simulateExpression=simulateExpression,
		  aOut=1L,
		  sigma2Out=1L,
		  simulateSigma2=1L,
		  oldCliqueInput=as.integer(oldClique),
		  oldComponentsInput=as.integer(unlist(oldComponents)),
		  nClique=as.integer(nClique),
		  nOldClique=as.integer(nOldClique),
		  nTotalInClique=as.integer(nTotalInClique))
##		  clique=as.integer(clique),
##		  nClique=as.integer(nClique),
##		  oldClique=as.integer(oldClique),
##		  oldComponents=as.integer(unlist(oldComponents)), ## not sure if this is the right order
##		  nOldClique=as.integer(nOldClique), ## length nClique, number of elements of the seperator i
##		  nNewClique=as.integer(nNewInClique), ## length nClique, the number of elements of P_i\S_i
##		  nTotalInClique=as.integer(nTotalInClique))
		  ##OmegaInput=Omega)
}
