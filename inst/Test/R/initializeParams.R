initializeParams <- function(G=3592L, nIt=1L, Q=3L, S=c(100L, 50L, 75L), seed=123L, x, ...){
	simulateExpression <- if(missing(x)) TRUE else stop("nothing here yet")
	## #########################################
	## Not sure this is the right position
	clique <- rep(0,G)
	rU <- runif(G)
	for (g in 2:G){
		if (rU[g]<0.95){
			clique[g] <- clique[g-1]}
		else clique[g] <- clique[g-1]+1
	}
	nNewInClique <- tabulate(clique+1)
	oldComponents <- vector("list", length(nNewInClique))
	oldClique <- nNewInClique
	oldClique[1] <- -1
	nTotalInClique <- nNewInClique
	rU <- runif(length(nNewInClique)-1)
	rU2 <- runif(length(nNewInClique)-1)
	for (c in 2:length(nNewInClique)){
		oldClique[c] =  round(rU[c-1]*(c-2))
		nOld = round(rU2[c-1]*rU2[c-1]*(nTotalInClique[oldClique[c]+1]))
		nTotalInClique[c] = nTotalInClique[c] + nOld
		oldComponents[[c]] = rep(0,nOld)
		nSampled = 0
		while(nSampled < nOld) {
			oldComponents[[c]][nSampled+1] = round((runif(1)*(nTotalInClique[oldClique[c]+1]-1))+1)
			existAlready = 0
			if (nSampled>0){
				for (i in 1:(nSampled)){
					if (oldComponents[[c]][i]==oldComponents[[c]][nSampled+1]){
						existAlready = 1
					}
				}
			}
			if (existAlready==0) {nSampled=nSampled+1}
		}
	}
	nClique <- max(clique) + 1
	nOldClique <- rep(0,nClique)
	for (i in 2:nClique){
		nOldClique[i] = length(oldComponents[[i]])
	}
	## Omega
	OmegaM = diag(rep(0,G))
	Si = NULL
	Omega = NULL
	CliquesF <- vector("list", nClique)
	for (i in 1:(nClique)){
		Ci = which(clique==i-1)
		cli = oldClique[i]
		if (cli>=0){
			list1 = oldComponents[[i]]
			CiS = CliquesF[[cli+1]]
			Si = CiS[list1]
			CliquesF[[i]] = sort(c(Ci,Si))
		} else {
			CliquesF[[i]] = Ci
		}
		setC = c(Ci,Si)
		OmegaM[setC,setC] = 2
		Om = OmegaM[setC,setC]
		Omega = c(Omega,c(Om[lower.tri(Om,diag=TRUE)]))
	}
	Omega[is.na(Omega)] <- 2
	gamma2 <- 0.5^2
	tau2Rho <- rep(1.0, Q)
	psi <- as.integer(runif(sum(S)) > 0.5)
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
		  nTotalInClique=as.integer(nTotalInClique),
		  nNewClique=as.integer(nNewInClique),
		  OmegaInput=as.numeric(Omega))
}
