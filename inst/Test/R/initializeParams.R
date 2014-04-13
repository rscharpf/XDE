simulateCliques <- function(G){
	u <- runif(G)
	clique <- rep(0L,G)
	for(g in 2:G){
		if(u[g] < 0.95){
			clique[g] <- clique[g-1]
		} else clique[g] <- clique[g-1]+1
	}
	as.integer(clique)
}

cliqueParams <- function(clique){
	##
	## nClique: integer (number of cliques)
	##
	## oldClique: vector of integers with length nClique
	##            - first element of oldClique is arbitrary
	##            - element i specifies for which clique j < i  P_i intersect P_j is not the empty set
	##            - C++ assumes the cliques to be numbered from zero to nClique-1
	##
	## nOldClique: vector of integers with length nClique
	##            - element i specifies the number of elements in separator i
	##            - the first element should always be zero
	##
	## nNewClique: vector of integers with length nClique
	##            - element i specifies number of elements of P_i\S_i
	##
	## oldComponents: vector of integers with length equal to the
	##                sum of all elements in the vector nOldClique
	##              - the first |S_1| elements specify the elements in S_1
	##              - the next |S_2| elements specify the elements in S_2
	##              - ...
	##              - assumes a local numbering of genes inside each clique from zero to |P_i|
	##
	## Toy example:
	##G <- 10
	## P_0={0,1,2,3}, P_1 = {2,4,5,6}, P_2={4,5,7,8}, P_3={4,6,9}
	## S_0={}, S_1=P_1 intersect P_0 = {2}, S_2 = P_2 intersect P_1 = {4,5}, S_3=P_3 intersect P_1={4,6}
	## then
##	 nClique=4
##	 oldClique=c(-1,0,1,1)
##	 nOldClique=c(0,1,2,2)
##	 nNewInClique=c(4,3,2,1)
##	 oldComponents=c(2,1,2,1,3)
	nNewInClique <- tabulate(clique+1)
	oldComponents <- vector("list", length(nNewInClique))
	oldClique <- nNewInClique
	oldClique[1] <- -1
	nTotalInClique <- nNewInClique
	u1 <- runif(length(nNewInClique)-1)
	u2 <- runif(length(nNewInClique)-1)
	u <- u1*u2
	for (c in 2:length(nNewInClique)){
		oldClique[c] <- round(u1[c-1]*(c-2))
		nOld <- round(u[c-1]*(nTotalInClique[oldClique[c]+1]))
		nTotalInClique[c] <- nTotalInClique[c] + nOld
		##oldComponents[[c]] <- if(nOld > 0) rep(0, nOld) else 0
		nSampled <- 0
		while(nSampled < nOld) {
			oldComponents[[c]][nSampled+1] <- round((runif(1)*(nTotalInClique[oldClique[c]+1]-1))+1)
			existAlready <- 0
			if (nSampled>0){
				for (i in 1:(nSampled)){
					if (oldComponents[[c]][i]==oldComponents[[c]][nSampled+1]){
						existAlready <- 1
					}
				}
			}
			if (existAlready==0) nSampled <- nSampled+1
		}
	}
	nClique <- max(clique) + 1
	nOldClique <- rep(0,nClique)
	for (i in 2:(nClique)){
		nOldClique[i] = length(oldComponents[[i]])
	}
	list(oldClique=as.integer(oldClique),
	     oldComponents=oldComponents,
	     nClique=as.integer(nClique),
	     nOldClique=as.integer(nOldClique),
	     nTotalInClique=as.integer(nTotalInClique),
	     nNewInClique=as.integer(nNewInClique))
}

simulateOmega <- function(clique, clique.params, G){
	nClique <- clique.params[["nClique"]]
	oldComponents <- clique.params[["oldComponents"]]
	oldClique <- clique.params[["oldClique"]]
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
	Omega
}


initializeParams <- function(G=3592L, nIt=1L, Q=3L,
			     S=c(100L, 50L, 75L),
			     clique,
			     clique.params,
			     psi,
			     seed=123L,
			     a=rep(0.5, Q),
			     b=rep(0.5, Q),
			     l=rep(1.0, Q),
			     t=rep(0.5*0.5, Q),
			     xi=rep(0.5, Q),
			     lambda=rep(0.9, Q),
			     theta=rep(0.1^2, Q),
			     x=as.numeric(0),
			     gamma2=0.5^2,
			     tau2Rho=rep(1.0,Q),
			     writeout=rep(1L, 20),
			     ...){
	simulateExpression <- if(length(x)<2) TRUE else FALSE
	x <- numeric(G*sum(S))
	if(missing(clique))
		clique <- simulateCliques(G)
	if(missing(clique.params))
		clique.params <- cliqueParams(clique)
	Omega <- simulateOmega(clique,clique.params, G)
	if(missing(psi)){
		message("simulating clinical variable")
		psi <- as.integer(runif(sum(S)) > 0.5)
	}
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
		  aIn=as.numeric(a),
		  bIn=as.numeric(b),
		  lIn=as.numeric(l),
		  tIn=as.numeric(t),
		  lambdaIn=as.numeric(lambda),
		  thetaIn=as.numeric(theta),
		  xiIn=as.numeric(xi),
		  simulateExpression=as.integer(simulateExpression),
		  x=as.numeric(x),
		  writeout=as.integer(writeout),
		  simulateSigma2=1L,
		  oldCliqueInput=clique.params[["oldClique"]],
		  oldComponentsInput=as.integer(unlist(clique.params[["oldComponents"]])),
		  nClique=as.integer(clique.params[["nClique"]]),
		  nOldClique=as.integer(clique.params[["nOldClique"]]),
		  nTotalInClique=as.integer(clique.params[["nTotalInClique"]]),
		  clique=as.integer(clique),
		  nNewInCliqueInput=as.integer(clique.params[["nNewInClique"]]))
}
