library(Test)
##trace(initializeParams, browser)
initializeParams(G=50L, nIt=1L, Q=3L, S=c(10L, 10L, 10L), seed=1L)

	## #################################################################################
	##str(oldClique)
	##str(unlist(oldComponents))
	##str(nOld)
	##str(nTotalInClique)
	## Omega
##	OmegaM = diag(rep(0,G))
##	Si = NULL
##	Omega = NULL
##	CliquesF <- vector("list", nClique)
##	for (i in 1:(nClique)){
##		Ci = which(clique==i-1)
##		cli = oldClique[i]
##		if (cli>=0){
##			list1 = oldComponents[[i]]
##			CiS = CliquesF[[cli+1]]
##			Si = CiS[list1]
##			CliquesF[[i]] = sort(c(Ci,Si))
##		}
##		else {
##			CliquesF[[i]] = Ci
##		}
##		setC = c(Ci,Si)
##		OmegaM[setC,setC] = 2
##		Om = OmegaM[setC,setC]
##		Omega = c(Omega,c(Om[lower.tri(Om,diag=TRUE)]))
##	}
##	Omega[is.na(Omega)] <- 2

## add arguments for number genes, etc. so that these are not hardcoded
##
## have updated parameters written to files
##
## Can only pass 65 arguments
## TODO
##
## - change initializeParams such that G, Q, psi (clinical
##   parameters), etc. are inputs
##
## - have an option that selects which parameters are
##   written to file
##
## - modify such that mcmc samples are appended to files according to
##  previous option (see how sigma2 and phi are currently written to
##  file, and adopt this for other parameters)
##
## - modify such that only selected iterations are appended to file
##
## - modify such that data is passed to the function and not simulated
## - (or have an option that allows passing data)
##
## - many of the specific update types are already available with R wrappers
##   (see R/Rupdates.R)
##   -> initializeParams could be changed to only provide initial parameters
##   -> using .C, we would need to specify the correct size for each parameter as input
##      to initializeParams function
##




