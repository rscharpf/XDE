## mainTest_v2.cpp contains contains code Haakon used for testing v2
##
## Rinterface_v2.cpp contains all of the entry points from R to update each of the parameters.
##


##---------------------------------------------------------------------------
##
## initialize fixed hyerparameters
##
##---------------------------------------------------------------------------
library(XDE)
hyper.params <- getHyperparameters(G=3592, Q=3, S=c(100,50,75))

##---------------------------------------------------------------------------
## initialize clinical variables
##---------------------------------------------------------------------------
S <- hyper.params[["S"]]
psi <- c(rbinom(S[1], 1, 0.5),
	 rbinom(S[2], 1, 0.5),
	 rbinom(S[3], 1, 0.5))

##---------------------------------------------------------------------------
##
## initialize parameters to be simulated
##
##---------------------------------------------------------------------------
trace(getParameters, browser)
tmp <- getParameters(hyper.params)

##
## simulate data from model
##

##---------------------------------------------------------------------------
##
##  run Metropolis-Hastings updates
##
##---------------------------------------------------------------------------


##---------------------------------------------------------------------------
##
##  initialize parameters to be simulated (to random values)
##
##---------------------------------------------------------------------------


##---------------------------------------------------------------------------
##
##  parameter updates
##
##---------------------------------------------------------------------------

## update potential
##...
##

## begin mcmc
## update nu
seed <- 123L
nu <- .C("updateANu", seed=seed,
	 nTry=Q
	 nAccept=0
	 epsilonANu=0.1,
	 a=a,
	 nu=nu,
	 Q=Q
	 G=G,
	 S=S,
	 x=x,
	 psi=psi,
	 delta=delta,
	 Delta=Delta,
	 gamma2=gamma2,
	 rho=rho,
	 sigma2=sigma2,
	 phi=phi,
	 tau2Rho=tau2Rho,
	 pA0=pA0,
	 pA1=pA1,
	 alphaA=alphaA,
	 betaA=betaA)







