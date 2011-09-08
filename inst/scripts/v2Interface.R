## mainTest_v2.cpp contains contains code Haakon used for testing v2
##
## Rinterface_v2.cpp contains all of the entry points from R to update each of the parameters.
##


##---------------------------------------------------------------------------
##
## initialize fixed hyerparameters
##
##---------------------------------------------------------------------------
G <- 1000            ## number of genes
Q <- 3               ## number of studies
S <- c(100, 50, 75)  ## study sizes

##
betaA <- alphaA <- 1
pA0 <- pA1 <- 0.1

##
betaB <- alphaB <- 1
pA0 <- pB0 <- 0.1

nuR <- 1+Q
nuRho <- 1+Q

betaXi <- alphaXi <- 1

c2Max <- 10

##---------------------------------------------------------------------------
## initialize clinical variables
##---------------------------------------------------------------------------
psi <- list(rbinom(S[1], 1, 0.5),
	    rbinom(S[2], 1, 0.5),
	    rbinom(S[3], 1, 0.5))

## psi is ordered by study, so
psi <- unlist(psi)

##---------------------------------------------------------------------------
##
## initialize parameters to be simulated
##
##---------------------------------------------------------------------------
gamma2 <- 0.5*0.5
c2 <- 0.5*0.5
tau2Rho <- rep(1, Q)
tau2R <- rep(1,Q)
a <- b <- rep(NA, Q)
l <- rep(1, Q)
t <- rep(0.5^2, Q)
sigma2 <- matrix(NA, G, Q)
## simulate variances from gamma
param2 <- l/t
param1 <- l*param2
for(q in seq(length=Q)){
	u <- runif(1)
	if(u < pA0){
		a[q] <- 0
	} else {
		if(u < pA0+pA1){
			a[q] <- 1
		} else {
			a[q] <-runif(1)
		}
	}
	u <- runif(1)
	if(u < pB0){
		b[q] <- 0
	} else {
		if(u < pB0+pB1){
			b[q] <- 1
		} else {
			b[q] <-runif(1)
		}
	}
	## ?
##	a[1] <- 0
##	a[2] <- 0.5
##	a[3] <- 1
##	b[1] <- 1
##	b[2] <- 0
##	b[3] <- 0.5
	## in R, mean is shape * scale and variance is shape*scale^2
	## in Haakon's formulation, mean = l and variance = t
	## => shape=l^2/t and scale =t/l
	##  (this means that param2 is the rate and param1 is the shape
	sigma2[, q] <- rgamma(G, shape=param1[q], rate=param2[q])
	## the platform index changes the fastest, then gene
	##
	## lambda
	##
	## theta
	##
	## phi
	##
	## r
	##
	## nu
	##
	## delta
	##
	## Delta
}
sigma2 <- as.numeric(t(sigma2))
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







