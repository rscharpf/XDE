library(Test)
set.seed(1)
source("Test/R/initializeParams.R")
S <- c(100L, 50L, 75L)
G <- 100L
x <- rnorm(G*sum(S))
cliques <- simulateCliques(G=100L)
clique.params <- cliqueParams(cliques)
initializeParams(nIt=2L, G=100L,
		 x=x,
		 clique=cliques,
		 clique.params=clique.params,
		 a=c(0.0, 0.5, 1.0),
		 b=c(0.0, 0.5, 1.0),
		 l=rep(1.0,3),
		 t=rep(0.5^2,3),
		 writeout=rep(1L, 20))
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





