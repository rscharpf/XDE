library(Test)
## add arguments for number genes, etc. so that these are not hardcoded
##
## have updated parameters written to files
##
G <- 3592
#########################################
# Not sure this is the right position
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

nClique <- max(clique) + 1
nOldClique <- rep(0,nClique)
for (i in 2:nClique){
  nOldClique[i] = length(oldComponents[[i]])
}
#################################################################################
str(oldClique)
str(unlist(oldComponents))
str(nOld)
str(nTotalInClique)
## Can only pass 65 arguments
tmp <- .C("initializeParams",
	  nIt=1L,
	  seedR=123L,
	  G=3592L,
	  Q=3L,
	  S=c(100L, 50L, 75L),
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
	  sigma2=rep(0.0, 3592*3), ## G * Q
	  aOut=1L,
	  sigma2Out=1L,
	  simulateSigma2=1L,
	  clique=as.integer(clique),
	  nClique=as.integer(nClique),
	  oldClique=as.integer(oldClique),
	  oldComponents=as.integer(unlist(oldComponents)), ## not sure if this is the right order
	  nOldClique=as.integer(nOldClique), ## length nClique, number of elements of the seperator i
	  nNewClique=as.integer(nNewInClique), ## length nClique, the number of elements of P_i\S_i
	  OmegaInput=as.numeric(nTotalInClique))

tmp <- list(
	  nIt=1L,
	  seedR=123L,
	  G=3592L,
	  Q=3L,
	  S=c(100L, 50L, 75L),
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
	  sigma2=rep(0.0, 3592*3), ## G * Q
	  aOut=1L,
	  sigma2Out=1L,
	  simulateSigma2=1L,
	  clique=as.integer(clique),
	  nClique=as.integer(nClique),
	  oldClique=as.integer(oldClique),
	  oldComponents=as.integer(unlist(oldComponents)), ## not sure if this is the right order
	  nOldClique=as.integer(nOldClique), ## length nClique, number of elements of the seperator i
	  nNewClique=as.integer(nNewInClique), ## length nClique, the number of elements of P_i\S_i
	  OmegaInput=as.numeric(nTotalInClique))
#	  clique,
#	  nNewInClique,
#	  nTotalInClique)

## only last iteration written to file...
phi <- as.matrix(read.table("phi.txt", sep=" ", header=FALSE,
			    colClasses=c("NULL", "NULL", rep("numeric", 3))))
## need to figure out whether this is row-major or column-major order
phi.array <- array(phi, dim=c(G, 3, 3))
sigma2 <- as.matrix(read.table("sigma2.txt", sep=" ", header=FALSE,
			    colClasses=c("NULL", "numeric")))
sigma2 <- matrix(sigma2, nrow=G)
## expression data (?)
x <- as.matrix(read.table("x.txt", sep=" ", header=FALSE,
	       colClasses=c("NULL", "NULL", "NULL", "numeric")))
## just the first study?
x <- matrix(x, 3592, 100)
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




