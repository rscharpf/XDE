library(Test)
## add arguments for number genes, etc. so that these are not hardcoded
##
## have updated parameters written to files
##
G <- 3592
tmp <- .C("initializeParams",
	  nIt=1L,
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
	  aOut=1L,
	  sigma2Out=1L)
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




