updateMissing <- function(orig, current){
	## the parameter list after an Rupdate type
	## must have all the parameters required for
	## a different update type.
	nms <- setdiff(names(orig), names(current))
	for(i in seq_along(nms)){
		eval(substitute(current$NAME_ARG, list(NAME_ARG=nms[i])))
	}
	return(current)
}

rupdateANu <- function(object,
		       nTry=length(object),
		       nAccept=0L,
		       epsilon=0.2,
		       dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    a=object[["a"]],
			    nu=object[["nu"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    gamma2=object[["gamma2"]],
			    rho=object[["rho"]],
			    sigma2=object[["sigma2"]],
			    phi=object[["phi"]],
			    tau2Rho=object[["tau2Rho"]],
			    pA0=object[["pA0"]],
			    pA1=object[["pA1"]],
			    alphaA=object[["alphaA"]],
			    betaA=object[["betaA"]])
		return(res)
	}
	res <- .C("updateANu",
		   seed=object[["seed"]],
		   nTry=as.integer(nTry),
		   nAccept=as.integer(nAccept),
		   epsilon=as.numeric(epsilon),
		   a=object[["a"]],
		   nu=object[["nu"]],
		   Q=object[["Q"]],
		   G=object[["G"]],
		   S=object[["S"]],
		   x=exprs(object),
		   psi=object[["phenodata"]],
		   delta=object[["delta"]],
		   Delta=object[["Delta"]],
		   gamma2=object[["gamma2"]],
		   rho=object[["rho"]],
		   sigma2=object[["sigma2"]],
		   phi=object[["phi"]],
		   tau2Rho=object[["tau2Rho"]],
		   pA0=object[["pA0"]],
		   pA1=object[["pA1"]],
		   alphaA=object[["alphaA"]],
		   betaA=object[["betaA"]])
	## update a and nu
	slot(object, "seed") <- res[["seed"]]
	slot(object, "a") <- res[["a"]]
	slot(object, "nu") <- res[["nu"]]
	return(object)
}

rupdateBDDelta <- function(object,
			   nTry=length(object),
			   nAccept=0L,
			   epsilon=0.2,
			   dryrun=FALSE){
	if(dryrun){
		obj <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    tau2Rho=object[["tau2Rho"]],
			    b=object[["b"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    nu=object[["nu"]],
			    delta=object[["delta"]],
			    c2=object[["c2"]],
			    r=object[["r"]],
			    sigma2=object[["sigma2"]],
			    phi=object[["phi"]],
			    tau2R=object[["tau2R"]],
			    pB0=object[["pB0"]],
			    pB1=object[["pB1"]],
			    alphaB=object[["alphaB"]],
			    betaB=object[["betaB"]])
		return(obj)
	}
	res <- .C("updateBDDelta",
		  seed=object[["seed"]],
		  nTry=as.integer(nTry),
		  nAccept=as.integer(nAccept),
		  epsilon=as.numeric(epsilon),
		  b=object[["b"]],
		  Delta=object[["Delta"]],
		  Q=object[["Q"]],
		  G=object[["G"]],
		  S=object[["S"]],
		  x=exprs(object),
		  psi=object[["phenodata"]],
		  nu=object[["nu"]],
		  delta=object[["delta"]],
		  c2=object[["c2"]],
		  r=object[["r"]],
		  sigma2=object[["sigma2"]],
		  phi=object[["phi"]],
		  tau2R=object[["tau2R"]],
		  pB0=object[["pB0"]],
		  pB1=object[["pB1"]],
		  alphaB=object[["alphaB"]],
		  betaB=object[["betaB"]])
	slot(object, "seed") <- res[["seed"]]
	slot(object, "b") <- res[["b"]]
	slot(object, "Delta") <- res[["Delta"]]
	return(object)
}




