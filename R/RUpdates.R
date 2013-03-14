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
		       nTry=10L,
		       nAccept=0L,
		       epsilon=0.2,
		       dryrun=FALSE){
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
			   nTry=10L,
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

	res <- .C("updateBDDelta_MI",
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
	return(object) ##nipun: changed updateBDDelta to updateBDDelta_MI
}

rupdateTau2RhoNu <- function(object,
			     nTry=10L,
			     nAccept=0L,
			     epsilon=0.2,
			     dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    tau2Rho=object[["tau2Rho"]],
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
			    a=object[["a"]])
		return(res)
	}
	res <- .C("updateTau2RhoNu",
		  seed=object[["seed"]],
		  nTry=nTry,
		  nAccept=nAccept,
		  epsilon=epsilon,
		  tau2Rho=object[["tau2Rho"]],
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
		  a=object[["a"]])
	## update tau2Rho and nu
	slot(object, "seed") <- res[["seed"]]
	slot(object, "tau2Rho") <- res[["tau2Rho"]]
	slot(object, "nu") <- res[["nu"]]
	return(object)
}


rupdateTau2RDDelta <- function(object,
			       nTry=10L,
			       nAccept=0L,
			       epsilon=0.2,
			       dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    tau2R=object[["tau2R"]],
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
			    b=object[["b"]])
		return(res)
	}
	res <- .C("updateTau2RDDelta_MI",
		  seed=object[["seed"]],
		  nTry=nTry,
		  nAccept=nAccept,
		  epsilon=epsilon,
		  tau2R=object[["tau2R"]],
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
		  b=object[["b"]])
	## update tau2R and Delta
	slot(object, "seed") <- res[["seed"]]
	slot(object, "tau2R") <- res[["tau2R"]]
	slot(object, "Delta") <- res[["Delta"]]
	return(object)
}

rupdateNu <- function(object,
		      nAccept=0L,
		      dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nAccept=nAccept,
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
			    a=object[["a"]])
		return(res)
	}
	res <- .C("updateNu",
		  seed=object[["seed"]],
		  nAccept=nAccept,
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
		  a=object[["a"]])
	slot(object, "seed") <- res[["seed"]]
	slot(object, "nu") <- res[["nu"]]
	return(object)
}

rupdateDDelta <- function(object,
			  nAccept=0L,
			  dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nAccept=nAccept,
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
			    b=object[["b"]])
		return(res)
	}
	res <- .C("updateDDelta_MI",
		  seed=object[["seed"]],
		  nAccept=nAccept,
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
			    b=object[["b"]])
	slot(object, "seed") <- res[["seed"]]
	slot(object, "Delta") <- res[["Delta"]]
	return(object) ##nipun: changed updateDDelta to updateDDelta_MI
}

rupdateC2 <- function(object,
		      nTry=10L,
		      nAccept=0L,
		      dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    c2=object[["c2"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    r=object[["r"]],
			    sigma2=object[["sigma2"]],
			    tau2R=object[["tau2R"]],
			    b=object[["b"]],
			    c2Max=object[["c2Max"]])
		return(res)
	}
	res <- .C("updateC2_MI",
		  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    c2=object[["c2"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    r=object[["r"]],
			    sigma2=object[["sigma2"]],
			    tau2R=object[["tau2R"]],
			    b=object[["b"]],
			    c2Max=object[["c2Max"]])
	slot(object, "seed") <- res[["seed"]]
	slot(object, "c2") <- res[["c2"]]
	return(object) ##nipun: changed updateC2 to updateC2_MI
}

rupdateC2DDelta <- function(object,
			    nTry=10L,
			    nAccept=0L,
			    epsilon=0.2,
			    dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    c2=object[["c2"]],
			    Delta=object[["Delta"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    nu=object[["nu"]],
			    delta=object[["delta"]],
			    r=object[["r"]],
			    sigma2=object[["sigma2"]],
			    phi=object[["phi"]],
			    tau2R=object[["tau2R"]],
			    b=object[["b"]],
			    c2Max=object[["c2Max"]])
		return(res)
	}
	res <- .C("updateC2DDelta_MI",
		  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    c2=object[["c2"]],
			    Delta=object[["Delta"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    nu=object[["nu"]],
			    delta=object[["delta"]],
			    r=object[["r"]],
			    sigma2=object[["sigma2"]],
			    phi=object[["phi"]],
			    tau2R=object[["tau2R"]],
			    b=object[["b"]],
			    c2Max=object[["c2Max"]])
	slot(object, "seed") <- res[["seed"]]
	slot(object, "c2") <- res[["c2"]]
	slot(object, "Delta") <- res[["Delta"]]
	return(object) ##nipun: changed updateC2DDelta to updateC2DDelta_MI
}


rupdateGamma2 <- function(object,
			  nAccept=0L,
			  dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nAccept=nAccept,
			    gamma2=object[["gamma2"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    nu=object[["nu"]],
			    rho=object[["rho"]],
			    sigma2=object[["sigma2"]],
			    tau2Rho=object[["tau2Rho"]],
			    a=object[["a"]])
		return(res)
	}
	res <- .C("updateGamma2",
		  seed=object[["seed"]],
			    nAccept=nAccept,
			    gamma2=object[["gamma2"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    nu=object[["nu"]],
			    rho=object[["rho"]],
			    sigma2=object[["sigma2"]],
			    tau2Rho=object[["tau2Rho"]],
			    a=object[["a"]])
	slot(object, "seed") <- res[["seed"]]
	slot(object, "gamma2") <- res[["gamma2"]]
	return(object)
}

rupdateGamma2Nu <- function(object,
			    nTry=10L,
			    nAccept=0L,
			    epsilon=0.2,
			    dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    gamma2=object[["gamma2"]],
			    nu=object[["nu"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    rho=object[["rho"]],
			    sigma2=object[["sigma2"]],
			    phi=object[["phi"]],
			    tau2Rho=object[["tau2Rho"]],
			    a=object[["a"]])
		return(res)
	}
		res <- .C("updateGamma2Nu_MI",
			  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    gamma2=object[["gamma2"]],
			    nu=object[["nu"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    rho=object[["rho"]],
			    sigma2=object[["sigma2"]],
			    phi=object[["phi"]],
			    tau2Rho=object[["tau2Rho"]],
			    a=object[["a"]])
	slot(object, "seed") <- res[["seed"]]
	slot(object, "gamma2") <- res[["gamma2"]]
	slot(object, "nu") <- res[["nu"]]
	return(object) ##nipun: changed updateGamma2Nu to updateGamma2Nu_MI
}

rupdateRDDelta <- function(object,
			   nTry=10L,
			   nAccept=0L,
			   epsilon=0.2,
			   dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    r=object[["r"]],
			    Delta=object[["Delta"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    nu=object[["nu"]],
			    delta=object[["delta"]],
			    c2=object[["c2"]],
			    sigma2=object[["sigma2"]],
			    phi=object[["phi"]],
			    tau2R=object[["tau2R"]],
			    b=object[["b"]],
			    nuR=object[["nuR"]])
		return(res)
	}
		res <- .C("updateRDDelta_MI",
			  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    r=object[["r"]],
			    Delta=object[["Delta"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    nu=object[["nu"]],
			    delta=object[["delta"]],
			    c2=object[["c2"]],
			    sigma2=object[["sigma2"]],
			    phi=object[["phi"]],
			    tau2R=object[["tau2R"]],
			    b=object[["b"]],
			    nuR=object[["nuR"]])
	slot(object, "seed") <- res[["seed"]]
	slot(object, "r") <- res[["r"]]
	slot(object, "Delta") <- res[["Delta"]]
	return(object)
}

rupdateRC2 <- function(object,
		       nTry=10L,
		       nAccept=0L,
		       epsilon=0.2,
		       dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    r=object[["r"]],
			    c2=object[["c2"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    sigma2=object[["sigma2"]],
			    tau2R=object[["tau2R"]],
			    b=object[["b"]],
			    nuR=object[["nuR"]],
			    c2Max=object[["c2Max"]])
		return(res)
	}
		res <- .C("updateRC2_MI",
			  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    r=object[["r"]],
			    c2=object[["c2"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    sigma2=object[["sigma2"]],
			    tau2R=object[["tau2R"]],
			    b=object[["b"]],
			    nuR=object[["nuR"]],
			    c2Max=object[["c2Max"]])
	slot(object, "seed") <- res[["seed"]]
	slot(object, "r") <- res[["r"]]
	slot(object, "c2") <- res[["c2"]]
	return(object)
}

rupdateRhoNu <- function(object,
			 nTry=10L,
			 nAccept=0L,
			 epsilon=0.2,
			 dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    rho=object[["rho"]],
			    nu=object[["nu"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    gamma2=object[["gamma2"]],
			    sigma2=object[["sigma2"]],
			    phi=object[["phi"]],
			    tau2Rho=object[["tau2Rho"]],
			    a=object[["a"]],
			    nuRho=object[["nuRho"]])
		return(res)
	}
		res <- .C("updateRhoNu",
			  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    rho=object[["rho"]],
			    nu=object[["nu"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    gamma2=object[["gamma2"]],
			    sigma2=object[["sigma2"]],
			    phi=object[["phi"]],
			    tau2Rho=object[["tau2Rho"]],
			    a=object[["a"]],
			    nuRho=object[["nuRho"]])
	slot(object, "seed") <- res[["seed"]]
	slot(object, "rho") <- res[["rho"]]
	slot(object, "nu") <- res[["nu"]]
	return(object)
}


rupdateRhoGamma2 <- function(object,
			 nTry=10L,
			 nAccept=0L,
			 epsilon=0.2,
			 dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    rho=object[["rho"]],
			    gamma2=object[["gamma2"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    nu=object[["nu"]],
			    sigma2=object[["sigma2"]],
			    tau2Rho=object[["tau2Rho"]],
			    a=object[["a"]],
			    nuRho=object[["nuRho"]])
		return(res)
	}
		res <- .C("updateRhoGamma2",
			  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    rho=object[["rho"]],
			    gamma2=object[["gamma2"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    nu=object[["nu"]],
			    sigma2=object[["sigma2"]],
			    tau2Rho=object[["tau2Rho"]],
			    a=object[["a"]],
			    nuRho=object[["nuRho"]])
	slot(object, "seed") <- res[["seed"]]
	slot(object, "rho") <- res[["rho"]]
	slot(object, "gamma2") <- res[["gamma2"]]
	return(object)
}

rupdateSigma2 <- function(object,
			  nTry=10L,
			  nAccept=0L,
			  epsilon=0.2,
			  dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    sigma2=object[["sigma2"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    nu=object[["nu"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    c2=object[["c2"]],
			    gamma2=object[["gamma2"]],
			    r=object[["r"]],
			    rho=object[["rho"]],
			    phi=object[["phi"]],
			    t=object[["t"]],
			    l=object[["l"]],
			    tau2R=object[["tau2R"]],
			    tau2Rho=object[["tau2Rho"]],
			    a=object[["a"]],
			    b=object[["b"]])
		return(res)
	}
		res <- .C("updateSigma2_MI",
			  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    sigma2=object[["sigma2"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    nu=object[["nu"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    c2=object[["c2"]],
			    gamma2=object[["gamma2"]],
			    r=object[["r"]],
			    rho=object[["rho"]],
			    phi=object[["phi"]],
			    t=object[["t"]],
			    l=object[["l"]],
			    tau2R=object[["tau2R"]],
			    tau2Rho=object[["tau2Rho"]],
			    a=object[["a"]],
			    b=object[["b"]])
	slot(object, "seed") <- res[["seed"]]
	slot(object, "sigma2") <- res[["sigma2"]]
	return(object)
}

rupdatePhi <- function(object,
			  nTry=10L,
			  nAccept=0L,
			  epsilon=0.2,
			  dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    phi=object[["phi"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    nu=object[["nu"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    sigma2=object[["sigma2"]],
			    theta=object[["theta"]],
			    lambda=object[["lambda"]])
		return(res)
	}
		res <- .C("updatePhi",
			  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    phi=object[["phi"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    nu=object[["nu"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    sigma2=object[["sigma2"]],
			    theta=object[["theta"]],
			    lambda=object[["lambda"]])
	slot(object, "seed") <- res[["seed"]]
	slot(object, "phi") <- res[["phi"]]
	return(object)
}

rupdateTheta <- function(object,
			 nTry=10L,
			 nAccept=0L,
			 epsilon=0.2,
			 dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    theta=object[["theta"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    phi=object[["phi"]],
			    lambda=object[["lambda"]])
		return(res)
	}
		res <- .C("updateTheta",
			  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    theta=object[["theta"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    phi=object[["phi"]],
			    lambda=object[["lambda"]])
	slot(object, "seed") <- res[["seed"]]
	slot(object, "theta") <- res[["theta"]]
	return(object)
}

rupdateLambda <- function(object,
			  nTry=10L,
			  nAccept=0L,
			  epsilon=0.2,
			  dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    lambda=object[["lambda"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    phi=object[["phi"]],
			    theta=object[["theta"]])
		return(res)
	}
		res <- .C("updateLambda",
			  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    lambda=object[["lambda"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    phi=object[["phi"]],
			    theta=object[["theta"]])
	slot(object, "seed") <- res[["seed"]]
	slot(object, "lambda") <- res[["lambda"]]
	return(object)
}

rupdateT <- function(object,
		     nTry=10L,
		     nAccept=0L,
		     epsilon=0.2,
		     dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    t=object[["t"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    sigma2=object[["sigma2"]],
			    l=object[["l"]])
		return(res)
	}
		res <- .C("updateT",
			  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    t=object[["t"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    sigma2=object[["sigma2"]],
			    l=object[["l"]])
	slot(object, "seed") <- res[["seed"]]
	slot(object, "t") <- res[["t"]]
	return(object)
}

rupdateL <- function(object,
		     nTry=10L,
		     nAccept=0L,
		     epsilon=0.2,
		     dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    l=object[["l"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    sigma2=object[["sigma2"]],
			    t=object[["t"]])
		return(res)
	}
		res <- .C("updateL",
			  seed=object[["seed"]],
			  nTry=nTry,
			  nAccept=nAccept,
			  epsilon=epsilon,
			  l=object[["l"]],
			  Q=object[["Q"]],
			  G=object[["G"]],
			  sigma2=object[["sigma2"]],
			  t=object[["t"]])
	slot(object, "seed") <- res[["seed"]]
	slot(object, "l") <- res[["l"]]
	return(object)
}

rupdateXi <- function(object,
		      nAccept=0L,
		      dryrun=FALSE,
		      one.delta=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nAccept=nAccept,
			    xi=object[["xi"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    delta=object[["delta"]],
			    t=object[["t"]],
			    alphaXi=object[["alphaXi"]],
			    betaXi=object[["betaXi"]])
		return(res)
	}
	if(one.delta){
		res <- .C("updateXi_MB",
			  seed=object[["seed"]],
			    nAccept=nAccept,
			    xi=object[["xi"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    delta=object[["delta"]],
			    t=object[["t"]],
			    alphaXi=object[["alphaXi"]],
			    betaXi=object[["betaXi"]])		##nipun: changed updateXi_onedelta to updateXi_MB

	} else {
		res <- .C("updateXi_MA",
			  seed=object[["seed"]],
			    nAccept=nAccept,
			    xi=object[["xi"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    delta=object[["delta"]],
			    t=object[["t"]],
			    alphaXi=object[["alphaXi"]],
			    betaXi=object[["betaXi"]])  ##nipun: changed updateXi to updateXi_MA
	}
	slot(object, "seed") <- res[["seed"]]
	## why do we need to do this?
	slot(object, "delta") <- res[["delta"]]
	return(object)
}

rupdateDeltaDDelta <- function(object,
			       nTry=10L,
			       nAccept=0L,
			       dryrun=FALSE,
			       one.delta=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    nu=object[["nu"]],
			    c2=object[["c2"]],
			    r=object[["r"]],
			    sigma2=object[["sigma2"]],
			    phi=object[["phi"]],
			    tau2R=object[["tau2R"]],
			    xi=object[["xi"]],
			    b=object[["b"]])
		return(res)
	}
	if(one.delta){
		res <- .C("updateDeltaDDelta_MBI",
			  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    nu=object[["nu"]],
			    c2=object[["c2"]],
			    r=object[["r"]],
			    sigma2=object[["sigma2"]],
			    phi=object[["phi"]],
			    tau2R=object[["tau2R"]],
			    xi=object[["xi"]],
			    b=object[["b"]]) ##nipun: changed updateDeltaDDelta_onedelta to updateDeltaDDelta_MBI
	} else {
		res <- .C("updateDeltaDDelta_MAI",
			  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    nu=object[["nu"]],
			    c2=object[["c2"]],
			    r=object[["r"]],
			    sigma2=object[["sigma2"]],
			    phi=object[["phi"]],
			    tau2R=object[["tau2R"]],
			    xi=object[["xi"]],
			    b=object[["b"]]) ##nipun: changed updateDeltaDDelta to updateDeltaDDelta_MAI
	}
	slot(object, "seed") <- res[["seed"]]
	slot(object, "delta") <- res[["delta"]]
	slot(object, "Delta") <- res[["Delta"]]
	return(object)
}


##nipun: added rupdateDeltaDDelta_MRF2
rupdateDeltaDDelta_MRF2 <- function(object,
                               nTry=10L,
                               nAccept=0L,
                               dryrun=FALSE,
                               one.delta){
  if(dryrun){
    res <- list(seed=object[["seed"]],
                nTry=nTry,
                nAccept=nAccept,
                delta=object[["delta"]],
                Delta=object[["Delta"]],
                Q=object[["Q"]],
                G=object[["G"]],
                S=object[["S"]],
                x=exprs(object),
                psi=object[["phenodata"]],
                nu=object[["nu"]],
                c2=object[["c2"]],
                r=object[["r"]],
                sigma2=object[["sigma2"]],
                phi=object[["phi"]],
                tau2R=object[["tau2R"]],
                xi=object[["xi"]],
                b=object[["b"]])
    return(res)
  }
  if(one.delta){
    res <- .C("updateDeltaDDelta_MDI",
              seed=object[["seed"]],
              nTry=nTry,
              nAccept=nAccept,
              delta=object[["delta"]],
              Delta=object[["Delta"]],
              Q=object[["Q"]],
              G=object[["G"]],
              S=object[["S"]],
              x=exprs(object),
              psi=object[["phenodata"]],
              nu=object[["nu"]],
              c2=object[["c2"]],
              r=object[["r"]],
              sigma2=object[["sigma2"]],
              phi=object[["phi"]],
              tau2R=object[["tau2R"]],
              b=object[["b"]],
              nNeighbour=object[["nNeighbour"]],
              neighbour=object[["neighbour"]],
              alpha = object[["alpha"]],
              beta=object[["beta"]])   ##nipun: changed updateDeltaDDelta_MRF2_onedelta to updateDeltaDDelta_MDI
  } else {
    res <- .C("updateDeltaDDelta_MCI",
              seed=object[["seed"]],
              nTry=nTry,
              nAccept=nAccept,
              delta=object[["delta"]],
              Delta=object[["Delta"]],
              Q=object[["Q"]],
              G=object[["G"]],
              S=object[["S"]],
              x=exprs(object),
              psi=object[["phenodata"]],
              nu=object[["nu"]],
              c2=object[["c2"]],
              r=object[["r"]],
              sigma2=object[["sigma2"]],
              phi=object[["phi"]],
              tau2R=object[["tau2R"]],
              b=object[["b"]],
              nNeighbour=object[["nNeighbour"]],
              neighbour=object[["neighbour"]],
              alpha = object[["alpha"]], ##working
              beta=object[["beta"]],
              betag=object[["betag"]]) ##nipun: changed updateDeltaDDelta_MRF2 to updateDeltaDDelta_MCI
  }
  slot(object, "seed") <- res[["seed"]]
  slot(object, "delta") <- res[["delta"]]
  slot(object, "Delta") <- res[["Delta"]]
  return(object)
}


rupdateLSigma2 <- function(object,
			   nTry=10L,
			   nAccept=0L,
			   epsilon=0.2,
			   dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    l=object[["l"]],
			    sigma2=object[["sigma2"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    nu=object[["nu"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    c2=object[["c2"]],
			    gamma2=object[["gamma2"]],
			    r=object[["r"]],
			    rho=object[["rho"]],
			    phi=object[["phi"]],
			    t=object[["t"]],
			    tau2R=object[["tau2R"]],
			    tau2Rho=object[["tau2Rho"]],
			    a=object[["a"]],
			    b=object[["b"]])
		return(res)
	}
		res <- .C("updateLSigma2_MI",
			  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    l=object[["l"]],
			    sigma2=object[["sigma2"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    nu=object[["nu"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    c2=object[["c2"]],
			    gamma2=object[["gamma2"]],
			    r=object[["r"]],
			    rho=object[["rho"]],
			    phi=object[["phi"]],
			    t=object[["t"]],
			    tau2R=object[["tau2R"]],
			    tau2Rho=object[["tau2Rho"]],
			    a=object[["a"]],
			    b=object[["b"]])
	slot(object, "seed") <- res[["seed"]]
	slot(object, "l") <- res[["l"]]
	slot(object, "sigma2") <- res[["sigma2"]]
	return(object) ##nipun: Changed updateLSigma2 to updateLSigma2_MI
}

rupdateTSigma2 <- function(object,
			   nTry=10L,
			   nAccept=0L,
			   epsilon=0.2,
			   dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    t=object[["t"]],
			    sigma2=object[["sigma2"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    nu=object[["nu"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    c2=object[["c2"]],
			    gamma2=object[["gamma2"]],
			    r=object[["r"]],
			    rho=object[["rho"]],
			    phi=object[["phi"]],
			    l=object[["l"]],
			    tau2R=object[["tau2R"]],
			    tau2Rho=object[["tau2Rho"]],
			    a=object[["a"]],
			    b=object[["b"]])
		return(res)
	}
		res <- .C("updateTSigma2_MI",
			  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    t=object[["t"]],
			    sigma2=object[["sigma2"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    nu=object[["nu"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    c2=object[["c2"]],
			    gamma2=object[["gamma2"]],
			    r=object[["r"]],
			    rho=object[["rho"]],
			    phi=object[["phi"]],
			    l=object[["l"]],
			    tau2R=object[["tau2R"]],
			    tau2Rho=object[["tau2Rho"]],
			    a=object[["a"]],
			    b=object[["b"]])   ##nipun: Changed updateTSigma2 to updateTSigma2_MI
	slot(object, "seed") <- res[["seed"]]
	slot(object, "t") <- res[["t"]]
	slot(object, "sigma2") <- res[["sigma2"]]
	return(object)
}

rupdateLambdaPhi <- function(object,
			     nTry=10L,
			     nAccept=0L,
			     epsilon=0.2,
			     dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    lambda=object[["lambda"]],
			    phi=object[["phi"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
  			    x=exprs(object),
			    psi=object[["phenodata"]],
			    nu=object[["nu"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    sigma2=object[["sigma2"]],
			    theta=object[["theta"]])
		return(res)
	}
		res <- .C("updateLambdaPhi",
			  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    lambda=object[["lambda"]],
			    phi=object[["phi"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    nu=object[["nu"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    sigma2=object[["sigma2"]],
			    theta=object[["theta"]])
	slot(object, "seed") <- res[["seed"]]
	slot(object, "lambda") <- res[["lambda"]]
	slot(object, "phi") <- res[["phi"]]
	return(object)
}

rupdateThetaPhi <- function(object,
			     nTry=10L,
			     nAccept=0L,
			     epsilon=0.2,
			     dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    theta=object[["theta"]],
			    phi=object[["phi"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    nu=object[["nu"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    sigma2=object[["sigma2"]],
			    lambda=object[["lambda"]])
		return(res)
	}
		res <- .C("updateThetaPhi",
			  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    theta=object[["theta"]],
			    phi=object[["phi"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    nu=object[["nu"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    sigma2=object[["sigma2"]],
			    lambda=object[["lambda"]])
	slot(object, "seed") <- res[["seed"]]
	slot(object, "theta") <- res[["theta"]]
	slot(object, "phi") <- res[["phi"]]
	return(object)
}

##rupdateDeltaDDelta_MRF <- function(object,
##				   nTry=10L,
##				   nAccept=0L,
##			   dryrun=FALSE,
##				   one.delta=FALSE){
##	if(dryrun){
##	res <- list(seed=object[["seed"]],
##		    nTry=nTry,
##		    nAccept=nAccept,
##		    delta=object[["delta"]],
##		    Delta=object[["Delta"]],
##		    Q=object[["Q"]],
##		    G=object[["G"]],
##		    S=object[["S"]],
##		    x=exprs(object),
##		    psi=object[["phenodata"]],
##		    nu=object[["nu"]],
##		    c2=object[["c2"]],
##		    r=object[["r"]],
##		    sigma2=object[["sigma2"]],
##		    phi=object[["phi"]],
##		    tau2R=object[["tau2R"]],
##		    b=object[["b"]],
##		    nNeighbour=object[["nNeighbour"]],
##		    neighbour=object[["neighbour"]],
##		    alpha=object[["alpha"]],
##		    beta=object[["beta"]],
##		    betag=object[["betag"]])
##	if(one.delta) res <- res[-length(res)]
##	return(res)
##}
##if(one.delta){
##	res <- .C("updateDeltaDDelta_MRF_onedelta",
##		  seed=object[["seed"]],
##		    nTry=nTry,
##		    nAccept=nAccept,
##		    delta=object[["delta"]],
##		    Delta=object[["Delta"]],
##		    Q=object[["Q"]],
##		    G=object[["G"]],
##		    S=object[["S"]],
##		    x=exprs(object),
##		    psi=object[["phenodata"]],
##		    nu=object[["nu"]],
##		    c2=object[["c2"]],
##		    r=object[["r"]],
##		    sigma2=object[["sigma2"]],
##		    phi=object[["phi"]],
##		    tau2R=object[["tau2R"]],
##		    b=object[["b"]],
##		    nNeighbour=object[["nNeighbour"]],
##		    neighbour=object[["neighbour"]],
##		    alpha=object[["alpha"]],
##		    beta=object[["beta"]])##nipun:updateDeltaDDelta_MRF_onedelta missing
##} else {
##	res <- .C("updateDeltaDDelta_MRF",
##		  seed=object[["seed"]],
##		    nTry=nTry,
##		    nAccept=nAccept,
##		    delta=object[["delta"]],
##		    Delta=object[["Delta"]],
##		    Q=object[["Q"]],
##		    G=object[["G"]],
##		    S=object[["S"]],
##		    x=exprs(object),
##		    psi=object[["phenodata"]],
##		    nu=object[["nu"]],
##		    c2=object[["c2"]],
##		    r=object[["r"]],
##		    sigma2=object[["sigma2"]],
##		    phi=object[["phi"]],
##		    tau2R=object[["tau2R"]],
##		    b=object[["b"]],
##		    nNeighbour=object[["nNeighbour"]],
##		    neighbour=object[["neighbour"]],
##		    alpha=object[["alpha"]],
##		    beta=object[["beta"]],
##		  betag=object[["betag"]]) ##nipun:updateDeltaDDelta_MRF missing
##}
##slot(object, "seed") <- res[["seed"]]
##slot(object, "delta") <- res[["delta"]]
##slot(object, "Delta") <- res[["Delta"]]
##return(object)
##}





##rupdateAlpha_MRF <- function(object,
##			     nTry=10L,
##		     nAccept=0L,
##		     epsilon=0.2,
##		     dryrun=FALSE,
##		     one.delta=FALSE){
##if(dryrun){
##	res <- list(seed=object[["seed"]],
##		    nTry=nTry,
##		    nAccept=nAccept,
##		    epsilon=epsilon,
##		    alpha=object[["alpha"]],
##		    Q=object[["Q"]],
##		    G=object[["G"]],
##		    delta=object[["delta"]],
##			    nNeighbour=object[["nNeighbour"]],
##			    neighbour=object[["neighbour"]],
##		    beta=object[["beta"]],
##		    betag=object[["betag"]])
##	if(one.delta) res <- res[-length(res)]
##	return(res)
##}
##if(one.delta){
##	res <- .C("updateAlpha_MRF_onedelta",
##		  seed=object[["seed"]],
##		    nTry=nTry,
##		    nAccept=nAccept,
##		    epsilon=epsilon,
##		    alpha=object[["alpha"]],
##		    Q=object[["Q"]],
##		    G=object[["G"]],
##		    delta=object[["delta"]],
##		    nNeighbour=object[["nNeighbour"]],
##		    neighbour=object[["neighbour"]],
##		    beta=object[["beta"]])
##} else {
##	res <- .C("updateAlpha_MRF",
##		  seed=object[["seed"]],
##		    nTry=nTry,
##		    nAccept=nAccept,
##		    epsilon=epsilon,
##		    alpha=object[["alpha"]],
##		    Q=object[["Q"]],
##		    G=object[["G"]],
##		    delta=object[["delta"]],
##		    nNeighbour=object[["nNeighbour"]],
##		    neighbour=object[["neighbour"]],
##		    beta=object[["beta"]],
##		    betag=object[["betag"]])
##}
##slot(object, "seed") <- res[["seed"]]
##slot(object, "alpha") <- res[["alpha"]]
##return(object)
##}
##
##rupdateBeta_MRF <- function(object,
##		     nTry=10L,
##		     nAccept=0L,
##		     epsilon=0.2,
##		     dryrun=FALSE,
##		     one.delta=FALSE){
##if(dryrun){
##	res <- list(seed=object[["seed"]],
##		    nTry=nTry,
##		    nAccept=nAccept,
##		    epsilon=epsilon,
##		    beta=object[["beta"]],
##		    Q=object[["Q"]],
##		    G=object[["G"]],
##		    delta=object[["delta"]],
##		    nNeighbour=object[["nNeighbour"]],
##		    neighbour=object[["neighbour"]],
##		    alpha=object[["alpha"]],
##		    betag=object[["betag"]])
##	if(one.delta) res <- res[-length(res)]
##	return(res)
##}
##if(one.delta){
##	res <- .C("updateBeta_MRF_onedelta",
##		  seed=object[["seed"]],
##		    nTry=nTry,
##		    nAccept=nAccept,
##		    epsilon=epsilon,
##		    beta=object[["beta"]],
##		    Q=object[["Q"]],
##		    G=object[["G"]],
##		    delta=object[["delta"]],
##		    nNeighbour=object[["nNeighbour"]],
##		    neighbour=object[["neighbour"]],
##		    alpha=object[["alpha"]])
##} else {
##	res <- .C("updateBeta_MRF",
##		  seed=object[["seed"]],
##		    nTry=nTry,
##		    nAccept=nAccept,
##		    epsilon=epsilon,
##		    beta=object[["beta"]],
##		    Q=object[["Q"]],
##		    G=object[["G"]],
##		    delta=object[["delta"]],
##		    nNeighbour=object[["nNeighbour"]],
##		    neighbour=object[["neighbour"]],
##		    alpha=object[["alpha"]],
##		    betag=object[["betag"]])
##}
##slot(object, "seed") <- res[["seed"]]
##slot(object, "beta") <- res[["beta"]]
##return(object)
##}
##

##rupdateBetag_MRF <- function(object,
##		     nTry=10L,
##		     nAccept=0L,
##		     epsilon=0.2,
##		     dryrun=FALSE){
##if(dryrun){
##	res <- list(seed=object[["seed"]],
##		    nTry=nTry,
##		    nAccept=nAccept,
##		    epsilon=epsilon,
##		    betag=object[["betag"]],
##		    Q=object[["Q"]],
##		    G=object[["G"]],
##		    delta=object[["delta"]],
##		    nNeighbour=object[["nNeighbour"]],
##		    neighbour=object[["neighbour"]],
##		    alpha=object[["alpha"]],
##		    beta=object[["beta"]])
##	return(res)
##}
##	res <- .C("updateBetag_MRF",
##		  seed=object[["seed"]],
##		    nTry=nTry,
##		    nAccept=nAccept,
##		    epsilon=epsilon,
##		    betag=object[["betag"]],
##		    Q=object[["Q"]],
##		    G=object[["G"]],
##		    delta=object[["delta"]],
##		    nNeighbour=object[["nNeighbour"]],
##		    neighbour=object[["neighbour"]],
##		    alpha=object[["alpha"]],
##		    beta=object[["beta"]])
##slot(object, "seed") <- res[["seed"]]
##slot(object, "betag") <- res[["betag"]]
##	return(object)
##}

##nipun: wrappers for the remaining functions

rupdateAlpha <- function(object,
                     nTry=10L,
                     nAccept=0L,
                     epsilon=0.2,
                     dryrun=FALSE){
  #if(dryrun){
  # res <- list(seed=object[["seed"]],
  #              nTry=nTry,
  #              nAccept=nAccept,
  #              epsilon=epsilon,
  #              alpha=object[["alpha"]],
  #              Q=object[["Q"]],
  #              G=object[["G"]],
  #              delta=object[["delta"]],
  #              nNeighbour=object[["nNeighbour"]],
  #              neighbour=object[["neighbour"]],
  #              beta=object[["beta"]],
  #              betag=object[["betag"]])
  #  return(res)
  #}
  res <- .C("updateAlpha_MC",
            seed=object[["seed"]],
            nTry=nTry,
            nAccept=nAccept,
            epsilon=epsilon,
            alpha=object[["alpha"]],
            Q=object[["Q"]],
            G=object[["G"]],
            delta=object[["delta"]],
            nNeighbour=object[["nNeighbour"]],
            neighbour=object[["neighbour"]],
            beta=object[["beta"]],
            betag=object[["betag"]])
  slot(object, "seed") <- res[["seed"]]
  slot(object, "alpha") <- res[["alpha"]]
  return(object)
}

rupdateAlpha_onedelta <- function(object,
                         nTry=10L,
                         nAccept=0L,
                         epsilon=0.2,
                         dryrun=FALSE){
  if(dryrun){
    res <- list(seed=object[["seed"]],
                nTry=nTry,
                nAccept=nAccept,
                epsilon=epsilon,
                alpha=object[["alpha"]],
                Q=object[["Q"]],
                G=object[["G"]],
                delta=object[["delta"]],
                nNeighbour=object[["nNeighbour"]],
                neighbour=object[["neighbour"]],
                beta=object[["beta"]])
    return(res)
  }
  res <- .C("updateAlpha_MD",
            seed=object[["seed"]],
            nTry=nTry,
            nAccept=nAccept,
            epsilon=epsilon,
            alpha=object[["alpha"]],
            Q=object[["Q"]],
            G=object[["G"]],
            delta=object[["delta"]],
            nNeighbour=object[["nNeighbour"]],
            neighbour=object[["neighbour"]],
            beta=object[["beta"]])
  slot(object, "seed") <- res[["seed"]]
  slot(object, "alpha") <- res[["alpha"]]
  return(object)
}

rupdateBeta <- function(object,
                         nTry=10L,
                         nAccept=0L,
                         epsilon=0.2,
                         dryrun=FALSE){
  if(dryrun){
    res <- list(seed=object[["seed"]],
                nTry=nTry,
                nAccept=nAccept,
                epsilon=epsilon,
                beta=object[["beta"]],
                Q=object[["Q"]],
                G=object[["G"]],
                delta=object[["delta"]],
                nNeighbour=object[["nNeighbour"]],
                neighbour=object[["neighbour"]],
                alpha=object[["alpha"]],
                betag=object[["betag"]])
    return(res)
  }
  res <- .C("updateBeta_MC",
            seed=object[["seed"]],
            nTry=nTry,
            nAccept=nAccept,
            epsilon=epsilon,
            beta=object[["beta"]],
            Q=object[["Q"]],
            G=object[["G"]],
            delta=object[["delta"]],
            nNeighbour=object[["nNeighbour"]],
            neighbour=object[["neighbour"]],
            alpha=object[["alpha"]],
            betag=object[["betag"]])
  slot(object, "seed") <- res[["seed"]]
  slot(object, "beta") <- res[["beta"]]
  return(object)
}

rupdateBeta_onedelta <- function(object,
                        nTry=10L,
                        nAccept=0L,
                        epsilon=0.2,
                        dryrun=FALSE){
  if(dryrun){
    res <- list(seed=object[["seed"]],
                nTry=nTry,
                nAccept=nAccept,
                epsilon=epsilon,
                beta=object[["beta"]],
                Q=object[["Q"]],
                G=object[["G"]],
                delta=object[["delta"]],
                nNeighbour=object[["nNeighbour"]],
                neighbour=object[["neighbour"]],
                alpha=object[["alpha"]])
    return(res)
  }
  res <- .C("updateBeta_MD",
            seed=object[["seed"]],
            nTry=nTry,
            nAccept=nAccept,
            epsilon=epsilon,
            beta=object[["beta"]],
            Q=object[["Q"]],
            G=object[["G"]],
            delta=object[["delta"]],
            nNeighbour=object[["nNeighbour"]],
            neighbour=object[["neighbour"]],
            alpha=object[["alpha"]])
  slot(object, "seed") <- res[["seed"]]
  slot(object, "beta") <- res[["beta"]]
  return(object)
}

rupdateBetag <- function(object,
                     nTry=10L,
                     nAccept=0L,
                     epsilon=0.2,
                     dryrun=FALSE){
  if(dryrun){
    res <- list(seed=object[["seed"]],
                nTry=nTry,
                nAccept=nAccept,
                epsilon=epsilon,
                betag=object[["betag"]],
                Q=object[["Q"]],
                G=object[["G"]],
                delta=object[["delta"]],
                nNeighbour=object[["nNeighbour"]],
                neighbour=object[["neighbour"]],
                alpha=object[["alpha"]],
                beta=object[["beta"]])
    return(res)
  }
  res <- .C("updateBetag_MC",
            seed=object[["seed"]],
            nTry=nTry,
            nAccept=nAccept,
            epsilon=epsilon,
            betag=object[["betag"]],
            Q=object[["Q"]],
            G=object[["G"]],
            delta=object[["delta"]],
            nNeighbour=object[["nNeighbour"]],
            neighbour=object[["neighbour"]],
            alpha=object[["alpha"]],
            beta=object[["beta"]])
  slot(object, "seed") <- res[["seed"]]
  slot(object, "betag") <- res[["betag"]]
  return(object)
}


rupdateDelta <- function(object,
                        nTry=10L,
                        nAccept=0L,
                        dryrun=FALSE,
                        one.delta=FALSE){
  if(dryrun){
    res <- list(seed=object[["seed"]],
                nTry=nTry,
                nAccept=nAccept,
                delta=object[["delta"]],
                Q=object[["Q"]],
                G=object[["G"]],
                S=object[["S"]],
                x=exprs(object),
                psi=object[["phenodata"]],
                nu=object[["nu"]],
                Delta=object[["Delta"]],
                r=object[["r"]],
                sigma2=object[["sigma2"]],
                phi=object[["phi"]],
                xi=object[["xi"]],
                b=object[["b"]])

    if(one.delta) res <- res[-length(res)]
    return(res)
  }
  if(one.delta){
    res <- .C("updateDelta_MBII",
              seed=object[["seed"]],
              nTry=nTry,
              nAccept=nAccept,
              delta=object[["delta"]],
              Q=object[["Q"]],
              G=object[["G"]],
              S=object[["S"]],
              x=exprs(object),
              psi=object[["phenodata"]],
              nu=object[["nu"]],
              Delta=object[["Delta"]],
              r=object[["r"]],
              sigma2=object[["sigma2"]],
              phi=object[["phi"]],
              xi=object[["xi"]],
              b=object[["b"]])
  } else {
    res <- .C("updateDelta_MAII",
              seed=object[["seed"]],
              nTry=nTry,
              nAccept=nAccept,
              delta=object[["delta"]],
              Q=object[["Q"]],
              G=object[["G"]],
              S=object[["S"]],
              x=exprs(object),
              psi=object[["phenodata"]],
              nu=object[["nu"]],
              Delta=object[["Delta"]],
              r=object[["r"]],
              sigma2=object[["sigma2"]],
              phi=object[["phi"]],
              xi=object[["xi"]],
              b=object[["b"]])
  }
  slot(object, "seed") <- res[["seed"]]
  slot(object, "delta") <- res[["delta"]]
  return(object)
}

rupdateDelta_MRF2 <- function(object,
                                   nTry=10L,
                                   nAccept=0L,
                                   dryrun=FALSE,
                                   one.delta){
  if(dryrun){
    res <- list(seed=object[["seed"]],
                nTry=nTry,
                nAccept=nAccept,
                delta=object[["delta"]],
                Q=object[["Q"]],
                G=object[["G"]],
                S=object[["S"]],
                x=exprs(object),
                psi=object[["phenodata"]],
                nu=object[["nu"]],
                Delta=object[["Delta"]],
                r=object[["r"]],
                sigma2=object[["sigma2"]],
                phi=object[["phi"]],
                b=object[["b"]],
                nNeighbour=object[["nNeighbour"]],
                neighbour=object[["neighbour"]],
                alpha=object[["alpha"]],
                beta=object[["beta"]],
                betag=object[["betag"]])
    if(one.delta) res <- res[-length(res)]
    return(res)
  }
  if(one.delta){
    res <- .C("updateDelta_MDII",
              seed=object[["seed"]],
              nTry=nTry,
              nAccept=nAccept,
              delta=object[["delta"]],
              Q=object[["Q"]],
              G=object[["G"]],
              S=object[["S"]],
              x=exprs(object),
              psi=object[["phenodata"]],
              nu=object[["nu"]],
              Delta=object[["Delta"]],
              r=object[["r"]],
              sigma2=object[["sigma2"]],
              phi=object[["phi"]],
              b=object[["b"]],
              nNeighbour=object[["nNeighbour"]],
              neighbour=object[["neighbour"]],
              alpha=object[["alpha"]],
              beta=object[["beta"]])
   }else {
    res <- .C("updateDelta_MCII",
              seed=object[["seed"]],
              nTry=nTry,
              nAccept=nAccept,
              delta=object[["delta"]],
              Q=object[["Q"]],
              G=object[["G"]],
              S=object[["S"]],
              x=exprs(object),
              psi=object[["phenodata"]],
              nu=object[["nu"]],
              Delta=object[["Delta"]],
              r=object[["r"]],
              sigma2=object[["sigma2"]],
              phi=object[["phi"]],
              b=object[["b"]],
              nNeighbour=object[["nNeighbour"]],
              neighbour=object[["neighbour"]],
              alpha=object[["alpha"]],
              beta=object[["beta"]],
              betag=object[["betag"]])
  }
  slot(object, "seed") <- res[["seed"]]
  slot(object, "delta") <- res[["delta"]]
  return(object)
}

rupdateSigma2_HyperInverseWishart <- function(object,
                     nTry=10L,
                     nAccept=0L,
                     epsilon=0.2,
                     dryrun=FALSE){
  if(dryrun){
    res <- list(seed=object[["seed"]],
                nTry=nTry,
                nAccept=nAccept,
                epsilon=epsilon,
                sigma2=object[["sigma2"]],
                Q=object[["Q"]],
                G=object[["G"]],
                S=object[["S"]],
                x=exprs(object),
                psi=object[["phenodata"]],
                nu=object[["nu"]],
                delta=object[["delta"]],
                Delta=object[["Delta"]],
                gamma2=object[["gamma2"]],
                r=object[["r"]],
                rho=object[["rho"]],
                phi=object[["phi"]],
                t=object[["t"]],
                l=object[["l"]],
                tau2R=object[["tau2R"]],
                tau2Rho=object[["tau2Rho"]],
                a=object[["a"]],
                b=object[["b"]],
                Omega=object[["Omega"]],
                nClique=object[["nClique"]],
                oldClique=object[["oldClique"]],
                nOldComponents=object[["nOldComponents"]],
                nNewComponents=object[["nNewComponents"]],
                oldComponents=object[["oldComponents"]])
    return(res)
  }
  res <- .C("updateSigma2_MII",
            seed=object[["seed"]],
            nTry=nTry,
            nAccept=nAccept,
            epsilon=epsilon,
            sigma2=object[["sigma2"]],
            Q=object[["Q"]],
            G=object[["G"]],
            S=object[["S"]],
            x=exprs(object),
            psi=object[["phenodata"]],
            nu=object[["nu"]],
            delta=object[["delta"]],
            Delta=object[["Delta"]],
            gamma2=object[["gamma2"]],
            r=object[["r"]],
            rho=object[["rho"]],
            phi=object[["phi"]],
            t=object[["t"]],
            l=object[["l"]],
            tau2R=object[["tau2R"]],
            tau2Rho=object[["tau2Rho"]],
            a=object[["a"]],
            b=object[["b"]],
            Omega=object[["Omega"]],
            nClique=object[["nClique"]],
            oldClique=object[["oldClique"]],
            nOldComponents=object[["nOldComponents"]],
            nNewComponents=object[["nNewComponents"]],
            oldComponents=object[["oldComponents"]])

  slot(object, "seed") <- res[["seed"]]
  slot(object, "sigma2") <- res[["sigma2"]]
  return(object)
}

rupdateLSigma2_HyperInverseWishart <- function(object,
                     nTry=10L,
                     nAccept=0L,
                     epsilon=0.2,
                     dryrun=FALSE){
  if(dryrun){
    res <- list(seed=object[["seed"]],
                nTry=nTry,
                nAccept=nAccept,
                epsilon=epsilon,
                l=object[["l"]],
                sigma2=object[["sigma2"]],
                Q=object[["Q"]],
                G=object[["G"]],
                S=object[["S"]],
                x=exprs(object),
                psi=object[["phenodata"]],
                nu=object[["nu"]],
                delta=object[["delta"]],
                Delta=object[["Delta"]],
                gamma2=object[["gamma2"]],
                r=object[["r"]],
                rho=object[["rho"]],
                phi=object[["phi"]],
                t=object[["t"]],
                tau2R=object[["tau2R"]],
                tau2Rho=object[["tau2Rho"]],
                a=object[["a"]],
                b=object[["b"]],
                Omega=object[["Omega"]],
                S=object[["S"]],
                nClique=object[["nClique"]],
                oldClique=object[["oldClique"]],
                nOldComponents=object[["nOldComponents"]],
                nNewComponents=object[["nNewComponents"]],
                oldComponents=object[["oldComponents"]])
    return(res)
  }
  res <- .C("updateLSigma2_MII",
            seed=object[["seed"]],
            nTry=nTry,
            nAccept=nAccept,
            epsilon=epsilon,
            l=object[["l"]],
            sigma2=object[["sigma2"]],
            Q=object[["Q"]],
            G=object[["G"]],
            S=object[["S"]],
            x=exprs(object),
            psi=object[["phenodata"]],
            nu=object[["nu"]],
            delta=object[["delta"]],
            Delta=object[["Delta"]],
            gamma2=object[["gamma2"]],
            r=object[["r"]],
            rho=object[["rho"]],
            phi=object[["phi"]],
            t=object[["t"]],
            tau2R=object[["tau2R"]],
            tau2Rho=object[["tau2Rho"]],
            a=object[["a"]],
            b=object[["b"]],
            Omega=object[["Omega"]],
            S=object[["S"]],
            nClique=object[["nClique"]],
            oldClique=object[["oldClique"]],
            nOldComponents=object[["nOldComponents"]],
            nNewComponents=object[["nNewComponents"]],
            oldComponents=object[["oldComponents"]])
  slot(object, "seed") <- res[["seed"]]
  slot(object, "l") <- res[["l"]]
  return(object)
}

rupdateTSigma2_HyperInverseWishart <- function(object,
                     nTry=10L,
                     nAccept=0L,
                     epsilon=0.2,
                     dryrun=FALSE){
  if(dryrun){
    res <- list(seed=object[["seed"]],
                nTry=nTry,
                nAccept=nAccept,
                epsilon=epsilon,
                t=object[["t"]],
                sigma2=object[["sigma2"]],
                Q=object[["Q"]],
                G=object[["G"]],
                S=object[["S"]],
                x=exprs(object),
                psi=object[["phenodata"]],
                nu=object[["nu"]],
                delta=object[["delta"]],
                Delta=object[["Delta"]],
                gamma2=object[["gamma2"]],
                r=object[["r"]],
                rho=object[["rho"]],
                phi=object[["phi"]],
                l=object[["l"]],
                tau2R=object[["tau2R"]],
                tau2Rho=object[["tau2Rho"]],
                a=object[["a"]],
                b=object[["b"]],
                Omega=object[["Omega"]],
                nClique=object[["nClique"]],
                oldClique=object[["oldClique"]],
                nOldComponents=object[["nOldComponents"]],
                nNewComponents=object[["nNewComponents"]],
                oldComponents=object[["oldComponents"]])
    return(res)
  }
  res <- .C("updateTSigma2_MII",
            seed=object[["seed"]],
            nTry=nTry,
            nAccept=nAccept,
            epsilon=epsilon,
            t=object[["t"]],
            sigma2=object[["sigma2"]],
            Q=object[["Q"]],
            G=object[["G"]],
            S=object[["S"]],
            x=exprs(object),
            psi=object[["phenodata"]],
            nu=object[["nu"]],
            delta=object[["delta"]],
            Delta=object[["Delta"]],
            gamma2=object[["gamma2"]],
            r=object[["r"]],
            rho=object[["rho"]],
            phi=object[["phi"]],
            l=object[["l"]],
            tau2R=object[["tau2R"]],
            tau2Rho=object[["tau2Rho"]],
            a=object[["a"]],
            b=object[["b"]],
            Omega=object[["Omega"]],
            nClique=object[["nClique"]],
            oldClique=object[["oldClique"]],
            nOldComponents=object[["nOldComponents"]],
            nNewComponents=object[["nNewComponents"]],
            oldComponents=object[["oldComponents"]])

  slot(object, "seed") <- res[["seed"]]
  slot(object, "t") <- res[["t"]]
  return(object)
}

rupdateOmega_HyperInverseWishart <- function(object,
                     nAccept=0L,
                     dryrun=FALSE){
  if(dryrun){
    res <- list(seed=object[["seed"]],
                nAccept=nAccept,
                Omega=object[["Omega"]],
                Q=object[["Q"]],
                G=object[["G"]],
                Delta=object[["Delta"]],
                r=object[["r"]],
                sigma2=object[["sigma2"]],
                tau2R=object[["tau2R"]],
                b=object[["b"]],
                zeta=object[["zeta"]],
                D=object[["D"]],
                nClique=object[["nClique"]],
                oldClique=object[["oldClique"]],
                nOldComponents=object[["nOldComponents"]],
                nNewComponents=object[["nNewComponents"]],
                oldComponents=object[["oldComponents"]])
    return(res)
  }
  res <- .C("updateOmega_MII",
            seed=object[["seed"]],
            nAccept=nAccept,
            Omega=object[["Omega"]],
            Q=object[["Q"]],
            G=object[["G"]],
            Delta=object[["Delta"]],
            r=object[["r"]],
            sigma2=object[["sigma2"]],
            tau2R=object[["tau2R"]],
            b=object[["b"]],
            zeta=object[["zeta"]],
            D=object[["D"]],
            nClique=object[["nClique"]],
            oldClique=object[["oldClique"]],
            nOldComponents=object[["nOldComponents"]],
            nNewComponents=object[["nNewComponents"]],
            oldComponents=object[["oldComponents"]])
  slot(object, "seed") <- res[["seed"]]
  slot(object, "Omega") <- res[["Omega"]]
  return(object)
}


rupdateDDeltaStar_HyperInverseWishart <- function(object,
                     nAccept=0L,
                     dryrun=FALSE){
  if(dryrun){
    res <- list(seed=object[["seed"]],
                nAccept=nAccept,
                Delta=object[["Delta"]],
                Q=object[["Q"]],
                G=object[["G"]],
                S=object[["S"]],
                x=exprs(object),
                psi=object[["phenodata"]],
                nu=object[["nu"]],
                delta=object[["delta"]],
                r=object[["r"]],
                sigma2=object[["sigma2"]],
                phi=object[["phi"]],
                tau2R=object[["tau2R"]],
                b=object[["b"]],
                Omega=object[["Omega"]],
                nClique=object[["nClique"]],
                oldClique=object[["oldClique"]],
                nOldComponents=object[["nOldComponents"]],
                nNewComponents=object[["nNewComponents"]],
                oldComponents=object[["oldComponents"]])
    return(res)
  }
  res <- .C("updateDeltaStar_MII",
            seed=object[["seed"]],
            nAccept=nAccept,
            Delta=object[["Delta"]],
            Q=object[["Q"]],
            G=object[["G"]],
            S=object[["S"]],
            x=exprs(object),
            psi=object[["phenodata"]],
            nu=object[["nu"]],
            delta=object[["delta"]],
            r=object[["r"]],
            sigma2=object[["sigma2"]],
            phi=object[["phi"]],
            tau2R=object[["tau2R"]],
            b=object[["b"]],
            Omega=object[["Omega"]],
            nClique=object[["nClique"]],
            oldClique=object[["oldClique"]],
            nOldComponents=object[["nOldComponents"]],
            nNewComponents=object[["nNewComponents"]],
            oldComponents=object[["oldComponents"]])

  slot(object, "seed") <- res[["seed"]]
  slot(object, "Delta") <- res[["Delta"]]
  return(object)
}


rupdateTau2RDDeltaStar_HyperInverseWishart <- function(object,
                     nTry=10L,
                     nAccept=0L,
                     epsilon=0.2,
                     dryrun=FALSE){
	phenodata <- object[["phenodata"]]
  if(dryrun){
    res <- list(seed=object[["seed"]],
                nTry=nTry,
                nAccept=nAccept,
                epsilon=epsilon,
                tau2R=object[["tau2R"]],
                Delta=object[["Delta"]],
                Q=object[["Q"]],
                G=object[["G"]],
                S=object[["S"]],
                x=exprs(object),
                psi=object[["phenodata"]],
                nu=object[["nu"]],
                delta=object[["delta"]],
                r=object[["r"]],
                sigma2=object[["sigma2"]],
                phi=object[["phi"]],
                b=object[["b"]],
                Omega=object[["Omega"]],
                nClique=object[["nClique"]],
                oldClique=object[["oldClique"]],
                nOldComponents=object[["nOldComponents"]],
                nNewComponents=object[["nNewComponents"]],
                oldComponents=object[["oldComponents"]])
    return(res)
  }
  res <- .C("updateTau2RDDeltaStar_MII",
            seed=object[["seed"]],
            nTry=nTry,
            nAccept=nAccept,
            epsilon=epsilon,
            tau2R=object[["tau2R"]],
            Delta=object[["Delta"]],
            Q=object[["Q"]],
            G=object[["G"]],
            S=object[["S"]],
            x=exprs(object),
            psi=object[["phenodata"]],
            nu=object[["nu"]],
            delta=object[["delta"]],
            r=object[["r"]],
            sigma2=object[["sigma2"]],
            phi=object[["phi"]],
            b=object[["b"]],
            Omega=object[["Omega"]],
            nClique=object[["nClique"]],
            oldClique=object[["oldClique"]],
            nOldComponents=object[["nOldComponents"]],
            nNewComponents=object[["nNewComponents"]],
            oldComponents=object[["oldComponents"]])
  slot(object, "seed") <- res[["seed"]]
  slot(object, "tau2R") <- res[["tau2R"]]
  return(object)
}


rupdateBDDeltaStar_HyperInverseWishart <- function(object,
                     nTry=10L,
                     nAccept=0L,
                     epsilon=0.2,
                     dryrun=FALSE){
  if(dryrun){
    res <- list(seed=object[["seed"]],
                nTry=nTry,
                nAccept=nAccept,
                epsilon=epsilon,
                b=object[["b"]],
                Delta=object[["Delta"]],
                Q=object[["Q"]],
                G=object[["G"]],
                S=object[["S"]],
                x=exprs(object),
                psi=object[["phenodata"]],
                nu=object[["nu"]],
                delta=object[["delta"]],
                r=object[["r"]],
                sigma2=object[["sigma2"]],
                phi=object[["phi"]],
                tau2R=object[["tau2R"]],
                pB0=object[["pB0"]],
                pB1=object[["pB1"]],
                alphaB=object[["alphaB"]],
                betaB=object[["betaB"]],
                Omega=object[["Omega"]],
                nClique=object[["nClique"]],
                oldClique=object[["oldClique"]],
                nOldComponents=object[["nOldComponents"]],
                nNewComponents=object[["nNewComponents"]],
                oldComponents=object[["oldComponents"]])
    return(res)
  }
  res <- .C("updateBDDeltaStar_MII",
            seed=object[["seed"]],
            nTry=nTry,
            nAccept=nAccept,
            epsilon=epsilon,
            b=object[["b"]],
            Delta=object[["Delta"]],
            Q=object[["Q"]],
            G=object[["G"]],
            S=object[["S"]],
            x=exprs(object),
            psi=phenodata,
            nu=object[["nu"]],
            delta=object[["delta"]],
            r=object[["r"]],
            sigma2=object[["sigma2"]],
            phi=object[["phi"]],
            tau2R=object[["tau2R"]],
            pB0=object[["pB0"]],
            pB1=object[["pB1"]],
            alphaB=object[["alphaB"]],
            betaB=object[["betaB"]],
            Omega=object[["Omega"]],
            nClique=object[["nClique"]],
            oldClique=object[["oldClique"]],
            nOldComponents=object[["nOldComponents"]],
            nNewComponents=object[["nNewComponents"]],
            oldComponents=object[["oldComponents"]])
  slot(object, "seed") <- res[["seed"]]
  slot(object, "b") <- res[["b"]]
  return(object)
}


rupdateRDDeltaStar_HyperInverseWishart <- function(object,
                     nTry=10L,
                     nAccept=0L,
                     epsilon=0.2,
                     dryrun=FALSE){
	phenodata <- object[["phenodata"]]
  if(dryrun){
    res <- list(seed=object[["seed"]],
                nTry=nTry,
                nAccept=nAccept,
                epsilon=epsilon,
                r=object[["r"]],
                Delta=object[["Delta"]],
                Q=object[["Q"]],
                G=object[["G"]],
                S=object[["S"]],
                x=exprs(object),
                psi=object[["phenodata"]],
                nu=object[["nu"]],
                delta=object[["delta"]],
                sigma2=object[["sigma2"]],
                phi=object[["phi"]],
                tau2R=object[["tau2R"]],
                b=object[["b"]],
                nuR=object[["nuR"]],
                Omega=object[["Omega"]],
                nClique=object[["nClique"]],
                oldClique=object[["oldClique"]],
                nOldComponents=object[["nOldComponents"]],
                nNewComponents=object[["nNewComponents"]],
                oldComponents=object[["oldComponents"]])
    return(res)
  }
  res <- .C("updateRDDeltaStar_MII",
            seed=object[["seed"]],
            nTry=nTry,
            nAccept=nAccept,
            epsilon=epsilon,
            r=object[["r"]],
            Delta=object[["Delta"]],
            Q=object[["Q"]],
            G=object[["G"]],
            S=object[["S"]],
            x=exprs(object),
            psi=object[["phenodata"]],
            nu=object[["nu"]],
            delta=object[["delta"]],
            sigma2=object[["sigma2"]],
            phi=object[["phi"]],
            tau2R=object[["tau2R"]],
            b=object[["b"]],
            nuR=object[["nuR"]],
            Omega=object[["Omega"]],
            nClique=object[["nClique"]],
            oldClique=object[["oldClique"]],
            nOldComponents=object[["nOldComponents"]],
            nNewComponents=object[["nNewComponents"]],
            oldComponents=object[["oldComponents"]])
  slot(object, "seed") <- res[["seed"]]
  slot(object, "r") <- res[["r"]]
  return(object)
}



##nipun: common to A, B and MI
modelA_B_MIupdates <- function(object,
			       nupdates = 1500,
			       nTry=10L,
			       nAccept=0L,
			       epsilon=0.2,
			       dryrun=FALSE,
                               one.delta=FALSE){
	G <- object@G; Q <- object@Q
	delta <- array(NA, dim=c(G, Q, nupdates))
	## Q parameters
	lambda <- theta <- xi <- t <- l <- r <- rho <- tau2R <- tau2Rho <- b <- a <- matrix(NA, Q, nupdates)
	## single parameter
	gamma2 <- c2 <- rep(NA, nupdates)
	## hyperparameters
	## alphaXi, betaXi,...
	phenodata <- object@phenodata
	d.i <- object[["delta"]]
	DD.i <- object[["Delta"]]
	S <- object[["S"]]
	for(i in seq_len(nupdates)){
		res <- .C("updateANu",
			  seed=object[["seed"]],
			  nTry=as.integer(nTry),
			  nAccept=as.integer(nAccept),
			  epsilon=as.numeric(epsilon),
			  a=object[["a"]],
			  nu=object[["nu"]],
			  Q=Q,
			  G=G,
			  S=S,
			  x=x,
			  psi=phenodata,
			  delta=d.i,
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
		slot(object, "seed") <- res[["seed"]]
		slot(object, "nu") <- res[["nu"]]
		a[, i] <- slot(object, "a") <- res[["a"]]


		res <- .C("updateTau2RhoNu",
			  seed=object[["seed"]],
			  nTry=nTry,
			  nAccept=nAccept,
			  epsilon=epsilon,
			  tau2Rho=object[["tau2Rho"]],
			  nu=object[["nu"]],
			  Q=object[["Q"]],
			  G=object[["G"]],
			  S=object[["S"]],
			  x=x,
			  psi=phenodata,
			  delta=object[["delta"]],
			  Delta=object[["Delta"]],
			  gamma2=object[["gamma2"]],
			  rho=object[["rho"]],
			  sigma2=object[["sigma2"]],
			  phi=object[["phi"]],
			  a=object[["a"]])
		## update tau2Rho and nu
		slot(object, "seed") <- res[["seed"]]
		tau2Rho[, i] <- slot(object, "tau2Rho") <- res[["tau2Rho"]]
		slot(object, "nu") <- res[["nu"]]


		res <- .C("updateNu",
			  seed=object[["seed"]],
			  nAccept=nAccept,
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
			  a=object[["a"]])
		slot(object, "seed") <- res[["seed"]]
		slot(object, "nu") <- res[["nu"]]

		res <- .C("updateGamma2",
			  seed=object[["seed"]],
			  nAccept=nAccept,
			  gamma2=object[["gamma2"]],
			  Q=object[["Q"]],
			  G=object[["G"]],
			  nu=object[["nu"]],
			  rho=object[["rho"]],
			  sigma2=object[["sigma2"]],
			  tau2Rho=object[["tau2Rho"]],
			  a=object[["a"]])
		slot(object, "seed") <- res[["seed"]]
		gamma2[i] <- slot(object, "gamma2") <- res[["gamma2"]]

		res <- .C("updateRhoNu",
			  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    rho=object[["rho"]],
			    nu=object[["nu"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    gamma2=object[["gamma2"]],
			    sigma2=object[["sigma2"]],
			    phi=object[["phi"]],
			    tau2Rho=object[["tau2Rho"]],
			    a=object[["a"]],
			    nuRho=object[["nuRho"]])
		slot(object, "seed") <- res[["seed"]]
		rho[, i] <- slot(object, "rho") <- res[["rho"]]
		slot(object, "nu") <- res[["nu"]]

		res <- .C("updateRhoGamma2",
			  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    rho=object[["rho"]],
			    gamma2=object[["gamma2"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    nu=object[["nu"]],
			    sigma2=object[["sigma2"]],
			    tau2Rho=object[["tau2Rho"]],
			    a=object[["a"]],
			    nuRho=object[["nuRho"]])
		slot(object, "seed") <- res[["seed"]]
		rho[, i] <- slot(object, "rho") <- res[["rho"]]
		gamma2[i] <- slot(object, "gamma2") <- res[["gamma2"]]

		res <- .C("updatePhi",
			  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    phi=object[["phi"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    nu=object[["nu"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    sigma2=object[["sigma2"]],
			    theta=object[["theta"]],
			    lambda=object[["lambda"]])
		slot(object, "seed") <- res[["seed"]]
		slot(object, "phi") <- res[["phi"]]

		res <- .C("updateTheta",
			  seed=object[["seed"]],
			  nTry=nTry,
			  nAccept=nAccept,
			  epsilon=epsilon,
			  theta=object[["theta"]],
			  Q=object[["Q"]],
			  G=object[["G"]],
			  phi=object[["phi"]],
			  lambda=object[["lambda"]])
		slot(object, "seed") <- res[["seed"]]
		theta[,i] <- slot(object, "theta") <- res[["theta"]]

		res <- .C("updateLambda",
			  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    lambda=object[["lambda"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    phi=object[["phi"]],
			    theta=object[["theta"]])
		slot(object, "seed") <- res[["seed"]]
		lambda[,i] <- slot(object, "lambda") <- res[["lambda"]]

		res <- .C("updateT",
			  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    t=object[["t"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    sigma2=object[["sigma2"]],
			    l=object[["l"]])
		slot(object, "seed") <- res[["seed"]]
		t[,i] <- slot(object, "t") <- res[["t"]]

		res <- .C("updateL",
			  seed=object[["seed"]],
			  nTry=nTry,
			  nAccept=nAccept,
			  epsilon=epsilon,
			  l=object[["l"]],
			  Q=Q,
			  G=G,
			  sigma2=object[["sigma2"]],
			  t=object[["t"]])
		slot(object, "seed") <- res[["seed"]]
		l[,i] <- slot(object, "l") <- res[["l"]]

		res <- .C("updateLambdaPhi",
			  seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    lambda=object[["lambda"]],
			    phi=object[["phi"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    nu=object[["nu"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    sigma2=object[["sigma2"]],
			    theta=object[["theta"]])
		slot(object, "seed") <- res[["seed"]]
		slot(object, "lambda") <- res[["lambda"]]
		slot(object, "phi") <- res[["phi"]]

		if(one.delta){
			##updateXi_MB
			res <- .C("updateXi_MB",
				  seed=object[["seed"]],
				  nAccept=nAccept,
				  xi=object[["xi"]],
				  Q=Q,
				  G=G,
				  delta=d.i,
				  t=object[["t"]],
				  alphaXi=object[["alphaXi"]],
				  betaXi=object[["betaXi"]])		##nipun: changed updateXi_onedelta to updateXi_MB
		} else {
			##updateXi_MA
			res <- .C("updateXi_MA",
				  seed=object[["seed"]],
				  nAccept=nAccept,
				  xi=object[["xi"]],
				  Q=Q,
				  G=G,
				  delta=d.i,
				  t=object[["t"]],
				  alphaXi=object[["alphaXi"]],
				  betaXi=object[["betaXi"]])  ##nipun: changed updateXi to updateXi_MA
		}
		slot(object, "seed") <- res[["seed"]]
		d.i <- res[["delta"]]
		xi[,i] <- slot(object, "xi") <- res[["xi"]]
		##updateBDDelta_MI
		res <- .C("updateBDDelta_MI",
			  seed=object[["seed"]],
			  nTry=as.integer(nTry),
			  nAccept=as.integer(nAccept),
			  epsilon=as.numeric(epsilon),
			  b=object[["b"]],
			  Delta=DD.i,
			  Q=Q,
			  G=G,
			  S=S,
			  x=x,
			  psi=phenodata,
			  nu=object[["nu"]],
			  delta=d.i,
			  c2=object[["c2"]],
			  r=object[["r"]],
			  sigma2=object[["sigma2"]],
			  phi=object[["phi"]],
			  tau2R=object[["tau2R"]],
			  pB0=object[["pB0"]],
			  pB1=object[["pB1"]],
			  alphaB=object[["alphaB"]],
			  betaB=object[["betaB"]])

		b.i <- res[["b"]]
		b[, i] <- b.i
		DD.i <- res[["Delta"]]
		slot(object, "seed") <- res[["seed"]]
		slot(object, "b") <- b.i
		slot(object, "Delta") <- DD.i

		##updateTau2RDDelta_MI
		res <- .C("updateTau2RDDelta_MI",
			  seed=object[["seed"]],
			  nTry=nTry,
			  nAccept=nAccept,
			  epsilon=epsilon,
			  tau2R=object[["tau2R"]],
			  Delta=DD.i,
			  Q=Q,
			  G=G,
			  S=S,
			  x=x,
			  psi=phenodata,
			  nu=object[["nu"]],
			  delta=d.i,
			  c2=object[["c2"]],
			  r=object[["r"]],
			  sigma2=object[["sigma2"]],
			  phi=object[["phi"]],
			  b=b.i)
		## update tau2R and Delta

		DD.i <- res[["Delta"]]
		tau2R.i <- res[["tau2R"]]
		tau2R[, i] <- tau2R.i
		slot(object, "Delta") <- DD.i
		slot(object, "seed") <- res[["seed"]]
		slot(object, "tau2R") <- tau2R.i

		##updateDDelta_MI
		res <- .C("updateDDelta_MI",
			  seed=object[["seed"]],
			  nAccept=nAccept,
			  Delta=DD.i,
			  Q=Q,
			  G=G,
			  S=object[["S"]],
			  x=x,
			  psi=phenodata,
			  nu=object[["nu"]],
			  delta=d.i,
			  c2=object[["c2"]],
			  r=object[["r"]],
			  sigma2=object[["sigma2"]],
			  phi=object[["phi"]],
			  tau2R=tau2R.i,
			  b=object[["b"]])
		slot(object, "seed") <- res[["seed"]]
		DD.i <- res[["Delta"]]
		slot(object, "Delta") <- DD.i
		##all.equal(d.i, res[["delta"]])

		##updateC2_MI
		res <- .C("updateC2_MI",
			  seed=object[["seed"]],
			  nTry=nTry,
			  nAccept=nAccept,
			  c2=object[["c2"]],
			  Q=Q,
			  G=G,
			  delta=d.i,
			  Delta=DD.i,
			  r=object[["r"]],
			  sigma2=object[["sigma2"]],
			  tau2R=tau2R.i,
			  b=b.i,
			  c2Max=object[["c2Max"]])
		slot(object, "seed") <- res[["seed"]]
		c2[i] <- slot(object, "c2") <- res[["c2"]]

		##updateC2DDelta_MI
		res <- .C("updateC2DDelta_MI",
			  seed=object[["seed"]],
			  nTry=nTry,
			  nAccept=nAccept,
			  epsilon=epsilon,
			  c2=object[["c2"]],
			  Delta=DD.i,
			  Q=Q,
			  G=G,
			  S=S,
			  x=x,
			  psi=phenodata,
			  nu=object[["nu"]],
			  delta=d.i,
			  r=object[["r"]],
			  sigma2=object[["sigma2"]],
			  phi=object[["phi"]],
			  tau2R=tau2R.i,
			  b=b.i,
			  c2Max=object[["c2Max"]])
		##all.equal(d.i, res[["delta"]])
		DD.i <- res[["Delta"]]
		slot(object, "seed") <- res[["seed"]]
		c2[i] <- slot(object, "c2") <- res[["c2"]]
		slot(object, "Delta") <- DD.i

		##updateGamma2Nu_MI
		res <- .C("updateGamma2Nu_MI",
			  seed=object[["seed"]],
			  nTry=nTry,
			  nAccept=nAccept,
			  epsilon=epsilon,
			  gamma2=object[["gamma2"]],
			  nu=object[["nu"]],
			  Q=Q,
			  G=G,
			  S=S,
			  x=x,
			  psi=phenodata,
			  delta=d.i,
			  Delta=DD.i,
			  rho=object[["rho"]],
			  sigma2=object[["sigma2"]],
			  phi=object[["phi"]],
			  tau2Rho=object[["tau2Rho"]],
			  a=object[["a"]])
		slot(object, "seed") <- res[["seed"]]
		gamma2[i] <- slot(object, "gamma2") <- res[["gamma2"]]
		slot(object,  "nu") <- res[["nu"]]

		##updateRDDelta_MI
		res <- .C("updateRDDelta_MI",
			  seed=object[["seed"]],
			  nTry=nTry,
			  nAccept=nAccept,
			  epsilon=epsilon,
			  r=object[["r"]],
			  Delta=DD.i,
			  Q=Q,
			  G=G,
			  S=S,
			  x=x,
			  psi=phenodata,
			  nu=object[["nu"]],
			  delta=d.i,
			  c2=object[["c2"]],
			  sigma2=object[["sigma2"]],
			  phi=object[["phi"]],
			  tau2R=tau2R.i,
			  b=b.i,
			  nuR=object[["nuR"]])
		slot(object, "seed") <- res[["seed"]]
		r[, i] <- slot(object, "r") <- res[["r"]]
		DD.i <- res[["Delta"]]
		slot(object, "Delta") <- DD.i


		##updateRC2_MI
		res <- .C("updateRC2_MI",
			  seed=object[["seed"]],
			  nTry=nTry,
			  nAccept=nAccept,
			  epsilon=epsilon,
			  r=object[["r"]],
			  c2=object[["c2"]],
			  Q=object[["Q"]],
			  G=object[["G"]],
			  delta=object[["delta"]],
			  Delta=object[["Delta"]],
			  sigma2=object[["sigma2"]],
			  tau2R=object[["tau2R"]],
			  b=object[["b"]],
			  nuR=object[["nuR"]],
			  c2Max=object[["c2Max"]])
		slot(object, "seed") <- res[["seed"]]
		r[, i] <- slot(object, "r") <- res[["r"]]
		c2[i] <- slot(object, "c2") <- res[["c2"]]

		##updateSigma2_MI
		res <- .C("updateSigma2_MI",
			  seed=object[["seed"]],
			  nTry=nTry,
			  nAccept=nAccept,
			  epsilon=epsilon,
			  sigma2=object[["sigma2"]],
			  Q=Q,
			  G=G,
			  S=S,
			  x=x,
			  psi=phenodata,
			  nu=object[["nu"]],
			  delta=d.i,
			  Delta=DD.i,
			  c2=object[["c2"]],
			  gamma2=object[["gamma2"]],
			  r=object[["r"]],
			  rho=object[["rho"]],
			  phi=object[["phi"]],
			  t=object[["t"]],
			  l=object[["l"]],
			  tau2R=object[["tau2R"]],
			  tau2Rho=object[["tau2Rho"]],
			  a=object[["a"]],
			  b=object[["b"]])
		slot(object, "seed") <- res[["seed"]]
		slot(object, "sigma2") <- res[["sigma2"]]

		##updateLSigma2_MI
		res <- .C("updateLSigma2_MI",
			  seed=object[["seed"]],
			  nTry=nTry,
			  nAccept=nAccept,
			  epsilon=epsilon,
			  l=object[["l"]],
			  sigma2=object[["sigma2"]],
			  Q=Q,
			  G=G,
			  S=S,
			  x=x,
			  psi=phenodata,
			  nu=object[["nu"]],
			  delta=d.i,
			  Delta=DD.i,
			  c2=object[["c2"]],
			  gamma2=object[["gamma2"]],
			  r=object[["r"]],
			  rho=object[["rho"]],
			  phi=object[["phi"]],
			  t=object[["t"]],
			  tau2R=object[["tau2R"]],
			  tau2Rho=object[["tau2Rho"]],
			  a=object[["a"]],
			  b=object[["b"]])
		slot(object, "seed") <- res[["seed"]]
		l[, i] <- slot(object, "l") <- res[["l"]]
		slot(object, "sigma2") <- res[["sigma2"]]

		##updateTSigma2_MI
		res <- .C("updateTSigma2_MI",
			  seed=object[["seed"]],
			  nTry=nTry,
			  nAccept=nAccept,
			  epsilon=epsilon,
			  t=object[["t"]],
			  sigma2=object[["sigma2"]],
			  Q=Q,
			  G=G,
			  S=S,
			  x=x,
			  psi=phenodata,
			  nu=object[["nu"]],
			  delta=object[["delta"]],
			  Delta=object[["Delta"]],
			  c2=object[["c2"]],
			  gamma2=object[["gamma2"]],
			  r=object[["r"]],
			  rho=object[["rho"]],
			  phi=object[["phi"]],
			  l=object[["l"]],
			  tau2R=object[["tau2R"]],
			  tau2Rho=object[["tau2Rho"]],
			  a=object[["a"]],
			  b=object[["b"]])   ##nipun: Changed updateTSigma2 to updateTSigma2_MI
		slot(object, "seed") <- res[["seed"]]
		t[, i] <- slot(object, "t") <- res[["t"]]
		slot(object, "sigma2") <- res[["sigma2"]]
	}##for
	list(object=object,
	     d=d,
	     t=t,
	     l=l,
	     lambda=lambda,
	     theta=theta,
	     c2=c2,
	     r=r,
	     rho=rho,
	     tau2R=tau2R,
	     tau2Rho=tau2Rho,
	     a=a,
	     b=b,
	     gamma2=gamma2,
	     c2=c2,
	     xi=xi)

	return(object)
}##modelA_B_MIupdates
