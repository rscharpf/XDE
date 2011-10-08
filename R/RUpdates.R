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
	res <- .C("updateTau2RDDelta",
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
	res <- .C("updateDDelta",
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
	return(object)
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
	res <- .C("updateC2",
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
	slot(object, "seed") <- object[["seed"]]
	slot(object, "c2") <- object[["c2"]]
	return(object)
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
	res <- .C("updateC2DDelta",
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
	return(object)
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
		res <- .C("updateGamma2Nu",
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
	return(object)
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
		res <- .C("updateRDDelta",
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
		res <- .C("updateRC2",
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
		res <- .C("updateSigma2",
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
		res <- .C("updateXi_onedelta",
			  seed=object[["seed"]],
			    nAccept=nAccept,
			    xi=object[["xi"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    delta=object[["delta"]],
			    t=object[["t"]],
			    alphaXi=object[["alphaXi"]],
			    betaXi=object[["betaXi"]])		##

	} else {
		res <- .C("updateXi",
			  seed=object[["seed"]],
			    nAccept=nAccept,
			    xi=object[["xi"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    delta=object[["delta"]],
			    t=object[["t"]],
			    alphaXi=object[["alphaXi"]],
			    betaXi=object[["betaXi"]])
	}
	slot(object, "seed") <- res[["seed"]]
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
		res <- .C("updateDeltaDDelta_onedelta",
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
			    b=object[["b"]])
	} else {
		res <- .C("updateDeltaDDelta",
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
			    b=object[["b"]])
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
		res <- .C("updateLSigma2",
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
	return(object)
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
		res <- .C("updateTSigma2",
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
			    b=object[["b"]])
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

rupdateDeltaDDelta_MRF <- function(object,
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
		res <- .C("updateDeltaDDelta_MRF_onedelta",
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
			    alpha=object[["alpha"]],
			    beta=object[["beta"]])
	} else {
		res <- .C("updateDeltaDDelta_MRF",
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
			    alpha=object[["alpha"]],
			    beta=object[["beta"]],
			  betag=object[["betag"]])
	}
	slot(object, "seed") <- res[["seed"]]
	slot(object, "delta") <- res[["delta"]]
	slot(object, "Delta") <- res[["Delta"]]
	return(object)
}


rupdateAlpha_MRF <- function(object,
			     nTry=10L,
			     nAccept=0L,
			     epsilon=0.2,
			     dryrun=FALSE,
			     one.delta=FALSE){
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
			    beta=object[["beta"]],
			    betag=object[["betag"]])
		if(one.delta) res <- res[-length(res)]
		return(res)
	}
	if(one.delta){
		res <- .C("updateAlpha_MRF_onedelta",
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
	} else {
		res <- .C("updateAlpha_MRF",
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
	}
	slot(object, "seed") <- res[["seed"]]
	slot(object, "alpha") <- res[["alpha"]]
	return(object)
}

rupdateBeta_MRF <- function(object,
			     nTry=10L,
			     nAccept=0L,
			     epsilon=0.2,
			     dryrun=FALSE,
			     one.delta=FALSE){
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
		if(one.delta) res <- res[-length(res)]
		return(res)
	}
	if(one.delta){
		res <- .C("updateBeta_MRF_onedelta",
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
	} else {
		res <- .C("updateBeta_MRF",
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
	}
	slot(object, "seed") <- res[["seed"]]
	slot(object, "beta") <- res[["beta"]]
	return(object)
}

rupdateBetag_MRF <- function(object,
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
		res <- .C("updateBetag_MRF",
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
