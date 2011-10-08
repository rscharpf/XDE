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
		      nTry=10L,
		      nAccept=0L,
		      epsilon=0.2,
		      dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
			    nu=object[["nu"]],
			    Q=object[["Q"]],
			    G=object[["G"]],
			    S=object[["S"]],
			    x=exprs(object),
			    psi=object[["phenodata"]],
			    delta=object[["delta"]],
			    Delta=object[["Delta"]],
			    nu=object[["nu"]],
			    delta=object[["delta"]],
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
		  nTry=nTry,
		  nAccept=nAccept,
		  epsilon=epsilon,
		  nu=object[["nu"]],
		  Q=object[["Q"]],
		  G=object[["G"]],
		  S=object[["S"]],
		  x=exprs(object),
		  psi=object[["phenodata"]],
		  delta=object[["delta"]],
		  Delta=object[["Delta"]],
		  nu=object[["nu"]],
		  delta=object[["delta"]],
		  gamma2=object[["gamma2"]],
		  rho=object[["rho"]],
		  sigma2=object[["sigma2"]],
		  phi=object[["phi"]],
		  tau2Rho=object[["tau2Rho"]],
		  a=object[["a"]])
	slot(object, "seed") <- res[["seed"]]
	slot(object, "nu") <- res[["nu"]]
	return(res)
}

rupdateDDelta <- function(object,
			  nTry=10L,
			  nAccept=0L,
			  epsilon=0.2,
			  dryrun=FALSE){
	if(dryrun){
		res <- list(seed=object[["seed"]],
			    nTry=nTry,
			    nAccept=nAccept,
			    epsilon=epsilon,
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
}

rupdateC2 <- function(object,
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
}


rupdateGamma2 <- function(object,
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
			    Q=object[["Q"]],
			    G=object[["G"]],
			    nu=object[["nu"]],
			    rho=object[["rho"]],
			    sigma2=object[["sigma2"]],
			    tau2Rho=object[["tau2Rho"]],
			    a=object[["a"]])
		return(res)
	}
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
			    nu=object[["nu"]],
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
}

rupdateXi <- function(object,
		      nTry=10L,
		      nAccept=0L,
		      epsilon=0.2,
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
		##

	} else {
		##

	}
}

rupdateDeltaDDelta <- function(object,
			       nTry=10L,
			       nAccept=0L,
			       epsilon=0.2,
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
		##

	} else {
		##

	}
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

	} else {
	}
}


rupdateAlpha_MRF <- function(object,
			     nTry=10L,
			     nAccept=0L,
			     epsilon=epsilon,
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

	} else {
	}
}

rupdateBeta_MRF <- function(object,
			     nTry=10L,
			     nAccept=0L,
			     epsilon=epsilon,
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

	} else {
	}
}

rupdateBetag_MRF <- function(object,
			     nTry=10L,
			     nAccept=0L,
			     epsilon=epsilon,
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
		if(one.delta) res <- res[-length(res)]
		return(res)
	}
	if(one.delta){

	} else {
	}
}
