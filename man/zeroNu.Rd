\name{zeroNu}
\alias{zeroNu}
\title{Option for not modeling Nu}

\description{
  Nu is the average expression value in each study.
}

\usage{
zeroNu(object, ...)
}

\arguments{
  \item{object}{object of class \code{ExpressionSetList}}
  \item{\dots}{ Not implemented}
}

\details{

  This function should be regarded as experimental.

  The nu parameter models the average expression value in each
  study. Modeling nu allows one to estimate differential expression
  across studies that may differ in location and scale (as often occurs
  when multiple platforms are used).  The price to pay for modeling nu
  are additional assumptions (the nu\'s are assumed Gaussian) and a more
  heavily parameterized model.

  The method zeroNu allows one to fit the Bayesian model without
  estimating nu:

  - each gene is centered at zero

  - initial values for the first MCMC are chosen on the basis of
  empirical starting values

  - the initial values for a and rho are set to zero.

  - the nu, a, gamma2, and rho parameters are not updated during MCMC

}

\value{

  object of class \code{XdeParameter}
  
}

\references{

    R. Scharpf et al. (2007), A Bayesian Model for Cross-Study Differential
  Gene Expression, Technical Report 158, Johns Hopkins University,
  Department of Biostatistics
  
}

\author{R. Scharpf}
\keyword{methods}

