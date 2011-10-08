\name{RUpdates}
\alias{rupdateANu}
\alias{rupdateBDDelta}
\alias{rupdateTau2RhoNu}
\alias{rupdateTau2RDDelta}
\alias{rupdateNu}
\alias{rupdateDDelta}
\alias{rupdateC2}
\alias{rupdateC2DDelta}
\alias{rupdateGamma2}
\alias{rupdateGamma2Nu}
\alias{rupdateRDDelta}
\alias{rupdateRC2}
\alias{rupdateRhoNu}
\alias{rupdateRhoGamma2}
\alias{rupdateSigma2}
\alias{rupdatePhi}
\alias{rupdateTheta}
\alias{rupdateLambda}
\alias{rupdateT}
\alias{rupdateL}
\alias{rupdateXi}
\alias{rupdateDeltaDDelta}
\alias{rupdateLSigma2}
\alias{rupdateTSigma2}
\alias{rupdateLambdaPhi}
\alias{rupdateThetaPhi}
\alias{rupdateDeltaDDelta_MRF}
\alias{rupdateAlpha_MRF}
\alias{rupdateBeta_MRF}
\alias{rupdateBetag_MRF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ R interface to MCMC updates in C }
\description{R interface to MCMC updates in C }
\usage{
rupdateANu(object, nTry = 10L, nAccept = 0L, epsilon = 0.2, dryrun = FALSE)
rupdateBDDelta(object, nTry = 10L, nAccept = 0L, epsilon = 0.2, dryrun = FALSE)
rupdateBDDelta(object, nTry = 10L, nAccept = 0L, epsilon = 0.2, dryrun = FALSE)
rupdateTau2RhoNu(object, nTry = 10L, nAccept = 0L, epsilon = 0.2, dryrun = FALSE)
rupdateTau2RDDelta(object, nTry = 10L, nAccept = 0L, epsilon = 0.2, dryrun = FALSE)
rupdateNu(object, nAccept = 0L, dryrun = FALSE)
rupdateDDelta(object, nAccept = 0L, dryrun = FALSE)
rupdateC2(object, nTry=10L, nAccept = 0L, dryrun = FALSE)
rupdateC2DDelta(object, nTry = 10L, nAccept = 0L, epsilon = 0.2, dryrun = FALSE)
rupdateGamma2(object, nAccept=0L, dryrun=FALSE)
rupdateGamma2Nu(object, nTry = 10L, nAccept = 0L, epsilon = 0.2, dryrun = FALSE)
rupdateRDDelta(object, nTry = 10L, nAccept = 0L, epsilon = 0.2, dryrun = FALSE)
rupdateRC2(object, nTry = 10L, nAccept = 0L, epsilon = 0.2, dryrun = FALSE)
rupdateRhoNu(object, nTry = 10L, nAccept = 0L, epsilon = 0.2, dryrun = FALSE)
rupdateRhoGamma2(object, nTry = 10L, nAccept = 0L, epsilon = 0.2, dryrun = FALSE)
rupdateSigma2(object, nTry = 10L, nAccept = 0L, epsilon = 0.2, dryrun = FALSE)
rupdatePhi(object, nTry = 10L, nAccept = 0L, epsilon = 0.2, dryrun = FALSE)
rupdateTheta(object, nTry = 10L, nAccept = 0L, epsilon = 0.2, dryrun = FALSE)
rupdateLambda(object, nTry = 10L, nAccept = 0L, epsilon = 0.2, dryrun = FALSE)
rupdateT(object, nTry = 10L, nAccept = 0L, epsilon = 0.2, dryrun = FALSE)
rupdateL(object, nTry = 10L, nAccept = 0L, epsilon = 0.2, dryrun = FALSE)
rupdateXi(object, nAccept = 0L, dryrun = FALSE, one.delta=FALSE)
rupdateDeltaDDelta(object, nTry = 10L, nAccept = 0L, dryrun = FALSE, one.delta=FALSE)
rupdateLSigma2(object, nTry = 10L, nAccept = 0L, epsilon = 0.2, dryrun = FALSE)
rupdateTSigma2(object, nTry = 10L, nAccept = 0L, epsilon = 0.2, dryrun = FALSE)
rupdateLambdaPhi(object, nTry = 10L, nAccept = 0L, epsilon = 0.2, dryrun = FALSE)
rupdateThetaPhi(object, nTry = 10L, nAccept = 0L, epsilon = 0.2, dryrun = FALSE)
rupdateDeltaDDelta_MRF(object, nTry = 10L, nAccept = 0L, dryrun = FALSE, one.delta=FALSE)
rupdateAlpha_MRF(object, nTry = 10L, nAccept = 0L, epsilon=0.2, dryrun = FALSE, one.delta=FALSE)
rupdateBeta_MRF(object, nTry = 10L, nAccept = 0L, epsilon=0.2, dryrun = FALSE, one.delta=FALSE)
rupdateBetag_MRF(object, nTry = 10L, nAccept = 0L, epsilon=0.2, dryrun = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    Object of class \code{Parameters}}
  \item{nTry}{Number of MCMC samples}
  \item{nAccept}{Number of samples accepted.  Only relevant for Metropolis-Hastings updates}
  \item{epsilon}{Tuning parameter for the proposal. Only relevant for Metropolis-Hastings updates}
  \item{dryrun}{Return a list of the arguments that are to be passed to the low-level C function, without calling the C function}
  \item{one.delta}{Logical.  If true, we assume that a gene is differentially expressed in all studies or in none of the studies}
}

\details{

  The parameters that are updated are indicated by the function name.
  For example, \code{rupdateSigma2} updates the parameter \code{sigma2}.

}

\value{

  A object of class \code{Parameters}.  The returned object contains any
  updates from the MCMC, as well as a new seed.

}

\author{RS}

\seealso{\code{\linkS4class{Parameters}}}

\keyword{models}
\keyword{interface}
\keyword{iteration}

