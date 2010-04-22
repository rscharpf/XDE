\name{iterations}
\alias{iterations}
\alias{iterations<-}
\title{Number of MCMC iterations}
\description{
  Number of MCMC iterations
}
\usage{
iterations(object)
}
\arguments{
  \item{object}{ An object of class \code{XdeParameter} or
    \code{XdeMcmc}.}
}

\details{

  For an object of class \code{XdeParameter}, \code{iterations}
  specifies the total number of MCMC iterations. Note that by setting
  the \code{thin} parameter to a value greater than 1, the number of
  MCMC iterations will be greater than the number of saved MCMC
  iterations (saved iterations = iterations / thin).

  For an object of class \code{XdeMcmc} (a class that stores output from
  the MCMC), \code{iterations} specifies the number of iterations that
  were saved.

  The replacement method is only defined for the \code{XdeParameter}
  class.  The class \code{XdeMcmc} is meant to reflect the information
  in an already run chain, whereas \code{XdeParameter} is a class for
  parameterizing the Bayesian model that has not yet been fit.

}

\value{
  An integer
}
\author{R. Scharpf}
\seealso{\code{\link{XdeParameter-class}}, \code{\link{XdeMcmc-class}}}
\keyword{methods}
