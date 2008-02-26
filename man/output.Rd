\name{output}
\alias{output}
\alias{output<-}
\title{Options for storing results of the MCMC chains}
\description{
  A numeric vector indicating which chains to write to file and, for those
  parameters that are written to file, how often the chains should be
  written to file.  
}
\usage{
output(object)
}
\arguments{
  \item{object}{An object of class \code{XdeParameter} or \code{XdeMcmc} }
}
\details{
Replacement methods are only available for objects
  of class \code{XdeParameter}.  Accessor methods are available for
  objects of class \code{XdeParameter} and \code{XdeMcmc}.
}
\value{
  A named numerical vector.  The first element (thin) specifies how often to write
  chains to file.  For instance, if output[1]=2 the chains will be
  written to file every other iteration.  Elements 2 - 22 of the vector
  are indicators for whether to write the write the chains of the
  Bayesian parameters to file. 
}
\author{R. Scharpf}

\note{

  Parameters indexed by gene and study (Delta, Phi, Nu, and sigma2) grow
  very large quickly.

}

\seealso{\code{\link{burnin}}, \code{\link{XdeParameter-class}}, \code{\link{XdeMcmc-class}}}
\examples{
data(xmcmc)
output(xmcmc)
}
\keyword{methods}

