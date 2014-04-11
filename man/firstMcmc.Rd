\name{firstMcmc}
\alias{firstMcmc}
\alias{firstMcmc<-}
\title{Values for the first MCMC iteration}
\description{
  Accessor method for the values of the first MCMC iteration
}
\usage{
firstMcmc(object)
}
\arguments{
  \item{object}{ An object of class \code{XdeParameter} }
}
\value{
  Returns a list of the values to be used in the first iteration of the
  MCMC.
}
\author{R. Scharpf}
\seealso{\code{\link{XdeParameter-class}}, \code{\link{lastMcmc}}}
\examples{
data(expressionSetList)
params <- new("XdeParameter", phenotypeLabel="adenoVsquamous",
              esetList=expressionSetList)
str(firstMcmc(params))
}
\keyword{methods}

