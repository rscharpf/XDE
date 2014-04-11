\name{updates}
\alias{updates}
\alias{updates<-}
\title{Frequency of updating a parameter per MCMC iteration}
\description{
  
  Accessor and replacement methods for the class \code{XdeParameter} are
  available.  Specifying an update of integer N for a
  Metropolis-Hastings parameter means that N values are proposed for
  that parameter for each MCMC iteration.
}
\usage{
updates(object)
}

\arguments{
  \item{object}{ An object of class \code{XdeParameter} }
}
\details{
  See the XdeParameterClass vignette
}
\value{
  A numerical vector
}
\author{R. Scharpf}
\keyword{methods}

