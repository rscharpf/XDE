\name{thin}
\alias{thin}
\alias{thin<-}
\title{How often to write MCMC iterations to file}
\description{
A value greater than one means that not every MCMC iteration is written
to file.
}
\usage{
thin(x, ...)
}
\arguments{
  \item{x}{An object of class \code{XdeParameter}}
  \item{\dots}{not implemented}
}
\details{
  \code{thin} is an accessor for the first element in the vector
  returned by the method \code{output}.

  The replacement method replaces the first element in the \code{output}
  vector.

}

\value{
An integer.
}
\author{R. Scharpf}

\seealso{\code{\link{output}}}
\keyword{methods}
