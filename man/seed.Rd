\name{seed}
\alias{seed}
\alias{seed<-}
\title{Seed for the MCMC}
\description{
  Setting a seed is useful for reproducing MCMC chains
}
\usage{
seed(object)
}
\arguments{
  \item{object}{ An object of \code{XdeParameter} or \code{XdeMcmc} }
}
\details{
  
  The seed stored in the slot of an object of class \code{XdeParameter}
  and an object of class \code{XdeMcmc} are useful in different ways.
  For the \code{XdeParameter} class, the seed indicates what seed was
  used to initialize an MCMC chain.  By contrast, an object of class
  \code{XdeMcmc} contains a seed that would be useful for running
  additional iterations -- the seed here is guaranteed to be different
  from the seed that was used to initiate the MCMC. 
}
\value{
  An integer
}
\author{R. Scharpf}
\keyword{methods}

