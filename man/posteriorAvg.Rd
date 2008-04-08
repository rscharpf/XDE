\name{posteriorAvg}
\alias{posteriorAvg}
\alias{posteriorAvg<-}
\title{Accessor and replacement methods for posterior averages of
  differential expression}
\description{

  Accessor and replacement methods for objects of class \code{XdeMcmc}
  for posterior averages of differential expression

}
\usage{
posteriorAvg(object)
posteriorAvg(object) <- value
}
\arguments{
  \item{object}{Object of class \code{XdeMcmc}}
  \item{value}{A matrix of dimension G x 3, where G is the number of
    genes and 3 are different ways of quantifying differential
    expression in the context of multiple studies (concordant,
    discordant, or the union).}
}
  
\value{

  A matrix of dimension G x 3, where G is the number of genes and 3 are
    different ways of quantifying differential expression in the context
    of multiple studies (concordant, discordant, or the union).
  
}
\author{RS}
\seealso{\code{\link{calculatePosteriorAvg}}}
\keyword{methods}

