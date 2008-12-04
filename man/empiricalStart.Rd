\name{empiricalStart}
\alias{empiricalStart}
\title{Empirical starting values for the MCMC}
\description{
  Empirical starting values for the MCMC are based on data in objects of
  class \code{ExpressionSetList} 
}
\usage{
empiricalStart(object, zeroNu = FALSE, phenotypeLabel)
}
\arguments{
  \item{object}{An object of class \code{ExpressionSetList}}
  \item{zeroNu}{Logical: if TRUE, the nu in the Bayesian model are not
    modeled -- set to zero and not updated in the MCMC.  Setting zeroNu
    to TRUE should be regarded as experimental}
  \item{phenotypeLabel}{character: binary phenotype.  phenotypeLabel
    must be in the varLabels of each ExpressionSet object}
}

\value{
  An object of class \code{XdeParameter}
}
\author{R. Scharpf}
\seealso{ \code{\link{zeroNu}}, \code{\link{XdeParameter-class}}, \code{\link{ExpressionSetList-class}}}
\keyword{methods}
