\name{standardizeSamples}
\alias{standardizeSamples}
\title{Centers the genes at zero and standardizes the samples to have
  variance 1}

\description{
  
  For each study (element) in an \code{ExpressionSetList} object, this
  function centers the genes to have mean zero (rows) and scales the
  variance of the samples to 1.
  
}

\usage{
standardizeSamples(object, ...)
}

\arguments{
  \item{object}{Object of class \code{ExpressionSetList}}
  \item{\dots}{ Additional arguments not implemented}
}

\value{
  An object of class \code{ExpressionSetList}
}
\author{R. Scharpf}
\note{

  Requires genefilter package
  
}
\keyword{methods}

