\name{posteriorAvg}
\alias{posteriorAvg}
\title{Accessor for posterior averages}

\description{ Accessor for posterior averages of the probability of
  differential expression, concordant differential expression, and
  discordant differential expression.  }

\usage{
posteriorAvg(object)
}
\arguments{
  \item{object}{An object of class \code{XdeMcmc}}
}

\value{
  
  A matrix. Rows are genes and columns are posterior averages for
  differential expression ("diffExpressed"), concordant differential
  expression ("concordant"), and discordant differential expression
  ("discordant").
  
}

\references{

  R. Scharpf et al., A Bayesian Model for Cross-Study Differential
  Gene Expression, Technical Report 158, Johns Hopkins University,
  Department of Biostatistics, 2007
  
}

\author{R. Scharpf}
\seealso{\code{\link{XdeMcmc-class}}}

\examples{
data(xmcmc)
posteriorAvg(xmcmc)[1:5, ]
hist(posteriorAvg(xmcmc)[, "concordant"], breaks=50, main="")
}
\keyword{methods}

