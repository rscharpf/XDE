\name{empiricalStart}
\alias{empiricalStart}
\title{Empirical starting values for the MCMC}
\description{
  Empirical starting values for the MCMC are based on data in objects of
  class \code{ExpressionSetList} 
}
\usage{
empiricalStart(object, zeroNu = FALSE, phenotypeLabel, one.delta=FALSE, T_THRESH=4)
}
\arguments{
  \item{object}{An object of class \code{ExpressionSetList}}
  \item{zeroNu}{Logical: if TRUE, the nu in the Bayesian model are not
    modeled -- set to zero and not updated in the MCMC.  Setting zeroNu
    to TRUE should be regarded as experimental}
  \item{phenotypeLabel}{character: binary phenotype.  phenotypeLabel
    must be in the varLabels of each ExpressionSet object}
  
  \item{one.delta}{delta in the Bayesian model is a gene-specific
    indicator for differential expression.  If one.delta is FALSE, we
    assume that a gene can be differentially expressed in a subset of
    studies.  When TRUE, we assume that a gene is differentially
    expressed in all studies or in none.}
  
  \item{T_THRESH}{A threshold of t-statistics (calculated row-wise for
    each study) for determining starting values of the differential
    expression indicator, delta.}
}

\value{

  A list containing starting values for the MCMC that are derived from
  empirical estimates of the data.  
  
}
\author{R. Scharpf}

\seealso{ \code{\link{zeroNu}}, \code{\link{XdeParameter-class}},
  \code{\link{ExpressionSetList-class}}}

\examples{
library(XDE)
data(expressionSetList)
eList <- studyCenter(expressionSetList)
empirical <- empiricalStart(eList, phenotypeLabel="adenoVsquamous", T_THRESH=3)
##By default, initial values for the MCMC are sampled from the prior
##when initializing an object of class XdeParamater
params <- new("XdeParameter", esetList=eList,
	      phenotypeLabel="adenoVsquamous", one.delta=FALSE, burnin=TRUE)
##The initial values can be replaced by empirical values as follows:
firstMcmc(params) <- empirical
}


\keyword{methods}
