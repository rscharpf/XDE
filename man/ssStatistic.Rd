\name{ssStatistic}
\alias{ssStatistic}
\title{Calculate single study estimates of effect size}
\description{
  Calculate single study estimates of effect size for lists of \code{ExpressionSets}
}
\usage{
ssStatistic(statistic = c("t", "sam", "z")[1], phenotypeLabel, esetList, ...)
}
\arguments{
  \item{statistic}{Character string indicating Welch t-statistic (t), SAM (sam), or a z-statistic (z)}
  \item{phenotypeLabel}{Character string indicating the name of the binary covariate}
  \item{esetList}{An object of class \code{ExpressionSetList}}
  \item{\dots}{Not implemented. Potentially additional arguments to the
  above methods that are implemented in other packages}
}

\details{

  This function is a wrapper that provides an estimate of effect size
  for each study (element) in an \code{ExpressionSetList} object.

  For Welch t-statistic, this function is a wrapper for mt.teststat in the
  multtest package.
  
  For SAM, this function is a wrapper for the sam function in the
  siggenes package.

  The "z" statistic is a standardized unbiased estimate of effect size
  (Hedges and Olkin, 1985) -- implementation is in the zScores
  function in the R package GeneMeta.

See the complete references below.                  
  
}

\value{
  
  A matrix: rows are genes and columns are studies
  
}

\references{

J.K. Choi, U. Yu, S. Kim, and O.J. Yoo (2003), Combining multiple
  microarray studies and modeling interstudy variation, Bioinformatics,
  19(1) I84-I90.

Y. Ge, S. Dudoit & T. P. Speed (2003), Resampling-based multiple testing
for microarray data hypothesis Test 12(1) : 1-44 (with discussions on
44-77).

L. Lusa R. Gentleman, and M. Ruschhaupt, GeneMeta: MetaAnalysis for High
Throughput Experiments

L.V. Hedges and I. Olkin, Statistical Methods for Meta-analysis
(1985), Academic Press

Tusher, Tibshirani and Chu (2001), Significance analysis of microarrays
applied to the ionizing radiation response, PNAS 2001 98: 5116-5121,
(Apr 24).

}
\author{R. Scharpf}
\examples{

data(expressionSetList)
if(require(siggenes)){
  sam <- ssStatistic("sam", esetList=expressionSetList, phenotypeLabel="adenoVsquamous")
}  
}
\keyword{methods}

