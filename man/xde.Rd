\name{xde}
\alias{xde}
\title{Fit the Bayesian hierarchical model for cross-study differential
  gene expression}

\description{

  Fits the Bayesian hierarchical model for cross-study differential
  gene expression.

}

\usage{

  xde(paramsMcmc, esetList, outputMcmc, batchSize=NULL, NCONC=2,
  center=TRUE, ...)

}


\arguments{
  \item{paramsMcmc}{Object of class \code{XdeParameter}}
  \item{esetList}{Object of class \code{ExpressionSetList}}
  \item{outputMcmc}{Object of class \code{XdeMcmc} (optional)}
  \item{batchSize}{Integer or NULL.  The number of iterations written to
  log files before summarizing the chain and then removing. Experimental.}
  \item{NCONC}{The number of studies for which a gene must be
  differentially expressed in the same direction to be considered as
  concordantly differentially expressed.}
\item{center}{Logical.  If TRUE, each study is centered to have mean
  zero.}
\item{\dots}{Additional arguments passed to \code{xdeFit}.}
}

\details{

  Details for fitting the Bayesian model are discussed elsewhere (see
  citation below and XdeParameterClass vignette)

  If an integer is specified for the batchSize, summary statistics for
  the log-files are calculated for every batchSize iterations.  The log
  files are then removed and the next iteration will start a new log
  file.  This allows one to do many iterations without creating enormous
  log files.  This is only reasonable to do if one has already assessed
  convergence.

}

\value{
  Object of class \code{XdeMcmc}
}

\references{

  R. Scharpf et al., A Bayesian Model for Cross-Study Differential Gene
  Expression, JASA 2009, p1295--1310.

}

\author{R. Scharpf}
\note{

  See the vignettes for XdeParameterClass and XDE.
}

\seealso{\code{\link{XdeMcmc-class}}, \code{\link{XdeParameter-class}}, \code{\link{ExpressionSetList-class}}}

\examples{
\dontrun{
  data(expressionSetList)
  xparam <- new("XdeParameter", phenotypeLabel="adenoVsquamous", esetList=expressionSetList)
  iterations(xparam) <- 10
  fit <- xde(xparam, esetList=expressionSetList)
}
}
\keyword{models}
\keyword{multivariate}

