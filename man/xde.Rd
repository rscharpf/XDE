\name{xde}
\alias{xde}
\title{Fit the Bayesian hierarchical model for cross-study differential
  gene expression}

\description{
  
  Fits the Bayesian hierarchical model for cross-study differential
  gene expression. 
  
}

\usage{
  
  xde(paramsMcmc, esetList, outputMcmc)

}


\arguments{
  \item{paramsMcmc}{Object of class \code{XdeParameter}}
  \item{esetList}{Object of class \code{ExpressionSetList}}
  \item{outputMcmc}{Object of class \code{XdeMcmc} (optional)}
}

\details{

  Details for fitting the Bayesian model are discussed elsewhere (see
  citation below and XdeParameterClass vignette)
  
}

\value{
  Object of class \code{XdeMcmc}
}

\references{

  R. Scharpf et al., A Bayesian Model for Cross-Study Differential Gene
  Expression, Technical Report 158, Johns Hopkins University, Department
  of Biostatistics, 2007

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

