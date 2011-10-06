\name{Parameters-utils}
\alias{HyperParams}
\alias{Params}
\title{
  Initialize both hyperparameters and parameters.
}
\description{
  Initialize hyperparameters and parameters from a
  \code{ExpressionSetList} object.
}

\usage{
HyperParams(object, G, Q, S, seed, clinicalCovariate, ...)
Params(object, one.delta = FALSE, MRF = FALSE)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{object}{For \code{HyperParams}, object must be an
  \code{ExpressionSetList}.  For \code{Params}, 'object' must be a \code{HyperParams} object.}
  \item{G}{ integer.  The number of genes.  Ignored if
    \code{object} is not missing.}
  \item{Q}{ integer.  The number of studies.  Ignored if
    \code{object} is not missing.}
  \item{S}{ vector of integers.  The number of samples in each study.  Ignored if
    \code{object} is not missing.}
  \item{clinicalCovariate}{character.  Name of binary clinical covariate.}
  \item{seed}{Integer.  Seed for reproducibility.}
  \item{\dots}{Ignored}
  \item{one.delta}{Logical.  If TRUE, we assume that a gene is
  differentiall expressed in all studies or none of the studies.}
  \item{MRF}{Logical.  If TRUE, we use a Markov Random Field prior on the deltas. }
}
\details{

}
\value{
  Return \code{HyperParams} or \code{Params} object, respectively.
}
\author{
  R. Scharpf
}

\seealso{
   \code{\linkS4class{HyperParams}}, \code{\linkS4class{Params}}
}
\examples{
data(expressionSetList)
hyper.params <- HyperParams(expressionSetList, clinicalCovariate="adenoVsquamous")
params <- Params(hyper.params, one.delta=FALSE, MRF=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{list}
\keyword{methods}

