\name{getParameters}
\alias{getParameters}
\alias{getHyperParameters}
\title{
  Initialize both hyperparameters and parameters.
}
\description{
  Initialize hyperparameters and parameters from a
  \code{ExpressionSetList} object.
}

\usage{
getHyperparameters(object, G, Q, S, seed, ...)
getParameters(hyper.parameters, one.delta = FALSE, MRF = FALSE)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{object}{ A \code{ExpressionSetList}.}
  \item{G}{ integer.  The number of genes.  Ignored if
    \code{object} is not missing.}
  \item{Q}{ integer.  The number of studies.  Ignored if
    \code{object} is not missing.}
  \item{S}{ vector of integers.  The number of samples in each study.  Ignored if
    \code{object} is not missing.}
  \item{seed}{Integer.  Seed for reproducibility.}
  \item{\dots}{Ignored}
  \item{hyper.parameters}{ A \code{list} of hyperparameters. }
  \item{one.delta}{Logical.  If TRUE, we assume that a gene is
  differentiall expressed in all studies or none of the studies.}
  \item{MRF}{Logical.  If TRUE, we use a Markov Random Field prior on the deltas. }
}
\details{

}
\value{
  A list.
}
\author{
  R. Scharpf
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(expressionSetList)
ohyper.params <- getHyperparameters(expressionSetList)
params <- getParameters(hyper.params, one.delta=FALSE, MRF=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{list}

