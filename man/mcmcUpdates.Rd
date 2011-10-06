\name{mcmcUpdates}
\alias{updateANu}
\title{
  Updates for Bayesian model parameters.
}
\description{
  Updates for Bayesian model parameters.
}
\usage{
updateANu(object, hyper.params, params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A \code{ExpressionSetList}.
}
  \item{hyper.params}{
    \code{list}. The hyperparameters.
}
  \item{params}{
    \code{list}. The parameters.
}
}
\details{

  \code{updateANu} updates nu and the power term \code{a}.  We assume
  that nu is multivariate normal.  nu is updated from its full
  conditional by a Gibbs step.  'a' is updated by a truncated random
  walk proposal.


}
\value{
  list
}
\references{
  Scharpf et al, 2009 JASA 104, p 1295--1310
}
\author{
R. Scharpf
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(expressionSetList)
hyper.params <- XDE:::getHyperparameters(object=expressionSetList)
params <- getParameters(hyper.params)
res <- updateANu(expressionSetList,
		 hyper.params,
		 params)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}

