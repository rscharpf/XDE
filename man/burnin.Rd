\name{burnin}
\alias{burnin}
\alias{burnin<-}
\title{Indicator for running a MCMC burnin}
\description{
  When \code{TRUE}, log files from MCMC chains are not written to file.
  When \code{FALSE}, log files are written for every parameter by default.
}
\usage{
burnin(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object of class \code{XdeParameter} }
}

\value{
  logical
}
\author{R. Scharpf}
\seealso{\code{\link{XdeParameter-class}}}
\examples{
\dontrun{
data(expressionSetList)
params <- new("XdeParameter", phenotypeLabel="adenoVsquamous",
              esetList=expressionSetList)

##the replacement method for burnin is called for its side effect of
##providing default values of storing MCMC chains
output(params)[2:22]
burnin(params) <- FALSE
output(params)[2:22]
burnin(params) <- TRUE
output(params)[2:22]
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

