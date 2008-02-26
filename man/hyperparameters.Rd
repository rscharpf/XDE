\name{hyperparameters}
\alias{hyperparameters}
\alias{hyperparameters<-}
\title{Accessor for hyperparameters of the Bayesian model}
\description{
  Accessor and replacement methods for hyperparameters of the Bayesian
  model are provided
}
\usage{
hyperparameters(object)
}
\arguments{
  \item{object}{An object of class \code{XdeParameter}}
}
\details{
  See the XdeParameterClass vignette for a more detailed discussion.
  The default values provided when initializing an object of class
  \code{XdeParameter} works well in most instances.
}
\value{
  A numerical vector
}
\references{

  R. Scharpf et al., A Bayesian Model for Cross-Study Differential Gene
  Expression, Technical Report 158, Johns Hopkins University, Department
  of Biostatistics, 2007

}
\author{R. Scharpf}
\examples{

data(expressionSetList)
xlist <- new("XdeParameter", esetList=expressionSetList, phenotypeLabel="adenoVsquamous")
hyperparameters(xlist)

}
\keyword{methods}
