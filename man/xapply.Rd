\name{xapply}
\alias{xapply}
\title{An lapply for ExpressionSetList}
\description{
  This is an lapply for ExpressionSetList instances
}
\usage{
xapply(X, FUN, ...)
}
\arguments{
  \item{X}{\code{ExpressionSetList} instance}
  \item{FUN}{function to apply to each element in \code{X}}
  \item{\dots}{Additional arguments to \code{FUN}}
}
\value{
  An object of class \code{ExpressionSetList}
}
\author{RS}
\seealso{\code{\link{ExpressionSetList-class}}}
\keyword{methods}
