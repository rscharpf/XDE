\name{studyCenter}
\alias{studyCenter}
\title{Center the expression values in a study to zero}
\description{
  Centers each study in a list so that the average expression value of
  each stuy is zero
}
\usage{
studyCenter(object)
}
\arguments{
  \item{object}{ An object of class \code{ExpressionSetList}}
}
\value{
  An object of class \code{ExpressionSetList}
}
\author{R. Scharpf}
\seealso{\code{\link{geneCenter}}, \code{\link{ExpressionSetList-class}}}
\examples{
  data(expressionSetList)
  centered <- studyCenter(expressionSetList)
  lapply(centered, function(object) round(mean(exprs(object)), 4))
}
\keyword{methods}

