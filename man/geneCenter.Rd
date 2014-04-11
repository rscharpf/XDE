\name{geneCenter}
\alias{geneCenter}
\title{Center the expression values for each gene in a study to zero}
\description{
  Mean centers the genes for each study in a list
}
\usage{
geneCenter(object)
}
\arguments{
  \item{object}{ Object of class \code{ExpressionSetList}}
}
\value{
  Object of class \code{ExpressionSetList}
}
\author{R. Scharpf}
\seealso{ \code{\link{studyCenter}}, \code{\link{ExpressionSetList-class}}}

\examples{
  data(expressionSetList)
  centered <- geneCenter(expressionSetList)
}
\keyword{methods}

