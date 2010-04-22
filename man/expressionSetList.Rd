\name{expressionSetList}
\alias{expressionSetList}
\docType{data}
\title{Example of ExpressionSetList}
\description{
  Object of class \code{ExpressionSetList} containing three
  studies. Each element in the list is an \code{ExpressionSet}
}
\usage{data(expressionSetList)}

\details{

  Parmigiani et al. (2004) performed a cross-study analysis of three
  lung cancer studies.  The studies used in this analysis were merged by
  UniGene identifiers to obtain a set of 3,171 gene.  The R experiment
  data package \code{lungExpression} that was developed to facilitate
  the reproducibility of this analysis contains the three studies as
  ExpressionSets.  Here, we take a random sample of 500 features from
  one study (the "stanford" study), and split this study into three
  artificial studies that each contain 4 squamous carcinomas and 3
  adenocarcinomas.  The three artificial studies are then used to create
  an instance of the \code{ExpressionSetList} class.

  See Garber et al. (2001) for the raw data and description of the
  \code{stanford} study.
     
}
\source{

  The experiment data package lungExpression (\url{www.bioconductor.org})
  
}
\references{

    Parmigiani et al. (2004) A cross-study comparison of gene expression
  studies for the molecular classification of lung cancer, Clin Cancer
  Res, 10(9): 2922-2927

  Garber et al. (2001) Diversity of gene expression in adenocarcinoma of
  the lung, PNAS, 98:13784-13789

}

\examples{
data(expressionSetList)
}
\keyword{datasets}
