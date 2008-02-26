\name{ExpressionSetList-class}
\docType{class}
\alias{ExpressionSetList-class}
\alias{.integrativeCorrelationFilter,ExpressionSetList-method}
\alias{[,ExpressionSetList-method}
\alias{coerce,list,ExpressionSetList-method}
\alias{featureNames,ExpressionSetList-method}
\alias{geneCenter,ExpressionSetList-method}
\alias{nSamples}
\alias{nSamples,ExpressionSetList-method}
\alias{nrow,ExpressionSetList-method}
\alias{pData,ExpressionSetList-method}
\alias{pca,ExpressionSetList-method}
\alias{standardizeSamples,ExpressionSetList-method}
\alias{studyCenter,ExpressionSetList-method}
\alias{zeroNu,ExpressionSetList-method}

\title{A class for containing a list of ExpressionSets}

\description{Each element in the list must be a valid
  \code{ExpressionSet}.  The \code{featureNames} must be identical for
  each \code{ExpressionSet}. }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ExpressionSetList", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"}}
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{AssayData}"}, by class "list", distance 2.
}
\section{Methods}{
  \describe{
    
    \item{.integrativeCorrelationFilter}{\code{signature(object =
  "ExpressionSetList")} Experimental function for filtering an arbitrary
  list of ExpressionSets by integrative correlation.  Genes are excluded
  that do not exceed the fdr threshold in at least 1 of the studies.}

\item{"["}{\code{signature(x = "ExpressionSetList")} Subsets each
  \code{ExpressionSet} element in the list.}

    
\item{coerce}{\code{signature(from = "list", to = "ExpressionSetList")}
  Coerces a list of \code{ExpressionSet} objects to an object of class
  \code{ExpressionSetList}. The \code{validityMethod} for the
  \code{ExpressionSetList} class will return an error if the
  \code{featureNames} for each \code{ExpressionSet} are not identical.}

\item{featureNames}{\code{signature(object = "ExpressionSetList")}
  Accessor for the \code{featureNames}}

\item{geneCenter}{\code{signature(object = "ExpressionSetList")} See \code{\link{geneCenter}}}

\item{nSamples}{\code{signature(x = "ExpressionSetList")} Numerical
  vector giving the number of samples in each \code{ExpressionSet}}

\item{nrow}{\code{signature(x = "ExpressionSetList")} Numerical: number
  of features or genes}

\item{pData}{\code{signature(object = "ExpressionSetList")} returns a
  list of \code{data.frames}.  The elements of the list correspond to
  the studies in the \code{ExpressionSetList} object.}

\item{.pca}{\code{signature(object = "ExpressionSetList")} Runs
  principal components to generate cross-study summary scores of
  differential expression in multiple studies.  Not meant to be called
  directly.  See \code{\link{xsScores}}.}

\item{standardizeSamples}{\code{signature(object = "ExpressionSetList")}
  See \code{\link{standardizeSamples}}
}

\item{studyCenter}{\code{signature(object = "ExpressionSetList")}
  See \code{\link{studyCenter}}
}

\item{zeroNu}{\code{signature(object = "ExpressionSetList")}
  See \code{`\link{zeroNu}}. }
}
}

\author{R. Scharpf}

\seealso{
  \code{\link{XdeMcmc-class}}, \code{\link{XdeParameter-class}}
}
\examples{
showClass("ExpressionSetList")
data(expressionSetList)
}
\keyword{classes}
