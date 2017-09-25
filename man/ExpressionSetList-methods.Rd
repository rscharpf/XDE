\name{ExpressionSetList-methods}
%\docType{methods}
\alias{ExpressionSetList-methods}
%\alias{[,ExpressionSetList-method}
%\alias{coerce,list,ExpressionSetList-method}
%\alias{dim,ExpressionSetList-method}
%\alias{lapply,ExpressionSetList-method}
\alias{[,ExpressionSetList,ANY,ANY,ANY-method}
\alias{phenotype}
\title{Methods for ExpressionSetList}
\description{
  Methods for objects of class \code{ExpressionSetList}.
}

\usage{
  phenotype(object, varLabel)
}



\arguments{

  \item{object}{ A \code{ExpressionSetList}.}

  \item{varLabel}{\code{character}.  Name of the clinical variable.}

}


\value{

  \code{phenotype} returns a matrix of the clinical variable where each
  column is a study.  We require that the clinical variable have the
  same name in each study (each element of the \code{ExpressionSetList}
  object) and that the clinical variable is binary with values 1 or 0.


}

\keyword{methods}
\keyword{manip}
