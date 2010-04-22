\name{pairs-methods}
\docType{methods}
\alias{pairs-methods}
\alias{pairs,matrix-method}
\alias{pairs,data.frame-method}
\title{pairs function for high-throughput data}

\description{

  A convenient wrapper for pairs that uses smoothScatter to plot the
  density of the points and displays the spearman correlation
  coefficient of the pairwise scatterplots.

}
\section{Methods}{
\describe{

\item{x = "matrix"}{ Typically a matrix of effect size estimates
  obtained in each study.  Rows are genes, columns are studies.}

\item{x = "data.frame"}{Typically a \code{data.frame} of effect size
  estimates obtained in each study.  Rows are genes, columns are studies.}
}
}
\keyword{methods}
\keyword{hplot}
