\name{HyperParams-class}
\Rdversion{1.1}
\docType{class}
\alias{HyperParams-class}
\alias{Params,HyperParams-method}
\title{Class \code{"HyperParams"}}
\description{
  Class for storing hyperparameters.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form HyperParams(object, ...)
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{AssayData}"}, by class "list", distance 2.
}
\section{Methods}{
No methods defined with class "HyperParams" in the signature.
}

\section{Details}{

  Elements stored in this class are fixed. See the \code{Params} class
  for storing parameters that are updated from MCMC sampling.

}

\author{
R. Scharpf
}


\seealso{
  \code{\linkS4class{Params}}
}
\examples{
showClass("HyperParams")
}
\keyword{classes}
