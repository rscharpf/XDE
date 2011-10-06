\name{Params-class}
\Rdversion{1.1}
\docType{class}
\alias{Params-class}

\title{Class \code{"Params"}}
\description{
  Class for storing parameters for the MCMC chain.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{Params(hyper.params,...)}.

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
No methods defined with class "Params" in the signature.
}
\author{
R. Scharpf
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\linkS4class{HyperParams}}
}
\examples{
showClass("Params")
}
\keyword{classes}
