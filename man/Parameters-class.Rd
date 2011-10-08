\name{Parameters-class}
\Rdversion{1.1}
\docType{class}
\alias{Parameters-class}
\alias{[[<-,Parameters-method}
\alias{[[,Parameters-method}
\alias{$<-,Parameters-method}
\alias{$,Parameters-method}
\alias{coerce,XdeParameter,Parameters-method}
\alias{show,Parameters-method}

\title{Container for XDE parameters}
\description{Container for XDE parameters}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Parameters", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{seed}:}{Object of class \code{"integer"} ~~ }
    \item{\code{data}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{phenodata}:}{Object of class \code{"integer"} ~~ }
    \item{\code{G}:}{Object of class \code{"integer"} ~~ }
    \item{\code{Q}:}{Object of class \code{"integer"} ~~ }
    \item{\code{S}:}{Object of class \code{"integer"} ~~ }
    \item{\code{alphaA}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{alphaB}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{betaA}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{betaB}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{pA0}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{pA1}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{pB0}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{pB1}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nuR}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nuRho}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{alphaXi}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{betaXi}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{c2Max}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{alphaEta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{betaEta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{pOmega0}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{lambdaOmega}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{lambdaKappa}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{gamma2}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{c2}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{tau2Rho}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{tau2R}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{a}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{b}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{l}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{t}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{lambda}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{theta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{phi}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{sigma2}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{r}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{rho}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nu}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{delta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Delta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{xi}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{"[[<-"}{\code{signature(x = "Parameters")}: ... }
    \item{"[["}{\code{signature(x = "Parameters")}: ... }
    \item{"$<-"}{\code{signature(x = "Parameters")}: ... }
    \item{$}{\code{signature(x = "Parameters")}: ... }
    \item{coerce}{\code{signature(from = "XdeParameter", to = "Parameters")}: ... }
    \item{show}{\code{signature(object = "Parameters")}: ... }
	 }
}

\examples{
showClass("Parameters")
}
\keyword{classes}
