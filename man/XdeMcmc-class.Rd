\name{XdeMcmc-class}
\docType{class}
\alias{XdeMcmc-class}
\alias{$}
\alias{$,XdeMcmc-method}
\alias{.standardizedDelta,XdeMcmc-method}
\alias{bayesianEffectSize}
\alias{bayesianEffectSize,XdeMcmc-method}
\alias{bayesianEffectSize<-}
\alias{bayesianEffectSize<-,XdeMcmc,matrix-method}
\alias{calculateBayesianEffectSize}
\alias{calculateBayesianEffectSize,XdeMcmc-method}
\alias{calculatePosteriorAvg,XdeMcmc-method}
\alias{directory}
\alias{directory,XdeMcmc-method}
\alias{featureNames,XdeMcmc-method}
\alias{initialize,XdeMcmc-method}
\alias{iterations,XdeMcmc-method}
\alias{lastMcmc,XdeMcmc-method}
\alias{nrow,XdeMcmc-method}
\alias{output,XdeMcmc-method}
\alias{plot,XdeMcmc-method}
\alias{posteriorAvg,XdeMcmc-method}
\alias{posteriorAvg<-,XdeMcmc,matrix-method}
\alias{seed,XdeMcmc-method}
\alias{show,XdeMcmc-method}
\alias{studyNames,XdeMcmc-method}

\title{Class for storing output from the Bayesian model}
\description{Stores output, including the last iteration of the MCMC.}

\section{Objects from the Class}{
  
  Objects can be created by calls of the form \code{new("XdeMcmc",
  studyNames, featureNames, iterations, seed, output, directory,
  lastMcmc, posteriorAvg, bayesianEffectSize)}.

}

\section{Slots}{
  \describe{
    \item{\code{studyNames}:}{Object of class \code{"character"} }
    \item{\code{featureNames}:}{Object of class \code{"character"} }
    \item{\code{iterations}:}{Object of class \code{"numeric"} }
    \item{\code{directory}:}{Object of class \code{"character"} }
    \item{\code{seed}:}{Object of class \code{"integer"} }
    \item{\code{output}:}{Object of class \code{"numeric"} }
    \item{\code{lastMcmc}:}{Object of class \code{"environment"} }
    \item{\code{posteriorAvg}:}{Object of class \code{"NULLorMatrix"} }
    \item{\code{bayesianEffectSize}:}{Object of class \code{"NULLorMatrix"} }
  }
}
\section{Methods}{
  \describe{
    \item{\code{$}}{\code{signature(x = "XdeMcmc")} }
    \item{.standardizedDelta}{\code{signature(object = "XdeMcmc")}}
    \item{bayesianEffectSize}{\code{signature(object = "XdeMcmc")} }
    \item{bayesianEffectSize<-}{\code{signature(object = "XdeMcmc",
	value = "matrix")} }
    \item{calculatePosteriorAvg}{\code{signature(object = "XdeMcmc")}:
  See \code{\link{calculatePosteriorAvg}}}
    \item{directory}{\code{signature(object = "XdeMcmc")} }
    \item{featureNames}{\code{signature(object = "XdeMcmc")} }
    \item{initialize}{\code{signature(.Object = "XdeMcmc")} }
    \item{iterations}{\code{signature(object = "XdeMcmc")} }
    \item{lastMcmc}{\code{signature(object = "XdeMcmc")} }
    \item{nrow}{\code{signature(x = "XdeMcmc")} }
    \item{output}{\code{signature(object = "XdeMcmc")} }
    \item{plot}{\code{signature(x = "XdeMcmc")} }
    \item{posteriorAvg}{\code{signature(object = "XdeMcmc")} }
    \item{seed}{\code{signature(object = "XdeMcmc")} }
    \item{show}{\code{signature(object = "XdeMcmc")} }
    \item{studyNames}{\code{signature(object = "XdeMcmc")} }
  }
}
\author{R. Scharpf}

\seealso{
  The class for storing the data: \code{\link{ExpressionSetList-class}}
  and the class that contains default options for fitting the Bayesian
  model: \code{\link{XdeParameter-class}}
}
      
\examples{
##See XDE vignette:
\dontrun{
openVignette(package="XDE")
}
}
\keyword{classes}
