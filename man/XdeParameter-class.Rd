\name{XdeParameter-class}
\docType{class}
\alias{XdeParameter-class}
\alias{burnin,XdeParameter-method}
\alias{burnin<-,XdeParameter,logical-method}
\alias{coerce,XdeParameter,Params-method}
\alias{directory,XdeParameter-method}
\alias{directory<-}
\alias{directory<-,XdeParameter-method}
\alias{firstMcmc,XdeParameter-method}
\alias{firstMcmc<-,XdeParameter,environment-method}
\alias{firstMcmc<-,XdeParameter,list-method}
\alias{hyperparameters,XdeParameter-method}
\alias{hyperparameters<-,XdeParameter-method}
\alias{initialize,XdeParameter-method}
\alias{iterations,XdeParameter-method}
\alias{iterations<-,XdeParameter,numeric-method}
\alias{iterations<-,XdeParameter,integer-method}
\alias{output,XdeParameter-method}
\alias{output<-,XdeParameter-method}
\alias{phenotypeLabel}
\alias{phenotypeLabel<-}
\alias{phenotypeLabel,XdeParameter-method}
\alias{phenotypeLabel<-,XdeParameter,character-method}
\alias{savedIterations}
\alias{savedIterations,XdeParameter-method}
\alias{seed,XdeParameter-method}
\alias{seed<-,XdeParameter,integer-method}
\alias{seed<-,XdeParameter,numeric-method}
\alias{show,XdeParameter-method}
\alias{showIterations}
\alias{showIterations<-}
\alias{showIterations,XdeParameter-method}
\alias{showIterations<-,XdeParameter-method}
\alias{studyNames}
\alias{studyNames<-}
\alias{studyNames,XdeParameter-method}
\alias{studyNames<-,XdeParameter-method}
\alias{thin,XdeParameter-method}
\alias{thin<-,XdeParameter,numeric-method}
\alias{tuning,XdeParameter-method}
\alias{tuning<-,XdeParameter-method}
\alias{updates,XdeParameter-method}
\alias{updates<-,XdeParameter-method}

\title{Container class for storing options of the Bayesian hierarchical model}
\description{This class contains initial values for the first iteration
  of the MCMC, options for saving MCMC chains, options for changing the
  tuning parameters of the Metropolis-Hastings algorithm, options for
  changing hyperparameters from their defaults, etc.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("XdeParameter",
  esetList, updates, tuning, hyperparameters, output, iterations,
  burnin, seed, randomSeed, genes, studies, firstMcmc,
  specifiedInitialValues, directory, phenotypeLabel, seed, showIterations,
  verbose, studyNames, one.delta)}.
}
\section{Slots}{
	 \describe{
    \item{\code{updates}:}{Object of class \code{numeric}. The frequency
    of  updates for each iteration of the chain.}
    \item{\code{tuning}:}{Object of class \code{numeric}.  Tuning
      parameters for the Metropolis-Hastings proposals }
    \item{\code{hyperparameters}:}{Object of class
      \code{numeric}. Hyperparameters for the Bayesian hierarchical model}
    \item{\code{output}:}{Object of class \code{numeric}. Indicator for
      whether to save the MCMC chain to file.  If the value is zero, the
    chain is not saved.}
    \item{\code{iterations}:}{Object of class \code{numeric}.  The total
    number of MCMC iterations.}
    \item{\code{burnin}:}{Object of class \code{logical}.  If set to
      FALSE, by default none of the chains will be saved (called for its
      side-effect of setting the output to zero for each parameter).}

    \item{\code{notes}:}{Object of class \code{character}. }
    \item{\code{firstMcmc}:}{Object of class \code{environment}. Values
      for the first iteration of the MCMC }
    \item{\code{seed}:}{Object of class \code{integer}.  Seed used for
    simulating random numbers.}
    \item{\code{showIterations}:}{Object of class \code{logical}.
      Whether to show the MCMC iteration when fitting the model}
    \item{\code{specifiedInitialValues}:}{Object of class
      \code{logical}.  If TRUE (the default), the values stored in
      \code{firstMcmc} will be used for the first iteration of the MCMC. }
    \item{\code{directory}:}{Object of class \code{character}. Specifies
    where to write the log files}
    \item{\code{phenotypeLabel}:}{Object of class \code{character}. The
      name of the binary covariate used for differential expression }
    \item{\code{verbose}:}{Object of class \code{logical} }
    \item{\code{studyNames}:}{Object of class \code{character}. Names of
      the datasets}
    \item{\code{one.delta}:}{Logical.  If TRUE, a gene is assumed to be
      differentially in all studies or none of the studies.}
  }
}
\section{Methods}{
  \describe{
    \item{burnin}{\code{signature(object = "XdeParameter")} logical.
      See \code{\link{burnin}}}
    \item{burnin<-}{\code{signature(object = "XdeParameter", value =
	"logical")} logical. See \code{\link{burnin}}}
    \item{directory}{\code{signature(object = "XdeParameter")} character
    string giving the path or relative path to store log files from the
    MCMC chain}
    \item{directory<-}{\code{signature(object = "XdeParameter")}Path to
      store log files. }
    \item{firstMcmc}{\code{signature(object = "XdeParameter")} See \code{\link{firstMcmc}}}
    \item{firstMcmc<-}{\code{signature(object = "XdeParameter", value = "environment")} }
    \item{firstMcmc<-}{\code{signature(object = "XdeParameter", value = "list")} }
    \item{hyperparameters}{\code{signature(object = "XdeParameter")} See
    the XdeParameterClass vignette}
    \item{hyperparameters<-}{\code{signature(object = "XdeParameter")} See
    the XdeParameterClass vignette}
    \item{initialize}{\code{signature(.Object = "XdeParameter")} Method
      for initializing an instance of the class.  The default values
      provided work well in most cases.}
    \item{iterations}{\code{signature(object = "XdeParameter")} Accessor
      for the total number of MCMC iterations to run}
    \item{iterations<-}{\code{signature(object = "XdeParameter", value =
	numeric)} The replacement method is useful for setting a
      different number of iterations.}
    \item{iterations<-}{\code{signature(object = "XdeParameter", value = "integer")} }
    \item{output}{\code{signature(object = "XdeParameter")} See also
      \code{\link{output}}. This method is also defined for class \code{XdeMcmc}}
    \item{output<-}{\code{signature(object = "XdeParameter")} See also \code{\link{output}} }
    \item{phenotypeLabel}{\code{signature(object = "XdeParameter")}
      The name of a binary covariate present in each study}
    \item{phenotypeLabel<-}{\code{signature(object = "XdeParameter", value = "character")} }
    \item{savedIterations}{\code{signature(object = "XdeParameter")} The
    number of MCMC iterations written to file.  It is the value of the
    total number of iterations divided by the thinning parameter.  See
    also \code{\link{output}}}
    \item{seed}{\code{signature(object = "XdeParameter")} See
      \code{\link{seed}}}
    \item{seed<-}{\code{signature(object = "XdeParameter",
    value="integer")} Replacement method. See also \code{\link{seed}}.}
    \item{show}{\code{signature(object = "XdeParameter")} Produces a
      short summary of objects that are instances of the
      \code{XdeParameter} class}
    \item{showIterations}{\code{signature(object = "XdeParameter")} logical}
    \item{showIterations<-}{\code{signature(object = "XdeParameter")} }
    \item{studyNames}{\code{signature(object = "XdeParameter")} Names of
    the high-throughput gene expression studies}
    \item{studyNames<-}{\code{signature(object = "XdeParameter")} }
    \item{thin}{\code{signature(x = "XdeParameter")} See
      \code{\link{output}} and \code{\link{thin}}}
    \item{thin<-}{\code{signature(x = "XdeParameter", value = numeric)}
      See \code{\link{thin}}}
    \item{tuning}{\code{signature(object = "XdeParameter")} See also \code{\link{tuning}}}
    \item{tuning<-}{\code{signature(object = "XdeParameter")} }
    \item{updates}{\code{signature(object = "XdeParameter")} See also \code{\link{updates}}}
    \item{updates<-}{\code{signature(object = "XdeParameter")} }
  }
}
\references{R. Scharpf}
\author{R. Scharpf}

\seealso{
  \code{\link{ExpressionSetList-class}}
}
\examples{
showClass("XdeParameter")
##See the XdeParameterClass vignette
}
\keyword{classes}
