\name{lastMcmc}
\alias{lastMcmc}
\alias{lastMcmc<-}
\title{MCMC values for the last iteration}
\description{
  MCMC values for the last iteration.  Useful if more iterations are needed. 
}

\usage{
  lastMcmc(object)
}

\arguments{
  \item{object}{Object of class \code{XdeMcmc}}
}

\value{
  An environment.
}

\author{R. Scharpf}

\seealso{\code{\link{firstMcmc}}}

\examples{
  \dontrun{
    data(expressionSetList)
    xparam <- new("XdeParameter", phenotypeLabel="adenoVsquamous",
                  esetList=expressionSetList)
    iterations(xparam) <- 10  
    fit <- xde(xparam, esetList=expressionSetList)
    ##Do more iterations and use a different seed
    firstMcmc(xparam) <- lastMcmc(fit)
    seed(xparam) <- 97814
    fit2 <- xde(xparam, esetList=expressionSetList)

    ##Or
    fit2 <- xde(xparam, esetList=expressionSetList, outputMcmc=fit)
  }
}
\keyword{methods}

