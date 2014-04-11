\name{xmcmc}
\alias{xmcmc}
\docType{data}
\title{Object of class XdeMcmc}
\description{
  An object of class XdeMcmc is created by fitting the Bayesian
  hierarchical model to the \code{expressionSetList} example data.
}
\usage{data(xmcmc)}

\details{
  The \code{xmcmc} data example was obtained as described in the
  XDE vignette.
}

\examples{

data(xmcmc)
xmcmc

##ordinarily, one should not need to change the directory in an object
##of class XdeMcmc -- therefore, a replacment method is not defined
pathToLogFiles <- system.file("logFiles", package="XDE")
xmcmc@directory <- pathToLogFiles

##The $ operator can be used to extract chains. For instance, here we
##extract the c2 chain
c2 <- xmcmc$c2
plot.ts(c2)
}
\keyword{datasets}
