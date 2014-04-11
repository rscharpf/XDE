\name{calculatePosteriorAvg}
\alias{calculatePosteriorAvg}
\title{Calculate the posterior average for indicators of concordant and
  discordant differential expression}

\description{
  
  This function calculates the posterior average for indicators of concordant and
  discordant differential expression from the saved log files. See details.  
  
}

\usage{
calculatePosteriorAvg(object, NCONC=2, NDIFF=1, burnin=0)
}

\arguments{
  \item{object}{ Object of class \code{XdeMcmc}}
  \item{NCONC}{Integer: number of studies for which the gene must be
    differentially expressed (in the same direction) to be classified as
  concordant differential expression}
  \item{NDIFF}{Integer: number of studies for which a gene must be
    up- or down-regulated to be classified as differentially expressed.
    It is the union of concordant and discordant differential
    expression.}
  \item{burnin}{Integer: number of MCMC iterations for the burnin. Posterior means
  are computed from the MCMC samples following burnin.}
}
\details{

 For each iteration,

1.  calculate the sign of delta * Delta
 
2.  For each gene, compute the number of positive signs (P) and the
number of negative signs (N) (a G x 2 matrix, where G is the number of
genes in common across all studies).  P + N is <= S, where S is the
number of studies.

3.  for a given gene, the discordant indicator is simply when P * N
is nonzero.

4.  The concordant indicator requires P * N = 0 AND P + N >=
NCONC, where NCONC is specified by the user.

5.  differential expression is simply   | P | + | N |  >= NDIFF.  By
default, NDIFF is 1 but can be user-specified.

The posterior average is then computed from the mean over all MCMC iterations.

}

\value{
  
  A G x 3 matrix.
  
}

\author{RS}
\seealso{\code{\link{posteriorAvg}}}
\keyword{manip}

