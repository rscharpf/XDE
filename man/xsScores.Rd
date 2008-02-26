\name{xsScores}
\alias{xsScores}
\title{Alternative cross-study scores of differential expression}

\description{
  Alternative cross-study scores of differential expression
}

\usage{
xsScores(statistic, N)
}

\arguments{
  \item{statistic}{a matrix of study-specific estimates of effect size.
  Rows are genes and columns are studies.}
  \item{N}{numerical vector: the number of samples in each study (the
    length should be the number of columns in \code{statistic})}
}

\value{
  numeric
}

\references{

  J.K. Choi, U. Yu, S. Kim, and O.J. Yoo (2003), Combining multiple
  microarray studies and modeling interstudy variation, Bioinformatics,
  19(1) I84-I90.  

  E. Garrett-Mayer, G. Parmigiani, X. Zhong, L. Cope, and E. Gabrielson
	(2007), Cross-study validation and combined analysis of gene
	expression microarray data, Biostatistics, September

  R. Scharpf et al., A Bayesian Model for Cross-Study Differential Gene
  Expression, Technical Report 158, Johns Hopkins University, Department
  of Biostatistics, 2007

}

\author{R. Scharpf}

\seealso{the GeneMeta package, \code{\link{ssStatistic}}}

\examples{
  data(expressionSetList)
  t <- ssStatistic(statistic="t", phenotypeLabel="adenoVsquamous", esetList=expressionSetList)
  tScores <- xsScores(t, N=nSamples(expressionSetList))
}
  
\keyword{htest}
\keyword{models}

