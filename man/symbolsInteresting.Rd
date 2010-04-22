\name{symbolsInteresting}
\alias{symbolsInteresting}
\title{Useful for changing the look of pairs plots to emphasize
  concordant or discordant genes}
\description{

  This function can be used to order genes in a matrix by the rank of a
  statistic and provide different plotting symbols and colors for genes
  that exceed a certain threshold of the ranking statistic.
  
}

\usage{
symbolsInteresting(rankingStatistic, percentile = 0.9, colors = c("grey50", "royalblue"), symbols = c(".", "o"), size = c(3, 1), background = c("white", "grey70"))
}
\arguments{
  \item{rankingStatistic}{ Any numerical vector }
  \item{percentile}{A percentile of the rankingStatistic -- above which
    a gene would be classified as 'interesting'}
  \item{colors}{character string of length 2: a color for genes not exceeding the
    percentile and a color for genes exceeding the thresold}
  \item{symbols}{two plotting symbols (numeric or character): symbol for
  genes not exceeding percentile and symbol for genes exceeding percentile}
  \item{size}{numeric vector of length 2: size of plotting symbol for genes not exceeding
    percentile and size of plotting symbol for genes exceeding percentile }
  \item{background}{character vector of length 2: background color of
    plotting symbols for gene not exceeding percentile and for genes
    exceeding the percentile}
}
\value{
  \item{order}{the order of the rankingStatistic}
  \item{pch}{plotting symbols (same length as rankingStatistic)}
  \item{col}{color of plotting symbols (same length as
    rankingStatistic)}
  \item{bg}{background color of plotting symbols (same length as
    rankingStatistic)}
  \item{cex}{size of plotting symbols (same length as
    rankingStatistic)}  
}
\author{R. Scharpf}
\examples{

  data(expressionSetList)
  data(xmcmc)
  pathToLogFiles <- system.file("logFiles", package="XDE")
  load(file.path(pathToLogFiles, "BES.rda"))
  load(file.path(pathToLogFiles, "postAvg.rda"))
  op.conc <- symbolsInteresting(rankingStatistic=postAvg[, "concordant"])
  graphics:::pairs(BES[op.conc$order, ], pch=op.conc$pch, col=op.conc$col, 
                   bg=op.conc$bg, upper.panel=NULL, cex=op.conc$cex)  

}
\keyword{dplot}

