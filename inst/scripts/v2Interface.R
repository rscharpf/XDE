## mainTest_v2.cpp contains contains code Haakon used for testing v2
##
## Rinterface_v2.cpp contains all of the entry points from R to update each of the parameters.
##


##---------------------------------------------------------------------------
##
## initialize fixed hyerparameters
##
##---------------------------------------------------------------------------
library(XDE)
data(expressionSetList)
##hyper.params <- XDE:::getHyperparameters(G=3592, Q=3, S=c(100,50,75))
hyper.params <- HyperParams(object=expressionSetList, clinicalCovariate="adenoVsquamous")
params <- Params(hyper.params)

## note that some values have changed, but nAccept is 'zero'
## Also, the seed has not changed
trace(rupdateANu, browser)
res <- rupdateANu(expressionSetList,
		  hyper.params,
		  params,
		  nTry=1000L,
		  epsilon=0.05)

trace(updateBDDelta, browser)
res <- rupdateBDDelta(expressionSetList,
		     hyper.params,
		     params,
		     nTry=100000L,
		     epsilon=.2)
















