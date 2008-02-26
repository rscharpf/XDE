library(XDE)
data(expressionSetList)
xlist <- expressionSetList
params <- new("XdeParameter", esetList=xlist, phenotypeLabel="adenoVsquamous")
iterations(params) <- 2000
##iterations(params) <- 10
thin(params) <- 2
burnin(params) <- FALSE

##do not save any of the parameters indexed by g and p
output(params)[c("potential", "acceptance", "nu", ##"DDelta",
                 "probDelta", "phi")] <- 0
##                 "delta", "sigma2"] <- 0
output(params)
directory(params) <- "/home/bst/student/rscharpf/projects/sandbox/xdeBranch/inst/logFiles"
if(file.exists(directory(params))){
  system.time(xmcmc <- xde(params, xlist))
  xmcmc@lastMcmc <- environment() ##save memory
  bayesianEffectSize(xmcmc) <- calculateBayesianEffectSize(xmcmc)
  unlink("../inst/logFiles/Delta.log")
  unlink("../inst/logFiles/delta.log")
  unlink("../inst/logFiles/sigma2.log")  
  save(xmcmc, file="~/projects/sandbox/xdeBranch/data/xmcmc.RData", compress=TRUE)
}

