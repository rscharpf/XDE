library(Test)
##
## break mainTestHyper_v2.cpp:200
set.seed(1)
trace(initializeParams, browser)
initializeParams(G=100)
## add arguments for number genes, etc. so that these are not hardcoded
##
## have updated parameters written to files
##
## Can only pass 65 arguments
## TODO
##
## - change initializeParams such that G, Q, psi (clinical
##   parameters), etc. are inputs
##
## - have an option that selects which parameters are
##   written to file
##
## - modify such that mcmc samples are appended to files according to
##  previous option (see how sigma2 and phi are currently written to
##  file, and adopt this for other parameters)
##
## - modify such that only selected iterations are appended to file
##
## - modify such that data is passed to the function and not simulated
## - (or have an option that allows passing data)
##
## - many of the specific update types are already available with R wrappers
##   (see R/Rupdates.R)
##   -> initializeParams could be changed to only provide initial parameters
##   -> using .C, we would need to specify the correct size for each parameter as input
##      to initializeParams function
##




