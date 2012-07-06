library(XDE)
data(expressionSetList)
## some initial values
params <- XDE:::Parameters(object=expressionSetList,
			   clinicalCovariate="adenoVsquamous")
dms <- dims(params)
params2 <- XDE:::startingValues(params,
				expressionSetList,
				"adenoVsquamous")
orig <- params2
dms2 <- dims(params2)
one.delta <- FALSE
all.equal(dms, dms2)
##for(i in 1:50){
	params2 <- rupdateANu(object=params2,
			      nTry=5L,
			      epsilon=0.1,
			      dryrun=FALSE)

	params2 <- rupdateBDDelta(object=params2,
				  nTry=5L,
				  epsilon=0.1,
				  dryrun=FALSE) 
                                      
	params2 <- rupdateTau2RhoNu(object=params2,
				    nTry=5L,
				    epsilon=0.1,
				    dryrun=FALSE)
                                      
	params2 <- rupdateTau2RDDelta(object=params2,
				      nTry=5L,
				      epsilon=0.1,
				      dryrun=FALSE) 
                                      
	params2 <- rupdateNu(object=params2,
			     dryrun=FALSE)
	## NaN's appear in nu after the update...?
	
  params2 <- rupdateDDelta(object=params2,
				 dryrun=FALSE) 

	params2 <- rupdateC2(object=params2,
			     nTry=5L,
			     dryrun=FALSE)

	params2 <- rupdateC2DDelta(object=params2,
				   nTry=5L,
				   dryrun=FALSE) 

	params2 <- rupdateGamma2(object=params2,
				 dryrun=FALSE)

	params2 <- rupdateGamma2Nu(object=params2,
				   nTry=5L,
				   epsilon=0.1,
				   dryrun=FALSE) 

	params2 <- rupdateRDDelta(object=params2,
				 nTry=5L,
				 epsilon=0.1,
				 dryrun=FALSE) 

	params2 <- rupdateRC2(object=params2,
				 nTry=5L,
				 epsilon=0.1,
				 dryrun=FALSE) 

	params2 <- rupdateRhoNu(object=params2,
				nTry=5L,
				epsilon=0.1,
				dryrun=FALSE)

	params2 <- rupdateRhoGamma2(object=params2,
				    nTry=5L,
				    epsilon=0.1,
				    dryrun=FALSE)

	params2 <- rupdateSigma2(object=params2,
				    nTry=5L,
				    epsilon=0.1,
				    dryrun=FALSE) 

	params2 <- rupdatePhi(object=params2,
				 nTry=5L,
				 epsilon=0.1,
				 dryrun=FALSE)

	params2 <- rupdateTheta(object=params2,
				 nTry=5L,
				 epsilon=0.1,
				 dryrun=FALSE)

	params2 <- rupdateLambda(object=params2,
				 nTry=5L,
				 epsilon=0.1,
				 dryrun=FALSE)
                                      
	params2 <- rupdateT(object=params2,
				 nTry=5L,
				 epsilon=0.1,
				 dryrun=FALSE)

	params2 <- rupdateL(object=params2,
				 nTry=5L,
				 epsilon=0.1,
				 dryrun=FALSE)

	if(!one.delta){
		params2 <- rupdateXi(object=params2,
				     dryrun=FALSE,
				     one.delta=FALSE) 
		params2 <- rupdateDeltaDDelta(object=params2,
					      nTry=5L,
					      dryrun=FALSE,
					      one.delta=FALSE) 
    
    ##nipun: for  new wrappers{
		params2 <- rupdateDeltaDDelta_MRF2(object=params2,
		                              nTry=5L,
		                              dryrun=FALSE,
		                              one.delta=FALSE)
    
    params2 <- rupdateDelta(object=params2,
                            nTry=5L,
                            dryrun=FALSE,
                            one.delta=FALSE)
    
    params2 <- rupdateDelta_MRF2(object=params2,
                                 nTry=5L,
                                 dryrun=FALSE,
                                 one.delta=FALSE)
    ##nipun: }
	}
	if(one.delta){
		## must make sure that the deltas
		## for a gene are the same
		params2 <- rupdateXi(object=params2,
				     dryrun=FALSE,
				     one.delta=TRUE) 
		params2 <- rupdateDeltaDDelta(object=params2,
					      nTry=5L,
					      dryrun=FALSE,
					      one.delta=TRUE) 
    
		##nipun: for new wrappers{
		params2 <- rupdateDeltaDDelta_MRF2(object=params2,
		                                   nTry=5L,
		                                   dryrun=FALSE,
		                                   one.delta=TRUE)
    
    params2 <- rupdateDelta(object=params2,
                            nTry=5L,
                            dryrun=FALSE,
                            one.delta=TRUE)
    
    params2 <- rupdateDelta_MRF2(object=params2,
                                 nTry=5L,
                                 dryrun=FALSE,
                                 one.delta=TRUE)
    ##nipun: }
	}
                                      
	params2 <- rupdateLSigma2(object=params2,
				  nTry=5L,
				  epsilon=0.1,
				  dryrun=FALSE) 

	params2 <- rupdateTSigma2(object=params2,
				      nTry=5L,
				      epsilon=0.1,
				      dryrun=FALSE) 

	params2 <- rupdateLambdaPhi(object=params2,
				      nTry=5L,
				      epsilon=0.1,
				      dryrun=FALSE)

	params2 <- rupdateThetaPhi(object=params2,
				      nTry=5L,
				      epsilon=0.1,
				      dryrun=FALSE)
                                      
	##nipun: for new wrappers{			                 
  params2 <- rupdateAlpha(object=params2,
				      nTry=5L,
				      epsilon=0.1,
				      dryrun=FALSE)
                                      
  params2 <- rupdateAlpha_onedelta(object=params2,
              nTry=5L,
              epsilon=0.1,
              dryrun=FALSE)                                      
	
  params2 <- rupdateBeta(object=params2,
                         nTry=5L,
                         epsilon=0.1,
                         dryrun=FALSE)
                                      
  params2 <- rupdateBeta_onedelta(object=params2,
                                  nTry=5L,
                                  epsilon=0.1,
                                  dryrun=FALSE)
                                    
  params2 <- rupdateBetag(object=params2,
                          nTry=5L,
                          epsilon=0.1,
                          dryrun=FALSE)     
                                      
  params2 <- rupdateSigma2_HyperInverseWishart(object=params2,
                                               nTry=5L,
                                               epsilon=0.1,
                                               dryrun=FALSE)
                                      
  params2 <- rupdateLSigma2_HyperInverseWishart(object=params2,
                                              nTry=5L,
                                              epsilon=0.1,
                                              dryrun=FALSE) 
                                      
  params2 <- rupdateTSigma2_HyperInverseWishart(object=params2,
                                                nTry=5L,
                                                epsilon=0.1,
                                                dryrun=FALSE)
                                      
  params2 <- rupdateOmega_HyperInverseWishart(object=params2,
                                              nTry=5L,
                                              epsilon=0.1,
                                              dryrun=FALSE)
                                      
  params2 <- rupdateDDeltaStar_HyperInverseWishart(object=params2,
                                                   nTry=5L,
                                                   epsilon=0.1,
                                                   dryrun=FALSE)
                                      
  params2 <- rupdateTau2RDDeltaStar_HyperInverseWishart(object=params2,
                                                      nTry=5L,
                                                      epsilon=0.1,
                                                      dryrun=FALSE)  
                                      
  params2 <- rupdateBDDeltaStar_HyperInverseWishart(object=params2,
                                                    nTry=5L,
                                                    epsilon=0.1,
                                                    dryrun=FALSE)
                                      
  params2 <- rupdateRDDeltaStar_HyperInverseWishart(object=params2,
                                                    nTry=5L,
                                                    epsilon=0.1,
                                                    dryrun=FALSE)                              
  ##nipun: }
                                      
	if(FALSE){
		## add appropriate slots to Parameters object
    
    params2 <- rupdateDeltaDDelta_MRF(object=params2,
						  nTry=5L,
					  epsilon=0.1,
						  dryrun=FALSE) 
						    

		params2 <- rupdateAlpha_MRF(object=params2,
					    nTry=5L,
					    epsilon=0.1,
					    dryrun=FALSE) 
					      
                                                                    
		params2 <- rupdateBeta_MRF(object=params2,
					   nTry=5L,
					   epsilon=0.1,
					   dryrun=FALSE)
                                                                    
		params2 <- rupdateBetag_MRF(object=params2,
					    nTry=5L,
					    epsilon=0.1,
					    dryrun=FALSE)
	}
##}//end for(i in 1:50) 
















##
	#rupdate <- function(object,
	#		    type=c("ANu",
	#		           "BDDelta",
	#		           "Tau2Rho",
	#		           "Tau2DDelta"),
	#		    nTry=10L,
	#		    epsilon=0.1,
	#		    dryrun=FALSE){
	#	update.type <- paste("rupdate", type, sep="")
	#	.C(update.type, ...)
