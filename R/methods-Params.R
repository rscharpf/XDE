setMethod("vectorize", signature(object="Params"), function(object){
	## gene index changes fastest, then study index
	params$nu <- vectorize(params$nu)
	params$delta <- vectorize(params$delta)
	params$Delta <- vectorize(params$Delta)
	params$phi <- vectorize(params$phi)
	rho <- params$rho
	rho <- rho[upper.tri(rho)]
	params$rho <- rho
	params$sigma2 <- vectorize(params$sigma2)
	return(params)
})
