THISPKG <- "XDE"

.onAttach <- function(libname, pkgname) {
	version <- packageDescription("XDE", fields="Version")
	packageStartupMessage("Welcome to XDE version ", version)
}

.onUnload <- function( libpath ){
	library.dynam.unload(THISPKG, libpath)
}
