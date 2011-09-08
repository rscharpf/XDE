.onAttach <- function(libname, pkgname) {
	version <- packageDescription("XDE", field="Version")
	packageStartupMessage("Welcome to XDE version ", version)
}
