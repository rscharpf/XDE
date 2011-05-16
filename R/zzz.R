.onAttach <- function(libname, pkgname) {
	version <- packageDescription("XDE", field="Version")
	message("Welcome to XDE version ", version)
}
