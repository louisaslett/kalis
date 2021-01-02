.onAttach <- function(libname, pkgname){
  packageStartupMessage(.Call(CCall_ComputeStatus))
}
