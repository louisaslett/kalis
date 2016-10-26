SetGenWorkingDir <- function(working.dir) {
  if(is.na(file.info(working.dir)$isdir) || !file.info(working.dir)$isdir) {
    stop(working.dir, " is not a directory.")
  }
  assign("working.dir", working.dir, envir = pkgCache)
}
GetGenWorkingDir <- function() {
  get("working.dir", envir = pkgCache)
}
