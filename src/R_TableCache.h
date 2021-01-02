#ifndef R_TABLECACHE_H
#define R_TABLECACHE_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP CopyFBTable(SEXP Rto, SEXP Rfrom);

#endif
