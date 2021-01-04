#include "R_TableMaker.h"

#include "R_Kalis.h"
#include "Cache.h"



SEXP MakeForwardTable(SEXP Rfrom_r, SEXP Rto_r, SEXP Rpars_sha) {
  if(num_inds == 0 || hap_size == 0) {
    REprintf("Error: Cache is empty.\n");
    KALIS_RETURN;
  }
  if((size_t) INTEGER(Rfrom_r)[0] > num_inds || (size_t) INTEGER(Rto_r)[0] > num_inds || INTEGER(Rfrom_r)[0] > INTEGER(Rto_r)[0]) {
    REprintf("Error: Invalid from/to receipient arguments.\n");
    KALIS_RETURN;
  }

  SEXP Ralpha   = PROTECT(Rf_allocMatrix(REALSXP, num_inds, INTEGER(Rto_r)[0]-INTEGER(Rfrom_r)[0]+1));
  SEXP Ralpha_f = PROTECT(Rf_allocVector(REALSXP, INTEGER(Rto_r)[0]-INTEGER(Rfrom_r)[0]+1));
  SEXP Rl       = PROTECT(Rf_allocVector(INTSXP, 1));
  INTEGER(Rl)[0] = 2147483647;

  const char *names[] = {"alpha", "alpha.f", "l", "from_recipient", "to_recipient", "pars.sha256", ""};
  SEXP fwd = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(fwd, 0, Ralpha);
  SET_VECTOR_ELT(fwd, 1, Ralpha_f);
  SET_VECTOR_ELT(fwd, 2, Rl);
  SET_VECTOR_ELT(fwd, 3, Rfrom_r);
  SET_VECTOR_ELT(fwd, 4, Rto_r);
  SET_VECTOR_ELT(fwd, 5, Rpars_sha);

  SEXP class = PROTECT(Rf_mkString("kalisForwardTable"));
  Rf_setAttrib(fwd, R_ClassSymbol, class);

  UNPROTECT(5);
  return(fwd);
}



SEXP MakeBackwardTable(SEXP Rfrom_r, SEXP Rto_r, SEXP Rpars_sha) {
  if(num_inds == 0 || hap_size == 0) {
    REprintf("Error: Cache is empty.\n");
    KALIS_RETURN;
  }
  if((size_t) INTEGER(Rfrom_r)[0] > num_inds || (size_t) INTEGER(Rto_r)[0] > num_inds || INTEGER(Rfrom_r)[0] > INTEGER(Rto_r)[0]) {
    REprintf("Error: Invalid from/to receipient arguments.\n");
    KALIS_RETURN;
  }

  SEXP Rbeta   = PROTECT(Rf_allocMatrix(REALSXP, num_inds, INTEGER(Rto_r)[0]-INTEGER(Rfrom_r)[0]+1));
  SEXP Rbeta_g = PROTECT(Rf_allocVector(REALSXP, INTEGER(Rto_r)[0]-INTEGER(Rfrom_r)[0]+1));
  SEXP Rl      = PROTECT(Rf_allocVector(INTSXP, 1));
  INTEGER(Rl)[0] = 2147483647;
  SEXP Rbt     = PROTECT(Rf_allocVector(LGLSXP, 1));
  LOGICAL(Rbt)[0] = FALSE;

  const char *names[] = {"beta", "beta.g", "l", "from_recipient", "to_recipient", "pars.sha256", "beta.theta", ""};
  SEXP bck = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(bck, 0, Rbeta);
  SET_VECTOR_ELT(bck, 1, Rbeta_g);
  SET_VECTOR_ELT(bck, 2, Rl);
  SET_VECTOR_ELT(bck, 3, Rfrom_r);
  SET_VECTOR_ELT(bck, 4, Rto_r);
  SET_VECTOR_ELT(bck, 5, Rpars_sha);
  SET_VECTOR_ELT(bck, 6, Rbt);

  SEXP class = PROTECT(Rf_mkString("kalisBackwardTable"));
  Rf_setAttrib(bck, R_ClassSymbol, class);

  UNPROTECT(6);
  return(bck);
}
