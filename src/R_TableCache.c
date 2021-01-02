#include "R_TableCache.h"

#include <string.h> // memcpy

#include "R_Kalis.h"
#include "Cache.h"

SEXP CopyFBTable(SEXP Rto, SEXP Rfrom) {
  // Quite minimal checks here, mostly to avoid C level errors, most checks done pre-call in R
  if(LENGTH(Rto) != LENGTH(Rfrom)) {
    REprintf("Error: cannot copy between different kalis object types (eg forward <-> backward).\n");
    KALIS_RETURN;
  }

  // FIRST: everything common to both forward and backward
  // Load values and do sanity check these objects can be copied between
  double *tbl_to, *tbl_from;
  KALIS_GET_TABLE(tbl_to, Rto);
  KALIS_GET_TABLE(tbl_from, Rfrom);
  if(Rf_nrows(VECTOR_ELT(Rto, KALIS_IDX_TABLE)) != Rf_nrows(VECTOR_ELT(Rfrom, KALIS_IDX_TABLE)) ||
     Rf_ncols(VECTOR_ELT(Rto, KALIS_IDX_TABLE)) != Rf_ncols(VECTOR_ELT(Rfrom, KALIS_IDX_TABLE))) {
    REprintf("Error: cannot copy between forward objects with differing sizes.\n");
    KALIS_RETURN;
  }
  if(tbl_to == tbl_from) {
    REprintf("Error: underlying forward objects point to same memory.\n");
    KALIS_RETURN;
  }

  double *tbl_aux_to, *tbl_aux_from;
  KALIS_GET_TABLE_AUX(tbl_aux_to, Rto);
  KALIS_GET_TABLE_AUX(tbl_aux_from, Rfrom);
  // Size of above implicitly match if tables do
  if(tbl_aux_to == tbl_aux_from) {
    REprintf("Error: underlying forward objects point to same memory.\n");
    KALIS_RETURN;
  }

  int *l_to, *l_from;
  KALIS_GET_TABLE_POS(l_to, Rto);
  KALIS_GET_TABLE_POS(l_from, Rfrom);

  int *from_rec_to, *from_rec_from;
  KALIS_GET_TABLE_FROM(from_rec_to, Rto);
  KALIS_GET_TABLE_FROM(from_rec_from, Rfrom);

  int *to_rec_to, *to_rec_from;
  KALIS_GET_TABLE_TO(to_rec_to, Rto);
  KALIS_GET_TABLE_TO(to_rec_from, Rfrom);
  // Both from and to recipient values will be acceptably spaced for available memory if tables are correctly sized per above

  const size_t tbl_size = XLENGTH(VECTOR_ELT(Rto, KALIS_IDX_TABLE)) * sizeof(double);
  const size_t tbl_aux_size = XLENGTH(VECTOR_ELT(Rto, KALIS_IDX_TABLE_AUX)) * sizeof(double);

  // Perform copy
  memcpy(tbl_to, tbl_from, tbl_size);
  memcpy(tbl_aux_to, tbl_aux_from, tbl_aux_size);
  *l_to = *l_from;
  *from_rec_to = *from_rec_from;
  *to_rec_to = *to_rec_from;

  // SECOND: see if this is a backward and do the extra variable if so
  if(LENGTH(Rto) == 7) {
    int *betatheta_to, *betatheta_from;
    KALIS_GET_TABLE_BETATHETA(betatheta_to, Rto);
    KALIS_GET_TABLE_BETATHETA(betatheta_from, Rfrom);
    *betatheta_to = *betatheta_from;
  }

  KALIS_RETURN;
}
