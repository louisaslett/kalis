#ifndef R_KALIS_H
#define R_KALIS_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include "Cache.h"



// Note that minimal error checking is performed to ensure onward code which
//   relies on manipulating only pointers does not segfault.
// Standard error checks are left to R's C interface on access to SEXPs.



// Define the list index position of elements of a forward/backward objects
#define KALIS_IDX_TABLE 0
#define KALIS_IDX_TABLE_AUX 1
#define KALIS_IDX_TABLE_POS 2
#define KALIS_IDX_TABLE_FROM 3
#define KALIS_IDX_TABLE_TO 4
#define KALIS_IDX_TABLE_BETATHETA 6



// Accessor macros to retrieve C versions of forward/backward objects
#define KALIS_GET_TABLE(DBLPTR, TBL)                                           \
if(num_inds == 0 || hap_size == 0) {                                           \
  REprintf("Error: Cache is empty.\n");                                        \
  KALIS_RETURN;                                                                \
}                                                                              \
if((size_t) Rf_nrows(VECTOR_ELT(TBL, KALIS_IDX_TABLE)) != num_inds) {          \
  REprintf("Error: Incorrect number of rows in table.\n");                     \
  KALIS_RETURN;                                                                \
}                                                                              \
/* nrows == to_recipient-from_recipient+1 && nrows <= L */                     \
if(Rf_ncols(VECTOR_ELT(TBL, KALIS_IDX_TABLE)) !=                               \
     Rf_asInteger(VECTOR_ELT(TBL, KALIS_IDX_TABLE_TO)) -                       \
       Rf_asInteger(VECTOR_ELT(TBL, KALIS_IDX_TABLE_FROM)) + 1 ||              \
   (size_t) Rf_ncols(VECTOR_ELT(TBL, KALIS_IDX_TABLE)) > num_inds) {           \
  REprintf("Error: Incorrect number of columns in table.\n");                  \
  KALIS_RETURN;                                                                \
}                                                                              \
DBLPTR = REAL(VECTOR_ELT(TBL, KALIS_IDX_TABLE));



#define KALIS_GET_TABLE_AUX(DBLPTR, TBL)                                       \
if(num_inds == 0 || hap_size == 0) {                                           \
  REprintf("Error: Cache is empty.\n");                                        \
  KALIS_RETURN;                                                                \
}                                                                              \
/* length == to_recipient-from_recipient+1 && length <= L */                   \
if(XLENGTH(VECTOR_ELT(TBL, KALIS_IDX_TABLE_AUX)) !=                            \
     Rf_asInteger(VECTOR_ELT(TBL, KALIS_IDX_TABLE_TO)) -                       \
       Rf_asInteger(VECTOR_ELT(TBL, KALIS_IDX_TABLE_FROM)) + 1 ||              \
   (size_t) Rf_ncols(VECTOR_ELT(TBL, KALIS_IDX_TABLE_AUX)) > num_inds) {       \
  REprintf("Error: Incorrect length of auxiliary table vector.\n");            \
  KALIS_RETURN;                                                                \
}                                                                              \
DBLPTR = REAL(VECTOR_ELT(TBL, KALIS_IDX_TABLE_AUX));



#define KALIS_GET_TABLE_POS(INTPTR, TBL)                                       \
if(XLENGTH(VECTOR_ELT(TBL, KALIS_IDX_TABLE_POS)) > 1) {                        \
  REprintf("Error: Table position is not scalar.\n");                          \
  KALIS_RETURN;                                                                \
}                                                                              \
INTPTR = INTEGER(VECTOR_ELT(TBL, KALIS_IDX_TABLE_POS));



#define KALIS_GET_TABLE_FROM(INTPTR, TBL)                                      \
if(XLENGTH(VECTOR_ELT(TBL, KALIS_IDX_TABLE_FROM)) > 1) {                       \
  REprintf("Error: Table from recipient is not scalar.\n");                    \
  KALIS_RETURN;                                                                \
}                                                                              \
INTPTR = INTEGER(VECTOR_ELT(TBL, KALIS_IDX_TABLE_FROM));                       \
if(*INTPTR < 1 || (size_t) *INTPTR > num_inds) {                               \
  REprintf("Error: Table from recipient invalid.\n");                          \
  KALIS_RETURN;                                                                \
}



#define KALIS_GET_TABLE_TO(INTPTR, TBL)                                        \
if(XLENGTH(VECTOR_ELT(TBL, KALIS_IDX_TABLE_TO)) > 1) {                         \
  REprintf("Error: Table to recipient is not scalar.\n");                      \
  KALIS_RETURN;                                                                \
}                                                                              \
INTPTR = INTEGER(VECTOR_ELT(TBL, KALIS_IDX_TABLE_TO));                         \
if(*INTPTR < 1 || (size_t) *INTPTR > num_inds) {                               \
  REprintf("Error: Table to recipient invalid.\n");                            \
  KALIS_RETURN;                                                                \
}



#define KALIS_GET_TABLE_BETATHETA(INTPTR, TBL)                                 \
if(XLENGTH(VECTOR_ELT(TBL, KALIS_IDX_TABLE_BETATHETA)) > 1) {                  \
  REprintf("Error: Table beta-theta status is not scalar.\n");                 \
  KALIS_RETURN;                                                                \
}                                                                              \
INTPTR = LOGICAL(VECTOR_ELT(TBL, KALIS_IDX_TABLE_BETATHETA));



//#define KALIS_RETURN return(Rf_allocVector(NILSXP, 1));
#define KALIS_RETURN return(R_NilValue);

#endif
