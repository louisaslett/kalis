#include <R_ext/Rdynload.h>

#include "R_Backward.h"
#include "R_Cache.h"
#include "R_CalcTraces.h"
#include "R_ComputeStatus.h"
#include "R_Forward.h"
#include "R_MatAndMul.h"
#include "R_MatAndMulBtwVar.h"
#include "R_OptCkpt.h"
#include "R_TableCache.h"
#include "R_TableMaker.h"
#include "R_CladeMat.h"


#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef callMethods[] = {
  CALLDEF(ResetBackwardTable, 1),
  CALLDEF(Backward, 8),
  CALLDEF(ClearHaplotypeCache2, 0),
  CALLDEF(CacheHaplotypes_matrix_2, 4),
  CALLDEF(CacheHaplotypes_hdf5_2, 4),
  CALLDEF(CacheHaplotypes_hapgz_ncols, 1),
  CALLDEF(CacheHaplotypes_hapgz_nlines, 1),
  CALLDEF(CacheHaplotypes_hapgz_2, 5),
  CALLDEF(QueryCache2_ind, 1),
  CALLDEF(QueryCache2_loc, 1),
  CALLDEF(CalcTraces, 6),
  CALLDEF(ComputeStatus, 0),
  CALLDEF(VectorBitWidth, 0),
  CALLDEF(ResetForwardTable, 1),
  CALLDEF(Forward, 8),
  CALLDEF(MatAndMul, 9),
  CALLDEF(MatAndMulBtwVar, 11),
  CALLDEF(OptCkpt, 3),
  CALLDEF(CopyFBTable, 2),
  CALLDEF(MakeForwardTable, 3),
  CALLDEF(MakeBackwardTable, 3),
  CALLDEF(CladeMat, 7),
  CALLDEF(UpdateRealInPlace,3),
  { NULL, NULL, 0 }
};

void R_init_kalis(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}
