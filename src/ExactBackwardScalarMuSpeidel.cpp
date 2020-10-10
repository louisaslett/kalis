#include "Stencil.h"

#define EXACTBACKWARDNOEXP ExactBackward_speidel_scmu
#define KALIS_MU MU_SCALAR
#define KALIS_PI PI_MATRIX
#define KALIS_SPEIDEL

#include "ExactBackwardStencil.cpp"

#undef EXACTBACKWARDNOEXP
#undef KALIS_MU
#undef KALIS_PI
#undef KALIS_SPEIDEL
