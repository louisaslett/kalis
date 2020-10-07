#include "Stencil.h"

#define EXACTBACKWARDNOEXP ExactBackward_speidel
#define KALIS_MU MU_VECTOR
#define KALIS_PI PI_MATRIX
#define KALIS_SPEIDEL

#include "ExactBackwardStencil.cpp"

#undef EXACTBACKWARDNOEXP
#undef KALIS_MU
#undef KALIS_PI
#undef KALIS_SPEIDEL



#define EXACTBACKWARDNOEXP ExactBackward
#define KALIS_MU MU_VECTOR
#define KALIS_PI PI_MATRIX

#include "ExactBackwardStencil.cpp"

#undef EXACTBACKWARDNOEXP
#undef KALIS_MU
#undef KALIS_PI
