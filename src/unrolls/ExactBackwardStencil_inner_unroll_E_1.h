/* File generated by R/unroll-backward-inner-loop.R */

double *betaNow0     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (0*KALIS_DOUBLEVEC_SIZE);

KALIS_DOUBLE _beta0  = KALIS_LOADU_DOUBLE(betaNow0);

KALIS_DOUBLE _theta0 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(0*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(0*KALIS_DOUBLEVEC_SIZE))%32));

#if KALIS_MU == MU_SCALAR
_theta0              = KALIS_FMA_DOUBLE(_theta0, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
#elif KALIS_MU == MU_VECTOR
_theta0              = KALIS_FMA_DOUBLE(_theta0, _muTmp1a, _muTmp2a); // theta = H * (2*mu - 1) - mu + 1
#endif

_beta0               = KALIS_MUL_DOUBLE(_beta0, _theta0); // (theta*beta)

KALIS_DOUBLE _pi0    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (0*KALIS_DOUBLEVEC_SIZE));

_g                   = KALIS_FMA_DOUBLE(_pi0, _beta0, _g); // g += Pi * (theta*beta)

KALIS_STOREU_DOUBLE(betaNow0, _beta0);