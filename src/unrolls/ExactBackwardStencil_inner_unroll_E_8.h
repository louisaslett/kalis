/* File generated by R/unroll-backward-inner-loop.R */

double *betaNow0     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (0*KALIS_DOUBLEVEC_SIZE);
double *betaNow1     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (1*KALIS_DOUBLEVEC_SIZE);
double *betaNow2     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (2*KALIS_DOUBLEVEC_SIZE);
double *betaNow3     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (3*KALIS_DOUBLEVEC_SIZE);
double *betaNow4     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (4*KALIS_DOUBLEVEC_SIZE);
double *betaNow5     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (5*KALIS_DOUBLEVEC_SIZE);
double *betaNow6     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (6*KALIS_DOUBLEVEC_SIZE);
double *betaNow7     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (7*KALIS_DOUBLEVEC_SIZE);

KALIS_DOUBLE _beta0  = KALIS_LOADU_DOUBLE(betaNow0);
KALIS_DOUBLE _beta1  = KALIS_LOADU_DOUBLE(betaNow1);
KALIS_DOUBLE _beta2  = KALIS_LOADU_DOUBLE(betaNow2);
KALIS_DOUBLE _beta3  = KALIS_LOADU_DOUBLE(betaNow3);
KALIS_DOUBLE _beta4  = KALIS_LOADU_DOUBLE(betaNow4);
KALIS_DOUBLE _beta5  = KALIS_LOADU_DOUBLE(betaNow5);
KALIS_DOUBLE _beta6  = KALIS_LOADU_DOUBLE(betaNow6);
KALIS_DOUBLE _beta7  = KALIS_LOADU_DOUBLE(betaNow7);

KALIS_DOUBLE _theta0 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(0*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(0*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta1 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(1*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(1*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta2 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(2*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(2*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta3 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(3*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(3*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta4 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(4*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(4*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta5 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(5*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(5*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta6 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(6*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(6*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta7 = KALIS_SPREADBITSTO_DOUBLE((HA[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(7*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(7*KALIS_DOUBLEVEC_SIZE))%32));

#if KALIS_MU == MU_SCALAR
_theta0              = KALIS_FMA_DOUBLE(_theta0, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta1              = KALIS_FMA_DOUBLE(_theta1, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta2              = KALIS_FMA_DOUBLE(_theta2, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta3              = KALIS_FMA_DOUBLE(_theta3, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta4              = KALIS_FMA_DOUBLE(_theta4, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta5              = KALIS_FMA_DOUBLE(_theta5, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta6              = KALIS_FMA_DOUBLE(_theta6, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta7              = KALIS_FMA_DOUBLE(_theta7, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
#elif KALIS_MU == MU_VECTOR
_theta0              = KALIS_FMA_DOUBLE(_theta0, _muTmp1a, _muTmp2a); // theta = H * (2*mu - 1) - mu + 1
_theta1              = KALIS_FMA_DOUBLE(_theta1, _muTmp1a, _muTmp2a); // theta = H * (2*mu - 1) - mu + 1
_theta2              = KALIS_FMA_DOUBLE(_theta2, _muTmp1a, _muTmp2a); // theta = H * (2*mu - 1) - mu + 1
_theta3              = KALIS_FMA_DOUBLE(_theta3, _muTmp1a, _muTmp2a); // theta = H * (2*mu - 1) - mu + 1
_theta4              = KALIS_FMA_DOUBLE(_theta4, _muTmp1a, _muTmp2a); // theta = H * (2*mu - 1) - mu + 1
_theta5              = KALIS_FMA_DOUBLE(_theta5, _muTmp1a, _muTmp2a); // theta = H * (2*mu - 1) - mu + 1
_theta6              = KALIS_FMA_DOUBLE(_theta6, _muTmp1a, _muTmp2a); // theta = H * (2*mu - 1) - mu + 1
_theta7              = KALIS_FMA_DOUBLE(_theta7, _muTmp1a, _muTmp2a); // theta = H * (2*mu - 1) - mu + 1
#endif

_beta0               = KALIS_MUL_DOUBLE(_beta0, _theta0); // (theta*beta)
_beta1               = KALIS_MUL_DOUBLE(_beta1, _theta1); // (theta*beta)
_beta2               = KALIS_MUL_DOUBLE(_beta2, _theta2); // (theta*beta)
_beta3               = KALIS_MUL_DOUBLE(_beta3, _theta3); // (theta*beta)
_beta4               = KALIS_MUL_DOUBLE(_beta4, _theta4); // (theta*beta)
_beta5               = KALIS_MUL_DOUBLE(_beta5, _theta5); // (theta*beta)
_beta6               = KALIS_MUL_DOUBLE(_beta6, _theta6); // (theta*beta)
_beta7               = KALIS_MUL_DOUBLE(_beta7, _theta7); // (theta*beta)

KALIS_DOUBLE _pi0    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (0*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi1    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (1*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi2    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (2*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi3    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (3*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi4    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (4*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi5    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (5*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi6    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (6*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi7    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (7*KALIS_DOUBLEVEC_SIZE));

_g                   = KALIS_FMA_DOUBLE(_pi0, _beta0, _g); // g += Pi * (theta*beta)
_g                   = KALIS_FMA_DOUBLE(_pi1, _beta1, _g); // g += Pi * (theta*beta)
_g                   = KALIS_FMA_DOUBLE(_pi2, _beta2, _g); // g += Pi * (theta*beta)
_g                   = KALIS_FMA_DOUBLE(_pi3, _beta3, _g); // g += Pi * (theta*beta)
_g                   = KALIS_FMA_DOUBLE(_pi4, _beta4, _g); // g += Pi * (theta*beta)
_g                   = KALIS_FMA_DOUBLE(_pi5, _beta5, _g); // g += Pi * (theta*beta)
_g                   = KALIS_FMA_DOUBLE(_pi6, _beta6, _g); // g += Pi * (theta*beta)
_g                   = KALIS_FMA_DOUBLE(_pi7, _beta7, _g); // g += Pi * (theta*beta)

KALIS_STOREU_DOUBLE(betaNow0, _beta0);
KALIS_STOREU_DOUBLE(betaNow1, _beta1);
KALIS_STOREU_DOUBLE(betaNow2, _beta2);
KALIS_STOREU_DOUBLE(betaNow3, _beta3);
KALIS_STOREU_DOUBLE(betaNow4, _beta4);
KALIS_STOREU_DOUBLE(betaNow5, _beta5);
KALIS_STOREU_DOUBLE(betaNow6, _beta6);
KALIS_STOREU_DOUBLE(betaNow7, _beta7);
