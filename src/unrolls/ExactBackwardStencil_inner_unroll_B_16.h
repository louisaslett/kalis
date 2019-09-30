/* File generated by R/unroll-backward-inner-loop.R */

double *betaNow0     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (0*KALIS_DOUBLEVEC_SIZE);
double *betaNow1     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (1*KALIS_DOUBLEVEC_SIZE);
double *betaNow2     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (2*KALIS_DOUBLEVEC_SIZE);
double *betaNow3     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (3*KALIS_DOUBLEVEC_SIZE);
double *betaNow4     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (4*KALIS_DOUBLEVEC_SIZE);
double *betaNow5     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (5*KALIS_DOUBLEVEC_SIZE);
double *betaNow6     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (6*KALIS_DOUBLEVEC_SIZE);
double *betaNow7     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (7*KALIS_DOUBLEVEC_SIZE);
double *betaNow8     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (8*KALIS_DOUBLEVEC_SIZE);
double *betaNow9     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (9*KALIS_DOUBLEVEC_SIZE);
double *betaNow10     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (10*KALIS_DOUBLEVEC_SIZE);
double *betaNow11     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (11*KALIS_DOUBLEVEC_SIZE);
double *betaNow12     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (12*KALIS_DOUBLEVEC_SIZE);
double *betaNow13     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (13*KALIS_DOUBLEVEC_SIZE);
double *betaNow14     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (14*KALIS_DOUBLEVEC_SIZE);
double *betaNow15     = betaRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (15*KALIS_DOUBLEVEC_SIZE);

KALIS_DOUBLE _beta0  = KALIS_LOADU_DOUBLE(betaNow0);
KALIS_DOUBLE _beta1  = KALIS_LOADU_DOUBLE(betaNow1);
KALIS_DOUBLE _beta2  = KALIS_LOADU_DOUBLE(betaNow2);
KALIS_DOUBLE _beta3  = KALIS_LOADU_DOUBLE(betaNow3);
KALIS_DOUBLE _beta4  = KALIS_LOADU_DOUBLE(betaNow4);
KALIS_DOUBLE _beta5  = KALIS_LOADU_DOUBLE(betaNow5);
KALIS_DOUBLE _beta6  = KALIS_LOADU_DOUBLE(betaNow6);
KALIS_DOUBLE _beta7  = KALIS_LOADU_DOUBLE(betaNow7);
KALIS_DOUBLE _beta8  = KALIS_LOADU_DOUBLE(betaNow8);
KALIS_DOUBLE _beta9  = KALIS_LOADU_DOUBLE(betaNow9);
KALIS_DOUBLE _beta10  = KALIS_LOADU_DOUBLE(betaNow10);
KALIS_DOUBLE _beta11  = KALIS_LOADU_DOUBLE(betaNow11);
KALIS_DOUBLE _beta12  = KALIS_LOADU_DOUBLE(betaNow12);
KALIS_DOUBLE _beta13  = KALIS_LOADU_DOUBLE(betaNow13);
KALIS_DOUBLE _beta14  = KALIS_LOADU_DOUBLE(betaNow14);
KALIS_DOUBLE _beta15  = KALIS_LOADU_DOUBLE(betaNow15);

KALIS_DOUBLE _theta0 = KALIS_SPREADBITSTO_DOUBLE((HB[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(0*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(0*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta1 = KALIS_SPREADBITSTO_DOUBLE((HB[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(1*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(1*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta2 = KALIS_SPREADBITSTO_DOUBLE((HB[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(2*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(2*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta3 = KALIS_SPREADBITSTO_DOUBLE((HB[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(3*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(3*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta4 = KALIS_SPREADBITSTO_DOUBLE((HB[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(4*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(4*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta5 = KALIS_SPREADBITSTO_DOUBLE((HB[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(5*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(5*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta6 = KALIS_SPREADBITSTO_DOUBLE((HB[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(6*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(6*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta7 = KALIS_SPREADBITSTO_DOUBLE((HB[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(7*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(7*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta8 = KALIS_SPREADBITSTO_DOUBLE((HB[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(8*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(8*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta9 = KALIS_SPREADBITSTO_DOUBLE((HB[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(9*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(9*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta10 = KALIS_SPREADBITSTO_DOUBLE((HB[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(10*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(10*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta11 = KALIS_SPREADBITSTO_DOUBLE((HB[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(11*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(11*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta12 = KALIS_SPREADBITSTO_DOUBLE((HB[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(12*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(12*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta13 = KALIS_SPREADBITSTO_DOUBLE((HB[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(13*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(13*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta14 = KALIS_SPREADBITSTO_DOUBLE((HB[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(14*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(14*KALIS_DOUBLEVEC_SIZE))%32));
KALIS_DOUBLE _theta15 = KALIS_SPREADBITSTO_DOUBLE((HB[(donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(15*KALIS_DOUBLEVEC_SIZE))/32]) >> ((donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL)+(15*KALIS_DOUBLEVEC_SIZE))%32));

KALIS_DOUBLE _pi0    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (0*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi1    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (1*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi2    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (2*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi3    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (3*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi4    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (4*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi5    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (5*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi6    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (6*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi7    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (7*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi8    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (8*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi9    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (9*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi10    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (10*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi11    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (11*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi12    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (12*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi13    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (13*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi14    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (14*KALIS_DOUBLEVEC_SIZE));
KALIS_DOUBLE _pi15    = KALIS_LOADU_DOUBLE(PiRow + donoroff*(32*KALIS_INTVEC_SIZE) + donor*(KALIS_DOUBLEVEC_SIZE*KALIS_UNROLL) + (15*KALIS_DOUBLEVEC_SIZE));

_beta0               = KALIS_FMA_DOUBLE(_beta0, _omRhoDivG, _rho); // (rho + [theta*beta] * [(1-rho)/g])
_beta1               = KALIS_FMA_DOUBLE(_beta1, _omRhoDivG, _rho); // (rho + [theta*beta] * [(1-rho)/g])
_beta2               = KALIS_FMA_DOUBLE(_beta2, _omRhoDivG, _rho); // (rho + [theta*beta] * [(1-rho)/g])
_beta3               = KALIS_FMA_DOUBLE(_beta3, _omRhoDivG, _rho); // (rho + [theta*beta] * [(1-rho)/g])
_beta4               = KALIS_FMA_DOUBLE(_beta4, _omRhoDivG, _rho); // (rho + [theta*beta] * [(1-rho)/g])
_beta5               = KALIS_FMA_DOUBLE(_beta5, _omRhoDivG, _rho); // (rho + [theta*beta] * [(1-rho)/g])
_beta6               = KALIS_FMA_DOUBLE(_beta6, _omRhoDivG, _rho); // (rho + [theta*beta] * [(1-rho)/g])
_beta7               = KALIS_FMA_DOUBLE(_beta7, _omRhoDivG, _rho); // (rho + [theta*beta] * [(1-rho)/g])
_beta8               = KALIS_FMA_DOUBLE(_beta8, _omRhoDivG, _rho); // (rho + [theta*beta] * [(1-rho)/g])
_beta9               = KALIS_FMA_DOUBLE(_beta9, _omRhoDivG, _rho); // (rho + [theta*beta] * [(1-rho)/g])
_beta10               = KALIS_FMA_DOUBLE(_beta10, _omRhoDivG, _rho); // (rho + [theta*beta] * [(1-rho)/g])
_beta11               = KALIS_FMA_DOUBLE(_beta11, _omRhoDivG, _rho); // (rho + [theta*beta] * [(1-rho)/g])
_beta12               = KALIS_FMA_DOUBLE(_beta12, _omRhoDivG, _rho); // (rho + [theta*beta] * [(1-rho)/g])
_beta13               = KALIS_FMA_DOUBLE(_beta13, _omRhoDivG, _rho); // (rho + [theta*beta] * [(1-rho)/g])
_beta14               = KALIS_FMA_DOUBLE(_beta14, _omRhoDivG, _rho); // (rho + [theta*beta] * [(1-rho)/g])
_beta15               = KALIS_FMA_DOUBLE(_beta15, _omRhoDivG, _rho); // (rho + [theta*beta] * [(1-rho)/g])

#if KALIS_MU == MU_SCALAR
_theta0              = KALIS_FMA_DOUBLE(_theta0, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta1              = KALIS_FMA_DOUBLE(_theta1, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta2              = KALIS_FMA_DOUBLE(_theta2, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta3              = KALIS_FMA_DOUBLE(_theta3, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta4              = KALIS_FMA_DOUBLE(_theta4, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta5              = KALIS_FMA_DOUBLE(_theta5, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta6              = KALIS_FMA_DOUBLE(_theta6, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta7              = KALIS_FMA_DOUBLE(_theta7, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta8              = KALIS_FMA_DOUBLE(_theta8, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta9              = KALIS_FMA_DOUBLE(_theta9, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta10              = KALIS_FMA_DOUBLE(_theta10, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta11              = KALIS_FMA_DOUBLE(_theta11, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta12              = KALIS_FMA_DOUBLE(_theta12, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta13              = KALIS_FMA_DOUBLE(_theta13, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta14              = KALIS_FMA_DOUBLE(_theta14, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
_theta15              = KALIS_FMA_DOUBLE(_theta15, _muTmp1, _muTmp2); // theta = H * (2*mu - 1) - mu + 1
#elif KALIS_MU == MU_VECTOR
_theta0              = KALIS_FMA_DOUBLE(_theta0, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
_theta1              = KALIS_FMA_DOUBLE(_theta1, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
_theta2              = KALIS_FMA_DOUBLE(_theta2, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
_theta3              = KALIS_FMA_DOUBLE(_theta3, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
_theta4              = KALIS_FMA_DOUBLE(_theta4, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
_theta5              = KALIS_FMA_DOUBLE(_theta5, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
_theta6              = KALIS_FMA_DOUBLE(_theta6, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
_theta7              = KALIS_FMA_DOUBLE(_theta7, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
_theta8              = KALIS_FMA_DOUBLE(_theta8, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
_theta9              = KALIS_FMA_DOUBLE(_theta9, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
_theta10              = KALIS_FMA_DOUBLE(_theta10, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
_theta11              = KALIS_FMA_DOUBLE(_theta11, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
_theta12              = KALIS_FMA_DOUBLE(_theta12, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
_theta13              = KALIS_FMA_DOUBLE(_theta13, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
_theta14              = KALIS_FMA_DOUBLE(_theta14, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
_theta15              = KALIS_FMA_DOUBLE(_theta15, _muTmp1b, _muTmp2b); // theta = H * (2*mu - 1) - mu + 1
#endif

_beta0              = KALIS_MUL_DOUBLE(_beta0, _theta0); // (theta*beta)
_beta1              = KALIS_MUL_DOUBLE(_beta1, _theta1); // (theta*beta)
_beta2              = KALIS_MUL_DOUBLE(_beta2, _theta2); // (theta*beta)
_beta3              = KALIS_MUL_DOUBLE(_beta3, _theta3); // (theta*beta)
_beta4              = KALIS_MUL_DOUBLE(_beta4, _theta4); // (theta*beta)
_beta5              = KALIS_MUL_DOUBLE(_beta5, _theta5); // (theta*beta)
_beta6              = KALIS_MUL_DOUBLE(_beta6, _theta6); // (theta*beta)
_beta7              = KALIS_MUL_DOUBLE(_beta7, _theta7); // (theta*beta)
_beta8              = KALIS_MUL_DOUBLE(_beta8, _theta8); // (theta*beta)
_beta9              = KALIS_MUL_DOUBLE(_beta9, _theta9); // (theta*beta)
_beta10              = KALIS_MUL_DOUBLE(_beta10, _theta10); // (theta*beta)
_beta11              = KALIS_MUL_DOUBLE(_beta11, _theta11); // (theta*beta)
_beta12              = KALIS_MUL_DOUBLE(_beta12, _theta12); // (theta*beta)
_beta13              = KALIS_MUL_DOUBLE(_beta13, _theta13); // (theta*beta)
_beta14              = KALIS_MUL_DOUBLE(_beta14, _theta14); // (theta*beta)
_beta15              = KALIS_MUL_DOUBLE(_beta15, _theta15); // (theta*beta)

_g                   = KALIS_FMA_DOUBLE(_pi0, _beta0, _g); // g += Pi * (theta*beta)
_g                   = KALIS_FMA_DOUBLE(_pi1, _beta1, _g); // g += Pi * (theta*beta)
_g                   = KALIS_FMA_DOUBLE(_pi2, _beta2, _g); // g += Pi * (theta*beta)
_g                   = KALIS_FMA_DOUBLE(_pi3, _beta3, _g); // g += Pi * (theta*beta)
_g                   = KALIS_FMA_DOUBLE(_pi4, _beta4, _g); // g += Pi * (theta*beta)
_g                   = KALIS_FMA_DOUBLE(_pi5, _beta5, _g); // g += Pi * (theta*beta)
_g                   = KALIS_FMA_DOUBLE(_pi6, _beta6, _g); // g += Pi * (theta*beta)
_g                   = KALIS_FMA_DOUBLE(_pi7, _beta7, _g); // g += Pi * (theta*beta)
_g                   = KALIS_FMA_DOUBLE(_pi8, _beta8, _g); // g += Pi * (theta*beta)
_g                   = KALIS_FMA_DOUBLE(_pi9, _beta9, _g); // g += Pi * (theta*beta)
_g                   = KALIS_FMA_DOUBLE(_pi10, _beta10, _g); // g += Pi * (theta*beta)
_g                   = KALIS_FMA_DOUBLE(_pi11, _beta11, _g); // g += Pi * (theta*beta)
_g                   = KALIS_FMA_DOUBLE(_pi12, _beta12, _g); // g += Pi * (theta*beta)
_g                   = KALIS_FMA_DOUBLE(_pi13, _beta13, _g); // g += Pi * (theta*beta)
_g                   = KALIS_FMA_DOUBLE(_pi14, _beta14, _g); // g += Pi * (theta*beta)
_g                   = KALIS_FMA_DOUBLE(_pi15, _beta15, _g); // g += Pi * (theta*beta)

KALIS_STOREU_DOUBLE(betaNow0, _beta0);
KALIS_STOREU_DOUBLE(betaNow1, _beta1);
KALIS_STOREU_DOUBLE(betaNow2, _beta2);
KALIS_STOREU_DOUBLE(betaNow3, _beta3);
KALIS_STOREU_DOUBLE(betaNow4, _beta4);
KALIS_STOREU_DOUBLE(betaNow5, _beta5);
KALIS_STOREU_DOUBLE(betaNow6, _beta6);
KALIS_STOREU_DOUBLE(betaNow7, _beta7);
KALIS_STOREU_DOUBLE(betaNow8, _beta8);
KALIS_STOREU_DOUBLE(betaNow9, _beta9);
KALIS_STOREU_DOUBLE(betaNow10, _beta10);
KALIS_STOREU_DOUBLE(betaNow11, _beta11);
KALIS_STOREU_DOUBLE(betaNow12, _beta12);
KALIS_STOREU_DOUBLE(betaNow13, _beta13);
KALIS_STOREU_DOUBLE(betaNow14, _beta14);
KALIS_STOREU_DOUBLE(betaNow15, _beta15);