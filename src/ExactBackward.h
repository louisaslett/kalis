#ifndef EXACTBACKWARD_H
#define EXACTBACKWARD_H

#include <Rcpp.h>
using namespace Rcpp;

void ExactBackwardNoExpAVX3_cpp(NumericMatrix beta,
                                NumericVector beta_g,
                                const int beta_from_rec,
                                const int beta_t,
                                const int t,
                                const int from_rec,
                                const int to_rec,
                                const int L,
                                const int N,
                                NumericMatrix Pi,
                                NumericVector mu,
                                NumericVector rho);

void ParExactBackwardNoExpAVX3_cpp(NumericMatrix beta,
                                   NumericVector beta_g,
                                   const int beta_from_rec,
                                   const int beta_t,
                                   const int t,
                                   const int from_rec,
                                   const int to_rec,
                                   const int L,
                                   const int N,
                                   NumericMatrix Pi,
                                   NumericVector mu,
                                   NumericVector rho,
                                   const int nthreads);

void ExactBackwardNoExpAVX3_scmu_cpp(NumericMatrix beta,
                                     NumericVector beta_g,
                                     const int beta_from_rec,
                                     const int beta_t,
                                     const int t,
                                     const int from_rec,
                                     const int to_rec,
                                     const int L,
                                     const int N,
                                     NumericMatrix Pi,
                                     const double mu,
                                     NumericVector rho);

void ParExactBackwardNoExpAVX3_scmu_cpp(NumericMatrix beta,
                                        NumericVector beta_g,
                                        const int beta_from_rec,
                                        const int beta_t,
                                        const int t,
                                        const int from_rec,
                                        const int to_rec,
                                        const int L,
                                        const int N,
                                        NumericMatrix Pi,
                                        const double mu,
                                        NumericVector rho,
                                        const int nthreads);

void ExactBackwardNoExpAVX3_scPi_cpp(NumericMatrix beta,
                                     NumericVector beta_g,
                                     const int beta_from_rec,
                                     const int beta_t,
                                     const int t,
                                     const int from_rec,
                                     const int to_rec,
                                     const int L,
                                     const int N,
                                     const double Pi,
                                     NumericVector mu,
                                     NumericVector rho);

void ParExactBackwardNoExpAVX3_scPi_cpp(NumericMatrix beta,
                                        NumericVector beta_g,
                                        const int beta_from_rec,
                                        const int beta_t,
                                        const int t,
                                        const int from_rec,
                                        const int to_rec,
                                        const int L,
                                        const int N,
                                        const double Pi,
                                        NumericVector mu,
                                        NumericVector rho,
                                        const int nthreads);

void ExactBackwardNoExpAVX3_scmuPi_cpp(NumericMatrix beta,
                                       NumericVector beta_g,
                                       const int beta_from_rec,
                                       const int beta_t,
                                       const int t,
                                       const int from_rec,
                                       const int to_rec,
                                       const int L,
                                       const int N,
                                       const double Pi,
                                       const double mu,
                                       NumericVector rho);

void ParExactBackwardNoExpAVX3_scmuPi_cpp(NumericMatrix beta,
                                          NumericVector beta_g,
                                          const int beta_from_rec,
                                          const int beta_t,
                                          const int t,
                                          const int from_rec,
                                          const int to_rec,
                                          const int L,
                                          const int N,
                                          const double Pi,
                                          const double mu,
                                          NumericVector rho,
                                          const int nthreads);

#endif
