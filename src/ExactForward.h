#ifndef EXACTFORWARD_H
#define EXACTFORWARD_H

#include <Rcpp.h>
using namespace Rcpp;

void ExactForwardNaiveC_cpp(NumericMatrix alpha,
                            NumericVector alpha_f,
                            NumericVector alpha_f2,
                            const int alpha_from_rec,
                            const int alpha_t,
                            const int t,
                            const int from_rec,
                            const int to_rec,
                            const int L,
                            const int N,
                            NumericMatrix Pi,
                            NumericVector mu,
                            NumericVector rho);

void ParExactForwardNaiveC_cpp(NumericMatrix alpha,
                               NumericVector alpha_f,
                               NumericVector alpha_f2,
                               const int alpha_from_rec,
                               const int alpha_t,
                               const int t,
                               const int from_rec,
                               const int to_rec,
                               const int L,
                               const int N,
                               NumericMatrix Pi,
                               NumericVector mu,
                               NumericVector rho,
                               const int nthreads);

void ExactForwardNoExpAVX3_cpp(NumericMatrix alpha,
                               NumericVector alpha_f,
                               NumericVector alpha_f2,
                               const int alpha_from_rec,
                               const int alpha_t,
                               const int t,
                               const int from_rec,
                               const int to_rec,
                               const int L,
                               const int N,
                               NumericMatrix Pi,
                               NumericVector mu,
                               NumericVector rho);

void ParExactForwardNoExpAVX3_cpp(NumericMatrix alpha,
                                  NumericVector alpha_f,
                                  NumericVector alpha_f2,
                                  const int alpha_from_rec,
                                  const int alpha_t,
                                  const int t,
                                  const int from_rec,
                                  const int to_rec,
                                  const int L,
                                  const int N,
                                  NumericMatrix Pi,
                                  NumericVector mu,
                                  NumericVector rho,
                                  const int nthreads);

#endif
