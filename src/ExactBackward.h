#ifndef EXACTBACKWARD_H
#define EXACTBACKWARD_H

#include <Rcpp.h>
using namespace Rcpp;

void ExactBackwardNaiveC_cpp(NumericMatrix beta,
                             NumericVector beta_g,
                             NumericVector beta_g2,
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

void ParExactBackwardNaiveC_cpp(NumericMatrix beta,
                                NumericVector beta_g,
                                NumericVector beta_g2,
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

void ExactBackwardNoExpAVX3_cpp(NumericMatrix beta,
                                NumericVector beta_g,
                                NumericVector beta_g2,
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
                                   NumericVector beta_g2,
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

#endif