#ifndef EXACTBACKWARD_H
#define EXACTBACKWARD_H

#include <Rcpp.h>
using namespace Rcpp;

void ExactBackward_speidel_cpp(NumericMatrix beta,
                               NumericVector beta_g,
                               LogicalVector cur_beta_theta,
                               LogicalVector end_beta_theta,
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

void ExactBackward_cpp(NumericMatrix beta,
                       NumericVector beta_g,
                       LogicalVector cur_beta_theta,
                       LogicalVector end_beta_theta,
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

void ParExactBackward_speidel_cpp(NumericMatrix beta,
                                  NumericVector beta_g,
                                  LogicalVector cur_beta_theta,
                                  LogicalVector end_beta_theta,
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
                                  IntegerVector nthreads);

void ParExactBackward_cpp(NumericMatrix beta,
                          NumericVector beta_g,
                          LogicalVector cur_beta_theta,
                          LogicalVector end_beta_theta,
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
                          IntegerVector nthreads);

void ExactBackward_speidel_scmu_cpp(NumericMatrix beta,
                                    NumericVector beta_g,
                                    LogicalVector cur_beta_theta,
                                    LogicalVector end_beta_theta,
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

void ExactBackward_scmu_cpp(NumericMatrix beta,
                            NumericVector beta_g,
                            LogicalVector cur_beta_theta,
                            LogicalVector end_beta_theta,
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

void ParExactBackward_speidel_scmu_cpp(NumericMatrix beta,
                                       NumericVector beta_g,
                                       LogicalVector cur_beta_theta,
                                       LogicalVector end_beta_theta,
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
                                       IntegerVector nthreads);

void ParExactBackward_scmu_cpp(NumericMatrix beta,
                               NumericVector beta_g,
                               LogicalVector cur_beta_theta,
                               LogicalVector end_beta_theta,
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
                               IntegerVector nthreads);

void ExactBackward_speidel_scPi_cpp(NumericMatrix beta,
                                    NumericVector beta_g,
                                    LogicalVector cur_beta_theta,
                                    LogicalVector end_beta_theta,
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

void ExactBackward_scPi_cpp(NumericMatrix beta,
                            NumericVector beta_g,
                            LogicalVector cur_beta_theta,
                            LogicalVector end_beta_theta,
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

void ParExactBackward_speidel_scPi_cpp(NumericMatrix beta,
                                       NumericVector beta_g,
                                       LogicalVector cur_beta_theta,
                                       LogicalVector end_beta_theta,
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
                                       IntegerVector nthreads);

void ParExactBackward_scPi_cpp(NumericMatrix beta,
                               NumericVector beta_g,
                               LogicalVector cur_beta_theta,
                               LogicalVector end_beta_theta,
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
                               IntegerVector nthreads);

void ExactBackward_speidel_scmuPi_cpp(NumericMatrix beta,
                                      NumericVector beta_g,
                                      LogicalVector cur_beta_theta,
                                      LogicalVector end_beta_theta,
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

void ExactBackward_scmuPi_cpp(NumericMatrix beta,
                              NumericVector beta_g,
                              LogicalVector cur_beta_theta,
                              LogicalVector end_beta_theta,
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

void ParExactBackward_speidel_scmuPi_cpp(NumericMatrix beta,
                                         NumericVector beta_g,
                                         LogicalVector cur_beta_theta,
                                         LogicalVector end_beta_theta,
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
                                         IntegerVector nthreads);

void ParExactBackward_scmuPi_cpp(NumericMatrix beta,
                                 NumericVector beta_g,
                                 LogicalVector cur_beta_theta,
                                 LogicalVector end_beta_theta,
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
                                 IntegerVector nthreads);

#endif
