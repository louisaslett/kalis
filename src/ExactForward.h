#ifndef EXACTFORWARD_H
#define EXACTFORWARD_H

#include <stddef.h>

void ExactForward_scalarMu_scalarPi(double *const restrict alpha,
                                    double *const restrict alpha_f,
                                    const size_t alpha_from_rec,
                                    const size_t alpha_t,
                                    const size_t t,
                                    const size_t from_rec,
                                    const size_t to_rec,
                                    const size_t L,
                                    const size_t N,
                                    const double Pi,
                                    const double mu,
                                    const double *const restrict rho,
                                    const int *const restrict nthreads,
                                    const int nthreads_len);

void ExactForward_scalarPi(double *const restrict alpha,
                           double *const restrict alpha_f,
                           const size_t alpha_from_rec,
                           const size_t alpha_t,
                           const size_t t,
                           const size_t from_rec,
                           const size_t to_rec,
                           const size_t L,
                           const size_t N,
                           const double Pi,
                           const double *const restrict mu,
                           const double *const restrict rho,
                           const int *const restrict nthreads,
                           const int nthreads_len);

void ExactForward_scalarMu(double *const restrict alpha,
                           double *const restrict alpha_f,
                           const size_t alpha_from_rec,
                           const size_t alpha_t,
                           const size_t t,
                           const size_t from_rec,
                           const size_t to_rec,
                           const size_t L,
                           const size_t N,
                           const double *const restrict Pi,
                           const double mu,
                           const double *const restrict rho,
                           const int *const restrict nthreads,
                           const int nthreads_len);

void ExactForward(double *const restrict alpha,
                  double *const restrict alpha_f,
                  const size_t alpha_from_rec,
                  const size_t alpha_t,
                  const size_t t,
                  const size_t from_rec,
                  const size_t to_rec,
                  const size_t L,
                  const size_t N,
                  const double *const restrict Pi,
                  const double *const restrict mu,
                  const double *const restrict rho,
                  const int *const restrict nthreads,
                  const int nthreads_len);

void ExactForward_speidel_scalarMu_scalarPi(double *const restrict alpha,
                                            double *const restrict alpha_f,
                                            const size_t alpha_from_rec,
                                            const size_t alpha_t,
                                            const size_t t,
                                            const size_t from_rec,
                                            const size_t to_rec,
                                            const size_t L,
                                            const size_t N,
                                            const double Pi,
                                            const double mu,
                                            const double *const restrict rho,
                                            const int *const restrict nthreads,
                                            const int nthreads_len);

void ExactForward_speidel_scalarPi(double *const restrict alpha,
                                   double *const restrict alpha_f,
                                   const size_t alpha_from_rec,
                                   const size_t alpha_t,
                                   const size_t t,
                                   const size_t from_rec,
                                   const size_t to_rec,
                                   const size_t L,
                                   const size_t N,
                                   const double Pi,
                                   const double *const restrict mu,
                                   const double *const restrict rho,
                                   const int *const restrict nthreads,
                                   const int nthreads_len);

void ExactForward_speidel_scalarMu(double *const restrict alpha,
                                   double *const restrict alpha_f,
                                   const size_t alpha_from_rec,
                                   const size_t alpha_t,
                                   const size_t t,
                                   const size_t from_rec,
                                   const size_t to_rec,
                                   const size_t L,
                                   const size_t N,
                                   const double *const restrict Pi,
                                   const double mu,
                                   const double *const restrict rho,
                                   const int *const restrict nthreads,
                                   const int nthreads_len);

void ExactForward_speidel(double *const restrict alpha,
                          double *const restrict alpha_f,
                          const size_t alpha_from_rec,
                          const size_t alpha_t,
                          const size_t t,
                          const size_t from_rec,
                          const size_t to_rec,
                          const size_t L,
                          const size_t N,
                          const double *const restrict Pi,
                          const double *const restrict mu,
                          const double *const restrict rho,
                          const int *const restrict nthreads,
                          const int nthreads_len);

void ExactForward_1step_scalarMu_scalarPi(double *const restrict alpha,
                                          double *const restrict alpha_f,
                                          const size_t alpha_from_rec,
                                          const size_t alpha_t,
                                          const size_t t,
                                          const size_t from_rec,
                                          const size_t to_rec,
                                          const size_t L,
                                          const size_t N,
                                          const double Pi,
                                          const double mu,
                                          const double *const restrict rho,
                                          const int *const restrict nthreads,
                                          const int nthreads_len);

void ExactForward_1step_speidel_scalarMu_scalarPi(double *const restrict alpha,
                                                  double *const restrict alpha_f,
                                                  const size_t alpha_from_rec,
                                                  const size_t alpha_t,
                                                  const size_t t,
                                                  const size_t from_rec,
                                                  const size_t to_rec,
                                                  const size_t L,
                                                  const size_t N,
                                                  const double Pi,
                                                  const double mu,
                                                  const double *const restrict rho,
                                                  const int *const restrict nthreads,
                                                  const int nthreads_len);

#endif
