#include <Rcpp.h>
using namespace Rcpp;

#include <algorithm>
#include <thread>
#include <cmath>
#include <stdlib.h>


void CalcTraces_A(const double* __restrict__ M,
                  const double* __restrict__ tX,
                  const double* __restrict__ tQ,
                  const double* __restrict__ Z,
                  const double* __restrict__ J,
                  double* __restrict__ res,
                  size_t r,
                  size_t from_off,
                  size_t from,
                  size_t N,
                  size_t p) {
  double temp;
  for(size_t j = from; j < from+N; j++) {
    for(size_t i = 0; i < r; i++) {
      temp = M[i + j*r];
      for(size_t l = 0; l < p; l++){
        // temp += temp;
        temp += tX[l + i*p] * Z[l + (from_off+j)*p] - tQ[l + i*p] * J[l + (from_off+j)*p];
      }
      res[0] += temp * temp;
    }
  }
}


// [[Rcpp::export]]
NumericVector CalcTraces(NumericMatrix M,  // r x c
                         NumericMatrix tX, // p x r
                         NumericMatrix tQ, // p x r
                         NumericMatrix Z,  // p x r (will only use a subset of these rows)
                         NumericMatrix J,  // p x r (will only use a subset of these rows)
                         size_t from_recipient,
                         size_t nthreads) {

  size_t p = (size_t) tX.nrow();
  size_t r = (size_t) M.nrow();
  size_t c = (size_t) M.ncol();

  NumericVector res(1);

  res[0] = 0.0;

  size_t from_off = (from_recipient-1);

  if(nthreads > 1) {
    double *res_perth = (double*) calloc((nthreads+1), sizeof(double));
    if (res_perth == NULL) {
      printf("Failed allocating res_perth!\n");
      exit(1);
    }

    size_t num_perth = c/nthreads;
    size_t rag_end   = c%nthreads;

    std::vector<std::thread> threads;
    for(size_t i=0; i<nthreads; ++i) {
      threads.push_back(std::thread(
          CalcTraces_A,
          &(M[0]), &(tX[0]), &(tQ[0]), &(Z[0]), &(J[0]), res_perth + i, r, from_off, i*num_perth, num_perth, p));
    }
    // Tidy ragged end
    if(rag_end != 0) {
      CalcTraces_A(&(M[0]), &(tX[0]), &(tQ[0]), &(Z[0]), &(J[0]), res_perth + nthreads, r, from_off, nthreads*num_perth, rag_end, p);
    }
    for(auto& th : threads) {
      th.join();
    }

    for(int j = 0; j < nthreads+1; j++) {
      res[0] += res_perth[j];
    }

    free(res_perth);
  } else {
    CalcTraces_A(&(M[0]), &(tX[0]), &(tQ[0]), &(Z[0]), &(J[0]), &(res[0]), r, from_off, 0, c, p);
  }

  return res ;
}