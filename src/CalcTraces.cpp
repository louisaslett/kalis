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
                  double* __restrict__ res2,
                  double* diag,
                  size_t r,
                  size_t from_off,
                  size_t from,
                  size_t N,
                  size_t p) {
  double temp;
  double temp2;
  for(size_t j = from; j < from+N; j++) {
    for(size_t i = 0; i < r; i++) {
      temp = M[i + j*r];
      for(size_t l = 0; l < p; l++){
        // temp += temp;
        temp += tX[l + i*p] * Z[l + (from_off+j)*p] - tQ[l + i*p] * J[l + (from_off+j)*p];
      }
      temp2 = temp * temp;
      res2[0] += temp2; // part of the HS norm
      if(i==j){
        res[0] += temp; // add to the trace (for the expected value of the quadratic form)
        diag[j] = temp; // add to the diagonal component of the varaince (for the variance of the quadratic form)
      }
    }
  }
}


// [[Rcpp::export]]
List CalcTraces(NumericMatrix M,  // r x c
                         NumericMatrix tX, // p x r
                         NumericMatrix tQ, // p x r
                         NumericMatrix Z,  // p x r (will only use a subset of these rows)
                         NumericMatrix J,  // p x r (will only use a subset of these rows)
                         size_t from_recipient,
                         size_t nthreads) {

  size_t p = (size_t) tX.nrow();
  size_t r = (size_t) M.nrow();
  size_t c = (size_t) M.ncol();

  NumericVector trace(1);
  NumericVector hsnorm(1);
  NumericVector diag(r);

  size_t from_off = (from_recipient-1);

  if(nthreads > 1) {
    double *res_perth = (double*) calloc((nthreads+1), sizeof(double));
    double *res_perth2 = (double*) calloc((nthreads+1), sizeof(double));

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
          &(M[0]), &(tX[0]), &(tQ[0]), &(Z[0]), &(J[0]), res_perth + i, res_perth2 + i, &(diag[0]), r, from_off, i*num_perth, num_perth, p));
    }
    // Tidy ragged end
    if(rag_end != 0) {
      CalcTraces_A(&(M[0]), &(tX[0]), &(tQ[0]), &(Z[0]), &(J[0]), res_perth + nthreads, res_perth2 + nthreads, &(diag[0]), r, from_off, nthreads*num_perth, rag_end, p);
    }
    for(auto& th : threads) {
      th.join();
    }

    for(int j = 0; j < nthreads+1; j++) {
      trace[0] += res_perth[j];
      hsnorm[0] += res_perth2[j];
    }

    free(res_perth);
    free(res_perth2);
  } else {
    CalcTraces_A(&(M[0]), &(tX[0]), &(tQ[0]), &(Z[0]), &(J[0]), &(trace[0]), &(hsnorm[0]), &(diag[0]), r, from_off, 0, c, p);
  }

  List L = List::create(Named("trace") = trace , Named("hsnorm2") = hsnorm,
                        Named("diag") = diag);
  return L;
}