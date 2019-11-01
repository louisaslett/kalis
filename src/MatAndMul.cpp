#include <Rcpp.h>
using namespace Rcpp;

#include <algorithm>
#include <thread>
#include <cmath>
#include <stdlib.h>



void MatAndMul_C(const double* __restrict__ alpha,
                 const double* __restrict__ beta,
                 const double* __restrict__ x,
                 double* __restrict__ c_1,
                 double* __restrict__ res,
                 size_t j,
                 size_t r,
                 size_t from_off) {
  double z1;
  z1 = 0.0;

  for(size_t i = 0; i < r; i++) {
    z1 += c_1[i] = *(alpha++) * *(beta++);
  }

  if(z1 <= 0){
    // if the sum is zero, then we assume that all c_1 are 0, so we say that all of the distances are really far, so distance - min possible distance  = 0...so all of the distances go to 0.
    // this is consistent with the standardized distance version
    for(size_t i = 0; i < r; i++) {
      c_1[i] = 0.0;
    }
  } else {

    z1 = log(z1);

    for(size_t i = 0; i < r; i++) {
      if(i!=j) {
        if(c_1[i] == 0){
          c_1[i] = 708.396418532264078749 + z1;
        } else {
          c_1[i] = -log(c_1[i]) + z1; // this is the final value for M[i,j]
        }
        res[j+from_off] += c_1[i]*x[i];
        res[i]          += c_1[i]*x[j+from_off];
      } else {
        c_1[i] = 0.0;
      }
    }
  }
}

void MatAndMul_C_standardize(const double* __restrict__ alpha,
                             const double* __restrict__ beta,
                             const double* __restrict__ x,
                             double* __restrict__ c_1,
                             double* __restrict__ res,
                             size_t j,
                             size_t r,
                             size_t from_off) {
  double z1 = 0.0;
  double z2 = 0.0;
  double temp = 0.0;
  double rd = (double)(r - 1);

  for(size_t i = 0; i < r; i++) {

    temp = *(alpha++) * *(beta++);

    if(i!=j) {
      if(temp <= 0) {
        z1 += c_1[i] = 708.396418532264078749;
      } else {
        z1 += c_1[i] = -log(temp);
      }
      z2 += std::pow(c_1[i], 2.0);
    }
  }

  // note that if all of the c_1 are 0, then all of the distances are 708 distance - mean distance  = 0...so all of the standardized distances go to 0.
  // this is consistent with the unstandardized distance version

  z1 = z1 / rd; // Calculate Mean

  z2 = (z2 - std::pow(z1, 2.0) * rd) / (rd - 1.0);

  z2 = std::pow(z2, -0.5); // Calculate standard deviation

  z1 = -z1 * z2;

  for(size_t i = 0; i < r; i++) {

    if(i!=j) {
      c_1[i] = c_1[i] * z2 + z1; // this is the final value for M[i,j]
      res[j+from_off] += c_1[i]*x[i];
      res[i]          += c_1[i]*x[j+from_off];
    } else {
      c_1[i] = 0.0;
    }
  }

}


void MatAndMul_B(double* M,
                 const double* __restrict__ alpha,
                 const double* __restrict__ beta,
                 const double* __restrict__ x,
                 bool standardize,
                 double* __restrict__ res,
                 size_t r,
                 size_t from_off,
                 size_t from,
                 size_t N) {
  if( standardize == true ){
    for(size_t j = from; j < from+N; j++) {
      MatAndMul_C_standardize(alpha+j*r,
                              beta+j*r,
                              x,
                              M+j*r,
                              res,
                              j,
                              r,
                              from_off);
    }
  } else {
    for(size_t j = from; j < from+N; j++) {
      MatAndMul_C(alpha+j*r,
                  beta+j*r,
                  x,
                  M+j*r,
                  res,
                  j,
                  r,
                  from_off);
    }
  }
}

void MatAndMul_A(double* __restrict__ res,
                 double* __restrict__ M,
                 const double* __restrict__ alpha,
                 const double* __restrict__ beta,
                 const double* __restrict__ x,
                 bool standardize,
                 size_t from_off,
                 size_t nthreads,
                 size_t r,
                 size_t c) {

  from_off = (from_off-1); // got rid of /2 here

  if(nthreads > 1) {
    double *res_perth = (double*) calloc(r*(nthreads+1), sizeof(double));
    if (res_perth == NULL) {
      printf("Failed allocating res_perth!\n");
      exit(1);
    }

    size_t num_perth = c/nthreads;
    size_t rag_end   = c%nthreads;

    std::vector<std::thread> threads;
    for(size_t i=0; i<nthreads; ++i) {
      threads.push_back(std::thread(
          MatAndMul_B,
          M, alpha, beta, x, standardize, res_perth + i*r, r, from_off, i*num_perth, num_perth));
    }
    // Tidy ragged end
    if(rag_end != 0) {
      MatAndMul_B(M, alpha, beta, x, standardize, res_perth+nthreads*r, r, from_off, nthreads*num_perth, rag_end);
    }
    for(auto& th : threads) {
      th.join();
    }

    for(size_t i = 0; i < r; i++) {
      for(size_t j = 0; j < nthreads+1; j++) {
        res[i] += res_perth[i+j*r];
      }
      res[i] *= 0.5;
    }

    free(res_perth);
  } else {

    MatAndMul_B(M, alpha, beta, x, standardize, res, r, from_off, 0, c);

    for(size_t i = 0; i < r; i++) {
      res[i] *= 0.5;
    }

  }
}




// [[Rcpp::export]]
NumericVector MatAndMul(NumericMatrix M,
                        NumericMatrix alpha,
                        NumericMatrix beta,
                        NumericVector x,
                        LogicalVector standardize,
                        int from_recipient,
                        int nthreads) {

  size_t r = (size_t) alpha.nrow(); // we need to make p=r and c2 = c
  size_t c = (size_t) alpha.ncol();

  NumericVector res(r);

  if(alpha.nrow() != r || alpha.ncol() != c || beta.nrow() != r || beta.ncol() != c) {
    Rcout << "alpha, beta aren't right!\n";
    return(NumericVector(1));
  }
  if(x.length() != r) {
    Rcout << "x isn't right!\n";
    return(NumericVector(1));
  }
  if(M.nrow() != r || M.ncol() != c) {
    Rcout << "M isn't right!\n";
    return(NumericVector(1));
  }
  if(res.length() != r) {
    Rcout << "res isn't right!\n";
    return(NumericVector(1));
  }

  bool stdz = is_true(all(standardize));

  MatAndMul_A(&(res[0]),
              &(M[0]),
              &(alpha[0]),
              &(beta[0]),
              &(x[0]),
              stdz,
              from_recipient,
              nthreads,
              r,
              c);

  return(res);
}


void MatOnlyMul_A(const double* __restrict__ M,
                  const double* __restrict__ x,
                  double* __restrict__ res,
                  size_t nrow,
                  size_t from_off,
                  size_t from,
                  size_t N) {
  for(size_t j = from; j < from+N; j++) {
    for(size_t i = 0; i < nrow; i++) {
      res[j+from_off] += M[i+j*nrow]*x[i];
      res[i]          += M[i+j*nrow]*x[j+from_off];
    }
  }
}


// [[Rcpp::export]]
NumericVector MatOnlyMul(NumericMatrix M,
                         NumericVector x,
                         size_t from_recipient,
                         size_t nthreads) {

  size_t r = (size_t) M.nrow();
  size_t c = (size_t) M.ncol();

  NumericVector res(r);

  if(x.length() != r) {
    Rcout << "x isn't right!\n";
    return(NumericVector(1));
  }
  if(res.length() != r) {
    Rcout << "res isn't right!\n";
    return(NumericVector(1));
  }

  size_t from_off = (from_recipient-1);  //just got rid of the / 2 here

  if(nthreads > 1) {
    double *res_perth = (double*) calloc(r*(nthreads+1), sizeof(double));
    if (res_perth == NULL) {
      printf("Failed allocating res_perth!\n");
      exit(1);
    }


    size_t num_perth = c/nthreads;
    size_t rag_end   = c%nthreads;

    std::vector<std::thread> threads;
    for(size_t i=0; i<nthreads; ++i) {
      threads.push_back(std::thread(
          MatOnlyMul_A,
          &(M[0]), &(x[0]), res_perth + i*r, r, from_off, i*num_perth, num_perth));
    }
    // Tidy ragged end
    if(rag_end != 0) {
      MatOnlyMul_A(&(M[0]), &(x[0]), res_perth + nthreads*r, r, from_off, nthreads*num_perth, rag_end);
    }
    for(auto& th : threads) {
      th.join();
    }

    for(int i = 0; i < res.length(); i++) {
      for(int j = 0; j < nthreads+1; j++) {
        res(i) += res_perth[i+j*r];
      }
      res(i) *= 0.5;
    }

    free(res_perth);
  } else {
    MatOnlyMul_A(&(M[0]), &(x[0]), &(res[0]), r, from_off, 0, c);

    for(int i = 0; i < res.length(); i++) {
      res(i) *= 0.5;
    }
  }

  return(res);
}