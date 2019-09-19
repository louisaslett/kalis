#include <Rcpp.h>
using namespace Rcpp;

#include <algorithm>
#include <thread>
#include <stdlib.h>









void DedipAndMul_C(const double* __restrict__ alpha_c1,
                   const double* __restrict__ alpha_c2,
                   const double* __restrict__ beta_c1,
                   const double* __restrict__ beta_c2,
                   const double* __restrict__ x,
                   double* __restrict__ c_1,
                   double* __restrict__ c_2,
                   double* __restrict__ res,
                   size_t j,
                   size_t r,
                   size_t from_off) {
  double c11, c21, c12, c22, z1, z2;
  z1 = 0.0; z2 = 0.0;
  for(size_t i = 0; i < r/2; i++) {
    z1 += c11 = *(alpha_c1++) * *(beta_c1++);
    z2 += c12 = *(alpha_c2++) * *(beta_c2++);
    z1 += c21 = *(alpha_c1++) * *(beta_c1++);
    z2 += c22 = *(alpha_c2++) * *(beta_c2++);

    c_1[i] = fmax(c11, c21);
    c_2[i] = fmax(c12, c22);
  }

  for(size_t i = 0; i < r/2; i++) {
    c_1[i] = -log(fmax(c_1[i]/z1, c_2[i]/z2)); // this is the final value for M[i,j]

    res[j+from_off] += c_1[i]*x[i];
    res[i]          += c_1[i]*x[j+from_off];
  }
}

void DedipAndMul_B(double* M,
                   const double* __restrict__ alpha,
                   const double* __restrict__ beta,
                   const double* __restrict__ x,
                   double* __restrict__ res,
                   double* __restrict__ tmp,
                   size_t r,
                   size_t from_off,
                   size_t from,
                   size_t N) {
  for(size_t j = from; j < from+N; j++) {
    DedipAndMul_C(alpha+2*j*r,
                  alpha+(2*j+1)*r,
                  beta+2*j*r,
                  beta+(2*j+1)*r,
                  x,
                  M+j*(r/2),
                  tmp,
                  res,
                  j,
                  r,
                  from_off);
  }
}

void DedipAndMul_A(double* __restrict__ res,
                   double* __restrict__ M,
                   const double* __restrict__ alpha,
                   const double* __restrict__ beta,
                   const double* __restrict__ x,
                   size_t from_off,
                   size_t nthreads,
                   size_t r,
                   size_t c,
                   size_t c2,
                   size_t p) {

  from_off = (from_off-1)/2;

  if(nthreads > 1) {
    double *res_perth = (double*) calloc(p*(nthreads+1), sizeof(double));
    if (res_perth == NULL) {
      printf("Failed allocating res_perth!\n");
      exit(1);
    }
    double *tmp       = (double*) calloc(p*(nthreads+1), sizeof(double));
    if (tmp == NULL) {
      printf("Failed allocating tmp!\n");
      exit(1);
    }

    size_t num_perth = c2/nthreads;
    size_t rag_end   = c2%nthreads;

    std::vector<std::thread> threads;
    for(size_t i=0; i<nthreads; ++i) {
      threads.push_back(std::thread(
          DedipAndMul_B,
          M, alpha, beta, x, res_perth + i*p, tmp + i*p, r, from_off, i*num_perth, num_perth));
    }
    // Tidy ragged end
    if(rag_end != 0) {
      DedipAndMul_B(M, alpha, beta, x, res_perth+nthreads*p, tmp+nthreads*p, r, from_off, nthreads*num_perth, rag_end);
    }
    for(auto& th : threads) {
      th.join();
    }

    for(size_t i = 0; i < p; i++) {
      for(size_t j = 0; j < nthreads+1; j++) {
        res[i] += res_perth[i+j*p];
      }
      res[i] *= 0.5;
    }

    free(res_perth);
    free(tmp);
  } else {
    double *tmp = (double*) calloc(p, sizeof(double));
    if (tmp == NULL) {
      printf("Failed allocating!\n");
      exit(1);
    }

    DedipAndMul_B(M, alpha, beta, x, res, tmp, r, from_off, 0, c2);

    for(size_t i = 0; i < p; i++) {
     res[i] *= 0.5;
    }

    free(tmp);
  }
}




// [[Rcpp::export]]
NumericVector DedipAndMul(NumericMatrix M,
                          NumericMatrix alpha,
                          NumericMatrix beta,
                          NumericVector x,
                          int from_recipient,
                          int nthreads) {

  size_t r = (size_t) alpha.nrow();
  size_t c = (size_t) alpha.ncol();
  size_t p = (size_t) r/2;
  size_t c2 = (size_t) c/2;
  //Rcout << "r = " << r << ", c = " << c << ", p = " << p << ", c2 = " << c2 << "\n";

  NumericVector res(p);

  if(alpha.nrow() != r || alpha.ncol() != c || beta.nrow() != r || beta.ncol() != c) {
    Rcout << "alpha, beta aren't right!\n";
    return(NumericVector(1));
  }
  if(x.length() != p) {
    Rcout << "x isn't right!\n";
    return(NumericVector(1));
  }
  if(M.nrow() != p || M.ncol() != c2) {
    Rcout << "M isn't right!\n";
    return(NumericVector(1));
  }
  if(res.length() != p) {
    Rcout << "res isn't right!\n";
    return(NumericVector(1));
  }

  DedipAndMul_A(&(res[0]),
                &(M[0]),
                &(alpha[0]),
                &(beta[0]),
                &(x[0]),
                from_recipient,
                nthreads,
                r,
                c,
                c2,
                p);

  return(res);
}


void DedipOnlyMul_A(const double* __restrict__ M,
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
NumericVector DedipOnlyMul(NumericMatrix M,
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

  size_t from_off = (from_recipient-1)/2;

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
          DedipOnlyMul_A,
          &(M[0]), &(x[0]), res_perth + i*r, r, from_off, i*num_perth, num_perth));
    }
    // Tidy ragged end
    if(rag_end != 0) {
      DedipOnlyMul_A(&(M[0]), &(x[0]), res_perth + nthreads*r, r, from_off, nthreads*num_perth, rag_end);
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
    DedipOnlyMul_A(&(M[0]), &(x[0]), &(res[0]), r, from_off, 0, c);

    for(int i = 0; i < res.length(); i++) {
      res(i) *= 0.5;
    }
  }

  return(res);
}

