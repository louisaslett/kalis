#include <Rcpp.h>
using namespace Rcpp;

#include <algorithm>
#include <thread>
#include <cmath>
#include <stdlib.h>



void DoubleMatAndMul_C(const double* __restrict__ alpha,
                       const double* __restrict__ beta,
                       const double* __restrict__ x,
                       double* __restrict__ c_1,
                       double* __restrict__ c_2,
                       double* __restrict__ res,
                       double* __restrict__ res2,
                       size_t j,
                       size_t r,
                       size_t from_off) {
  double z0 = 0.0;
  double z1 = 0.0;
  double z2 = 0.0;
  double temp = 0.0;
  double rd = (double)(r - 1);

  // Calculate column sums
  for(size_t i = 0; i < r; i++) {
    z0 += c_2[i] = *(alpha++) * *(beta++);
  }

  if(z0<=0) { // if all of the donors have "0" probability of being copied
    for(size_t i=0; i < r; i++){
      if(i!=j) {
        c_1[i] = 744.4400719213812180897; // final value for M_raw[i,j]
        c_2[i] = 0.0; // final value for M_std[i,j]

        res[j+from_off] += c_1[i]*x[i];
        res[i]          += c_1[i]*x[j+from_off];

        res2[j+from_off] += 0.0;
        res2[i]          += 0.0;

      } else {
        c_2[i] = 0.0;
        c_1[i] = 0.0;
      }
    }

  } else {

    for(size_t i=0; i < r; i++){
      if(i!=j) {

        c_2[i] = c_2[i] / z0;

        if(c_2[i]==0) {
          c_1[i] = 744.4400719213812180897;  // -log of smallest representable double
        } else {
          c_1[i] = -log(c_2[i]);
        }

        z1 += c_1[i];
        z2 += std::pow(c_1[i], 2.0);

      }
    }

    z1 = z1 / rd; // Calculate Mean

    z2 = (z2 - std::pow(z1, 2.0) * rd) / (rd - 1.0); // Unbiased Estimate of the Variance

    if(z2 <= 0){ // if there is no variance in the distances

      for(size_t i = 0; i < r; i++) {
        if(i!=j) {
          c_2[i] = 0.0;

          res[j+from_off] += c_1[i]*x[i];
          res[i]          += c_1[i]*x[j+from_off];

          res2[j+from_off] += 0.0;
          res2[i]          += 0.0;

        } else {
          c_1[i] = 0.0;
          c_2[i] = 0.0;
        }
      }

    } else {

      z2 = std::pow(z2, -0.5); // Calculate 1 / standard deviation

      z1 = -z1 * z2;

      for(size_t i = 0; i < r; i++) {
        if(i!=j) {
          c_2[i] = c_1[i] * z2 + z1; // this is the final value for M[i,j]

          res[j+from_off] += c_1[i]*x[i];
          res[i]          += c_1[i]*x[j+from_off];

          res2[j+from_off] += c_2[i]*x[i];
          res2[i]          += c_2[i]*x[j+from_off];

        } else {
          c_1[i] = 0.0;
          c_2[i] = 0.0;
        }
      }
    }
  }
}


void DoubleMatAndMul_B(double* M,
                       double* M2,
                       const double* __restrict__ alpha,
                       const double* __restrict__ beta,
                       const double* __restrict__ x,
                       double* __restrict__ res,
                       double* __restrict__ res2,
                       size_t r,
                       size_t from_off,
                       size_t from,
                       size_t N) {
  for(size_t j = from; j < from+N; j++) {
    DoubleMatAndMul_C(alpha+j*r,
                      beta+j*r,
                      x,
                      M+j*r,
                      M2+j*r,
                      res,
                      res2,
                      j,
                      r,
                      from_off);
  }
}

void DoubleMatAndMul_A(double* __restrict__ res,
                       double* __restrict__ res2,
                       double* __restrict__ M,
                       double* __restrict__ M2,
                       const double* __restrict__ alpha,
                       const double* __restrict__ beta,
                       const double* __restrict__ x,
                       size_t from_off,
                       size_t nthreads,
                       size_t r,
                       size_t c) {

  from_off = (from_off-1); // got rid of /2 here

  if(nthreads > 1) {
    double *res_perth = (double*) calloc(r*(nthreads+1), sizeof(double));
    double *res_perth2 = (double*) calloc(r*(nthreads+1), sizeof(double));
    if (res_perth == NULL) {
      printf("Failed allocating res_perth!\n");
      exit(1);
    }

    size_t num_perth = c/nthreads;
    size_t rag_end   = c%nthreads;

    std::vector<std::thread> threads;
    for(size_t i=0; i<nthreads; ++i) {
      threads.push_back(std::thread(
          DoubleMatAndMul_B,
          M, M2, alpha, beta, x, res_perth + i*r, res_perth2 + i*r, r, from_off, i*num_perth, num_perth));
    }
    // Tidy ragged end
    if(rag_end != 0) {
      DoubleMatAndMul_B(M, M2, alpha, beta, x, res_perth + nthreads*r, res_perth2 + nthreads*r, r, from_off, nthreads*num_perth, rag_end);
    }
    for(auto& th : threads) {
      th.join();
    }

    for(size_t i = 0; i < r; i++) {
      for(size_t j = 0; j < nthreads+1; j++) {
        res[i] += res_perth[i+j*r];
        res2[i] += res_perth2[i+j*r];

      }
      res[i] *= 0.5;
      res2[i] *= 0.5;
    }

    free(res_perth);
    free(res_perth2);

  } else {

    DoubleMatAndMul_B(M, M2, alpha, beta, x, res, res2, r, from_off, 0, c);

    for(size_t i = 0; i < r; i++) {
      res[i] *= 0.5;
      res2[i] *= 0.5;
    }

  }
}




// [[Rcpp::export]]
List DoubleMatAndMul(NumericMatrix M,
                     NumericMatrix M2,
                     List fwd,
                     List bck,
                     NumericVector x,
                     int from_recipient,
                     int nthreads) {

  NumericMatrix alpha = fwd["alpha"];
  NumericMatrix beta = bck["beta"];

  size_t r = (size_t) alpha.nrow(); // we need to make p=r and c2 = c
  size_t c = (size_t) alpha.ncol();

  NumericVector res(r);
  NumericVector res2(r);

  if(r < 3) {
    Rcout << "each matrix must have at least three rows to perform standardization \n";
    List L = List::create(Named("raw") = NumericVector(1) , Named("std") = NumericVector(1));
    return(L);
  }
  if(alpha.nrow() != r || alpha.ncol() != c || beta.nrow() != r || beta.ncol() != c) {
    Rcout << "alpha, beta aren't right!\n";
    List L = List::create(Named("raw") = NumericVector(1) , Named("std") = NumericVector(1));
    return(L);
  }
  if(x.length() != r) {
    Rcout << "x isn't right!\n";
    List L = List::create(Named("raw") = NumericVector(1) , Named("std") = NumericVector(1));
    return(L);
    }
  if(M.nrow() != r || M.ncol() != c) {
    Rcout << "M isn't right!\n";
    List L = List::create(Named("raw") = NumericVector(1) , Named("std") = NumericVector(1));
    return(L);
    }
  if(M2.nrow() != r || M2.ncol() != c) {
    Rcout << "M isn't right!\n";
    List L = List::create(Named("raw") = NumericVector(1) , Named("std") = NumericVector(1));
    return(L);
    }

  if(&(M[0]) == &(M2[0])) {
    Rcout << "M and M2 have the same memory address; they must be created separately in R (eg: do not use M2 <- M as a shortcut) \n";
    List L = List::create(Named("raw") = NumericVector(1) , Named("std") = NumericVector(1));
    return(L);
  }


  if(res.length() != r) {
    Rcout << "res isn't right!\n";
    List L = List::create(Named("raw") = NumericVector(1) , Named("std") = NumericVector(1));
    return(L);
    }

  DoubleMatAndMul_A(&(res[0]),
                    &(res2[0]),
                    &(M[0]),
                    &(M2[0]),
                    &(alpha[0]),
                    &(beta[0]),
                    &(x[0]),
                    from_recipient,
                    nthreads,
                    r,
                    c);

  List L = List::create(Named("raw") = res , Named("std") = res2);

  return L;
}
