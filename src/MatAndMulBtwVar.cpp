#include <Rcpp.h>
using namespace Rcpp;

#include <algorithm>
#include <thread>
#include <cmath>
#include <stdlib.h>




void MatAndMulBtwVar_C_probs(const double* __restrict__ alpha,
                             const double* __restrict__ f,
                             const double* __restrict__ beta,
                             const double* __restrict__ g,
                             const double* __restrict__ x,
                             double* __restrict__ c_1,
                             double* __restrict__ res,
                             bool useunif,
                             size_t j,
                             size_t r,
                             size_t from_off,
                             double fwd_rho,
                             double bck_rho) {
  double z0 = 0.0;
  double rd = (double)(r - 1);
  double ccc = 0.0;

  if(useunif){
    ccc = 1/rd;
  }

  if(*g != 0 & *f != 0){
    double beta_coeff_c1 = (1 - bck_rho) / (*g);
    double alpha_coeff_c1 = (1 - fwd_rho) / (*f);

    double alpha_offset = fwd_rho / rd;

    // We will assume in this accumulation that alpha * beta is always 0 along the matrix diagonal
    for(size_t i = 0; i < r; i++) {
      z0 += c_1[i] = (alpha_coeff_c1 * *(alpha++) + alpha_offset) * (beta_coeff_c1 * *(beta++) + bck_rho) ;
    }
  }

  if(z0 <= 0){
    // if the sum is zero (all of the alpha*beta entries = 0), then all of the distances are beyond the resolution
    // of numerical precision, so we set all of the distances to the -log of the smallest respresentable double
    for(size_t i = 0; i < r; i++) {
      if(i==j){
        c_1[i] = 0.0;
      } else {
        c_1[i] = ccc;
      }
    }

  } else { // z0 > 0

    for(size_t i = 0; i < r; i++) {
      if(i==j) {
        c_1[i] = 0.0;
      } else {

        c_1[i] = c_1[i] / z0;

        res[j+from_off] += c_1[i]*x[i];
        res[i]          += c_1[i]*x[j+from_off];
      }
    }
  }
}

void MatAndMulBtwVar_C_raw_dist(const double* __restrict__ alpha,
                       const double* __restrict__ f,
                       const double* __restrict__ beta,
                       const double* __restrict__ g,
                       const double* __restrict__ x,
                       double* __restrict__ c_1,
                       double* __restrict__ res,
                       size_t j,
                       size_t r,
                       size_t from_off,
                       double fwd_rho,
                       double bck_rho) {
  double z0 = 0.0;
  double rd = (double)(r - 1);

  if(*g != 0 & *f != 0){
    double beta_coeff_c1 = (1 - bck_rho) / (*g);
    double alpha_coeff_c1 = (1 - fwd_rho) / (*f);

    double alpha_offset = fwd_rho / rd;

    // We will assume in this accumulation that alpha * beta is always 0 along the matrix diagonal
    for(size_t i = 0; i < r; i++) {
      z0 += c_1[i] = (alpha_coeff_c1 * *(alpha++) + alpha_offset) * (beta_coeff_c1 * *(beta++) + bck_rho) ;
    }
  }

  if(z0 <= 0){
    // if the sum is zero (all of the alpha*beta entries = 0), then all of the distances are beyond the resolution
    // of numerical precision, so we set all of the distances to the -log of the smallest respresentable double
    for(size_t i = 0; i < r; i++) {
      if(i==j){
        c_1[i] = 0.0;
      } else {
        c_1[i] = 744.4400719213812180897;
      }
    }

  } else { // z0 > 0

    for(size_t i = 0; i < r; i++) {
      if(i==j) {
        c_1[i] = 0.0;
      } else {

        c_1[i] = c_1[i] / z0;

        if(c_1[i] == 0){
          c_1[i] = 744.4400719213812180897;  // -log of smallest representable double
        } else {
          c_1[i] = -log(c_1[i]); // this is the final value for M[i,j]
        }
        res[j+from_off] += c_1[i]*x[i];
        res[i]          += c_1[i]*x[j+from_off];
      }
    }
  }
}



void MatAndMulBtwVar_C_standardize_dist(const double* __restrict__ alpha,
                                   const double* __restrict__ f,
                                   const double* __restrict__ beta,
                                   const double* __restrict__ g,
                                   const double* __restrict__ x,
                                   double* __restrict__ c_1,
                                   double* __restrict__ res,
                                   size_t j,
                                   size_t r,
                                   size_t from_off,
                                   double fwd_rho,
                                   double bck_rho) {
  double z0 = 0.0;
  double z1 = 0.0;
  double z2 = 0.0;
  double rd = (double)(r - 1);


  if(*g != 0 & *f != 0){

    double beta_coeff_c1 = (1 - bck_rho) / (*g);
    double alpha_coeff_c1 = (1 - fwd_rho) / (*f);
    double alpha_offset = fwd_rho / rd;


    // We will assume in this accumulation that alpha * beta is always 0 along the matrix diagonal
    for(size_t i = 0; i < r; i++) {
      z0 += c_1[i] = (alpha_coeff_c1 * *(alpha++) + alpha_offset) * (beta_coeff_c1 * *(beta++) + bck_rho) ;
    }
  }

  // if all of the donors have "0" probability of being copied, so all c_1[i] = 0, we exit here...otherwise...
  if(z0<=0){

    for(size_t i = 0; i < r; i++) {
      c_1[i] = 0.0;
    }

  } else { // z0 > 0

    for(size_t i=0; i < r; i++){
      if(i==j) {
        continue;
      } else {
        c_1[i] = c_1[i] / z0;

        if(c_1[i]==0) {
          c_1[i] = 744.4400719213812180897;  // -log of smallest representable double
        } else {
          c_1[i] = -log(c_1[i]);
        }

        z1 += c_1[i];
        z2 += std::pow(c_1[i],2.0);
      }
    }

    z1 = z1 / rd; // Calculate Mean

    z2 = (z2 - std::pow(z1, 2.0) * rd) / (rd - 1.0); // Unbiased Estimate of the Variance

    if(z2 <= 0){ // if there is no variance in the distances, we just set them all to 0
      for(size_t i = 0; i < r; i++) {
        c_1[i] = 0.0;
      }
    } else {

      z2 = std::pow(z2, -0.5); // Calculate standard deviation

      z1 = -z1 * z2;

      for(size_t i = 0; i < r; i++) {
        if(i==j) {
          c_1[i] = 0.0;
        } else {
          c_1[i] = c_1[i] * z2 + z1; // this is the final value for M[i,j]
          res[j+from_off] += c_1[i]*x[i];
          res[i]          += c_1[i]*x[j+from_off];
        }
      }
    }
  }
}



void MatAndMulBtwVar_B(double* M,
                       const double* __restrict__ alpha,
                       const double* __restrict__ f,
                       const double* __restrict__ beta,
                       const double* __restrict__ g,
                       const double* __restrict__ x,
                       bool standardize,
                       bool probs,
                       bool useunif,
                       double* __restrict__ res,
                       size_t r,
                       size_t from_off,
                       size_t from,
                       size_t N,
                       double fwd_rho,
                       double bck_rho) {

  if(probs){

    for(size_t j = from; j < from+N; j++) {
      MatAndMulBtwVar_C_probs(alpha+j*r,
                              f+j,
                              beta+j*r,
                              g+j,
                              x,
                              M+j*r,
                              res,
                              useunif,
                              j,
                              r,
                              from_off,
                              fwd_rho,
                              bck_rho);
    }

  } else {


    if( standardize == true ){
      for(size_t j = from; j < from+N; j++) {
        MatAndMulBtwVar_C_standardize_dist(alpha+j*r,
                                      f+j,
                                      beta+j*r,
                                      g+j,
                                      x,
                                      M+j*r,
                                      res,
                                      j,
                                      r,
                                      from_off,
                                      fwd_rho,
                                      bck_rho);
      }
    } else {
      for(size_t j = from; j < from+N; j++) {
        MatAndMulBtwVar_C_raw_dist(alpha+j*r,
                          f+j,
                          beta+j*r,
                          g+j,
                          x,
                          M+j*r,
                          res,
                          j,
                          r,
                          from_off,
                          fwd_rho,
                          bck_rho);
      }
    }
  }
}



void MatAndMulBtwVar_A(double* __restrict__ res,
                       double* __restrict__ M,
                       const double* __restrict__ alpha,
                       const double* __restrict__ f,
                       const double* __restrict__ beta,
                       const double* __restrict__ g,
                       const double* __restrict__ x,
                       bool standardize,
                       bool probs,
                       bool useunif,
                       size_t from_off,
                       size_t nthreads,
                       size_t r,
                       size_t c,
                       double fwd_rho,
                       double bck_rho) {

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
          MatAndMulBtwVar_B,
          M, alpha, f, beta, g, x, standardize, probs, useunif, res_perth + i*r, r, from_off, i*num_perth, num_perth, fwd_rho, bck_rho));
    }
    // Tidy ragged end
    if(rag_end != 0) {
      MatAndMulBtwVar_B(M, alpha, f, beta, g, x, standardize, probs, useunif, res_perth + nthreads*r, r, from_off, nthreads*num_perth, rag_end, fwd_rho, bck_rho);
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

    MatAndMulBtwVar_B(M, alpha, f, beta, g, x, standardize, probs, useunif, res, r, from_off, 0, c, fwd_rho, bck_rho);

    for(size_t i = 0; i < r; i++) {
      res[i] *= 0.5;
    }

  }
}




// [[Rcpp::export]]
NumericVector MatAndMulBtwVar(NumericMatrix M,
                              List fwd,
                              List bck,
                              NumericVector x,
                              LogicalVector standardize,
                              LogicalVector calcprobs,
                              LogicalVector unifonunderflow,
                              double fwd_rho,
                              double bck_rho,
                              int from_recipient,
                              int nthreads) {

  M.attr("class") = "kalisDistanceMatrix";

  NumericMatrix alpha = fwd["alpha"];
  NumericMatrix beta = bck["beta"];

  NumericVector f = fwd["alpha.f"];
  NumericVector g = bck["beta.g"];

  size_t r = (size_t) alpha.nrow(); // we need to make p=r and c2 = c
  size_t c = (size_t) alpha.ncol();

  NumericVector res(r);

  if(from_recipient % 2 == 0) {
    Rcout << "from_recipient must be odd. \n";
    return(NumericVector(1));
  }

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
  bool probs = is_true(all(calcprobs));
  bool useunif = is_true(all(unifonunderflow));

  if(stdz == true & r < 3) {
    Rcout << "each matrix must have at least three rows to perform standardization \n";
    return(NumericVector(1));
  }


  MatAndMulBtwVar_A(&(res[0]),
                    &(M[0]),
                    &(alpha[0]),
                    &(f[0]),
                    &(beta[0]),
                    &(g[0]),
                    &(x[0]),
                    stdz,
                    probs,
                    useunif,
                    from_recipient,
                    nthreads,
                    r,
                    c,
                    fwd_rho,
                    bck_rho);

  return(res);
}