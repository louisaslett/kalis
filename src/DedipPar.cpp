#include <Rcpp.h>
using namespace Rcpp;

#include <algorithm>
#include <thread>
#include <stdlib.h>


void DedipPar_C(const double* __restrict__ alpha_c1,
                const double* __restrict__ alpha_c2,
                const double* __restrict__ beta_c1,
                const double* __restrict__ beta_c2,
                const double* __restrict__ x,
                double* __restrict__ M_min_min,
                double* __restrict__ M_min_mean,
                double* __restrict__ M_min_max,
                double* __restrict__ M_min2nd,
                double* __restrict__ res_min_min,
                double* __restrict__ res_min_mean,
                double* __restrict__ res_min_max,
                double* __restrict__ res_min2nd,
                size_t j,
                size_t r,
                size_t from_off) {

  double quad[4];
  double m1, m2, s1, s2;
  double rd = (double)(r - 1);
  m1 = 0.0; m2 = 0.0;
  s1 = 0.0; s2 = 0.0;

  for(size_t i = 0; i < r/2; i++) {
    // each column is an independent HMM
    // quad is arranged
    // quad[0]  quad[2]
    // quad[1]  quad[3]

    quad[0] = *(alpha_c1++) * *(beta_c1++);
    quad[1] = *(alpha_c1++) * *(beta_c1++);
    quad[2] = *(alpha_c2++) * *(beta_c2++);
    quad[3] = *(alpha_c2++) * *(beta_c2++);


    // We always do the off-diagonal elements because we want to count distances between sister chromosomes even when i=j
    if(quad[1] <= 0) {
      m1 += M_min_mean[i] = 708.396418532264078749;
    } else {
      m1 += M_min_mean[i] = -log(quad[1]);
    }
    s1 += std::pow(M_min_mean[i], 2.0);

    if(quad[2] <= 0) {
      m2 += M_min_max[i] = 708.396418532264078749;
    } else {
      m2 += M_min_max[i] = -log(quad[2]);
    }
    s2 += std::pow(M_min_max[i], 2.0);


    if(i!=j) {

      if(quad[0] <= 0) {
        m1 += M_min_min[i] = 708.396418532264078749;
      } else {
        m1 += M_min_min[i] = -log(quad[0]);
      }
      s1 += std::pow(M_min_min[i], 2.0);

      if(quad[3] <= 0) {
        m2 += M_min2nd[i] = 708.396418532264078749;
      } else {
        m2 += M_min2nd[i] = -log(quad[3]);
      }
      s2 += std::pow(M_min2nd[i], 2.0);

    }

  }

  m1 = m1 / rd; // Calculate Mean
  m2 = m2 / rd;

  s1 = (s1 - std::pow(m1, 2.0) * rd) / (rd - 1.0);
  s2 = (s2 - std::pow(m2, 2.0) * rd) / (rd - 1.0);

  s1 = std::pow(s1, -0.5); // Calculate standard deviation
  s2 = std::pow(s2, -0.5);

  m1 = -m1 * s1;
  m2 = -m2 * s2;


  for(size_t i = 0; i < r/2; i++) {
    if(i!=j) {

      quad[0] = M_min_min[i] * s1 + m1;
      quad[1] = M_min_mean[i] * s1 + m1;
      quad[2] = M_min_max[i] * s2 + m2;
      quad[3] = M_min2nd[i] * s2 + m2;

      M_min_mean[i] = 0.5 * fmin(quad[0] + quad[3], quad[1] + quad[2]);
      M_min_max[i] = fmin(fmax(quad[0], quad[3]), fmax(quad[1], quad[2]));

      std::sort(quad,quad+4); //recall that this does the ascending sort in place

      M_min_min[i] = quad[0];
      M_min2nd[i]  = quad[1];

      res_min_min[j+from_off]   += M_min_min[i]*x[i];
      res_min_min[i]            += M_min_min[i]*x[j+from_off];
      res_min_mean[j+from_off]  += M_min_mean[i]*x[i];
      res_min_mean[i]           += M_min_mean[i]*x[j+from_off];
      res_min_max[j+from_off]   += M_min_max[i]*x[i];
      res_min_max[i]            += M_min_max[i]*x[j+from_off];
      res_min2nd[j+from_off]    += M_min2nd[i]*x[i];
      res_min2nd[i]             += M_min2nd[i]*x[j+from_off];

    } else {

      M_min_min[i] = 0.0;
      M_min_mean[i] = 0.0;
      M_min_max[i] = 0.0;
      M_min2nd[i] = 0.0;

    }
  }
}

void DedipPar_B(double* M_min_min,
                double* M_min_mean,
                double* M_min_max,
                double* M_min2nd,
                const double* __restrict__ alpha,
                const double* __restrict__ beta,
                const double* __restrict__ x,
                double* __restrict__ res_min_min,
                double* __restrict__ res_min_mean,
                double* __restrict__ res_min_max,
                double* __restrict__ res_min2nd,
                size_t r,
                size_t from_off,
                size_t from,
                size_t N) {
  for(size_t j = from; j < from+N; j++) {
    DedipPar_C(alpha+2*j*r,
               alpha+(2*j+1)*r,
               beta+2*j*r,
               beta+(2*j+1)*r,
               x,
               M_min_min + j*(r/2),
               M_min_mean + j*(r/2),
               M_min_max + j*(r/2),
               M_min2nd + j*(r/2),
               res_min_min,
               res_min_mean,
               res_min_max,
               res_min2nd,
               j,
               r,
               from_off);
  }
}

void DedipPar_A(double* __restrict__ res_min_min,
                double* __restrict__ res_min_mean,
                double* __restrict__ res_min_max,
                double* __restrict__ res_min2nd,
                double* __restrict__ M_min_min,
                double* __restrict__ M_min_mean,
                double* __restrict__ M_min_max,
                double* __restrict__ M_min2nd,
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
    double *res_perth_min_min = (double*) calloc(p*(nthreads+1), sizeof(double));
    double *res_perth_min_mean = (double*) calloc(p*(nthreads+1), sizeof(double));
    double *res_perth_min_max = (double*) calloc(p*(nthreads+1), sizeof(double));
    double *res_perth_min2nd = (double*) calloc(p*(nthreads+1), sizeof(double));

    if (res_perth_min_min == NULL) {
      printf("Failed allocating res_perth_min_min!\n");
      exit(1);
    }

    size_t num_perth = c2/nthreads;
    size_t rag_end   = c2%nthreads;

    std::vector<std::thread> threads;
    for(size_t i=0; i<nthreads; ++i) {
      threads.push_back(std::thread(
          DedipPar_B,
          M_min_min, M_min_mean, M_min_max, M_min2nd, alpha, beta, x, res_perth_min_min + i*p, res_perth_min_mean + i*p,
          res_perth_min_max + i*p, res_perth_min2nd + i*p, r, from_off, i*num_perth, num_perth));
    }
    // Tidy ragged end
    if(rag_end != 0) {
      DedipPar_B(M_min_min, M_min_mean, M_min_max, M_min2nd, alpha, beta, x, res_perth_min_min + nthreads*p, res_perth_min_mean + nthreads*p,
                 res_perth_min_max + nthreads*p, res_perth_min2nd + nthreads*p, r, from_off, nthreads*num_perth, rag_end);
    }
    for(auto& th : threads) {
      th.join();
    }

    for(size_t i = 0; i < p; i++) {
      for(size_t j = 0; j < nthreads+1; j++) {
        res_min_min[i]  += res_perth_min_min[i+j*p];
        res_min_mean[i] += res_perth_min_mean[i+j*p];
        res_min_max[i]  += res_perth_min_max[i+j*p];
        res_min2nd[i]   += res_perth_min2nd[i+j*p];
      }
      res_min_min[i]  *= 0.5;
      res_min_mean[i] *= 0.5;
      res_min_max[i]  *= 0.5;
      res_min2nd[i]   *= 0.5;
    }

    free(res_perth_min_min);
    free(res_perth_min_mean);
    free(res_perth_min_max);
    free(res_perth_min2nd);

  } else {

    DedipPar_B(M_min_min, M_min_mean, M_min_max, M_min2nd, alpha, beta, x, res_min_min, res_min_mean,
               res_min_max, res_min2nd, r, from_off, 0, c2);

    for(size_t i = 0; i < p; i++) {
      res_min_min[i]  *= 0.5;
      res_min_mean[i] *= 0.5;
      res_min_max[i]  *= 0.5;
      res_min2nd[i]   *= 0.5;
    }
  }
}




// [[Rcpp::export]]
List DedipPar(List M,
              List fwd,
              List bck,
              NumericVector x,
              int from_recipient,
              int nthreads) {

  NumericMatrix alpha = fwd["alpha"];
  NumericMatrix beta = bck["beta"];

  size_t r = (size_t) alpha.nrow();
  size_t c = (size_t) alpha.ncol();
  size_t p = (size_t) r/2;
  size_t c2 = (size_t) c/2;
  //Rcout << "r = " << r << ", c = " << c << ", p = " << p << ", c2 = " << c2 << "\n";

  // Unpack M
  NumericMatrix M_min_min   = M[0];
  NumericMatrix M_min_mean  = M[1];
  NumericMatrix M_min_max   = M[2];
  NumericMatrix M_min2nd    = M[3];

  // Initialize Output
  NumericVector res_min_min(p); // formerly called phi_min
  NumericVector res_min_mean(p); // fomerly called phi_add
  NumericVector res_min_max(p); // formerly called phi_dom
  NumericVector res_min2nd(p); // still overall 2nd min

  if(from_recipient % 2 == 0) {
    Rcout << "from_recipient must be odd. \n";
    List L = List::create(Named("min_min") = NumericVector(1) , Named("min_mean") = NumericVector(1) ,
                          Named("min_max") = NumericVector(1) , Named("min2nd") = NumericVector(1));
    return(L);
  }

  if(alpha.nrow() != r || alpha.ncol() != c || beta.nrow() != r || beta.ncol() != c) {
    Rcout << "alpha, beta aren't right!\n";
    List L = List::create(Named("min_min") = NumericVector(1) , Named("min_mean") = NumericVector(1) ,
                          Named("min_max") = NumericVector(1) , Named("min2nd") = NumericVector(1));
    return(L);  }

  if(x.length() != p) {
    Rcout << "x isn't right!\n";
    List L = List::create(Named("min_min") = NumericVector(1) , Named("min_mean") = NumericVector(1) ,
                          Named("min_max") = NumericVector(1) , Named("min2nd") = NumericVector(1));
    return(L);  }

  if(r < 3) {
    Rcout << "each matrix must have at least three rows to perform standardization \n";
    List L = List::create(Named("min_min") = NumericVector(1) , Named("min_mean") = NumericVector(1) ,
                          Named("min_max") = NumericVector(1) , Named("min2nd") = NumericVector(1));
    return(L);
  }

  if(&(M_min_min[0]) == &(M_min_mean[0])) {
    Rcout << "All of the matrices in the provided list M must be created separately in R so that they have different memory addresses (eg: do not use M[[2]] <- M[[1]] as a shortcut) \n";
    List L = List::create(Named("min_min") = NumericVector(1) , Named("min_mean") = NumericVector(1) ,
                          Named("min_max") = NumericVector(1) , Named("min2nd") = NumericVector(1));
    return(L);
  }


  DedipPar_A(&(res_min_min[0]),
             &(res_min_mean[0]),
             &(res_min_max[0]),
             &(res_min2nd[0]),
             &(M_min_min[0]),
             &(M_min_mean[0]),
             &(M_min_max[0]),
             &(M_min2nd[0]),
             &(alpha[0]),
             &(beta[0]),
             &(x[0]),
             from_recipient,
             nthreads,
             r,
             c,
             c2,
             p);

  List L = List::create(Named("min_min") = res_min_min , Named("min_mean") = res_min_mean ,
                        Named("min_max") = res_min_max , Named("min2nd") = res_min2nd);

  return L;
}

