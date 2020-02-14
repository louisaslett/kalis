#include <Rcpp.h>
using namespace Rcpp;

#include <algorithm>
#include <thread>
#include <stdlib.h>


void DedipParBtwVarScalarPi_C(const double* __restrict__ alpha_c1,
                const double* __restrict__ alpha_c2,
                const double* __restrict__ f_c1,
                const double* __restrict__ f_c2,
                const double* __restrict__ beta_c1,
                const double* __restrict__ beta_c2,
                const double* __restrict__ g_c1,
                const double* __restrict__ g_c2,
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
                size_t from_off,
                double fwd_rho,
                double bck_rho) {

  double quad[4];
  double z1 = 0.0;
  double z2 = 0.0;
  double m1, m2, s1, s2;
  double rd = (double)(r - 1);
  m1 = 0.0; m2 = 0.0;
  s1 = 0.0; s2 = 0.0;

  double beta_coeff_c1 = (1 - bck_rho) / (*g_c1);
  double alpha_coeff_c1 = (1 - fwd_rho) / (*f_c1);

  double beta_coeff_c2 = (1 - bck_rho) / (*g_c2);
  double alpha_coeff_c2 = (1 - fwd_rho) / (*f_c2);

  double alpha_offset = fwd_rho / rd;

  for(size_t i = 0; i < r/2; i++) {
    // each column is an independent HMM
    // quad is arranged
    // M_min_min[i]  M_min_max[i]
    // M_min_mean[i]  M_min2nd[i]

    // We always do the off-diagonal elements because we want to count distances between sister chromosomes even when i=j
    z1 += M_min_min[i]  = (alpha_coeff_c1 * *(alpha_c1++) + alpha_offset) * (beta_coeff_c1 * *(beta_c1++) + bck_rho) ;
    z1 += M_min_mean[i] = (alpha_coeff_c1 * *(alpha_c1++) + alpha_offset) * (beta_coeff_c1 * *(beta_c1++) + bck_rho) ;
    z2 += M_min_max[i]  = (alpha_coeff_c2 * *(alpha_c2++) + alpha_offset) * (beta_coeff_c2 * *(beta_c2++) + bck_rho) ;
    z2 += M_min2nd[i]   = (alpha_coeff_c2 * *(alpha_c2++) + alpha_offset) * (beta_coeff_c2 * *(beta_c2++) + bck_rho) ;
  }

  // NEED TO ADD IN HANDLING HERE FOR THE CASE WHERE THE COLSUM IS 0!!!!!!!!

  for(size_t i = 0; i < r/2; i++) {

    if(z1 <= 0){
      M_min_mean[i] = 0.0;
    } else {

      M_min_mean[i] = M_min_mean[i] / z1;

      if(M_min_mean[i] == 0) {
        M_min_mean[i] = 744.4400719213812180897;  // -log of smallest representable double
      } else {
        M_min_mean[i] = -log(M_min_mean[i]);
      }

      m1 += M_min_mean[i];
      s1 += std::pow(M_min_mean[i], 2.0);
    }

    if(z2 <= 0){
      M_min_max[i] = 0.0;
    } else {

      M_min_max[i] = M_min_max[i] / z2;

      if(M_min_max[i] == 0) {
        M_min_max[i] = 744.4400719213812180897;  // -log of smallest representable double
      } else {
        M_min_max[i] = -log(M_min_max[i]);
      }

      m2 += M_min_max[i];
      s2 += std::pow(M_min_max[i], 2.0);
    }


    if(i!=j) {

      if(z1 <= 0){
        M_min_min[i] = 0.0;
      } else {

        M_min_min[i] = M_min_min[i] / z1;

        if(M_min_min[i] == 0) {
          M_min_min[i] = 744.4400719213812180897;  // -log of smallest representable double
        } else {
          M_min_min[i] = -log(M_min_min[i]);
        }

        m1 += M_min_min[i];
        s1 += std::pow(M_min_min[i], 2.0);
      }

      if(z2 <= 0){
        M_min_min[i] = 0.0;
      } else {
        M_min2nd[i] = M_min2nd[i] / z2;

        if(M_min2nd[i] == 0) {
          M_min2nd[i] = 744.4400719213812180897;  // -log of smallest representable double
        } else {
          M_min2nd[i] = -log(M_min2nd[i]);
        }

        m2 += M_min2nd[i];
        s2 += std::pow(M_min2nd[i], 2.0);
      }

    }

  }

  m1 = m1 / rd; // Calculate Mean
  m2 = m2 / rd;

  s1 = (s1 - std::pow(m1, 2.0) * rd) / (rd - 1.0);
  s2 = (s2 - std::pow(m2, 2.0) * rd) / (rd - 1.0);

  // The lines here below where we set s1 and s2 handle both the case where all of the raw distances are 0
  // in a column and the case where the varance is zero -- in both cases, all standardized distances remain 0.

  if(s1!=0) {
    s1 = std::pow(s1, -0.5); // Calculate standard deviation
  }

  if(s2!=0) {
    s2 = std::pow(s2, -0.5); // Calculate standard deviation
  }

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

void DedipParBtwVarScalarPi_B(double* M_min_min,
                double* M_min_mean,
                double* M_min_max,
                double* M_min2nd,
                const double* __restrict__ alpha,
                const double* __restrict__ f,
                const double* __restrict__ beta,
                const double* __restrict__ g,
                const double* __restrict__ x,
                double* __restrict__ res_min_min,
                double* __restrict__ res_min_mean,
                double* __restrict__ res_min_max,
                double* __restrict__ res_min2nd,
                size_t r,
                size_t from_off,
                size_t from,
                size_t N,
                double fwd_rho,
                double bck_rho) {
  for(size_t j = from; j < from+N; j++) {
    DedipParBtwVarScalarPi_C(alpha+2*j*r,
               alpha+(2*j+1)*r,
               f+(2*j),
               f+(2*j+1),
               beta+2*j*r,
               beta+(2*j+1)*r,
               g+(2*j),
               g+(2*j+1),
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
               from_off,
               fwd_rho,
               bck_rho);
  }
}

void DedipParBtwVarScalarPi_A(double* __restrict__ res_min_min,
                double* __restrict__ res_min_mean,
                double* __restrict__ res_min_max,
                double* __restrict__ res_min2nd,
                double* __restrict__ M_min_min,
                double* __restrict__ M_min_mean,
                double* __restrict__ M_min_max,
                double* __restrict__ M_min2nd,
                const double* __restrict__ alpha,
                const double* __restrict__ f,
                const double* __restrict__ beta,
                const double* __restrict__ g,
                const double* __restrict__ x,
                size_t from_off,
                size_t nthreads,
                size_t r,
                size_t c,
                size_t c2,
                size_t p,
                double fwd_rho,
                double bck_rho) {

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
          DedipParBtwVarScalarPi_B,
          M_min_min, M_min_mean, M_min_max, M_min2nd, alpha, f, beta, g, x, res_perth_min_min + i*p, res_perth_min_mean + i*p,
          res_perth_min_max + i*p, res_perth_min2nd + i*p, r, from_off, i*num_perth, num_perth, fwd_rho, bck_rho));
    }
    // Tidy ragged end
    if(rag_end != 0) {
      DedipParBtwVarScalarPi_B(M_min_min, M_min_mean, M_min_max, M_min2nd, alpha, f, beta, g, x, res_perth_min_min + nthreads*p, res_perth_min_mean + nthreads*p,
                 res_perth_min_max + nthreads*p, res_perth_min2nd + nthreads*p, r, from_off, nthreads*num_perth, rag_end, fwd_rho, bck_rho);
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

    DedipParBtwVarScalarPi_B(M_min_min, M_min_mean, M_min_max, M_min2nd, alpha, f, beta, g, x, res_min_min, res_min_mean,
               res_min_max, res_min2nd, r, from_off, 0, c2, fwd_rho, bck_rho);

    for(size_t i = 0; i < p; i++) {
      res_min_min[i]  *= 0.5;
      res_min_mean[i] *= 0.5;
      res_min_max[i]  *= 0.5;
      res_min2nd[i]   *= 0.5;
    }
  }
}




// [[Rcpp::export]]
List DedipParBtwVarScalarPi(List M,
                            List fwd,
                            List bck,
                            NumericVector x,
                            double fwd_rho,
                            double bck_rho,
                            int from_recipient,
                            int nthreads) {

  NumericMatrix alpha = fwd["alpha"];
  NumericMatrix beta = bck["beta"];

  NumericVector f = fwd["alpha.f"];
  NumericVector g = bck["beta.g"];

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


  DedipParBtwVarScalarPi_A(&(res_min_min[0]),
                           &(res_min_mean[0]),
                           &(res_min_max[0]),
                           &(res_min2nd[0]),
                           &(M_min_min[0]),
                           &(M_min_mean[0]),
                           &(M_min_max[0]),
                           &(M_min2nd[0]),
                           &(alpha[0]),
                           &(f[0]),
                           &(beta[0]),
                           &(g[0]),
                           &(x[0]),
                           from_recipient,
                           nthreads,
                           r,
                           c,
                           c2,
                           p,
                           fwd_rho,
                           bck_rho);

  List L = List::create(Named("min_min") = res_min_min , Named("min_mean") = res_min_mean ,
                        Named("min_max") = res_min_max , Named("min2nd") = res_min2nd);

  return L;
}

