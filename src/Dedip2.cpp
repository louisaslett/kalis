#include <Rcpp.h>
using namespace Rcpp;

#include <algorithm>

// [[Rcpp::export]]
NumericVector Dedip2_min(NumericMatrix fwd, NumericMatrix bck, NumericVector s) {
  int N = fwd.ncol();
  if(N%2!=0) {
    stop("Dimension not even.");
  }

  NumericVector res(((N/2+1)*N/2)/2);
  int k=0;
  NumericVector quad(4);
  for(int i=0; i<N/2; i++) {
    for(int j=i; j<N/2; j++) {
      if(i==j) {
        res[k++] = 0.0;
        continue;
      }
      quad[0] = ((fwd(2*i  ,2*j)   * bck(2*i,  2*j))   * (fwd(2*j  ,2*i)   * bck(2*j  ,2*i))) / (s[2*j] * s[2*i]) ;
      quad[1] = ((fwd(2*i  ,2*j+1) * bck(2*i,  2*j+1)) * (fwd(2*j+1,2*i)   * bck(2*j+1,2*i))) / (s[2*j+1] * s[2*i]);
      quad[2] = ((fwd(2*i+1,2*j)   * bck(2*i+1,2*j))   * (fwd(2*j  ,2*i+1) * bck(2*j  ,2*i+1))) / (s[2*j] * s[2*i+1]);
      quad[3] = ((fwd(2*i+1,2*j+1) * bck(2*i+1,2*j+1)) * (fwd(2*j+1,2*i+1) * bck(2*j+1,2*i+1))) / (s[2*j+1] * s[2*i+1]);
      res[k++] = -0.5*log(max(quad));
    }
  }

  return(res);
}

// [[Rcpp::export]]
NumericVector Dedip2_2nd_min(NumericMatrix fwd, NumericMatrix bck, NumericVector s) {
  int N = fwd.ncol();
  if(N%2!=0) {
    stop("Dimension not even.");
  }

  NumericVector res(((N/2+1)*N/2)/2);
  int k=0;
  NumericVector quad(4);
  for(int i=0; i<N/2; i++) {
    for(int j=i; j<N/2; j++) {
      if(i==j) {
        res[k++] = 0.0;
        continue;
      }
      quad[0] = ((fwd(2*i  ,2*j)   * bck(2*i,  2*j))   * (fwd(2*j  ,2*i)   * bck(2*j  ,2*i))) / (s[2*j] * s[2*i]) ;
      quad[1] = ((fwd(2*i  ,2*j+1) * bck(2*i,  2*j+1)) * (fwd(2*j+1,2*i)   * bck(2*j+1,2*i))) / (s[2*j+1] * s[2*i]);
      quad[2] = ((fwd(2*i+1,2*j)   * bck(2*i+1,2*j))   * (fwd(2*j  ,2*i+1) * bck(2*j  ,2*i+1))) / (s[2*j] * s[2*i+1]);
      quad[3] = ((fwd(2*i+1,2*j+1) * bck(2*i+1,2*j+1)) * (fwd(2*j+1,2*i+1) * bck(2*j+1,2*i+1))) / (s[2*j+1] * s[2*i+1]);
      res[k++] = -0.5*log(quad.sort()[2]);
    }
  }

  return(res);
}

// [[Rcpp::export]]
NumericVector Dedip2_max(NumericMatrix fwd, NumericMatrix bck, NumericVector s) {
  int N = fwd.ncol();
  if(N%2!=0) {
    stop("Dimension not even.");
  }

  NumericVector res(((N/2+1)*N/2)/2);
  int k=0;
  NumericVector quad(4);
  for(int i=0; i<N/2; i++) {
    for(int j=i; j<N/2; j++) {
      if(i==j) {
        res[k++] = 0.0;
        continue;
      }
      quad[0] = ((fwd(2*i  ,2*j)   * bck(2*i,  2*j))   * (fwd(2*j  ,2*i)   * bck(2*j  ,2*i))) / (s[2*j] * s[2*i]) ;
      quad[1] = ((fwd(2*i  ,2*j+1) * bck(2*i,  2*j+1)) * (fwd(2*j+1,2*i)   * bck(2*j+1,2*i))) / (s[2*j+1] * s[2*i]);
      quad[2] = ((fwd(2*i+1,2*j)   * bck(2*i+1,2*j))   * (fwd(2*j  ,2*i+1) * bck(2*j  ,2*i+1))) / (s[2*j] * s[2*i+1]);
      quad[3] = ((fwd(2*i+1,2*j+1) * bck(2*i+1,2*j+1)) * (fwd(2*j+1,2*i+1) * bck(2*j+1,2*i+1))) / (s[2*j+1] * s[2*i+1]);
      res[k++] = -0.5*log(min(quad));
    }
  }

  return(res);
}

// [[Rcpp::export]]
NumericVector Dedip2_dom(NumericMatrix fwd, NumericMatrix bck, NumericVector s) {
  int N = fwd.ncol();
  if(N%2!=0) {
    stop("Dimension not even.");
  }

  NumericVector res(((N/2+1)*N/2)/2);
  int k=0;
  NumericVector quad(4);
  for(int i=0; i<N/2; i++) {
    for(int j=i; j<N/2; j++) {
      if(i==j) {
        res[k++] = 0.0;
        continue;
      }
      quad[0] = ((fwd(2*i  ,2*j)   * bck(2*i,  2*j))   * (fwd(2*j  ,2*i)   * bck(2*j  ,2*i))) / (s[2*j] * s[2*i]) ;
      quad[1] = ((fwd(2*i  ,2*j+1) * bck(2*i,  2*j+1)) * (fwd(2*j+1,2*i)   * bck(2*j+1,2*i))) / (s[2*j+1] * s[2*i]);
      quad[2] = ((fwd(2*i+1,2*j)   * bck(2*i+1,2*j))   * (fwd(2*j  ,2*i+1) * bck(2*j  ,2*i+1))) / (s[2*j] * s[2*i+1]);
      quad[3] = ((fwd(2*i+1,2*j+1) * bck(2*i+1,2*j+1)) * (fwd(2*j+1,2*i+1) * bck(2*j+1,2*i+1))) / (s[2*j+1] * s[2*i+1]);
      res[k++] = -0.5*log(std::max(std::min(quad[0], quad[3]), std::min(quad[1], quad[2])));
    }
  }

  return(res);
}

// [[Rcpp::export]]
NumericVector Dedip2_add(NumericMatrix fwd, NumericMatrix bck, NumericVector s) {
  int N = fwd.ncol();
  if(N%2!=0) {
    stop("Dimension not even.");
  }

  NumericVector res(((N/2+1)*N/2)/2);
  int k=0;
  double tot1, tot2;
  for(int i=0; i<N/2; i++) {
    for(int j=i; j<N/2; j++) {
      if(i==j) {
        res[k++] = 0.0;
        continue;
      }
      tot1 = (fwd(2*i  ,2*j) * bck(2*i,  2*j)) * (fwd(2*j  ,2*i) * bck(2*j  ,2*i))
      * (fwd(2*i+1,2*j+1) * bck(2*i+1,2*j+1)) * (fwd(2*j+1,2*i+1) * bck(2*j+1,2*i+1));
      tot2 = (fwd(2*i  ,2*j+1) * bck(2*i,  2*j+1)) * (fwd(2*j+1,2*i) * bck(2*j+1,2*i))
      * (fwd(2*i+1,2*j) * bck(2*i+1,2*j)) * (fwd(2*j  ,2*i+1) * bck(2*j  ,2*i+1));
      res[k++] = -0.25*log(std::max(tot1, tot2) / ((s[2*i] * s[2*i+1])*(s[2*j]*s[2*j+1])));
    }
  }

  return(res);
}

// [[Rcpp::export]]
NumericVector Dedip2_mean(NumericMatrix fwd, NumericMatrix bck, NumericVector s) {
  int N = fwd.ncol();
  if(N%2!=0) {
    stop("Dimension not even.");
  }

  NumericVector res(((N/2+1)*N/2)/2);
  int k=0;
  double tot;
  double normconst;
  for(int i=0; i<N/2; i++) {
    for(int j=i; j<N/2; j++) {
      if(i==j) {
        res[k++] = 0.0;
        continue;
      }
      tot = (fwd(2*i  ,2*j) * bck(2*i,  2*j)) * (fwd(2*j  ,2*i) * bck(2*j  ,2*i))
      * (fwd(2*i  ,2*j+1) * bck(2*i,  2*j+1)) * (fwd(2*j+1,2*i) * bck(2*j+1,2*i))
      * (fwd(2*i+1,2*j) * bck(2*i+1,2*j)) * (fwd(2*j  ,2*i+1) * bck(2*j  ,2*i+1))
      * (fwd(2*i+1,2*j+1) * bck(2*i+1,2*j+1)) * (fwd(2*j+1,2*i+1) * bck(2*j+1,2*i+1));
      normconst = (s[2*i] * s[2*i+1])*(s[2*j]*s[2*j+1]);
      res[k++] = -0.125*log(tot / (normconst * normconst) );
    }
  }

  return(res);
}

// [[Rcpp::export]]
List Dedip2_all(NumericMatrix fwd, NumericMatrix bck, NumericVector s) {
  int N = fwd.ncol();
  if(N%2!=0) {
    stop("Dimension not even.");
  }

  NumericVector res_min(((N/2+1)*N/2)/2);
  NumericVector res_min2nd(((N/2+1)*N/2)/2);
  NumericVector res_dom(((N/2+1)*N/2)/2);
  NumericVector res_max(((N/2+1)*N/2)/2);
  NumericVector res_add(((N/2+1)*N/2)/2);
  NumericVector res_mean(((N/2+1)*N/2)/2);
  int k=0;
  NumericVector quad(4);
  for(int i=0; i<N/2; i++) {
    for(int j=i; j<N/2; j++) {
      if(i==j) {
        res_min[k] = res_max[k] = res_mean[k] = 0.0;
        k++;
        continue;
      }
      quad[0] = ((fwd(2*i  ,2*j)   * bck(2*i,  2*j))   * (fwd(2*j  ,2*i)   * bck(2*j  ,2*i))) / (s[2*j] * s[2*i]) ;
      quad[1] = ((fwd(2*i  ,2*j+1) * bck(2*i,  2*j+1)) * (fwd(2*j+1,2*i)   * bck(2*j+1,2*i))) / (s[2*j+1] * s[2*i]);
      quad[2] = ((fwd(2*i+1,2*j)   * bck(2*i+1,2*j))   * (fwd(2*j  ,2*i+1) * bck(2*j  ,2*i+1))) / (s[2*j] * s[2*i+1]);
      quad[3] = ((fwd(2*i+1,2*j+1) * bck(2*i+1,2*j+1)) * (fwd(2*j+1,2*i+1) * bck(2*j+1,2*i+1))) / (s[2*j+1] * s[2*i+1]);
      res_min[k]    = -0.5*log(max(quad));
      res_min2nd[k] = -0.5*log(quad.sort()[2]);
      res_dom[k]    = -0.5*log(std::max(std::min(quad[0], quad[3]), std::min(quad[1], quad[2])));
      res_max[k]    = -0.5*log(min(quad));
      res_add[k]    = -0.25*log(std::max(quad[0]*quad[3], quad[1]*quad[2]));
      res_mean[k]   = -0.125*log((quad[0]*quad[1])*(quad[2]*quad[3]));
      k++;
    }
  }

  List res;
  res["min"] = res_min;
  res["min2nd"] = res_min2nd;
  res["dom"] = res_dom;
  res["max"] = res_max;
  res["add"] = res_add;
  res["mean"] = res_mean;
  return(res);
}
