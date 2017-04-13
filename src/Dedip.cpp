#include <Rcpp.h>
using namespace Rcpp;

#include <algorithm>

// [[Rcpp::export]]
NumericVector Dedip_min(NumericMatrix fwd, NumericMatrix bck) {
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
      quad[0] = (fwd(2*i  ,2*j)   * bck(2*i,  2*j))   * (fwd(2*j  ,2*i)   * bck(2*j  ,2*i));
      quad[1] = (fwd(2*i  ,2*j+1) * bck(2*i,  2*j+1)) * (fwd(2*j+1,2*i)   * bck(2*j+1,2*i));
      quad[2] = (fwd(2*i+1,2*j)   * bck(2*i+1,2*j))   * (fwd(2*j  ,2*i+1) * bck(2*j  ,2*i+1));
      quad[3] = (fwd(2*i+1,2*j+1) * bck(2*i+1,2*j+1)) * (fwd(2*j+1,2*i+1) * bck(2*j+1,2*i+1));
      res[k++] = -0.5*log(max(quad));
    }
  }

  return(res);
}

// [[Rcpp::export]]
NumericVector Dedip_2nd_min(NumericMatrix fwd, NumericMatrix bck) {
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
      quad[0] = (fwd(2*i  ,2*j)   * bck(2*i,  2*j))   * (fwd(2*j  ,2*i)   * bck(2*j  ,2*i));
      quad[1] = (fwd(2*i  ,2*j+1) * bck(2*i,  2*j+1)) * (fwd(2*j+1,2*i)   * bck(2*j+1,2*i));
      quad[2] = (fwd(2*i+1,2*j)   * bck(2*i+1,2*j))   * (fwd(2*j  ,2*i+1) * bck(2*j  ,2*i+1));
      quad[3] = (fwd(2*i+1,2*j+1) * bck(2*i+1,2*j+1)) * (fwd(2*j+1,2*i+1) * bck(2*j+1,2*i+1));
      res[k++] = -0.5*log(quad.sort()[2]);
    }
  }

  return(res);
}

// [[Rcpp::export]]
NumericVector Dedip_max(NumericMatrix fwd, NumericMatrix bck) {
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
      quad[0] = (fwd(2*i  ,2*j)   * bck(2*i,  2*j))   * (fwd(2*j  ,2*i)   * bck(2*j  ,2*i));
      quad[1] = (fwd(2*i  ,2*j+1) * bck(2*i,  2*j+1)) * (fwd(2*j+1,2*i)   * bck(2*j+1,2*i));
      quad[2] = (fwd(2*i+1,2*j)   * bck(2*i+1,2*j))   * (fwd(2*j  ,2*i+1) * bck(2*j  ,2*i+1));
      quad[3] = (fwd(2*i+1,2*j+1) * bck(2*i+1,2*j+1)) * (fwd(2*j+1,2*i+1) * bck(2*j+1,2*i+1));
      res[k++] = -0.5*log(min(quad));
    }
  }

  return(res);
}

// [[Rcpp::export]]
NumericVector Dedip_dom(NumericMatrix fwd, NumericMatrix bck) {
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
      quad[0] = (fwd(2*i  ,2*j)   * bck(2*i,  2*j))   * (fwd(2*j  ,2*i)   * bck(2*j  ,2*i));
      quad[1] = (fwd(2*i  ,2*j+1) * bck(2*i,  2*j+1)) * (fwd(2*j+1,2*i)   * bck(2*j+1,2*i));
      quad[2] = (fwd(2*i+1,2*j)   * bck(2*i+1,2*j))   * (fwd(2*j  ,2*i+1) * bck(2*j  ,2*i+1));
      quad[3] = (fwd(2*i+1,2*j+1) * bck(2*i+1,2*j+1)) * (fwd(2*j+1,2*i+1) * bck(2*j+1,2*i+1));
      res[k++] = -0.5*log(std::max(std::min(quad[0], quad[3]), std::min(quad[1], quad[2])));
    }
  }

  return(res);
}

// [[Rcpp::export]]
NumericVector Dedip_mean(NumericMatrix fwd, NumericMatrix bck) {
  int N = fwd.ncol();
  if(N%2!=0) {
    stop("Dimension not even.");
  }

  NumericVector res(((N/2+1)*N/2)/2);
  int k=0;
  double tot;
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
      res[k++] = -0.125*log(tot);
    }
  }

  return(res);
}

// [[Rcpp::export]]
List Dedip_all(NumericMatrix fwd, NumericMatrix bck) {
  int N = fwd.ncol();
  if(N%2!=0) {
    stop("Dimension not even.");
  }

  NumericVector res_min(((N/2+1)*N/2)/2);
  NumericVector res_min2nd(((N/2+1)*N/2)/2);
  NumericVector res_dom(((N/2+1)*N/2)/2);
  NumericVector res_max(((N/2+1)*N/2)/2);
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
      quad[0] = (fwd(2*i  ,2*j)   * bck(2*i,  2*j))   * (fwd(2*j  ,2*i)   * bck(2*j  ,2*i));
      quad[1] = (fwd(2*i  ,2*j+1) * bck(2*i,  2*j+1)) * (fwd(2*j+1,2*i)   * bck(2*j+1,2*i));
      quad[2] = (fwd(2*i+1,2*j)   * bck(2*i+1,2*j))   * (fwd(2*j  ,2*i+1) * bck(2*j  ,2*i+1));
      quad[3] = (fwd(2*i+1,2*j+1) * bck(2*i+1,2*j+1)) * (fwd(2*j+1,2*i+1) * bck(2*j+1,2*i+1));
      res_min[k]    = -0.5*log(max(quad));
      res_min2nd[k] = -0.5*log(quad.sort()[2]);
      res_dom[k]    = -0.5*log(std::max(std::min(quad[0], quad[3]), std::min(quad[1], quad[2])));
      res_max[k]    = -0.5*log(min(quad));
      res_mean[k]   = -0.125*log((quad[0]*quad[1])*(quad[2]*quad[3]));
      k++;
    }
  }

  List res;
  res["min"] = res_min;
  res["min2nd"] = res_min2nd;
  res["dom"] = res_dom;
  res["max"] = res_max;
  res["mean"] = res_mean;
  return(res);
}
