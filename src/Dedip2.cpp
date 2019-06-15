#include <Rcpp.h>
using namespace Rcpp;

#include <algorithm>

// [[Rcpp::export]]
NumericVector Dedip2_min(NumericMatrix M) {
  int N = M.ncol();
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
      quad[0] = M(2*i, 2*j) + M(2*j, 2*i);
      quad[1] = M(2*i, 2*j+1) + M(2*j+1, 2*i);
      quad[2] = M(2*i+1, 2*j) + M(2*j, 2*i+1);
      quad[3] = M(2*i+1  ,2*j+1) + M(2*j+1, 2*i+1);
      res[k++] = 0.5*min(quad);
    }
  }

  return(res);
}

// [[Rcpp::export]]
NumericVector Dedip2_2nd_min(NumericMatrix M) {
  int N = M.ncol();
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
      quad[0] = M(2*i, 2*j) + M(2*j, 2*i);
      quad[1] = M(2*i, 2*j+1) + M(2*j+1, 2*i);
      quad[2] = M(2*i+1, 2*j) + M(2*j, 2*i+1);
      quad[3] = M(2*i+1  ,2*j+1) + M(2*j+1, 2*i+1);
      res[k++] = 0.5*quad.sort()[3];
    }
  }

  return(res);
}

// [[Rcpp::export]]
NumericVector Dedip2_dom(NumericMatrix M) {
  int N = M.ncol();
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
      quad[0] = M(2*i, 2*j) + M(2*j, 2*i);
      quad[1] = M(2*i, 2*j+1) + M(2*j+1, 2*i);
      quad[2] = M(2*i+1, 2*j) + M(2*j, 2*i+1);
      quad[3] = M(2*i+1  ,2*j+1) + M(2*j+1, 2*i+1);
      res[k++] = 0.5*std::min(std::max(quad[0], quad[3]), std::max(quad[1], quad[2]));
    }
  }

  return(res);
}



// [[Rcpp::export]]
List Dedip2_all(NumericMatrix M) {
  int N = M.ncol();
  if(N%2!=0) {
    stop("Dimension not even.");
  }

  NumericVector res_min(((N/2+1)*N/2)/2);
  NumericVector res_min2nd(((N/2+1)*N/2)/2);
  NumericVector res_dom(((N/2+1)*N/2)/2);

  int k=0;
  NumericVector quad(4);
  for(int i=0; i<N/2; i++) {
    for(int j=i; j<N/2; j++) {
      if(i==j) {
        res_min[k] = res_min2nd[k] = res_dom[k] = 0.0;
        k++;
        continue;
      }
      quad[0] = M(2*i, 2*j) + M(2*j, 2*i);
      quad[1] = M(2*i, 2*j+1) + M(2*j+1, 2*i);
      quad[2] = M(2*i+1, 2*j) + M(2*j, 2*i+1);
      quad[3] = M(2*i+1  ,2*j+1) + M(2*j+1, 2*i+1);
      res_min[k]    = 0.5*min(quad);
      res_min2nd[k] = 0.5*quad.sort()[3];
      res_dom[k]    = 0.5*std::min(std::max(quad[0], quad[3]), std::max(quad[1], quad[2]));
      k++;
    }
  }

  List res;
  res["min"] = res_min;
  res["min2nd"] = res_min2nd;
  res["dom"] = res_dom;
  return(res);
}