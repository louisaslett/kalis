// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

#include <fstream>
#include <algorithm>
// #include <bitset>

#include "Cache.h"
char *seq_data = NULL;
char **seq_ind = NULL;
int num_seqs = 0;
int seq_size = 0;

// [[Rcpp::export]]
int CacheAllSequences2(CharacterVector seqs, int bufsize) {
  num_seqs = seqs.size();

  seq_data = new char[bufsize];
  seq_ind = new char*[num_seqs];
  seq_ind[0] = seq_data;

  for(int i=0; i<num_seqs; i++) {
    std::ifstream input(seqs[i], std::ios::binary);
    input.read((char*)&seq_size, sizeof(seq_size));
    input.read(seq_ind[i], int(ceil(seq_size/8.0)));
    input.close();
    if(i<num_seqs-1) {
      seq_ind[i+1] = seq_ind[i] + int(ceil(seq_size/8.0));
    }
  }

  Rcout << "Cache loaded: " << num_seqs << " sequences of length " << seq_size << " (consuming " << (bufsize+num_seqs)/1073741824.0 << " GB RAM)" << std::endl;
  // std::bitset<8> x(seq_ind[3][4]);
  // Rcout << seqs[3];
  // Rcout << x;

  return(seq_size);
}

// [[Rcpp::export]]
RawVector QueryCache2(int idx) {
  RawVector res(int(ceil(seq_size/8.0)));

  if(idx > num_seqs) {
    return(res);
  }

  std::copy_n(seq_ind[idx], int(ceil(seq_size/8.0)), res.begin());

  return(res);
}

// [[Rcpp::export]]
void ClearSequenceCache2() {
  if(seq_data != NULL) {
    delete [] seq_data;
    seq_data = NULL;
  }
  if(seq_ind != NULL) {
    delete [] seq_ind;
    seq_ind = NULL;
  }
  num_seqs = 0;
  seq_size = 0;
}
