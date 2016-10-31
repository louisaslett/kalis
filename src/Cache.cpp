// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

#include <fstream>
#include <algorithm>
// #include <bitset>

#include "Cache.h"
uint32_t *seq_data = NULL;
uint32_t **seq_locus = NULL;
int32_t num_inds = 0;
int32_t seq_size = 0;

// [[Rcpp::export]]
int CacheAllSequences2(CharacterVector seqs, int bufsize) {
  num_inds = seqs.size();

  // Process first individual separately to learn sequence lengths
  std::ifstream input(seqs[0], std::ios::binary);
  input.read((char*)&seq_size, sizeof(seq_size));

  // Setup cache storage
  seq_data = new uint32_t[seq_size*int(ceil(num_inds/32.0))];
  seq_locus = new uint32_t*[seq_size];
  for(int l=0; l<seq_size; l++) {
    seq_locus[l] = seq_data + l*int(ceil(num_inds/32.0));
  }
  uint32_t seq_tmp[int(ceil(seq_size/32.0))];

  // Finish processing first individual separately
  input.read((char*) seq_tmp, int(ceil(seq_size/8.0)));
  for(int_fast32_t l=0; l<seq_size; l++) {
    seq_locus[l][0] ^= (-((seq_tmp[l/32] >> l%32) & 1) ^ seq_locus[l][0]) & (1 << 0);
  }
  input.close();

  // Process remainder
  int32_t this_seq_size;
  for(int_fast32_t i=1; i<num_inds; i++) {
    std::ifstream input(seqs[i], std::ios::binary);
    input.read((char*)&this_seq_size, sizeof(this_seq_size));
    if(this_seq_size != seq_size) {
      Rcout << "Error: individual " << seqs[i] << " has differing sequence length (" << this_seq_size << ") --- skipping\n";
    } else {
      input.read((char*) seq_tmp, int(ceil(seq_size/8.0)));
      for(int_fast32_t l=0; l<seq_size; l++) {
        seq_locus[l][i/32] ^= (-((seq_tmp[l/32] >> l%32) & 1) ^ seq_locus[l][i/32]) & (1 << i%32);
      }
    }
    input.close();
  }

  Rcout << "Cache loaded: " << num_inds << " sequences of length " << seq_size << " (consuming " << (seq_size*int(ceil(num_inds/32.0))*sizeof(uint32_t))/1073741824.0 << " GB RAM)" << std::endl;
  // std::bitset<8> x(seq_ind[3][4]);
  // Rcout << seqs[3];
  // Rcout << x;

  return(seq_size);
}

// [[Rcpp::export]]
IntegerVector QueryCache2(int idx) {
  IntegerVector res(int(ceil(seq_size/32.0)));

  if(idx > num_inds) {
    return(res);
  }

  // Extract individual from locus oriented layout
  uint32_t seq_tmp[int(ceil(seq_size/32.0))];
  for(int_fast32_t l=0; l<seq_size; l++) {
    seq_tmp[l/32] ^= (-((seq_locus[l][idx/32] >> idx%32) & 1) ^ seq_tmp[l/32]) & (1 << l%32);
  }

  // Copy and return
  std::copy_n(seq_tmp, int(ceil(seq_size/32.0)), res.begin());

  return(res);
}

// [[Rcpp::export]]
void ClearSequenceCache2() {
  if(seq_data != NULL) {
    delete[] seq_data;
    seq_data = NULL;
  }
  if(seq_locus != NULL) {
    delete[] seq_locus;
    seq_locus = NULL;
  }
  num_inds = 0;
  seq_size = 0;
}
