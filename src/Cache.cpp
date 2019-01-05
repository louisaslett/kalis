// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

#include <fstream>
#include <algorithm>
// #include <bitset>
#include <stdlib.h>

#include "Cache.h"
uint32_t *hap_data = NULL;
uint32_t **hap_locus = NULL;
int32_t num_inds = 0;
int32_t hap_size = 0;

// [[Rcpp::export]]
int CacheAllHaplotypes2(CharacterVector haps, int bufsize) {
  num_inds = haps.size();

  // Process first individual separately to learn haplotype lengths
  std::ifstream input(haps[0], std::ios::binary);
  input.read((char*)&hap_size, sizeof(hap_size));

  // Setup cache storage
  //hap_data = new uint32_t[hap_size*int(ceil(num_inds/32.0))];
  // NB need 32-byte aligned for AVX.  Note also that *each* new locus must be on a 32-byte boundary!
  //   int(ceil(num_inds/32.0)) ints are required per locus to store all individuals
  //   int(ceil((num_inds/32.0)/8.0)) __m256i's are required per locus to store all individuals
  if(posix_memalign((void**) &hap_data, 32, hap_size*int(ceil((num_inds/32.0)/8.0))*8*sizeof(uint32_t)) != 0) {
    Rcout << "Unable to assign enough aligned memory for cache storage\n";
    return(0);
  }
  hap_locus = new uint32_t*[hap_size];
  for(int l=0; l<hap_size; l++) {
    hap_locus[l] = hap_data + l*int(ceil((num_inds/32.0)/8.0))*8;
  }
  uint32_t hap_tmp[int(ceil(hap_size/32.0))];

  // Finish processing first individual separately
  input.read((char*) hap_tmp, int(ceil(hap_size/8.0)));
  for(int_fast32_t l=0; l<hap_size; l++) {
    hap_locus[l][0] ^= (-((hap_tmp[l/32] >> l%32) & 1) ^ hap_locus[l][0]) & (1 << 0);
  }
  input.close();

  // Process remainder
  int32_t this_hap_size;
  for(int_fast32_t i=1; i<num_inds; i++) {
    std::ifstream input(haps[i], std::ios::binary);
    input.read((char*)&this_hap_size, sizeof(this_hap_size));
    if(this_hap_size != hap_size) {
      Rcout << "Error: individual " << haps[i] << " has differing haplotype length (" << this_hap_size << ") --- skipping\n";
    } else {
      input.read((char*) hap_tmp, int(ceil(hap_size/8.0)));
      for(int_fast32_t l=0; l<hap_size; l++) {
        hap_locus[l][i/32] ^= (-((hap_tmp[l/32] >> l%32) & 1) ^ hap_locus[l][i/32]) & (1 << i%32);
      }
    }
    input.close();
  }

  Rcout << "Cache loaded: " << num_inds << " haplotypes of length " << hap_size << " (consuming " << (hap_size*int(ceil((num_inds/32.0)/8.0))*8*sizeof(uint32_t))/1073741824.0 << " GB RAM)" << std::endl;
  // std::bitset<8> x(hap_ind[3][4]);
  // Rcout << haps[3];
  // Rcout << x;

  return(hap_size);
}

// [[Rcpp::export]]
int CacheAllHaplotypesH52(Function nexthaps, int N, int L) {
  num_inds = N;
  hap_size = L;

  // Setup cache storage
  // NB need 32-byte aligned for AVX.  Note also that *each* new locus must be on a 32-byte boundary!
  //   int(ceil(num_inds/32.0)) ints are required per locus to store all individuals
  //   int(ceil((num_inds/32.0)/8.0)) __m256i's are required per locus to store all individuals
  if(posix_memalign((void**) &hap_data, 32, hap_size*int(ceil((num_inds/32.0)/8.0))*8*sizeof(uint32_t)) != 0) {
    Rcout << "Unable to assign enough aligned memory for cache storage\n";
    return(0);
  }
  hap_locus = new uint32_t*[hap_size];
  for(int l=0; l<hap_size; l++) {
    hap_locus[l] = hap_data + l*int(ceil((num_inds/32.0)/8.0))*8;
  }

  // Process haplotypes in chunks from input function
  IntegerMatrix nexthapmat;
  nexthapmat = nexthaps();
  int_fast32_t ind = 0;
  while(nexthapmat.nrow() > 0) {
    for(int_fast32_t i=0; i<nexthapmat.ncol(); i++) {
      for(int_fast32_t l=0; l<hap_size; l++) {
        hap_locus[l][ind/32] ^= (-(nexthapmat(l,i)) ^ hap_locus[l][ind/32]) & (1 << ind%32);
      }
      ind++;
    }
    nexthapmat = nexthaps();
  }

  Rcpp::Function msg("message");
  msg(std::string("Cache loaded: ") +
    std::to_string(num_inds) +
    std::string(" haplotypes of length ") +
    std::to_string(hap_size) +
    std::string(" (consuming ") +
    std::to_string((hap_size*int(ceil((num_inds/32.0)/8.0))*8*sizeof(uint32_t))/1073741824.0) +
    std::string(" GB RAM)"));
  // std::bitset<8> x(hap_ind[3][4]);
  // Rcout << haps[3];
  // Rcout << x;

  return(hap_size);
}

// [[Rcpp::export]]
IntegerVector QueryCache2_ind(int idx) {
  IntegerVector res(int(ceil(hap_size/32.0)));

  if(idx > num_inds) {
    return(res);
  }

  // Extract individual from locus oriented layout
  uint32_t hap_tmp[int(ceil(hap_size/32.0))];
  for(int_fast32_t l=0; l<hap_size; l++) {
    hap_tmp[l/32] ^= (-((hap_locus[l][idx/32] >> idx%32) & 1) ^ hap_tmp[l/32]) & (1 << l%32);
  }

  // Copy and return
  std::copy_n(hap_tmp, int(ceil(hap_size/32.0)), res.begin());

  return(res);
}

// [[Rcpp::export]]
IntegerVector QueryCache2_loc(int idx) {
  IntegerVector res(int(ceil(num_inds/32.0)));

  if(idx > hap_size) {
    return(res);
  }

  // Extract locus from locus oriented layout
  uint32_t hap_tmp[int(ceil(num_inds/32.0))];
  for(int_fast32_t i=0; i<ceil(num_inds/32.0); i++) {
    hap_tmp[i] = hap_locus[idx][i];
  }

  // Copy and return
  std::copy_n(hap_tmp, int(ceil(num_inds/32.0)), res.begin());

  return(res);
}

// [[Rcpp::export]]
void ClearHaplotypeCache2() {
  if(hap_data != NULL) {
    free(hap_data);
    hap_data = NULL;
  }
  if(hap_locus != NULL) {
    delete[] hap_locus;
    hap_locus = NULL;
  }
  num_inds = 0;
  hap_size = 0;
}
