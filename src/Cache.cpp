// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

#include <fstream>
#include <algorithm>
// #include <bitset>
#include <stdlib.h>

#include <zlib.h>


#include "Cache.h"
uint32_t *hap_data = NULL;
uint32_t **hap_locus = NULL;
int32_t num_inds = 0;
int32_t hap_size = 0;

// [[Rcpp::export]]
int CacheHaplotypes_matrix_2(IntegerMatrix x, int N, int L, int transpose) {
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
  int_fast32_t ind = 0;
  for(int_fast32_t i=0; i<num_inds; i++) {
    for(int_fast32_t l=0; l<hap_size; l++) {
      if(transpose) {
        hap_locus[l][ind/32] ^= (-(x(i,l)) ^ hap_locus[l][ind/32]) & (1 << ind%32);
      } else {
        hap_locus[l][ind/32] ^= (-(x(l,i)) ^ hap_locus[l][ind/32]) & (1 << ind%32);
      }
    }
    ind++;
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
int CacheHaplotypes_hdf5_2(Function nexthaps, int N, int L) {
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

// Count number of columns in first line of the file NOT including EOF or new
//   line
// NB: When using to determine buffer size, must add 2 (one for new line and one
//   for NULL termination)
// NB: When using to determine number of haplotypes, add 1 and divide by 2.
int CacheHaplotypes_hapgz_ncols_2(gzFile fd) {
  int linelength = 0;

  int next_ch = gzgetc(fd);
  while(next_ch != -1 && (char) next_ch != '\n') {
    linelength++;
    next_ch = gzgetc(fd);
  }
  return(linelength);
}

// [[Rcpp::export]]
int CacheHaplotypes_hapgz_ncols(std::string file) {
  gzFile fd = gzopen(file.c_str(), "rb");
  if(fd == NULL) {
    Rcout << "Error in gzopen\n" << strerror(errno) << "\n";
  }

  int linelength = CacheHaplotypes_hapgz_ncols_2(fd);

  gzclose(fd);

  return(linelength);
}

// NB: When using to determine the number of loci from:
//   legend.gz => subtract 1
//   hap.gz    => use as is
// [[Rcpp::export]]
int CacheHaplotypes_hapgz_nlines(std::string file) {
  gzFile fd = gzopen(file.c_str(), "rb");
  if(fd == NULL) {
    Rcout << "Error in gzopen\n" << strerror(errno) << "\n";
  }

  // This is optimal buffer size for .hap.gz files.  It's
  // suboptimal for .legend.gz, but those are small anyway
  int bufsize = CacheHaplotypes_hapgz_ncols_2(fd) + 2;
  char buf[bufsize];

  int nlines = 1, c = 0;

  char *line;
  while((line = gzgets(fd, buf, bufsize)) != NULL) {
    while(buf[c] != '\0') {
      if(buf[c++] == '\n')
        nlines++;
    }
    c = 0;
  }

  gzclose(fd);

  return(nlines);
}

// NB: N and L arguments are the sizes *in the file* ... the extracted size is
//   inferred from loci_idx and hap_idx
// [[Rcpp::export]]
int CacheHaplotypes_hapgz_2(std::string file,
                            IntegerVector loci_idx,
                            IntegerVector hap_idx,
                            int L, int N) {
  num_inds = hap_idx.length();
  hap_size = loci_idx.length();

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

  gzFile fd = gzopen(file.c_str(), "rb");
  if(fd == NULL) {
    Rcout << "Error in gzopen\n" << strerror(errno) << "\n";
  }

  int bufsize = CacheHaplotypes_hapgz_ncols_2(fd) + 2;
  gzrewind(fd);
  char buf[bufsize];

  int nlines = 1;

  int_fast32_t next_l = loci_idx(0)-1;
  int next_ll = 0;
  int_fast32_t next_i = hap_idx(0)-1;
  int next_ii = 0;
  char x;

  char *line;
  for(int_fast32_t l=0; l<L; l++) {
    line = gzgets(fd, buf, bufsize);
    if(line == NULL) {
      Rcout << "Error: only reached line " << l+1 << " ... there are not " << L << " loci in the file!\n";
      // TODO: CLEAN UP AND EXIT
    }

    // Do we even want this line?  If not continue, otherwise figure out the
    //   line after
    if(l != next_l) {
      continue;
    }

    for(int_fast32_t i=0; i<N; i++) {
      // Do we even want this col?  If not continue, otherwise figure out the
      //   col after
      if(i != next_i) {
        // Check that there is at least valid data at this position ...
        if(*line != '0' && *line != '1') {
          Rcout << "Error: line " << l << " contains an invalid character!\n";
          // TODO: CLEAN UP AND EXIT
        }
        if((i < N-1 && *(line+1) != ' ') || (i == N-1 && *(line+1) != '\n')) {
          Rcout << "Error: line " << l << " does not contain a space (or EOL) after haplotype << " << i+1 << "!\n";
          // TODO: CLEAN UP AND EXIT
        }
        // ... then move on if still on track
        line += 2;
        continue;
      }

      // Process line
      x = *(line++); // should be a 0/1
      if(x == '\0') {
        Rcout << "Error: line " << l << " is of the incorrect length!\n";
        // TODO: CLEAN UP AND EXIT
      }
      if(x == '1') {
        hap_locus[next_ll][next_ii/32] ^= (-1 ^ hap_locus[next_ll][next_ii/32]) & (1 << next_ii%32);
      } else if(x == '0') {
        hap_locus[next_ll][next_ii/32] ^= (-0 ^ hap_locus[next_ll][next_ii/32]) & (1 << next_ii%32);
      } else {
        Rcout << "Error: line " << l << " contains an invalid character!\n";
        // TODO: CLEAN UP AND EXIT
      }

      x = *(line++); // should be a space following the 0/1 unless last hap on the line
      if(i < N-1 && x != ' ') {
        Rcout << "Error: line " << l << " does not contain a space after haplotype << " << i+1 << "!\n";
        // TODO: CLEAN UP AND EXIT
      }
      if(i == N-1 && x != '\n') {
        Rcout << "Error: line " << l << " does not end at the right place!\n";
        // TODO: CLEAN UP AND EXIT
      }

      next_ii++;
      if(next_ii >= num_inds)
        break;
      next_i = hap_idx(next_ii)-1;
    }
    // Reset haps we're reading ready for next loop
    next_i = hap_idx(0)-1;
    next_ii = 0;

    next_ll++;
    if(next_ll >= hap_size)
      break;
    next_l = loci_idx(next_ll)-1;
  }
  //finishparsing:;

  gzclose(fd);

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
