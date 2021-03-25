#include "R_Cache.h"

#include <stdlib.h> // posix_memalign (C11)
#include <zlib.h>

#include "StencilVec.h"
#include "R_Kalis.h"
#include "Cache.h"
uint32_t *hap_data = NULL;
uint32_t **hap_locus = NULL;
size_t num_inds = 0;
size_t hap_size = 0;



SEXP ClearHaplotypeCache2() {
  if(hap_data != NULL) {
    free(hap_data);
    hap_data = NULL;
  }
  if(hap_locus != NULL) {
    Free(hap_locus);
    hap_locus = NULL;
  }
  num_inds = 0;
  hap_size = 0;

  KALIS_RETURN;
}



SEXP CacheHaplotypes_matrix_2(SEXP Rx, SEXP RN, SEXP RL, SEXP Rtranspose) {
  num_inds = (size_t) Rf_asInteger(RN);
  hap_size = (size_t) Rf_asInteger(RL);
  int *x = INTEGER(Rx);
  int transpose = Rf_asLogical(Rtranspose);
  if((!transpose && ((size_t) Rf_nrows(Rx) != hap_size || (size_t) Rf_ncols(Rx) != num_inds)) ||
      (transpose && ((size_t) Rf_nrows(Rx) != num_inds || (size_t) Rf_ncols(Rx) != hap_size))) {
    REprintf("Error: Invalid size input matrix.\n");
    ClearHaplotypeCache2();
    KALIS_RETURN;
  }

  // Setup cache storage
  // NB need for example 32-byte aligned for AVX.  Note also that *each* new locus must be on a 32-byte boundary!
  //   (size_t) ceil(num_inds/32.0) ints are required per locus to store all individuals
  //   (size_t) ceil((num_inds/32.0)/8.0) __m256i's are required per locus to store all individuals
  // So, in general,
  //   (size_t) ceil((((double) num_inds)/32.0)/((double) KALIS_INTVEC_SIZE)) of the integer vector units are required
  //   (size_t) ceil((((double) num_inds)/32.0)/((double) KALIS_INTVEC_SIZE))*KALIS_INTVEC_SIZE ints worth of space in total
  // The size is 4 bytes times #ints, ie align to 4*KALIS_INTVEC_SIZE
  // aligned_alloc is newer, but posix_memalign supported on older platforms (and older MacOS)
  size_t alignment = 4*KALIS_INTVEC_SIZE;
  while(alignment < sizeof(void*)) { // POSIX alignment must be at least sizeof(void*)
    alignment *= 2;
  }
  //REprintf("Alignment: %d\n", alignment);
  int err = posix_memalign((void**) &hap_data, alignment, hap_size*((size_t) ceil((((double) num_inds)/32.0)/((double) KALIS_INTVEC_SIZE))*KALIS_INTVEC_SIZE)*sizeof(uint32_t));
  if(err != 0) {
    REprintf("Error: Unable to assign aligned memory for cache storage (err # %d)\n", err);
    ClearHaplotypeCache2();
    KALIS_RETURN;
  }
  hap_locus = (uint32_t**) R_Calloc(hap_size, uint32_t*);
  for(size_t l=0; l<hap_size; l++) {
    hap_locus[l] = hap_data + l*((size_t) (ceil((num_inds/32.0)/8.0)))*8;
  }

  // Process haplotypes in chunks from input function
  size_t ind = 0;
  for(size_t i=0; i<num_inds; i++) {
    for(size_t l=0; l<hap_size; l++) {
      if(transpose) {
        hap_locus[l][ind/32] ^= (-(x[i+l*num_inds]) ^ hap_locus[l][ind/32]) & (1 << ind%32);
      } else {
        hap_locus[l][ind/32] ^= (-(x[l+i*hap_size]) ^ hap_locus[l][ind/32]) & (1 << ind%32);
      }
    }
    ind++;
  }

  KALIS_RETURN;
}



SEXP CacheHaplotypes_hdf5_2(SEXP Rnexthaps, SEXP Rnexthapsenv, SEXP RN, SEXP RL) {
  if(!Rf_isLanguage(Rnexthaps)) {
    REprintf("Error: 'Rnexthaps' must be a expression evaluating a closure.");
    KALIS_RETURN;
  }
  if(!Rf_isEnvironment(Rnexthapsenv)) {
    REprintf("Error: 'Rnexthapsenv' should be an environment");
    KALIS_RETURN;
  }
  num_inds = (size_t) Rf_asInteger(RN);
  hap_size = (size_t) Rf_asInteger(RL);

  // Setup cache storage
  // NB need for example 32-byte aligned for AVX.  Note also that *each* new locus must be on a 32-byte boundary!
  //   (size_t) ceil(num_inds/32.0) ints are required per locus to store all individuals
  //   (size_t) ceil((num_inds/32.0)/8.0) __m256i's are required per locus to store all individuals
  // So, in general,
  //   (size_t) ceil((((double) num_inds)/32.0)/((double) KALIS_INTVEC_SIZE)) of the integer vector units are required
  //   (size_t) ceil((((double) num_inds)/32.0)/((double) KALIS_INTVEC_SIZE))*KALIS_INTVEC_SIZE ints worth of space in total
  // The size is 4 bytes times #ints, ie align to 4*KALIS_INTVEC_SIZE
  // aligned_alloc is newer, but posix_memalign supported on older platforms (and older MacOS)
  size_t alignment = 4*KALIS_INTVEC_SIZE;
  while(alignment < sizeof(void*)) { // POSIX alignment must be at least sizeof(void*)
    alignment *= 2;
  }
  //REprintf("Alignment: %d\n", alignment);
  int err = posix_memalign((void**) &hap_data, alignment, hap_size*((size_t) ceil((((double) num_inds)/32.0)/((double) KALIS_INTVEC_SIZE))*KALIS_INTVEC_SIZE)*sizeof(uint32_t));
  if(err != 0) {
    REprintf("Error: Unable to assign aligned memory for cache storage (err # %d)\n", err);
    ClearHaplotypeCache2();
    KALIS_RETURN;
  }
  hap_locus = (uint32_t**) R_Calloc(hap_size, uint32_t*);
  for(size_t l=0; l<hap_size; l++) {
    hap_locus[l] = hap_data + l*((size_t) (ceil((num_inds/32.0)/8.0)))*8;
  }

  // Process haplotypes in chunks from input function
  SEXP Rnexthapmat;
  Rnexthapmat = PROTECT(Rf_eval(Rnexthaps, Rnexthapsenv));
  size_t ind = 0, nexthapmatncol;
  int *nexthapmatptr;
  while(Rf_nrows(Rnexthapmat) > 0) {
    nexthapmatptr = INTEGER(Rnexthapmat);
    nexthapmatncol = (size_t) Rf_ncols(Rnexthapmat);
    if((size_t) Rf_nrows(Rnexthapmat) != hap_size) {
      REprintf("Error: Incorrectly size matrix provided during haplotype loading\n");
      ClearHaplotypeCache2();
      UNPROTECT(1);
      KALIS_RETURN;
    }
    for(size_t i=0; i<nexthapmatncol; i++) {
      for(size_t l=0; l<hap_size; l++) {
        hap_locus[l][ind/32] ^= (-(nexthapmatptr[l+i*hap_size]) ^ hap_locus[l][ind/32]) & (1 << ind%32);
      }
      ind++;
    }
    UNPROTECT(1);
    Rnexthapmat = PROTECT(Rf_eval(Rnexthaps, Rnexthapsenv));
  }

  UNPROTECT(1);
  KALIS_RETURN;
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

SEXP CacheHaplotypes_hapgz_ncols(SEXP Rfile) {
  const char *file = CHAR(STRING_ELT(Rfile, 0));
  gzFile fd = gzopen(file, "rb");
  if(fd == NULL) {
    int err = 0;
    REprintf("Error: Cannot open gzip file (err #: %d) ... %s.\n", err, gzerror(fd, &err));
    KALIS_RETURN;
  }

  int linelength = CacheHaplotypes_hapgz_ncols_2(fd);

  gzclose(fd);

  return(Rf_ScalarInteger(linelength));
}

// NB: When using to determine the number of loci from:
//   legend.gz => subtract 1
//   hap.gz    => use as is
SEXP CacheHaplotypes_hapgz_nlines(SEXP Rfile) {
  const char *file = CHAR(STRING_ELT(Rfile, 0));
  gzFile fd = gzopen(file, "rb");
  if(fd == NULL) {
    int err = 0;
    REprintf("Error: Cannot open gzip file (err #: %d) ... %s.\n", err, gzerror(fd, &err));
    KALIS_RETURN;
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

  return(Rf_ScalarInteger(nlines));
}

// NB: N and L arguments are the sizes *in the file* ... the extracted size is
//   inferred from loci_idx and hap_idx
SEXP CacheHaplotypes_hapgz_2(SEXP Rfile, SEXP Rloci_idx, SEXP Rhap_idx, SEXP RL, SEXP RN) {
  const char *file = CHAR(STRING_ELT(Rfile, 0));
  int *loci_idx, *hap_idx;
  loci_idx = INTEGER(Rloci_idx);
  hap_idx = INTEGER(Rhap_idx);
  num_inds = (size_t) LENGTH(Rhap_idx);
  size_t num_inds_file = (size_t) Rf_asInteger(RN);
  hap_size = (size_t) LENGTH(Rloci_idx);
  size_t hap_size_file = (size_t) Rf_asInteger(RL);
  if(hap_size > hap_size_file) {
    REprintf("Mismatch between RL and length of Rloci_idx.\n");
    KALIS_RETURN;
  }
  if(num_inds > num_inds_file) {
    REprintf("Mismatch between RN and length of Rhap_idx.\n");
    KALIS_RETURN;
  }

  gzFile fd = gzopen(file, "rb");
  if(fd == NULL) {
    int err = 0;
    REprintf("Error: Cannot open gzip file (err #: %d) ... %s.\n", err, gzerror(fd, &err));
    KALIS_RETURN;
  }

  // Setup cache storage
  // NB need for example 32-byte aligned for AVX.  Note also that *each* new locus must be on a 32-byte boundary!
  //   (size_t) ceil(num_inds/32.0) ints are required per locus to store all individuals
  //   (size_t) ceil((num_inds/32.0)/8.0) __m256i's are required per locus to store all individuals
  // So, in general,
  //   (size_t) ceil((((double) num_inds)/32.0)/((double) KALIS_INTVEC_SIZE)) of the integer vector units are required
  //   (size_t) ceil((((double) num_inds)/32.0)/((double) KALIS_INTVEC_SIZE))*KALIS_INTVEC_SIZE ints worth of space in total
  // The size is 4 bytes times #ints, ie align to 4*KALIS_INTVEC_SIZE
  // aligned_alloc is newer, but posix_memalign supported on older platforms (and older MacOS)
  size_t alignment = 4*KALIS_INTVEC_SIZE;
  while(alignment < sizeof(void*)) { // POSIX alignment must be at least sizeof(void*)
    alignment *= 2;
  }
  //REprintf("Alignment: %d\n", alignment);
  int err = posix_memalign((void**) &hap_data, alignment, hap_size*((size_t) ceil((((double) num_inds)/32.0)/((double) KALIS_INTVEC_SIZE))*KALIS_INTVEC_SIZE)*sizeof(uint32_t));
  if(err != 0) {
    REprintf("Error: Unable to assign aligned memory for cache storage (err # %d)\n", err);
    ClearHaplotypeCache2();
    gzclose(fd);
    KALIS_RETURN;
  }
  hap_locus = (uint32_t**) R_Calloc(hap_size, uint32_t*);
  for(size_t l=0; l<hap_size; l++) {
    hap_locus[l] = hap_data + l*((size_t) (ceil((num_inds/32.0)/8.0)))*8;
  }

  int bufsize = CacheHaplotypes_hapgz_ncols_2(fd) + 2;
  gzrewind(fd);
  char buf[bufsize];

  size_t next_l = (size_t) loci_idx[0]-1;
  size_t next_ll = 0;
  size_t next_i = (size_t) hap_idx[0]-1;
  size_t next_ii = 0;
  char x;

  char *line;
  for(size_t l=0; l<hap_size_file; l++) {
    line = gzgets(fd, buf, bufsize);
    if(line == NULL) {
      REprintf("Error: only reached line %d ... there are not %d loci in the file!\n", l+1, hap_size_file);
      ClearHaplotypeCache2();
      gzclose(fd);
      KALIS_RETURN;
    }

    // Do we even want this line?  If not continue, otherwise figure out the
    //   line after
    if(l != next_l) {
      continue;
    }

    for(size_t i=0; i<num_inds_file; i++) {
      // Do we even want this col?  If not continue, otherwise figure out the
      //   col after
      if(i != next_i) {
        // Check that there is at least valid data at this position ...
        if(*line != '0' && *line != '1') {
          REprintf("Error: line %d contains an invalid character!\n", l);
          ClearHaplotypeCache2();
          gzclose(fd);
          KALIS_RETURN;
        }
        if((i < num_inds_file-1 && *(line+1) != ' ') || (i == num_inds_file-1 && *(line+1) != '\n')) {
          REprintf("Error: line %d does not contain a space (or EOL) after haplotype << %d!\n", l, i+1);
          ClearHaplotypeCache2();
          gzclose(fd);
          KALIS_RETURN;
        }
        // ... then move on if still on track
        line += 2;
        continue;
      }

      // Process line
      x = *(line++); // should be a 0/1
      if(x == '\0') {
        REprintf("Error: line %d is of the incorrect length!\n", l);
        ClearHaplotypeCache2();
        gzclose(fd);
        KALIS_RETURN;
      }
      if(x == '1') {
        hap_locus[next_ll][next_ii/32] ^= (-1 ^ hap_locus[next_ll][next_ii/32]) & (1 << next_ii%32);
      } else if(x == '0') {
        hap_locus[next_ll][next_ii/32] ^= (-0 ^ hap_locus[next_ll][next_ii/32]) & (1 << next_ii%32);
      } else {
        REprintf("Error: line %d contains an invalid character!\n", l);
        ClearHaplotypeCache2();
        gzclose(fd);
        KALIS_RETURN;
      }

      x = *(line++); // should be a space following the 0/1 unless last hap on the line
      if(i < num_inds_file-1 && x != ' ') {
        REprintf("Error: line %d does not contain a space after haplotype << %d!\n", l, i+1);
        ClearHaplotypeCache2();
        gzclose(fd);
        KALIS_RETURN;
      }
      if(i == num_inds_file-1 && x != '\n') {
        REprintf("Error: line %d does not end at the right place!\n", l);
        ClearHaplotypeCache2();
        gzclose(fd);
        KALIS_RETURN;
      }

      next_ii++;
      if(next_ii >= num_inds)
        break;
      next_i = (size_t) hap_idx[next_ii]-1;
    }
    // Reset haps we're reading ready for next loop
    next_i = (size_t) hap_idx[0]-1;
    next_ii = 0;

    next_ll++;
    if(next_ll >= hap_size)
      break;
    next_l = (size_t) loci_idx[next_ll]-1;
  }

  gzclose(fd);

  KALIS_RETURN;
}



SEXP QueryCache2_ind(SEXP Ridx) {
  size_t idx = (size_t) Rf_asInteger(Ridx);
  SEXP res = PROTECT(Rf_allocVector(INTSXP, (int) ceil(hap_size/32.0)));

  if(idx > num_inds) {
    REprintf("Error: request for haplotype index exceeding number of loaded haplotypes.\n");
    KALIS_RETURN;
  }

  // Extract individual from locus oriented layout
  int *hap_tmp = INTEGER(res);
  for(size_t l=0; l<hap_size; l++) {
    hap_tmp[l/32] ^= (-((hap_locus[l][idx/32] >> idx%32) & 1) ^ hap_tmp[l/32]) & (1 << l%32);
  }

  UNPROTECT(1);
  return(res);
}

SEXP QueryCache2_loc(SEXP Ridx) {
  size_t idx = (size_t) Rf_asInteger(Ridx);
  SEXP res = PROTECT(Rf_allocVector(INTSXP, (int) ceil(num_inds/32.0)));

  if(idx > hap_size) {
    REprintf("Error: request for locus index exceeding length of loaded haplotypes.\n");
    KALIS_RETURN;
  }

  // Extract locus from locus oriented layout
  int *hap_tmp = INTEGER(res);
  for(size_t i=0; i<ceil(num_inds/32.0); i++) {
    hap_tmp[i] = hap_locus[idx][i];
  }

  UNPROTECT(1);
  return(res);
}
