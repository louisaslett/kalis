#include "R_CladeMat.h"

#define _GNU_SOURCE
#include <pthread.h>
#include "R_Kalis.h"
#include "Cache.h"

#define min(X, Y)  ((X) < (Y) ? (X) : (Y))

typedef struct blob {
  struct blob* blob;
  size_t num_in_blob;
  double lower;
  double upper;
  struct blob* next;
  struct blob* prev;
  struct blob* morepop;
  struct blob* lesspop;
  double c3;
} blob;


void printblob(blob* bobtheblob) {
  Rprintf("Blob %p (%p): # = %d, range = (%lf, %lf), prv = %p, nxt = %p, lp = %p, mp = %p \n",
          bobtheblob,
          (bobtheblob==bobtheblob->blob?NULL:bobtheblob->blob),
          (int) bobtheblob->num_in_blob,
          bobtheblob->lower,
          bobtheblob->upper,
          bobtheblob->prev,
          bobtheblob->next,
          bobtheblob->lesspop,
          bobtheblob->morepop);
}

void blobby_pop_contest(blob* cur, blob** headpop, blob** tailpop) {
  if(cur->morepop != NULL && cur->morepop->num_in_blob < cur->num_in_blob) {
    blob *other;
    other = cur->morepop;

    if(other->morepop != NULL) {
      other->morepop->lesspop = cur;
    }
    cur->morepop = other->morepop;
    other->morepop = cur;
    if(cur->morepop == NULL) {
      *headpop = cur;
    }

    if(cur->lesspop != NULL) {
      cur->lesspop->morepop = other;
    }
    other->lesspop = cur->lesspop;
    cur->lesspop = other;
    if(other->lesspop == NULL) {
      *tailpop = other;
    }

    blobby_pop_contest(cur, headpop, tailpop);
  } else if(cur->lesspop != NULL && cur->lesspop->num_in_blob > cur->num_in_blob) {
    blob *other;
    other = cur->lesspop;

    if(cur->morepop != NULL) {
      cur->morepop->lesspop = other;
    }
    other->morepop = cur->morepop;
    cur->morepop = other;
    if(other->morepop == NULL) {
      *headpop = other;
    }

    if(other->lesspop != NULL) {
      other->lesspop->morepop = cur;
    }
    cur->lesspop = other->lesspop;
    other->lesspop = cur;
    if(cur->lesspop == NULL) {
      *tailpop = cur;
    }

    blobby_pop_contest(cur, headpop, tailpop);
  }
}


blob* hunttheblob(blob* blob) {
  if(blob->blob == blob)
    return(blob);
  return(hunttheblob(blob->blob));
}

double alphabetascaling(double x, double z0) {
  x = x * z0;
  if(x == 0.0) {
    return(744.4400719213812180897);
  } else {
    return(-log(x));
  }
}

blob* blobby_BB(const double* alpha, const double* beta, const size_t recipient, size_t n, blob* blobs, blob*** x_in_blob, double* n_clade, const double thres, const double maxd, const double unitdist, const int max1var) {
  double z0 = 0.0;
  for(size_t i = 0; i<n; i++) {
    z0 += alpha[i] * beta[i];
  }
  z0 = 1/z0;

  const double* a_ptr = alpha+1;
  const double* b_ptr = beta+1;

  double x, xnext;

  blob *head, *tail, *headpop, *tailpop;
  head = tail = headpop = tailpop = blobs;

  blobs[0].blob = blobs;
  blobs[0].num_in_blob = 1;
  if(recipient==0)
    x = 0.0;
  else
    x = alphabetascaling(*alpha * *beta, z0);
  if(recipient==1)
    xnext = 0.0;
  else
    xnext = alphabetascaling(*a_ptr * *b_ptr, z0);
  blobs[0].lower = x - thres;
  blobs[0].upper = x + thres;
  blobs[0].next = NULL;
  blobs[0].prev = NULL;
  blobs[0].morepop = NULL;
  blobs[0].lesspop = NULL;

  x_in_blob[0] = &(blobs[0].blob);

  // What element of blob array is the next new one?
  size_t next_new_blob = 1;

  const size_t recipientm1 = recipient - 1;
  blob* cur;

  //Rprintf("Adding obs %lf\n", x);

  for(size_t i = 1; i < n; i++) {
    x = xnext;
    if(i<n-1){
      if(recipientm1==i)
      {
        xnext = 0.0;
      }
      else
      {
        xnext = alphabetascaling(*(a_ptr+1) * *(b_ptr+1), z0);
      }
    }

    cur = headpop;

    // DETAILED DEBUG
    // Rprintf("COMPUTED BLOBS (hp = %p, tp = %p): \n", headpop, tailpop);
    // blob* prt = head;
    // while(prt != NULL) {
    //   printblob(prt);
    //   prt = prt->next;
    // }
    // Rprintf("Adding obs: %d = %lf\n", i, x);

    //Rprintf("Adding obs %lf\n", x);

    while(cur != NULL) {


      // this observation is between blobs
      if(cur->prev != NULL && x < cur->lower && x >= cur->prev->upper) { // do we want to check to the right here or just leave to that pointer iteration?
        blobs[next_new_blob].blob = blobs + next_new_blob;
        blobs[next_new_blob].num_in_blob = 1;
        blobs[next_new_blob].lower = x - thres;
        blobs[next_new_blob].upper = x + thres;
        blobs[next_new_blob].next = cur;
        blobs[next_new_blob].prev = cur->prev;
        blobs[next_new_blob].morepop = tailpop;
        blobs[next_new_blob].lesspop = NULL;
        cur->prev->next = blobs + next_new_blob;
        cur->prev = blobs + next_new_blob;
        tailpop->lesspop = blobs + next_new_blob;
        tailpop = blobs + next_new_blob;
        x_in_blob[i] = &(blobs[next_new_blob].blob);
        next_new_blob++;
        a_ptr++;
        b_ptr++;
        break;
      }

      // this observation is between and causes blobs to merge (to the right)
      if(cur->next != NULL && x >= cur->next->lower && x < cur->upper) {
        cur->next->blob = cur->blob;
        cur->num_in_blob += cur->next->num_in_blob + 1;
        if(tail == cur->next) {
          tail = cur->blob;
        } else {
          cur->next->next->prev = cur;
        }
        if(cur->next == headpop) {
          headpop = cur->next->lesspop;
        } else {
          cur->next->morepop->lesspop = cur->next->lesspop;
        }
        if(cur->next == tailpop) {
          tailpop = cur->next->morepop;
        } else {
          cur->next->lesspop->morepop = cur->next->morepop;
        }
        cur->upper = cur->next->upper;
        cur->next = cur->next->next;
        blobby_pop_contest(cur, &headpop, &tailpop);
        x_in_blob[i] = &(cur->blob);
        a_ptr++;
        b_ptr++;
        break;
      }

      // this observation is between and causes blobs to merge (to the left)
      if(cur->prev != NULL && x >= cur->lower && x < cur->prev->upper) {
        cur->prev->blob = cur->blob;
        cur->num_in_blob += cur->prev->num_in_blob + 1;
        if(head == cur->prev) {
          head = cur->blob;
        } else {
          cur->prev->prev->next = cur;
        }
        if(cur->prev == headpop) {
          headpop = cur->prev->lesspop;
        } else {
          cur->prev->morepop->lesspop = cur->prev->lesspop;
        }
        if(cur->prev == tailpop) {
          tailpop = cur->prev->morepop;
        } else {
          cur->prev->lesspop->morepop = cur->prev->morepop;
        }
        cur->lower = cur->prev->lower;
        cur->prev = cur->prev->prev;
        blobby_pop_contest(cur, &headpop, &tailpop);
        x_in_blob[i] = &(cur->blob);
        a_ptr++;
        b_ptr++;
        break;
      }

      // this observation is in this blob
      if(x >= cur->lower && x < cur->upper) {
        (cur->num_in_blob)++;
        if(x - thres < cur->lower) {
          cur->lower = x - thres;
        } else if(x + thres > cur->upper) {
          cur->upper = x + thres;
        }
        blobby_pop_contest(cur, &headpop, &tailpop);
        x_in_blob[i] = &(cur->blob);
        a_ptr++;
        b_ptr++;
        break;
      }

      cur = cur->lesspop;
    }
    // this observation is smaller than any previously seen blob
    if(cur == NULL && x < head->lower) {
      blobs[next_new_blob].blob = blobs + next_new_blob;
      blobs[next_new_blob].num_in_blob = 1;
      blobs[next_new_blob].lower = x - thres;
      blobs[next_new_blob].upper = x + thres;
      blobs[next_new_blob].next = head;
      blobs[next_new_blob].prev = NULL;
      blobs[next_new_blob].morepop = tailpop;
      blobs[next_new_blob].lesspop = NULL;
      head->prev = blobs + next_new_blob;
      head = blobs + next_new_blob;
      tailpop->lesspop = blobs + next_new_blob;
      tailpop = blobs + next_new_blob;
      x_in_blob[i] = &(blobs[next_new_blob].blob);
      next_new_blob++;
      a_ptr++;
      b_ptr++;
      continue;
    }
    // this observation is larger than any previously seen blob
    if(cur == NULL && x >= tail->upper) {
      // this observation is beyond the end of all blobs
      blobs[next_new_blob].blob = blobs + next_new_blob;
      blobs[next_new_blob].num_in_blob = 1;
      blobs[next_new_blob].lower = x - thres;
      blobs[next_new_blob].upper = x + thres;
      blobs[next_new_blob].next = NULL;
      blobs[next_new_blob].prev = tail;
      blobs[next_new_blob].morepop = tailpop;
      blobs[next_new_blob].lesspop = NULL;
      tail->next = blobs + next_new_blob;
      tail = blobs + next_new_blob;
      tailpop->lesspop = blobs + next_new_blob;
      tailpop = blobs + next_new_blob;
      x_in_blob[i] = &(blobs[next_new_blob].blob);
      next_new_blob++;
      a_ptr++;
      b_ptr++;
      continue;
    }
  }

  // Chase down pointers
  // DEBUG
  // Rprintf("ALL USED BLOBS (pre-parse):\n\n");
  // for(size_t i = 0; i<next_new_blob; i++) {
  //   printblob(blobs+i);
  // }

  for(size_t i = 0; i<next_new_blob; i++) {
    if(blobs[i].blob == blobs+i || (blobs[i].blob)->blob == blobs[i].blob)
      continue;

    blobs[i].blob = hunttheblob(blobs[i].blob);
  }

  // DEBUG
  // Rprintf("ALL USED BLOBS:\n\n");
  // for(size_t i = 0; i<next_new_blob; i++) {
  //   printblob(blobs+i);
  // }

  // Col 3 calc
  double j = n, diff = 0.0, temp_n_mut = 0.0;
  cur = tail;
  cur->c3 = 0.0;
  // int num_surv_blob = 0; // This does not seem to be used, except on line 337, which has no side effects
  *n_clade = 0.0;
  while(cur->prev != NULL) {
    // num_surv_blob++;
    j -= cur->num_in_blob;
    diff = cur->lower - cur->prev->upper + 2*thres;
    temp_n_mut = diff/unitdist;
    if(max1var){
      cur->prev->c3 = *n_clade += min(temp_n_mut,1.0)/j;
    } else {
      cur->prev->c3 = *n_clade += temp_n_mut/j;
    }
    cur = cur->prev;
  }
  // DEBUG
  // cur = head;
  // while(cur != NULL) {
  //   Rprintf("%lf\n", cur->c3);
  //   cur = cur->next;
  // }

  //DEBUG
  // Rprintf("\nCOMPUTED BLOBS:\n\n");
  // cur = head;
  // int i=0;
  // while(cur != NULL) {
  //   i++;
  //   printblob(cur);
  //   cur = cur->next;
  // }
  // Rprintf("%d blobs used in total, %d survive.\n", next_new_blob+1, num_surv_blob+1);

  return(head);
}

void blobby_B1(double* alpha1, double* beta1, size_t recipient, size_t n, double thres, const double unitdist, const int max1var) {
  thres *= unitdist;
  double maxd = 744.4400719213812180897;
  blob blobs[(int) (maxd/thres+2)];
  blob** x_in_blob[n];
  double* n_clade = NULL;
  blobby_BB(alpha1, beta1, recipient, n, blobs, x_in_blob, n_clade, thres, maxd, unitdist, max1var);
}


SEXP blobbyB1(SEXP ALPHA, SEXP BETA, SEXP FROMRECIPIENT, SEXP THRES, SEXP UNITDIST, SEXP MAX1VAR) {
  double* alpha = REAL(ALPHA);
  double* beta = REAL(BETA);
  double* thres = REAL(THRES);
  double* unitdist = REAL(UNITDIST);
  int* max1var = INTEGER(MAX1VAR);
  size_t cur_left_recipient = (size_t) *INTEGER(FROMRECIPIENT) - 1;

  blobby_B1(alpha, beta, cur_left_recipient, (size_t) Rf_length(ALPHA), *thres, *unitdist, *max1var);

  return(R_NilValue);
}



void blobby_B2(const double* alpha1, const double* beta1, const double* alpha2, const double* beta2, size_t cur_left_recipient,
               size_t n, double thres, const double unitdist, const int max1var,
               double* res,
               int** neigh1, int** neigh2, int* n_neigh1, int* n_neigh2,
               double* n_clade1, double* n_clade2,
               double* similarities1, double* similarities2,
               blob*** x_in_blob1, blob*** x_in_blob2,
               const double maxd, blob* blobs1, blob* blobs2) {
  thres *= unitdist;
  blob* head1 = blobby_BB(alpha1, beta1, cur_left_recipient, n, blobs1, x_in_blob1, n_clade1, thres, maxd, unitdist, max1var);
  blob* head2 = blobby_BB(alpha2, beta2, cur_left_recipient + 1, n, blobs2, x_in_blob2, n_clade2, thres, maxd, unitdist, max1var);

  int both1, both2;
  if(head1->num_in_blob > 1) {
    *n_neigh1 = head1->num_in_blob;
    both1 = 0;
    similarities1[0] = head1->c3;
    similarities1[1] = head1->c3;
    if(head1->next != NULL) {
      similarities1[2] = head1->next->c3;
    } else {
      similarities1[2] = 0.0;
    }
  } else {
    *n_neigh1 = 1 + head1->next->num_in_blob;
    both1 = 1;
    similarities1[0] = head1->c3;
    similarities1[1] = head1->next->c3;
    if(head1->next->next != NULL) {
      similarities1[2] = head1->next->next->c3;
    } else {
      similarities1[2] = 0.0;
    }
  }
  if(head2->num_in_blob > 1) {
    *n_neigh2 = head2->num_in_blob;
    both2 = 0;
    similarities2[0] = head2->c3;
    similarities2[1] = head2->c3;
    if(head2->next != NULL) {
      similarities2[2] = head2->next->c3;
    } else {
      similarities2[2] = 0.0;
    }
  } else {
    *n_neigh2 = 1 + head2->next->num_in_blob;
    both2 = 1;
    similarities2[0] = head2->c3;
    similarities2[1] = head2->next->c3;
    if(head2->next->next != NULL) {
      similarities2[2] = head2->next->next->c3;
    } else {
      similarities2[2] = 0.0;
    }
  }

  *neigh1 = malloc(sizeof(int)* *n_neigh1);
  *neigh2 = malloc(sizeof(int)* *n_neigh2);

  // Neighbours and dedip
  size_t j = 0, ni1 = 0, ni2 = 0;
  for(size_t i = 0; i<n/2; i++) {
    if((*(x_in_blob1[j]))->blob == head1 || (both1 && (*(x_in_blob1[j]))->blob == head1->next)) {
      (*neigh1)[ni1++] = j+1;
    }
    if((*(x_in_blob1[j+1]))->blob == head1 || (both1 && (*(x_in_blob1[j+1]))->blob == head1->next)) {
      (*neigh1)[ni1++] = j+2;
    }
    if((*(x_in_blob2[j]))->blob == head2 || (both2 && (*(x_in_blob2[j]))->blob == head2->next)) {
      (*neigh2)[ni2++] = j+1;
    }
    if((*(x_in_blob2[j+1]))->blob == head2 || (both2 && (*(x_in_blob2[j+1]))->blob == head2->next)) {
      (*neigh2)[ni2++] = j+2;
    }
    res[i] = (*(x_in_blob1[j]))->blob->c3 + (*(x_in_blob1[j+1]))->blob->c3 + (*(x_in_blob2[j]))->blob->c3 + (*(x_in_blob2[j+1]))->blob->c3;
    j += 2;
  }
}

SEXP blobbyB2(SEXP ALPHA, SEXP BETA, SEXP FROMRECIPIENT, SEXP THRES, SEXP UNITDIST, SEXP MAX1VAR, SEXP DEDIP) {
  // ALPHA and BETA must each be a R matrix with two columns

  double* alpha = REAL(ALPHA);
  double* beta = REAL(BETA);
  size_t from_recipient = (size_t) Rf_asInteger(FROMRECIPIENT) - 1;
  double thres = Rf_asReal(THRES);
  double unitdist = Rf_asReal(UNITDIST);
  int max1var = Rf_asInteger(MAX1VAR);
  double* dedip = REAL(DEDIP);
  size_t n = Rf_nrows(ALPHA);

  if(Rf_ncols(ALPHA) != Rf_ncols(BETA)) {
    Rf_error("All alphas/betas must have same number of columns");
  }
  if(Rf_nrows(ALPHA) != Rf_nrows(BETA)) {
    Rf_error("All alphas/betas must have same number of rows");
  }
  if(Rf_ncols(ALPHA) != 2 || Rf_ncols(BETA) != 2) {
    Rf_error("alpha and beta must both have two columns");
  }
  if(Rf_nrows(ALPHA)/2 != Rf_length(DEDIP)) {
    Rf_error("Length of DEDIP must equal nrows(alpha)/2.");
  }


  int* neigh1;
  int* neigh2;
  int n_neigh[2];

  double simi1[3], simi2[3];

  blob*** x_in_blob1 = malloc(sizeof(blob**)*n);
  blob*** x_in_blob2 = malloc(sizeof(blob**)*n);

  const double maxd = 744.4400719213812180897;
  blob* blobs1 = malloc(sizeof(blob)*((int) (maxd/thres+2)));
  blob* blobs2 = malloc(sizeof(blob)*((int) (maxd/thres+2)));

  double* n_clade1 = NULL;
  double* n_clade2 = NULL;

  blobby_B2(alpha, beta, alpha+n, beta+n, from_recipient, n, thres, unitdist, max1var, dedip,
            &neigh1, &neigh2, &n_neigh[0], &n_neigh[1], n_clade1, n_clade2,
                                                   simi1, simi2,
                                                   x_in_blob1, x_in_blob2,
                                                   maxd, blobs1, blobs2);

  free(x_in_blob1); free(x_in_blob2);
  free(blobs1); free(blobs2);

  return(R_NilValue);
}



struct blobby_core_args {
  const double* const restrict alpha;
  const double* const restrict beta;
  double* const restrict dedip;
  const size_t from_recipient;
  const size_t n;
  const double thres;
  const double unitdist;
  const int max1var;
  int** const neigh;
  int* const n_neigh;
  double* const n_clade;
  double* const similarities1;
  double* const similarities2;
  double* const similarities3;
};

struct blobby_args {
  struct blobby_core_args *core_args;
  size_t from;
  size_t N;
};

void* blobby_B(void *args) {
  struct blobby_args *b_args;
  b_args = (struct blobby_args *) args;
  const double* restrict alpha = b_args->core_args->alpha;
  const double* restrict beta = b_args->core_args->beta;
  double* restrict dedip = b_args->core_args->dedip;
  size_t from_recipient = b_args->core_args->from_recipient;
  const size_t n = b_args->core_args->n;
  const double thres = b_args->core_args->thres;
  const double unitdist = b_args->core_args->unitdist;
  const int max1var = b_args->core_args->max1var;
  int** neigh = b_args->core_args->neigh;
  int* n_neigh = b_args->core_args->n_neigh;
  double* n_clade = b_args->core_args->n_clade;
  double* similarities1 = b_args->core_args->similarities1;
  double* similarities2 = b_args->core_args->similarities2;
  double* similarities3 = b_args->core_args->similarities3;
  size_t from = b_args->from;
  size_t N = b_args->N;

  alpha += n*from;
  beta += n*from;
  dedip += (n/2)*(from/2);
  from_recipient += from;
  neigh += from;
  n_neigh += from;
  n_clade += from;
  similarities1 += from;
  similarities2 += from;
  similarities3 += from;

  double simi1[3], simi2[3];

  blob*** x_in_blob1 = malloc(sizeof(blob**)*n);
  blob*** x_in_blob2 = malloc(sizeof(blob**)*n);

  const double maxd = 744.4400719213812180897;
  blob* blobs1 = malloc(sizeof(blob)*((int) (maxd/thres+2)));
  blob* blobs2 = malloc(sizeof(blob)*((int) (maxd/thres+2)));

  for(size_t i = 0; i < N; i+=2) {
    blobby_B2(alpha, beta, alpha+n, beta+n, from_recipient, n, thres, unitdist, max1var, dedip,
              neigh, neigh+1, n_neigh, n_neigh+1, n_clade, n_clade+1,
              simi1, simi2,
              x_in_blob1, x_in_blob2,
              maxd, blobs1, blobs2);

    similarities1[i] = simi1[0];
    similarities2[i] = simi1[1];
    similarities3[i] = simi1[2];
    similarities1[i+1] = simi2[0];
    similarities2[i+1] = simi2[1];
    similarities3[i+1] = simi2[2];

    alpha += 2*n;
    beta += 2*n;
    dedip += n/2;
    from_recipient += 2;
    neigh += 2;
    n_neigh += 2;
    n_clade += 2;
  }

  free(x_in_blob1); free(x_in_blob2);
  free(blobs1); free(blobs2);

  return(NULL);
}

void blobby_A(const double* const restrict alpha,
              const double* const restrict beta,
              double* const restrict dedip,
              size_t from_recipient,
              size_t n,
              size_t p,
              double thres,
              double unitdist,
              int max1var,
              int** const neigh,
              int* const n_neigh,
              double* const n_clade,
              double* const similarities1,
              double* const similarities2,
              double* const similarities3,
              size_t nthreads) {


  struct blobby_core_args core_args = {
    .alpha = alpha,
    .beta = beta,
    .dedip = dedip,
    .from_recipient = from_recipient,
    .n = n,
    .thres = thres,
    .unitdist = unitdist,
    .max1var = max1var,
    .neigh = neigh,
    .n_neigh = n_neigh,
    .n_clade = n_clade,
    .similarities1 = similarities1,
    .similarities2 = similarities2,
    .similarities3 = similarities3
  };

  if(nthreads > 1) {

    pthread_t threads[nthreads];
    pthread_attr_t attr;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    size_t num_perth = (p/2)/nthreads;
    size_t rag_end   = (p/2)%nthreads;

    struct blobby_args args[nthreads+1];
    for(size_t i=0; i<nthreads; i++) {
      args[i].core_args = &core_args;
      args[i].from = i*2*num_perth;
      args[i].N = 2*num_perth;
    }

    for(size_t i=0; i<nthreads; ++i) {
      pthread_create(&threads[i], &attr, blobby_B, (void*) &args[i]);
      //blobby_B((void*) &args[i]);
    }
    Rprintf("%d threads created (rag end left over = %d)\n", (int) nthreads, (int) rag_end);

    // Tidy ragged end
    if(rag_end != 0) {
      args[nthreads].core_args = &core_args;
      args[nthreads].from = 2*nthreads*num_perth;
      args[nthreads].N = 2*rag_end;
      blobby_B((void*) &args[nthreads]);
    }

    for(size_t i=0; i<nthreads; i++) {
      pthread_join(threads[i], NULL);
    }
    pthread_attr_destroy(&attr);
  } else {

    //Rprintf("\n(not using threading)\n");
    struct blobby_args args;
    args.core_args = &core_args;
    args.from = 0;
    args.N = p;

    blobby_B((void*) &args);

  }
}



SEXP CladeMat(SEXP Rfwd,
              SEXP Rbck,
              SEXP RM,
              SEXP Runitdist,
              SEXP Rthresh,
              SEXP Rmax1var,
              SEXP Rnthreads) {

  // Extract table variables
  double *alpha, *beta;
  int *from_rec, *to_rec;
  KALIS_GET_TABLE(alpha, Rfwd);
  KALIS_GET_TABLE(beta, Rbck);
  KALIS_GET_TABLE_FROM(from_rec, Rbck);
  KALIS_GET_TABLE_TO(to_rec, Rbck);

  double unitdist = Rf_asReal(Runitdist);
  double thres = Rf_asReal(Rthresh);
  int max1var = Rf_asInteger(Rmax1var);
  double* dedip = REAL(RM);
  size_t n = num_inds; // from Cache.h
  size_t p = *to_rec - *from_rec + 1;
  size_t from_recipient = (size_t) *from_rec - 1;
  size_t nthreads = (size_t) Rf_asInteger(Rnthreads);


  SEXP SIMILARITIES1 = PROTECT(Rf_allocVector(REALSXP, p));
  double* similarities1 = REAL(SIMILARITIES1);
  SEXP SIMILARITIES2 = PROTECT(Rf_allocVector(REALSXP, p));
  double* similarities2 = REAL(SIMILARITIES2);
  SEXP SIMILARITIES3 = PROTECT(Rf_allocVector(REALSXP, p));
  double* similarities3 = REAL(SIMILARITIES3);
  SEXP RES_SIMI = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(RES_SIMI, 0, SIMILARITIES1);
  SET_VECTOR_ELT(RES_SIMI, 1, SIMILARITIES2);
  SET_VECTOR_ELT(RES_SIMI, 2, SIMILARITIES3);

  int** neigh;
  neigh = malloc(sizeof(int*)*p);
  if(neigh == NULL) {
    printf("Failed allocating neigh!\n");
    exit(1);
  }
  SEXP RES_NNEIGH = PROTECT(Rf_allocVector(INTSXP, p+1));
  int* n_neigh = INTEGER(RES_NNEIGH);
  *n_neigh = 1;

  SEXP RES_NCLADE = PROTECT(Rf_allocVector(REALSXP, p));
  double* n_clade = REAL(RES_NCLADE);

  blobby_A(alpha,
           beta,
           dedip,
           from_recipient,
           n,
           p,
           thres,
           unitdist,
           max1var,
           neigh,
           n_neigh+1,
           n_clade,
           similarities1,
           similarities2,
           similarities3,
           nthreads);

  int neigh_sz = 0;
  for(size_t i=0; i<p; i++) {
    neigh_sz += n_neigh[i+1];
  }

  double tot_clade = 0;
  for(size_t i=0; i<p; i++) {
    tot_clade += n_clade[i];
  }

  SEXP RES_TOT_CLADE = PROTECT(Rf_ScalarReal(tot_clade));

  SEXP RES_NEIGH = PROTECT(Rf_allocVector(INTSXP, neigh_sz));
  int* res_neigh = INTEGER(RES_NEIGH);
  int* tmp = res_neigh;
  for(size_t i=0; i<p; i++) {
    for(size_t j=0; j<n_neigh[i+1]; j++) {
      *(tmp++) = neigh[i][j];
    }
    free(neigh[i]);
  }
  free(neigh);

  for(size_t i=0; i<p; i++) {
    n_neigh[i+1] += n_neigh[i];
  }

  SEXP RES_NE = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(RES_NE, 0, RES_NNEIGH);
  SET_VECTOR_ELT(RES_NE, 1, RES_NEIGH);
  SEXP RES = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(RES, 0, RES_NE);
  SET_VECTOR_ELT(RES, 1, RES_SIMI);
  SET_VECTOR_ELT(RES, 2, RES_TOT_CLADE);


  UNPROTECT(10);
  return(RES);
}





SEXP UpdateRealInPlace(SEXP RM,
                       SEXP Ridx,
                       SEXP Rvec) {

  double* M = REAL(RM);
  int* idx = INTEGER(Ridx);
  double* vec = REAL(Rvec);
  int p = Rf_length(Ridx);

  for(size_t i=0; i<p; i++){
    M[*idx - 1] += *vec;
    idx++;
    vec++;
  }

  return(R_NilValue);
}
