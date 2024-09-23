#include "R_OptCkpt.h"

#include <limits.h> // DBL_MAX

#include "R_Kalis.h"

SEXP OptCkpt(SEXP Rcost_table, SEXP Rindex_table, SEXP Rpropagation_cost) {

  int cost_nr = Rf_nrows(Rcost_table);
  int index_nr = Rf_nrows(Rindex_table);

  if(cost_nr != index_nr + 1) {
    REprintf("Error: incorrect number of rows in cost/index table combination.\n");
    KALIS_RETURN;
  }
  if(Rf_ncols(Rcost_table) != Rf_ncols(Rindex_table) + 1) {
    REprintf("Error: incorrect number of cols in cost/index table combination.\n");
    KALIS_RETURN;
  }

  int k,n,i;

  double x = DBL_MAX;
  int xi = 0;
  double y = 0;

  int maxckpt = Rf_ncols(Rindex_table);
  int maxn = index_nr;

  double * p_c = REAL(Rcost_table);
  int * p_i = INTEGER(Rindex_table);
  double * s = REAL(Rpropagation_cost);


  for(k = 0; k < maxckpt; k++)
  {
    for(n = 0; n < maxn; n++)
    {
      x = DBL_MAX;
      xi = 0;
      y = 0;

      for(i = n; i > -1; i--)
      {
        y = p_c[i + cost_nr*(k+1)] + s[i] + p_c[n-i + cost_nr*k];
        if(y <= x) {
          xi = i;
          x = y;
        } else {
          break;
        }
      }
      p_i[n + index_nr * k] = xi+1;
      p_c[n+1 + cost_nr * (k+1)] = x;
    }
  }
  KALIS_RETURN;
}
