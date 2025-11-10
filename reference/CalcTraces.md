# Fast Calculation of Matrix Trace and Hilbert Schmidt Norm

Provides multithreaded calculation of trace and Hilbert Schmidt Norm of
a matrix \\PMP\\ (where \\P\\ is a projection matrix and \\M\\ is real
symmetric) without explicitly forming \\PMP\\.

## Usage

``` r
CalcTraces(
  M,
  tX,
  tQ,
  J,
  from_recipient = 1L,
  nthreads = min(parallel::detectCores(logical = FALSE), ncol(M))
)
```

## Arguments

- M:

  a real symmetric R matrix

- tX:

  `t((Q %*% (J%*%Q)) - (M %*% Q))`

- tQ:

  `t(Q)`

- J:

  `crossprod(Q, M)`

- from_recipient:

  haplotype index at which to start trace calculation — useful for
  distributed computation (experimental feature, more documentation to
  come)

- nthreads:

  the number of CPU cores to use. By default uses the `parallel` package
  to detect the number of physical cores.

## Value

A list containing three elements:

- `trace`:

  the trace, \\\mathrm{tr}(PMP)\\;

- `hsnorm2`:

  the *squared* Hilbert Schmidt Norm of \\PMP\\,
  \\\mathrm{tr}((PMP)'PMP)\\;

- `diag`:

  the diagonal of \\PMP\\.

## Details

\\P\\ here is assumed to have the form \\I-QQ'\\ for some matrix \\Q\\
of orthogonal columns.

## References

Christ, R.R., Wang, X., Aslett, L.J.M., Steinsaltz, D. and Hall, I.
(2024) "Clade Distillation for Genome-wide Association Studies", bioRxiv
2024.09.30.615852. Available at:
[doi:10.1101/2024.09.30.615852](https://doi.org/10.1101/2024.09.30.615852)
.

## Examples

``` r
# TODO
```
