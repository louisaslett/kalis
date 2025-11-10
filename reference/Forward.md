# Propagate a forward table

Propagates a `kalisForwardTable` to a downstream variant. The table is
updated in-place.

## Usage

``` r
Forward(
  fwd,
  pars,
  t = fwd$l + 1,
  nthreads = min(parallel::detectCores(logical = FALSE), fwd$to_recipient -
    fwd$from_recipient + 1)
)
```

## Arguments

- fwd:

  a `kalisForwardTable` object, as returned by
  [`MakeForwardTable()`](https://kalis.louisaslett.com/reference/MakeForwardTable.md).

- pars:

  a `kalisParameters` object, as returned by
  [`Parameters()`](https://kalis.louisaslett.com/reference/Parameters.md).

- t:

  a variant to move the forward table to. Must be greater than or equal
  to current variant of `fwd`. By default, it advances to the next
  variant downstream (or if an uninitialised table, to the first
  variant).

- nthreads:

  if a scalar, the number of CPU cores to use. If a vector, launch as
  many threads as the length of the vector and attempt to pin the
  threads to those CPU cores (requires system to support thread
  affinity). By default uses the `parallel` package to detect the number
  of physical cores.

## Value

There is nothing returned.

**NOTE:** for performance reasons, `fwd` is updated in-place.

## Details

`Forward` implements the forward algorithm to advance the Li and
Stephens rescaled hidden Markov model forward probabilities to a new
target variant. Naturally, this can only propagate a table to variants
downstream of its current position.

For mathematical details please see Section 2 of the kalis paper (Aslett
and Christ, 2024). Note that the precise formulation of the forward
equation is determined by whether the flag `use.spiedel` is set in the
parameters provided in `pars`.

## References

Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation
of the Li & Stephens model for local ancestry inference in R", *BMC
Bioinformatics*, **25**(1). Available at:
[doi:10.1186/s12859-024-05688-8](https://doi.org/10.1186/s12859-024-05688-8)
.

## See also

[`MakeForwardTable()`](https://kalis.louisaslett.com/reference/MakeForwardTable.md)
to generate a forward table;
[`Backward()`](https://kalis.louisaslett.com/reference/Backward.md) for
the analogous backward propagation function;
[`CopyTable()`](https://kalis.louisaslett.com/reference/CopyTable.md) to
create a copy of table.

## Examples

``` r
# Load the toy haplotype example and set toy parameters
data("SmallHaps")
data("SmallMap")

CacheHaplotypes(SmallHaps)
#> Warning: haplotypes already cached ... overwriting existing cache.

rho <- CalcRho(diff(SmallMap))
pars <- Parameters(rho)

# Create the forward table we want to propagate
fwd <- MakeForwardTable(pars)

# Calling Forward on this uninitialised table moves it to the first variant
Forward(fwd, pars)
fwd
#> Full Forward Table object for 300 haplotypes. 
#>   Current variant = 1 
#>   Memory consumed: 723.98 kB 

# And again moves it to the next variant (etcetera)
Forward(fwd, pars)
fwd
#> Full Forward Table object for 300 haplotypes. 
#>   Current variant = 2 
#>   Memory consumed: 723.98 kB 

# Or zoom to a particular variant
Forward(fwd, pars, 50)
fwd
#> Full Forward Table object for 300 haplotypes. 
#>   Current variant = 50 
#>   Memory consumed: 723.98 kB 
```
