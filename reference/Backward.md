# Propagate a backward table

Propagates a `kalisBackwardTable` to an upstream variant. The table is
updated in-place.

## Usage

``` r
Backward(
  bck,
  pars,
  t = bck$l - 1,
  nthreads = min(parallel::detectCores(logical = FALSE), bck$to_recipient -
    bck$from_recipient + 1),
  beta.theta = FALSE
)
```

## Arguments

- bck:

  a `kalisBackwardTable` object, as returned by
  [`MakeBackwardTable()`](https://kalis.louisaslett.com/reference/MakeBackwardTable.md).

- pars:

  a `kalisParameters` object, as returned by
  [`Parameters()`](https://kalis.louisaslett.com/reference/Parameters.md).

- t:

  a variant to move the backward table to. Must be less than or equal to
  current variant of `bck`. By default, it advances to the next variant
  upstream (or if an uninitialised table, to the last variant).

- nthreads:

  if a scalar, the number of CPU cores to use. If a vector, launch as
  many threads as the length of the vector and attempt to pin the
  threads to those CPU cores (requires system to support thread
  affinity). By default uses the `parallel` package to detect the number
  of physical cores.

- beta.theta:

  logical indicating whether the table should be returned in beta-theta
  space or in the standard space upon reaching the target variant `t`.
  See the Details section.

## Value

There is nothing returned.

**NOTE:** for performance reasons, `bck` is updated in-place.

## Details

`Backward` implements the backward algorithm to propagate the Li and
Stephens rescaled hidden Markov model backward probabilities to a new
target variant. Naturally, this can only propagate a table to variants
upstream of its current position.

For mathematical details please see Section 2 of the kalis paper (Aslett
and Christ, 2024). Note that the precise formulation of the backward
equation is determined by whether the flag `use.spiedel` is set in the
parameters provided in `pars`.

**Beta-theta space**

The rescaled HMM backward probabilities incorporate all of the haplotype
relatedness information downstream from but NOT including the target
variant, as is standard in the definition of HMM backward probabilities
– we refer to this as "beta space", or "rescaled probability space." In
contrast, the rescaled forward probabilities, naturally incorporate all
of the haplotype relatedness information upstream from AND including the
target variant. Setting `beta.theta = TRUE` incorporates the target
variant (analogous to the rescaled forward probabilities) – we refer to
this as "beta-theta space."

A backward table in beta-theta space (with `beta.theta = TRUE`) can be
propagated to an upstream variant without incorporating that variant,
thereby moving to beta space (`beta.theta = FALSE`), and vice versa.
However, while a backward table in beta space (`beta.theta = FALSE`) can
be updated to incorporate the current variant, a backward table that is
already in beta-theta space can not move to beta space without changing
variants – that would involve "forgetting" the current variant (see
Examples).

## References

Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation
of the Li & Stephens model for local ancestry inference in R", *BMC
Bioinformatics*, **25**(1). Available at:
[doi:10.1186/s12859-024-05688-8](https://doi.org/10.1186/s12859-024-05688-8)
.

## See also

[`MakeBackwardTable()`](https://kalis.louisaslett.com/reference/MakeBackwardTable.md)
to generate a backward table;
[`Forward()`](https://kalis.louisaslett.com/reference/Forward.md) for
the analogous forward propagation function;
[`CopyTable()`](https://kalis.louisaslett.com/reference/CopyTable.md) to
create a copy of table.

## Examples

``` r
# Load the toy haplotype example and set toy parameters
data("SmallHaps")
data("SmallMap")

CacheHaplotypes(SmallHaps)

rho <- CalcRho(diff(SmallMap))
pars <- Parameters(rho)

# Create the backward table we want to propagate
bck <- MakeBackwardTable(pars)

# Calling Backward on this uninitialised table moves it to the last variant
Backward(bck, pars)
bck
#> Full Backward Table object for 300 haplotypes, in rescaled probability space. 
#>   Current variant = 400 
#>   Memory consumed: 724.13 kB 

# And again moves it to the next variant (etcetera)
Backward(bck, pars)
bck
#> Full Backward Table object for 300 haplotypes, in rescaled probability space. 
#>   Current variant = 399 
#>   Memory consumed: 724.13 kB 

# Or zoom to a particular variant
Backward(bck, pars, 150)
bck
#> Full Backward Table object for 300 haplotypes, in rescaled probability space. 
#>   Current variant = 150 
#>   Memory consumed: 724.13 kB 

# Now moving to variant 125 AND specifying "beta space" (though this is the
# default, just being very clear)
Backward(bck, pars, 125, beta.theta = FALSE)
bck
#> Full Backward Table object for 300 haplotypes, in rescaled probability space. 
#>   Current variant = 125 
#>   Memory consumed: 724.13 kB 

# Now just moving to "beta-theta space" (can be done with or without
# propagating to a new variant)
Backward(bck, pars, 125, beta.theta = TRUE)
bck
#> Full Backward Table object for 300 haplotypes, in beta-theta space. 
#>   Current variant = 125 
#>   Memory consumed: 724.13 kB 

# However, attempting to move from "beta-theta space" back to "beta space"
# without propagating is not possible (see Details).
# The following will give an error (hence wrapped in try())
try(Backward(bck, pars, 125, beta.theta = FALSE))
#> Error in Backward(bck, pars, 125, beta.theta = FALSE) : 
#>   Cannot move from beta-theta space to beta space without moving at least one variant.
bck
#> Full Backward Table object for 300 haplotypes, in beta-theta space. 
#>   Current variant = 125 
#>   Memory consumed: 724.13 kB 
```
