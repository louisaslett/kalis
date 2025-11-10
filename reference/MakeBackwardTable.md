# Create a `kalisBackwardTable`

Allocates the memory for and initializes a backward table.

## Usage

``` r
MakeBackwardTable(pars, from_recipient = 1, to_recipient = Inf)
```

## Arguments

- pars:

  a `kalisParameters` object specifying the genetics parameters to be
  associated with this backward table. These parameters can be set up by
  using the
  [`Parameters()`](https://kalis.louisaslett.com/reference/Parameters.md)
  function.

- from_recipient:

  first recipient haplotype included if creating a partial backward
  table. By default all are included from the first recipient haplotype.
  Haplotypes are indexed from 1.

- to_recipient:

  last recipient haplotype included if creating a partial backward
  table. By default all are included upto the last recipient haplotype.
  Haplotypes are indexed from 1.

## Value

A specialized list of class `kalisBackwardTable`. The elements of the
backward table list are:

- `l`:

  denotes the current variant position. 2147483647 indicates a newly
  created backward table which has not yet been propagated to any
  variant.

- `beta`:

  is a matrix of rescaled backward probabilities under the Li and
  Stephens HMM. Each column of `beta` corresponds to an independent HMM
  such that \\\beta^l\_{ji}\\ is proportional to the probability of
  observing haplotype \\i\\ from variant \\l+1\\ up through variant
  \\L\\ given that haplotype \\j\\ is copied by haplotype \\i\\ at
  variant \\l\\.

- `beta.g`:

  is a vector containing scaling constants needed to continue
  propagating the HMM (TODO: add ref please see kalis paper for
  details).

- `beta.theta`:

  boolean indicator for whether the matrix `beta` is currently in
  so-called "beta-theta space" or not. See
  [`Backward()`](https://kalis.louisaslett.com/reference/Backward.md)
  for details about "beta-theta space" which is specified during the
  run.

A `kalisBackwardTable` also carries with it a checksum key for the
parameters `pars` it was provided. If one attempts to interact a
`kalisBackwardTable` with a `kalisForwardTable` with mismatched
parameters, an error will be thrown.

## Details

`MakeBackwardTable` initializes a `kalisBackwardTable` object that will
store the rescaled backward probabilities for each recipient haplotype.
Note: since all other haplotypes loaded by
[`CacheHaplotypes()`](https://kalis.louisaslett.com/reference/CacheHaplotypes.md)
into the kalis cache are taken as potential donor haplotypes, this table
will be of size (total haplotypes in cache)x(number of recipients). Also
note that this table object is linked to a given set of Li and Stephens
hidden Markov model parameters, created by
[`Parameters()`](https://kalis.louisaslett.com/reference/Parameters.md),
to ensure consistency for any given HMM run.

The returned `kalisBackwardTable` is ready to be propagated to a given
target variant with the function
[`Backward()`](https://kalis.louisaslett.com/reference/Backward.md).

Since each column corresponds to an independent Li and Stephens hidden
Markov model (ie for each recipient), it is possible to create a partial
backward table object which corresponds to a subset of recipients using
the `from_recipient` and `to_recipient` arguments.

## References

Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation
of the Li & Stephens model for local ancestry inference in R", *BMC
Bioinformatics*, **25**(1). Available at:
[doi:10.1186/s12859-024-05688-8](https://doi.org/10.1186/s12859-024-05688-8)
.

## See also

[`Backward()`](https://kalis.louisaslett.com/reference/Backward.md) to
propagate the newly created `kalisBackwardTable`;
[`MakeForwardTable()`](https://kalis.louisaslett.com/reference/MakeForwardTable.md)
to create a corresponding `kalisForwardTable`.

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
bck <- MakeBackwardTable(pars)
bck
#> Full Backward Table object for 300 haplotypes, in rescaled probability space. 
#>   Newly created table, currently uninitialised to any variant (ready for Backward function next).
#>   Memory consumed: 724.13 kB 

# Now ready to run the HMM backward recursions with Backward() ...
```
