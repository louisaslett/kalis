# Create a `kalisForwardTable`

Allocates the memory for and initializes a forward table.

## Usage

``` r
MakeForwardTable(pars, from_recipient = 1, to_recipient = Inf)
```

## Arguments

- pars:

  a `kalisParameters` object specifying the genetics parameters to be
  associated with this forward table. These parameters can be set up by
  using the
  [`Parameters()`](https://kalis.louisaslett.com/reference/Parameters.md)
  function.

- from_recipient:

  first recipient haplotype included if creating a partial forward
  table. By default all are included from the first recipient haplotype.
  Haplotypes are indexed from 1.

- to_recipient:

  last recipient haplotype included if creating a partial forward table.
  By default all are included upto the last recipient haplotype.
  Haplotypes are indexed from 1.

## Value

A specialized list of class `kalisForwardTable`. The elements of the
forward table list are:

- `l`:

  denotes the current variant position. 2147483647 indicates a newly
  created forward table which has not yet been propagated to any
  variant.

- `alpha`:

  is a matrix of rescaled forward probabilities under the Li and
  Stephens HMM. Each column of `alpha` corresponds to an independent HMM
  such that \\\alpha^l\_{ji}\\ is proportional to the probability that
  haplotype \\j\\ is copied by haplotype \\i\\ at variant \\l\\ and
  observing haplotype \\i\\ from variant 1 up through variant \\l\\.

- `alpha.f`:

  is a vector containing scaling constants needed to continue
  propagating the HMM (TODO: add ref please see kalis paper for
  details).

A `kalisForwardTable` also carries with it a checksum key for the
parameters `pars` it was provided. If one attempts to interact a
`kalisForwardTable` with a `kalisBackwardTable` with mismatched
parameters, an error will be thrown.

## Details

`MakeForwardTable` initializes a `kalisForwardTable` object that will
store the rescaled forward probabilities for each recipient haplotype.
Note: since all other haplotypes loaded by
[`CacheHaplotypes()`](https://kalis.louisaslett.com/reference/CacheHaplotypes.md)
into the kalis cache are taken as potential donor haplotypes, this table
will be of size (total haplotypes in cache)x(number of recipients). Also
note that this table object is linked to a given set of Li and Stephens
hidden Markov model parameters, created by
[`Parameters()`](https://kalis.louisaslett.com/reference/Parameters.md),
to ensure consistency for any given HMM run.

The returned `kalisForwardTable` is ready to be propagated to a given
target variant with the function
[`Forward()`](https://kalis.louisaslett.com/reference/Forward.md).

Since each column corresponds to an independent Li and Stephens hidden
Markov model (ie for each recipient), it is possible to create a partial
forward table object which corresponds to a subset of recipients using
the `from_recipient` and `to_recipient` arguments.

## References

Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation
of the Li & Stephens model for local ancestry inference in R", *BMC
Bioinformatics*, **25**(1). Available at:
[doi:10.1186/s12859-024-05688-8](https://doi.org/10.1186/s12859-024-05688-8)
.

## See also

[`Forward()`](https://kalis.louisaslett.com/reference/Forward.md) to
propagate the newly created `kalisForwardTable`;
[`MakeBackwardTable()`](https://kalis.louisaslett.com/reference/MakeBackwardTable.md)
to create a corresponding `kalisBackwardTable`.

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
fwd
#> Full Forward Table object for 300 haplotypes. 
#>   Newly created table, currently uninitialised to any variant (ready for Forward function next).
#>   Memory consumed: 723.98 kB 

# Now ready to run the HMM forward recursions with Forward() ...
```
