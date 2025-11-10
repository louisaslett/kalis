# Create cache for forward tables

Create an in-memory cache for forward tables to improve efficiency when
iterating in reverse along the haplotype sequences.

## Usage

``` r
CreateForwardTableCache(
  pars,
  size = 1,
  from_recipient = 1,
  to_recipient = Inf,
  max.tables = 0
)
```

## Arguments

- pars:

  a `kalisParameters` object, as returned by
  [`Parameters()`](https://kalis.louisaslett.com/reference/Parameters.md).

- size:

  the maximum amount of RAM (in GB) to devote to this cache. By default,
  1GB.

- from_recipient:

  first recipient haplotype if creating a partial forward table cache.
  By default includes from the first recipient haplotype.

- to_recipient:

  last recipient haplotype if creating a partial forward table cache. By
  default includes to the last recipient haplotype.

- max.tables:

  positive integer indicating the maximum number of tables to use in the
  cache. Both this and `size` will be honoured, so the number of tables
  may be smaller than this. By default, equals
  \\\lfloor\log_2(L)\rfloor\\.

## Value

A list of forward tables representing a cache and ready to be filled is
returned.

## Details

If the objective is to run the Li and Stephens hidden Markov model both
forwards and backwards to the same locus (and to do so for every
possible locus), then considerable efficiency can be achieved by first
performing a full scan forwards, filling a geometrically spaced cache
whilst doing so. Then, by working backwards, the backward propagation
moves one locus at a time and the forward propagation can move backwards
by moving forward from a recently cached local table.

Memory for a cache can be allocated using this function and should then
be filled with
[`FillTableCache()`](https://kalis.louisaslett.com/reference/FillTableCache.md).
To use the cache, then instead of using the
[`Forward()`](https://kalis.louisaslett.com/reference/Forward.md)
function, use
[`ForwardUsingTableCache()`](https://kalis.louisaslett.com/reference/ForwardUsingTableCache.md).

## References

Christ, R.R., Wang, X., Aslett, L.J.M., Steinsaltz, D. and Hall, I.
(2024) "Clade Distillation for Genome-wide Association Studies", bioRxiv
2024.09.30.615852. Available at:
[doi:10.1101/2024.09.30.615852](https://doi.org/10.1101/2024.09.30.615852)
.

## See also

[`MakeForwardTable()`](https://kalis.louisaslett.com/reference/MakeForwardTable.md)
to make a forward table;
[`FillTableCache()`](https://kalis.louisaslett.com/reference/FillTableCache.md)
to fill a cache;
[`ForwardUsingTableCache()`](https://kalis.louisaslett.com/reference/ForwardUsingTableCache.md)
to use a cache;
[`Forward()`](https://kalis.louisaslett.com/reference/Forward.md) for
forward function without using a cache.

Alternatively, see
[`ForwardIterator()`](https://kalis.louisaslett.com/reference/ForwardIterator.md)
to create an iterator which internally creates or uses a table cache.

## Examples

``` r
if (FALSE) { # \dontrun{
# This code assumes you have already:
#  i) cached the haplotypes using CacheHaplotypes function
#  ii) setup parameters in a variable called pars
#  iii) set the number of loci in a variable called L

# Allocate up to 10GB to a cache, with parameters already setup in pars ...
cache <- CreateForwardTableCache(pars, 10)
# ... and fill it
FillTableCache(cache, pars, nthreads = 8)

# Create forward and backward tables
fwd <- MakeForwardTable(pars)
bck <- MakeBackwardTable(pars)

# Then reach every locus faster by iterating backwards, using the cache to
# move the forward table into position faster
for(l in L:1) {
  Backward(bck, pars, l, nthreads = 8)
  ForwardUsingTableCache(fwd, pars, cache, l, nthreads = 8)
  # Do whatever work is required at
  # every locus here using fwd and bck
}
} # }
```
