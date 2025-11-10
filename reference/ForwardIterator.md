# Build an efficient iterator over loci

Function factory to create a `kalisForwardIterator` for propagating a
forward table iteratively over target variants using a table cache and
optimal checkpointing.

## Usage

``` r
ForwardIterator(
  pars,
  ram.ckpts = 1L,
  targets = 1:L(),
  base.fwd.table = NULL,
  disk.ckpts = 0,
  disk.dir = NULL,
  from_recipient = 1,
  to_recipient = Inf,
  lookup.tables = NULL,
  cache = NULL,
  save.cache = FALSE,
  force.unif = FALSE
)
```

## Arguments

- pars:

  a `kalisParameters` object, as returned by
  [`Parameters()`](https://kalis.louisaslett.com/reference/Parameters.md).

- ram.ckpts:

  an integer specifying the number of checkpoints to store in RAM.

- targets:

  a vector of variants to iterate over (starting with the most
  downstream target).

- base.fwd.table:

  a `kalisForwardTable` either at the most upstream target, or if the
  targets are evenly spaced, one interval upstream of the most upstream
  target. If `NULL` (the default), this is interpreted as the prior
  `Pi`, see
  [`Parameters()`](https://kalis.louisaslett.com/reference/Parameters.md).

- disk.ckpts:

  an integer specifying the number of checkpoints to store on disk.

- disk.dir:

  a path to a directory where a temporary folder may be made to store
  checkpoints on disk.

- from_recipient:

  first recipient haplotype included in the tables of the cache, if
  creating a partial forward table. By default all are included from the
  first recipient haplotype. Haplotypes are indexed from 1.

- to_recipient:

  last recipient haplotype included in the tables of the cache, if
  creating a partial forward table. By default all are included upto the
  last recipient haplotype. Haplotypes are indexed from 1.

- lookup.tables:

  an optional list as returned by
  [`CalcCheckpointTables()`](https://kalis.louisaslett.com/reference/CalcCheckpointTables.md).

- cache:

  a `kalisCheckpointTable` object, as returned by
  [`CreateForwardTableCache()`](https://kalis.louisaslett.com/reference/CreateForwardTableCache.md)
  or this function. By default `NULL`, which causes this function to
  create a new cache.

- save.cache:

  a logical. When `TRUE` does not reliquish the table cache upon
  exhaustion of the iterator. Defaults to `FALSE`.

- force.unif:

  a logical, if `TRUE` iterate over targets as if they were uniformly
  spaced. WARNING: DO NOT use this in conjunction with the targets
  method, still experimental. With `force.unif = TRUE`, the resulting
  iterator will appear to be targeting the first `length(targets)`
  variants with all methods, but in fact will be silently iterating over
  the original targets.

## Value

A function for iterating over the set of target variants. The returned
function has prototype:

`function(fwd, pars, t, nthreads = 1)`

which matches the standard
[`Forward()`](https://kalis.louisaslett.com/reference/Forward.md)
function, but which uses the table cache to speed up propagation to the
target variant. See
[`Forward()`](https://kalis.louisaslett.com/reference/Forward.md) for an
explanation of arguments.

## Details

See example.

## References

Christ, R.R., Wang, X., Aslett, L.J.M., Steinsaltz, D. and Hall, I.
(2024) "Clade Distillation for Genome-wide Association Studies", bioRxiv
2024.09.30.615852. Available at:
[doi:10.1101/2024.09.30.615852](https://doi.org/10.1101/2024.09.30.615852)
.

## See also

[`MakeForwardTable()`](https://kalis.louisaslett.com/reference/MakeForwardTable.md)
to create a `kalisForwardTable`;
[`CreateForwardTableCache()`](https://kalis.louisaslett.com/reference/CreateForwardTableCache.md)
to create a cache which can be used with this function.

## Examples

``` r
if (FALSE) { # \dontrun{
data("SmallHaps")
CacheHaplotypes(SmallHaps)
pars <- Parameters()
fwd <- MakeForwardTable(pars)
bck <- MakeBackwardTable(pars)
Iter <- ForwardIterator(2)
for(t in targets(Iter)){
  Iter(fwd,pars,t)
  Backward(bck,pars,t)
  print(paste("Mean Distance at locus",t,"is",mean(DistMat(fwd,bck))))
}
} # }
```
