# Retrieve haplotypes from memory cache

Retrieve haplotypes from the memory cache, converting the raw binary
into an integer 0/1 matrix for inspection and use in R.

## Usage

``` r
QueryCache(loci.idx = NULL, hap.idx = NULL)
```

## Arguments

- loci.idx:

  which variants to retrieve from the cache, specified as a (vector)
  index. This enables specifying variants by offset in the order they
  were loaded into the cache (from 1 to the number of variants).

- hap.idx:

  which haplotypes to retrieve from the cache, specified as a (vector)
  index. This enables specifying haplotypes by offset in the order they
  were loaded into the cache (from 1 to the number of haplotypes).

## Value

A matrix of 0/1 integers with `length(loci.idx)` rows and
`length(hap.idx)` columns, such that haplotypes appear in columns.

## Details

To achieve higher performance, kalis internally represents haplotypes in
an efficient raw binary format in memory which cannot be directly viewed
or manipulated in R. This function enables you to copy whole or partial
views of haplotypes/variants out of this low-level format and into a
standard R matrix of 0's and 1's.

## References

Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation
of the Li & Stephens model for local ancestry inference in R", *BMC
Bioinformatics*, **25**(1). Available at:
[doi:10.1186/s12859-024-05688-8](https://doi.org/10.1186/s12859-024-05688-8)
.

## See also

[`CacheHaplotypes()`](https://kalis.louisaslett.com/reference/CacheHaplotypes.md)
to fill the memory cache with haplotype data.

## Examples

``` r
# For the purposes of an example, fill the cache with random haplotypes ...
n.haps <- 100
n.vars <- 200
haps <- matrix(sample(0:1, n.haps*n.vars, replace = TRUE),
               nrow = n.vars, ncol = n.haps)
CacheHaplotypes(haps)
#> Warning: haplotypes already cached ... overwriting existing cache.

# ... and confirm we can read a chosen portion back.  Try to read back the
# 10th and 11th haplotypes from variants 50 to 150 inclusive
res <- QueryCache(50:150, 10:11)
all(res == haps[50:150, 10:11])
#> [1] TRUE
```
