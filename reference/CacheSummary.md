# Retrieve information about the haplotype cache

Retrieve information about the haplotype cache

## Usage

``` r
CacheSummary()

N()

L()
```

## Value

`CacheSummary()` prints information about the current state of the kalis
cache. Also invisibly returns a vector giving the dimensions of the
cached haplotype data (num variants, num haplotypes), or `NULL` if the
cache is empty.

`N()` returns the number of haplotypes currently in the kalis cache, or
`NULL` if the cache is empty.

`L()` returns the number of variants currently in the kalis cache, or
`NULL` if the cache is empty.

## References

Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation
of the Li & Stephens model for local ancestry inference in R", *BMC
Bioinformatics*, **25**(1). Available at:
[doi:10.1186/s12859-024-05688-8](https://doi.org/10.1186/s12859-024-05688-8)
.

## Examples

``` r
# First fill the cache with the toy data included in the package
data("SmallHaps")
CacheHaplotypes(SmallHaps)
#> Warning: haplotypes already cached ... overwriting existing cache.

# View full summary
CacheSummary()
#> Cache currently loaded with 300 haplotypes, each with 400 variants. 
#>   Memory consumed: 16 kB. 

# Also note the invisible return
dims <- CacheSummary()
#> Cache currently loaded with 300 haplotypes, each with 400 variants. 
#>   Memory consumed: 16 kB. 
dims
#> [1] 400 300

# Get just numbers of haplotypes and variants separately
N()
#> [1] 300
L()
#> [1] 400
```
