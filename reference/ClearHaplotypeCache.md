# Remove all cached haplotypes and free memory

Remove all the haplotypes that were cached by a previous caching call
and free up the memory that was allocated for future use.

## Usage

``` r
ClearHaplotypeCache()
```

## Value

Nothing is returned.

## Details

To achieve higher performance, kalis internally represents haplotypes in
an efficient raw binary format in memory which cannot be directly viewed
or manipulated in R, though you can extract a view from it using
[`QueryCache()`](https://kalis.louisaslett.com/reference/QueryCache.md).
In particular, this cache sits outside R's memory management and will
never be garbage collected (unless R is quit or the package is
unloaded). Therefore, this function is provided to enable freeing the
memory used by this cache.

## References

Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation
of the Li & Stephens model for local ancestry inference in R", *BMC
Bioinformatics*, **25**(1). Available at:
[doi:10.1186/s12859-024-05688-8](https://doi.org/10.1186/s12859-024-05688-8)
.

## See also

[`CacheHaplotypes()`](https://kalis.louisaslett.com/reference/CacheHaplotypes.md)
to create a haplotype cache;
[`QueryCache()`](https://kalis.louisaslett.com/reference/QueryCache.md)
to view the cache contents in an R matrix.

## Examples

``` r
# First fill the cache with the toy data included in the package
data("SmallHaps")
CacheHaplotypes(SmallHaps)
#> Warning: haplotypes already cached ... overwriting existing cache.

# Verify it is there
CacheSummary()
#> Cache currently loaded with 300 haplotypes, each with 400 variants. 
#>   Memory consumed: 16 kB. 

# Now clear
ClearHaplotypeCache()

# Verify it is gone
CacheSummary()
#> Cache currently empty.
```
