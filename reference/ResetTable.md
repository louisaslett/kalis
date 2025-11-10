# Reset Forward/Backward tables

Resets the provided forward/backward table so it can propagate from the
start/end of the Li and Stephens hidden Markov model.

## Usage

``` r
ResetTable(tbl)
```

## Arguments

- tbl:

  a `kalisForwardTable` or `kalisBackwardTable` object which is to be
  reset.

## Details

It is *much* faster to reset a forward/backward table rather than remove
and make a new one. This function marks a table as reset so that it will
be propagated as if freshly allocated.

## References

Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation
of the Li & Stephens model for local ancestry inference in R", *BMC
Bioinformatics*, **25**(1). Available at:
[doi:10.1186/s12859-024-05688-8](https://doi.org/10.1186/s12859-024-05688-8)
.

## See also

[`MakeForwardTable()`](https://kalis.louisaslett.com/reference/MakeForwardTable.md),
[`MakeBackwardTable()`](https://kalis.louisaslett.com/reference/MakeBackwardTable.md)
to create fresh tables.

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

# Move to variant 10
Forward(fwd, pars, 10)
fwd
#> Full Forward Table object for 300 haplotypes. 
#>   Current variant = 10 
#>   Memory consumed: 723.98 kB 

# This does **NOT** work as intended, cannot reverse from variant 10:
try(Forward(fwd, pars, 5))
#> Error in Forward(fwd, pars, 5) : 
#>   The forward table provided is for variant position 10 which is already past requested variant 5

# Could reset the table though:
ResetTable(fwd)
Forward(fwd, pars, 5)
fwd
#> Full Forward Table object for 300 haplotypes. 
#>   Current variant = 5 
#>   Memory consumed: 723.98 kB 
```
