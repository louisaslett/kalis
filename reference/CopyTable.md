# Copy Forward/Backward tables

Copies the contents of one forward/backward table into another.

## Usage

``` r
CopyTable(to, from)
```

## Arguments

- to:

  a `kalisForwardTable` or `kalisBackwardTable` object which is to be
  copied into.

- from:

  a `kalisForwardTable` or `kalisBackwardTable` object which is to be
  copied from.

## Details

The core code in kalis operates on forward and backward tables at a very
low level, both for speed (using low level CPU vector instructions) but
also to avoid unnecessary memory copies since these tables will tend to
be very large in serious genetics applications. As a result, if you
attempt to copy a table in the standard idomatic way in R:

    fwd2 <- fwd

then these two variables merely point to the *same* table: running the
forward algorithm on `fwd` would result in `fwd2` also changing.

This function is therefore designed to enable explicit copying of
tables.

## References

Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation
of the Li & Stephens model for local ancestry inference in R", *BMC
Bioinformatics*, **25**(1). Available at:
[doi:10.1186/s12859-024-05688-8](https://doi.org/10.1186/s12859-024-05688-8)
.

## See also

[`MakeForwardTable()`](https://kalis.louisaslett.com/reference/MakeForwardTable.md),
[`MakeBackwardTable()`](https://kalis.louisaslett.com/reference/MakeBackwardTable.md)
to create tables which can be copied into.

## Examples

``` r
# Load the toy haplotype example and set toy parameters
data("SmallHaps")
data("SmallMap")

CacheHaplotypes(SmallHaps)

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

# This does **NOT** work as intended:
fwd2 <- fwd
Forward(fwd, pars, 20)

# Both tables are now at variant 20
fwd
#> Full Forward Table object for 300 haplotypes. 
#>   Current variant = 20 
#>   Memory consumed: 723.98 kB 
fwd2
#> Full Forward Table object for 300 haplotypes. 
#>   Current variant = 20 
#>   Memory consumed: 723.98 kB 

# Instead, to copy we create another table and use this function
fwd2 <- MakeForwardTable(pars)
CopyTable(fwd2, fwd)

# Now propagate fwd again
Forward(fwd, pars, 30)

# But notice this fwd2 wasn't affected
fwd
#> Full Forward Table object for 300 haplotypes. 
#>   Current variant = 30 
#>   Memory consumed: 723.98 kB 
fwd2
#> Full Forward Table object for 300 haplotypes. 
#>   Current variant = 20 
#>   Memory consumed: 723.98 kB 
```
