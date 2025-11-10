# Small example phased haplotype dataset and recombination map

Simulated dataset with \\N=300\\ haplotypes at \\L=400\\ variants
generated using msprime (Kelleher et al., 2016), together with the
recombination map.

## Usage

``` r
data("SmallHaps")

data("SmallMap")
```

## Format

For `SmallHaps`, an object of class `matrix` of dimensions \\400 \times
300\\ where each column is a simulated haplotype. This is suitable for
passing directly to
[`CacheHaplotypes()`](https://kalis.louisaslett.com/reference/CacheHaplotypes.md).

For `SmallMap`, a vector of length 400 representing the recombination
map for the `SmallHaps` data. This can be used with
[`CalcRho()`](https://kalis.louisaslett.com/reference/CalcRho.md), by
converting to recombination distances using `diff(SmallMap)`.

## References

Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation
of the Li & Stephens model for local ancestry inference in R", *BMC
Bioinformatics*, **25**(1). Available at:
[doi:10.1186/s12859-024-05688-8](https://doi.org/10.1186/s12859-024-05688-8)
.

Kelleher, J., Etheridge, A.M. and McVean, G. (2016) "Efficient
coalescent simulation and genealogical analysis for large sample sizes",
*PLoS Computational Biology*, **12**(5). Available at:
[doi:10.1371/journal.pcbi.1004842](https://doi.org/10.1371/journal.pcbi.1004842)
.

## Examples

``` r
data("SmallHaps")

# \donttest{
# Plot Allele Frequencies
hist(rowMeans(SmallHaps),breaks=20)

# }

# Import into kalis cache directly ...
CacheHaplotypes(SmallHaps)
#> Warning: haplotypes already cached ... overwriting existing cache.

data("SmallMap")

# Find parameters
pars <- Parameters(CalcRho(diff(SmallMap)))
pars
#> Parameters object with:
#>   rho   = (0.030722976221863, 0.00188322071638501, 6.79809190046651e-06, ..., 1.56898508137998e-05, 1.7924679243824e-05, 1)
#>   mu    = 1e-08
#>   Pi    = 0.00334448160535117 
```
