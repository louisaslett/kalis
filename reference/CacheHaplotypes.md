# Load haplotypes into optimised package cache

Load haplotypes from hard drive or an R matrix into an optimised kalis
package memory cache (overwrites any previous load).

## Usage

``` r
CacheHaplotypes(
  haps,
  loci.idx = NULL,
  hap.idx = NULL,
  warn.singletons = TRUE,
  format = "auto",
  ...
)
```

## Arguments

- haps:

  can be the name of a file from which the haplotypes are to be read, or
  can be an R matrix containing only 0/1s. See Details section for
  supported file types.

- loci.idx:

  an optional vector of indices specifying the variants to load into the
  cache, indexed from 1.

- hap.idx:

  an optional vector of indices specifying the haplotypes to load into
  the cache, indexed from 1.

- warn.singletons:

  a logical, if `FALSE`, suppress warning that singletons (variants
  where there is only one 1 or only one 0) are present in the loaded
  `haps`.

- format:

  the file format that `haps` is stored in, or `"auto"` to detect the
  format based on the file extension. Recognised options are `"hapgz"`
  (format used by IMPUTE2 and SHAPEIT) or `"hdf5"` (custom). See Details
  section for more information, and for easy conversion from VCF/BCF and
  other formats see the Examples section.

- ...:

  format specific options for reading in `haps`. Supported optional
  arguments for each format are:

  1.  For `"hapgz"`

      - `legendgz.file` a string for faster loading: a `.legend.gz` file
        can be supplied and will be used to more efficiently determine
        the number of variants in the `.hap.gz` file

      - `L` an integer for faster loading: the number of variants in the
        `.hap.gz` file can be directly provided

      - `N` an integer for faster loading: the number of haplotypes in
        the `.hap.gz` file can be directly provided

  2.  For `"hdf5"`

      - `transpose` a logical, if `TRUE`, switch the interpretation of
        rows and columns in `haps`: hence switching the number of
        haplotypes and the number of variants (the HDF5 specification
        does not prescribe row/column interpretation, only defining the
        slowest changing dimension as 'first'). Defaults to `FALSE`.

      - `haps.path` a string giving the path to a 2-dimensional object
        in the HDF5 file specifying the haplotype matrix. Defaults to
        `/haps`

      - `hdf5.pkg` a string giving the HDF5 R package to use to load the
        file from disk. The packages `rhdf5` (BioConductor) and `hdf5r`
        (CRAN) are both supported. Default is to use `hdf5r` if both
        packages are available, with fallback to `rhdf5`. This should
        never need to be specified unless you have both packages but
        want to force the use of the `rhdf5` package.

  3.  R matrix

      - `transpose` a logical, if `TRUE`, switch the interpretation of
        rows and columns in `haps`: hence switching the number of
        haplotypes and the number of variants. Defaults to `FALSE`,
        meaning variants are taken to be in rows with haplotypes in
        columns (ie a num variants x num haplotypes matrix)

## Value

A vector giving the dimensions of the cached haplotype data is invisibly
returned (num variants, num haplotypes). It is highly recommended that
you run
[`CacheSummary()`](https://kalis.louisaslett.com/reference/CacheSummary.md)
after `CacheHaplotypes`, especially if you are uncertain about the
interpretation of rows and columns in `haps`. If
[`CacheSummary()`](https://kalis.louisaslett.com/reference/CacheSummary.md)
shows that the number of haplotypes and variants are reversed, try
calling `CacheHaplotypes` again with the extra argument
`transpose = TRUE`.

## Details

To achieve higher performance, kalis internally represents haplotypes in
an efficient raw binary format in memory. This function will load
haplotypes from a file or from a binary R matrix and convert this into
kalis' internal format ready for use by the other functions in this
package. Note that only one set of haplotypes can be cached at a time
and calling this function twice overwrites cache of haplotypes created
by the first function call.

Including singletons (variants where there is only one 1 or only one 0)
in the loaded haplotypes can lead to numerical instability and columns
of `NaN`s in the resulting forward and backward tables when `mu` (see
[`Parameters()`](https://kalis.louisaslett.com/reference/Parameters.md))
is small. Thus, kalis throws a warning when loaded haplotypes contain
singletons.

At present, hap.gz and hdf5 are supported natively, see the Examples
section below showing for how to convert from a VCF/BCF to hap.gz with
one `bcftools` command.

**hap.gz format**

This is the HAP/LEGEND/SAMPLE format used by IMPUTE2 and SHAPEIT. Only
the `.hap.gz` file is required for loading with `CacheHaplotypes`,
though the `.legend.gz` file can speed up reading the haplotypes. See
<http://samtools.github.io/bcftools/bcftools.html#convert> for more
details on this format.

**R matrix**

If supplying an R matrix, it must consist of only 0's or 1's. The
haplotypes should be stored in columns, with variants in rows. That is,
the dimensions should be:

(num rows)x(num cols) = (num variants)x(num haplotypes).

It is fine to delete this matrix from R after calling
`CacheHaplotypes()`.

**HDF5 format**

For HDF5 files, kalis expects a 2-dimensional object named `haps` at the
root level of the HDF5 file. Haplotypes should be stored in the slowest
changing dimension as defined in the HDF5 specification (note that
different languages treat this as rows or columns). If the haplotypes
are stored in the other dimension then simply set the argument
`transpose = TRUE`. If the user is unsure of the convention of the
language they used to create the HDF5 file, then the simplest approach
is to simply load the data specifying only the HDF5 file name and then
confirm that number of haplotypes and their length have not been
exchanged in the diagnostic output which kalis prints.

## References

Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation
of the Li & Stephens model for local ancestry inference in R", *BMC
Bioinformatics*, **25**(1). Available at:
[doi:10.1186/s12859-024-05688-8](https://doi.org/10.1186/s12859-024-05688-8)
.

## See also

[`CacheSummary()`](https://kalis.louisaslett.com/reference/CacheSummary.md)
for a list detailing the current cache status;
[`QueryCache()`](https://kalis.louisaslett.com/reference/QueryCache.md)
to copy the haplotypes in the `kalis` cache into an R matrix;
[`ClearHaplotypeCache()`](https://kalis.louisaslett.com/reference/ClearHaplotypeCache.md)
to remove the haplotypes from the cache and free the memory;
[`N()`](https://kalis.louisaslett.com/reference/CacheSummary.md) for the
number of haplotypes cached;
[`L()`](https://kalis.louisaslett.com/reference/CacheSummary.md) for the
number of variants cached.

## Examples

``` r
if (FALSE) { # \dontrun{
# If starting from a VCF/BCF first use bcftools to convert to
# HAP/SAMPLE/LEGEND format (bcftools can take in several starting formats)
# See http://samtools.github.io/bcftools/bcftools.html#convert
system("bcftools convert -h my.vcf.gz")
CacheHaplotypes("my.hap.gz")
CacheSummary()
} # }

# If starting directly from a hap.gz file on disk (HAP/LEGEND/SAMPLE format)
if (FALSE) { # \dontrun{
CacheHaplotypes("my.hap.gz")
} # }
# For example, to load the mini example built into the package:
CacheHaplotypes(system.file("small_example/small.hap.gz", package = "kalis"))
#> Warning: haplotypes already cached ... overwriting existing cache.
CacheSummary()
#> Cache currently loaded with 300 haplotypes, each with 400 variants. 
#>   Memory consumed: 16 kB. 


# If starting from an HDF5 file on disk
if (FALSE) { # \dontrun{
CacheHaplotypes("my.h5")
} # }
# For example, to load the mini example built into the package:
CacheHaplotypes(system.file("small_example/small.h5", package = "kalis"))
#> Warning: haplotypes already cached ... overwriting existing cache.
CacheSummary()
#> Cache currently loaded with 300 haplotypes, each with 400 variants. 
#>   Memory consumed: 16 kB. 


# If CacheSummary() indicates that the numbers of haplotypes and variants are
# the wrong way around, reload with argument transpose set to TRUE
if (FALSE) { # \dontrun{
CacheHaplotypes("myhaps.h5", transpose = TRUE)
CacheSummary()
} # }


# Alternatively, if you have an exotic file format that can be loaded in to R
# by other means, then a binary matrix can be supplied.  This example
# randomly simulates a binary matrix to illustrate.
n.haps <- 100
n.vars <- 200
haps <- matrix(sample(0:1, n.haps*n.vars, replace = TRUE),
               nrow = n.vars, ncol = n.haps)
CacheHaplotypes(haps)
#> Warning: haplotypes already cached ... overwriting existing cache.
# For example, to load the mini example built into the package:
data("SmallHaps")
CacheHaplotypes(SmallHaps)
#> Warning: haplotypes already cached ... overwriting existing cache.
```
