# I/O for haplotype matrices in HDF5 files

Reads/writes an R matrix of 0/1s to the HDF5 format which is used for
reading to the kalis optimised memory cache. If you're working with a
large haplotype dataset, we recommend that you convert it directly to
this HDF5 format (see vignette) rather than read it into R.

## Usage

``` r
WriteHaplotypes(
  hdf5.file,
  haps,
  hap.ids = NA,
  loci.ids = NA,
  haps.name = "/haps",
  hap.ids.name = "/hap.ids",
  loci.ids.name = "/loci.ids",
  append = FALSE
)

ReadHaplotypes(
  hdf5.file,
  loci.idx = NA,
  hap.idx = NA,
  loci.ids = NA,
  hap.ids = NA,
  haps.name = "/haps",
  loci.ids.name = "/loci.ids",
  hap.ids.name = "/hap.ids",
  transpose = FALSE
)
```

## Arguments

- hdf5.file:

  the name of the file which the haplotypes are to be written to.

- haps:

  a vector or a matrix where each column is a haplotype to be stored in
  the file `hdf5.file`.

- hap.ids:

  a character vector naming haplotypes when writing, or which haplotypes
  are to be read.

- loci.ids:

  a character vector naming variants when writing, or which variants are
  to be read.

- haps.name:

  a string providing the full path and object name where the haplotype
  matrix should be read/written.

- hap.ids.name:

  a string providing the full path and object name where the haplotype
  names (in `haps.ids`) should be read/written.

- loci.ids.name:

  a string providing the full path and object name where the variant
  names (in `loci.ids`) should be read/written.

- append:

  a logical indicating whether overwrite (default) or append to an
  existing `haps` dataset if it already exists in `hdf5.file`.

- loci.idx:

  an integer vector of the indices of which variants are to be read (for
  naming, use `hap.ids`).

- hap.idx:

  an integer vector of the indices of which haplotypes are to be read
  (for naming, use `hap.ids`).

- transpose:

  a logical indicating whether to transpose the logic of
  haplotypes/variants when reading.

## Value

`WriteHaplotypes` does not return anything.

`ReadHaplotypes` returns a binary matrix containing the haplotypes that
were specified in `ids`.

## Details

The primary method to load data into kalis' internal optimised cache is
from an HDF5 storage file. If the user has a collection of haplotypes
already represented as a matrix of 0's and 1's in R, this function can
be used to write to HDF5 in the format required to load into cache.

kalis expects a 2-dimensional object named `haps` at the root level of
the HDF5 file. Haplotypes should be stored in the slowest changing
dimension as defined in the HDF5 specification (note that different
languages treat this as rows or columns).

Note that if `hdf5.file` exists but does not contain a dataset named
`haps`, then `WriteHaplotypes` will simply create a `haps` dataset
within the existing file.

## References

Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation
of the Li & Stephens model for local ancestry inference in R", *BMC
Bioinformatics*, **25**(1). Available at:
[doi:10.1186/s12859-024-05688-8](https://doi.org/10.1186/s12859-024-05688-8)
.

## See also

[`CacheHaplotypes()`](https://kalis.louisaslett.com/reference/CacheHaplotypes.md)
to fill the kalis cache with haplotypes.

## Examples

``` r
# \donttest{
# Generate a random mini set of haplotypes to write
n.haps <- 20
n.vars <- 200
haps <- matrix(sample(0:1, n.haps*n.vars, replace = TRUE),
               nrow = n.vars, ncol = n.haps)

# ... write them to a file, giving alphabetic letters "A" through "T" as the
# haplotype names ...
WriteHaplotypes("~/myhaps.h5", haps, hap.ids = LETTERS[1:20])
#> Creating HDF5 file ...
#> Writing 20 haplotype(s) of size 200 ...

# ... and confirm we can read a chosen portion back.  Try to read back
# the 10th and 11th haplotypes by using their name (J and K are 10th and 11th
# letter of the alphabet)
h5 <- ReadHaplotypes("~/myhaps.h5", hap.ids = c("J","K"))
all(h5$haps == haps[, 10:11])
#> [1] TRUE

# Read from the .h5 file into the kalis cache and confirm that what we wrote
# out to the HDF5 file matches the original matrix we generated in R
CacheHaplotypes("~/myhaps.h5")
#> Warning: haplotypes already cached ... overwriting existing cache.
all(haps == QueryCache())
#> [1] TRUE
# }
```
