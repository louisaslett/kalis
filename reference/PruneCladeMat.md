# Prune called sprigs or singletons from inferred clade matrix

Use haplotype nearest neighbourhoods and other information returned by
[`CladeMat()`](https://kalis.louisaslett.com/reference/CladeMat.md) to
efficiently remove structure corresponding to singletons or called
sprigs from clade matrix 'M'

## Usage

``` r
PruneCladeMat(M, neigh, sprigs, prune = "singleton.info", from_recipient = 1L)
```

## Arguments

- M:

  clade matrix such as that updated by
  [`CladeMat()`](https://kalis.louisaslett.com/reference/CladeMat.md)

- neigh:

  a list of nearest neighborhoods as retuned by
  [`CladeMat()`](https://kalis.louisaslett.com/reference/CladeMat.md)

- sprigs:

  a sprigs object as returned by
  [`Sprigs()`](https://kalis.louisaslett.com/reference/Sprigs.md)

- prune:

  a character indicating the type of information to be removed from the
  [`CladeMat()`](https://kalis.louisaslett.com/reference/CladeMat.md).
  See Details.

- from_recipient:

  haplotype index at which to start trace calculation — useful for
  distributed computation (experimental feature, more documentation to
  come)

## Value

There is nothing returned.

**NOTE:** for performance reasons, `M` is updated in-place.

## References

Christ, R.R., Wang, X., Aslett, L.J.M., Steinsaltz, D. and Hall, I.
(2024) "Clade Distillation for Genome-wide Association Studies", bioRxiv
2024.09.30.615852. Available at:
[doi:10.1101/2024.09.30.615852](https://doi.org/10.1101/2024.09.30.615852)
.

## Examples

``` r
# TODO
```
