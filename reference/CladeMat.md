# Fast clade matrix construction

Constructs a clade matrix using forward and backward tables. The clade
matrix captures genetic relatedness information in the distances from
the Li & Stephens model that are not captured in the called clades.

## Usage

``` r
CladeMat(
  fwd,
  bck,
  M,
  unit.dist,
  thresh = 0.2,
  max1var = FALSE,
  nthreads = min(parallel::detectCores(logical = FALSE), fwd$to_recipient -
    fwd$from_recipient + 1)
)
```

## Arguments

- fwd:

  a `kalisForwardTable` object, as returned by
  [`MakeForwardTable()`](https://kalis.louisaslett.com/reference/MakeForwardTable.md)
  and propagated to a target variant by
  [`Forward()`](https://kalis.louisaslett.com/reference/Forward.md).
  This table must be at the same variant location as argument `bck`.

- bck:

  a `kalisBackwardTable` object, as returned by
  [`MakeBackwardTable()`](https://kalis.louisaslett.com/reference/MakeBackwardTable.md)
  and propagated to a target variant by
  [`Backward()`](https://kalis.louisaslett.com/reference/Backward.md).
  This table must be at the same variant location as argument `fwd`.

- M:

  a matrix with half the number of rows and columns as the corresponding
  forward/backward tables. This matrix is overwritten in place with the
  clade matrix result for performance reasons.

- unit.dist:

  the change in distance that is expected to correspond to a single
  mutation (typically \\-\log(\mu)\\) for the LS model)

- thresh:

  a regularization parameter: differences of distances must exceed this
  threshold (in `unit.dist` units) in order to cause the introduction of
  a probabilistic clade. Defaults to `0.2`.

- max1var:

  a logical regularization parameter. When `TRUE`, differences in
  distances exceeding 1 `unit.dist` are set to 1 (so that any edge in
  the latent ancestral tree with multiple mutations on them are treated
  as if only one mutation was on it).

- nthreads:

  the number of CPU cores to use. By default uses the `parallel` package
  to detect the number of physical cores.

## Value

A list, the first element contains a list of tied nearest neighbours
(one for each haplotype). Other elements of the returned list are for
internal use by
[`PruneCladeMat()`](https://kalis.louisaslett.com/reference/PruneCladeMat.md)
to allow for efficient removal of singletons and sprigs.

## Details

`CladeMat()` uses the forward and backward tables to construct the
corresponding clade matrix which can then be tested, for example using a
standard quadratic form score statistic.

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
