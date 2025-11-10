# Plotting function for a distance matrix

Clusters the given distance matrix and generates a heatmap to display
it.

## Usage

``` r
# S3 method for class 'kalisDistanceMatrix'
plot(x, cluster.method = "average", ...)
```

## Arguments

- x:

  a distance matrix, such as returned by
  [`DistMat()`](https://kalis.louisaslett.com/reference/DistMat.md).

- cluster.method:

  the agglomeration method to be used, which is passed to the
  [`fastcluster::hclust()`](https://rdrr.io/pkg/fastcluster/man/hclust.html)
  function. This must be (an unambiguous abbreviation of) one of
  "single", "complete", "average", "mcquitty", "ward.D", "ward.D2",
  "centroid" or "median".

- ...:

  further arguments to be passed on to the underlying
  [`lattice::levelplot()`](https://rdrr.io/pkg/lattice/man/levelplot.html)
  plotting function.

## Value

No return value, called for side effects.

## References

Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation
of the Li & Stephens model for local ancestry inference in R", *BMC
Bioinformatics*, **25**(1). Available at:
[doi:10.1186/s12859-024-05688-8](https://doi.org/10.1186/s12859-024-05688-8)
.

## Examples

``` r
# TODO
```
