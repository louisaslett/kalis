# Posterior marginal probabilities

Calculate the posterior marginal probabilities at a given variant using
forward and backward tables propagated to that position.

## Usage

``` r
PostProbs(
  fwd,
  bck,
  unif.on.underflow = FALSE,
  M = NULL,
  beta.theta.opts = NULL,
  nthreads = min(parallel::detectCores(logical = FALSE), fwd$to_recipient -
    fwd$from_recipient + 1)
)
```

## Arguments

- fwd:

  a forward table as returned by
  [`MakeForwardTable()`](https://kalis.louisaslett.com/reference/MakeForwardTable.md)
  and propagated to a target variant by
  [`Forward()`](https://kalis.louisaslett.com/reference/Forward.md).
  Must be at the same variant as `bck` (unless `bck` is in "beta-theta
  space" in which case if must be downstream ... see
  [`Backward()`](https://kalis.louisaslett.com/reference/Backward.md)
  for details).

- bck:

  a backward table as returned by
  [`MakeBackwardTable()`](https://kalis.louisaslett.com/reference/MakeBackwardTable.md)
  and propagated to a target variant by
  [`Backward()`](https://kalis.louisaslett.com/reference/Backward.md).
  Must be at the same variant as `fwd` (unless `bck` is in "beta-theta
  space" in which case if must be downstream ... see
  [`Backward()`](https://kalis.louisaslett.com/reference/Backward.md)
  for details).

- unif.on.underflow:

  a logical; if `TRUE`, then if all probabilities in a column underflow,
  then they will be set to \\\frac{1}{N-1}\\ instead of 0

- M:

  a pre-existing matrix into which to write the probabilities, can yield
  substantial speed up but requires special attention (see Details)

- beta.theta.opts:

  a list; see Details.

- nthreads:

  the number of CPU cores to use. By default uses the `parallel` package
  to detect the number of physical cores.

## Value

Matrix of posterior marginal probabilities. The \\(j,i)\\-th element of
of the returned matrix is the probability that \\j\\ is copied by \\i\\
at the current variant, \\l\\, of the two tables, given the haplotypes
observed (over the whole sequence).

## Details

The forward and backward tables must be at the same variant in order for
them to be combined to yield the posterior marginal probabilities at
variant \\l\\. The \\(j,i)\\-th element of of the returned matrix is the
probability that \\j\\ is copied by \\i\\ at the current variant, \\l\\,
of the two tables, given the haplotypes observed (over the whole
sequence).

Note that each column represents an independent HMM.

By convention, every diagonal element is zero.

**Notes on `beta.theta.opts`**

In order to obtain posterior marginal probability matrices between
variants `fwd$l` and `bck$l`, then `bck` must be in "beta-theta space",
see [`Backward()`](https://kalis.louisaslett.com/reference/Backward.md)
for details. This allows the forward and backward tables to be
transitioning both tables to some genomic position between `fwd$l` and
`bck$l`. The precise recombination distance by which each table is
propagated can be controlled by passing optional arguments in a list via
`beta.theta.opts`. The recombination distances used can be specified in
one of two ways.

1.  Manually. In this case, `beta.theta.opts` is a list containing two
    named elements:

    - `"rho.fwd"` \\\in (0,1)\\ specifies the transition probability
      `rho` for propagating the forward table.

    - `"rho.bck"` \\\in (0,1)\\ specifies the transition probability
      `rho` for propagating the backward table.

2.  Implicitly. In this case, `beta.theta.opts` is a list containing two
    named elements:

- `"pars"`: a `kalisParameters` object that implicitly defines the
  recombination distance \\\rho^\star\\ between `fwd$l` and `bck$l`

- `"bias"` \\\in (0,1)\\. The forward table is propagated a distance of
  `(bias)`\\\rho^\star\\ and the backward table is propagated a distance
  of `(1-bias)`\\\rho^\star\\.

**Performance notes**

If calculating many posterior probability matrices in succession,
providing a pre-existing matrix `M` that can be updated in-place can
dramatically increase speed by eliminating the time needed for memory
allocation. Be warned, since the matrix is updated in-place, if any
other variables point to the same memory address, they will also be
simultaneously overwritten. For example, writing

    M <- matrix(0, nrow(fwd$alpha), ncol(fwd$alpha))
    P <- M
    PostProbs(fwd, bck, M = M)

will update `M` and `P` simultaneously.

When provided, `M` must have dimensions matching that of `fwd$alpha`.
Typically, that is simply \\N \times N\\ for \\N\\ haplotypes. However,
if kalis is being run in a distributed manner, `M` will be a \\N \times
R\\ matrix where \\R\\ is the number of recipient haplotypes on the
current machine.

## References

Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation
of the Li & Stephens model for local ancestry inference in R", *BMC
Bioinformatics*, **25**(1). Available at:
[doi:10.1186/s12859-024-05688-8](https://doi.org/10.1186/s12859-024-05688-8)
.

## See also

[`DistMat()`](https://kalis.louisaslett.com/reference/DistMat.md) to
generate calculate \\d\_{ji}\\ distances directly;
[`Forward()`](https://kalis.louisaslett.com/reference/Forward.md) to
propagate a Forward table to a new variant;
[`Backward()`](https://kalis.louisaslett.com/reference/Backward.md) to
propagate a Backward table to a new variant

## Examples

``` r
# To get the posterior probabilities at, say, variants 100 of the toy data
# built into kalis
data("SmallHaps")
data("SmallMap")

CacheHaplotypes(SmallHaps)
#> Warning: haplotypes already cached ... overwriting existing cache.

rho <- CalcRho(diff(SmallMap))
pars <- Parameters(rho)

fwd <- MakeForwardTable(pars)
bck <- MakeBackwardTable(pars)

Forward(fwd, pars, 100)
Backward(bck, pars, 100)

p <- PostProbs(fwd, bck)
d <- DistMat(fwd, bck)
```
