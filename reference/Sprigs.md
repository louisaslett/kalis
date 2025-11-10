# Greedy sprig calling based on nearest neighbourhoods

Infer sprigs (very small clades) within a local phylogeny based on
one-directional nearest neighbourhoods assigned to each haplotype.

## Usage

``` r
Sprigs(x, old.sprigs = FALSE)
```

## Arguments

- x:

  list of integers where each entry is a nearest neighbourhood as
  returned by
  [`CladeMat()`](https://kalis.louisaslett.com/reference/CladeMat.md)

- old.sprigs:

  If `TRUE`, use an earlier (undocumented) version of maximal clique
  calling, defaults to `FALSE`.

## Value

a list containing:

- `assignments`:

  a vector of integers such that `$assignments[i]` gives the sprig to
  which haplotype `i` belongs. `NA` if a haplotype was not assigned to a
  sprig.

- `to.prune`:

  a vector of logicals, `$to.prune[i]==TRUE` when haplotype i has been
  assigned to a sprig

- `num.sprigs`:

  total number of sprigs calls, equal to
  `max($assignements, rm.na = TRUE)`

## Details

Call maximal cliques within a directed graph where edges correspond to
nearest neighbour relationships.

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
