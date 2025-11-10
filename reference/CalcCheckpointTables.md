# Calculate Checkpoint Tables

Calculate look up tables for solving optimal checkpointing problems with
dynamic programming.

## Usage

``` r
CalcCheckpointTables(propagation.cost, max.num.checkpoints, use.R = FALSE)
```

## Arguments

- propagation.cost:

  a non-negative vector such that `propagation.cost[i]` gives the
  relative amount of time or cost required to propagate `i` steps.

- max.num.checkpoints:

  the maximum number of checkpoints that should be considered when
  building the checkpoint table.

- use.R:

  a logical, when `TRUE` use base R rather than C implementation of
  table building. Defaults to `FALSE`.

## Value

A list containing:

- `cost`:

  the matrix \\F\\ in Christ et al. (2024)

- `index`:

  the matrix \\H\\ in Christ et al. (2024)

## References

Christ, R.R., Wang, X., Aslett, L.J.M., Steinsaltz, D. and Hall, I.
(2024) "Clade Distillation for Genome-wide Association Studies", bioRxiv
2024.09.30.615852. Available at:
[doi:10.1101/2024.09.30.615852](https://doi.org/10.1101/2024.09.30.615852)
.
