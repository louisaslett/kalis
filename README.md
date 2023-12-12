# kalis: High Performance Li &amp; Stephens Local Ancestry Inference

kalis provides a high performance implementation of the Li & Stephens model (doi: [10.1093/genetics/165.4.2213](https://doi.org/10.1093/genetics/165.4.2213)) for local ancestry inference (local referring to a region of the genome).
For a set of N phased haplotypes, kalis computes the posterior marginal probability of each haplotype copying every other haplotype by running N hidden Markov models in parallel.
This yields an N x N matrix that summarizes the recent local ancestry at each variant of interest.
The package provides functionality for specifying a recombination map, site-specific mutation rates, and differing prior copying probabilities for each recipient haplotype.
Extensive use is made of low level threading and CPU vector instructions.

## Installation

kalis will appear on CRAN in due course.
For now, the current development version can be installed using the [remotes](https://github.com/r-lib/remotes) package as follows:

```
install.packages("remotes")
remotes::install_github("louisaslett/kalis", build_vignettes = TRUE)
```

Note that kalis uses various low-level optimisations meaning that you should ensure the compiler is targeting your local CPU architecture.
The simplest way to do this is to pass configure variables setting the correct `CFLAGS` at install time:

```
remotes::install_github("louisaslett/kalis", build_vignettes = TRUE,
  configure.vars = c(kalis = "PKG_CFLAGS='-march=native -mtune=native -O3'"))
```

If you do not have the correct flags set, kalis will fall back to an implementation which does not use the special vector instruction set architecture of your CPU and will provide a warning when you load the package in your R session.
