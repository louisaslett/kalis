# Calculate recombination probabilities

Calculate the recombination probabilities, rho, from a given
recombination map for the haplotype data currently in the kalis cache.

## Usage

``` r
CalcRho(cM = 0, s = 1, gamma = 1, floor = TRUE)
```

## Arguments

- cM:

  a vector specifying the recombination distance between variants in
  centimorgans. Note element i of this vector should be the distance
  between variants `i` and `i+1` (not `i` and `i-1`), and thus length
  one less than the number of variants. This can be easily obtained by
  applying [`diff()`](https://rdrr.io/r/base/diff.html) to a
  recombination map 'CDF'. By default, recombination probabilities are
  zero, meaning that all variants are perfectly linked.

- s:

  a scalar multiplier on the recombination map (related to effective
  population size).

- gamma:

  a scalar power to which the Morgan distances are raised.

- floor:

  if `TRUE` (default) any recombination probabilities below machine
  precision (`1e-16`) will be zeroed out. If `FALSE` raw recombination
  probabilities will be preserved.

## Value

A vector of recombination probabilities which can be used as the `rho`
argument to the
[`Parameters()`](https://kalis.louisaslett.com/reference/Parameters.md)
function when creating a parameter set.

## Details

This is a utility function to calculate the recombination probabilities
rho (the Li and Stephens hidden Markov model transition probabilities)
from a recombination map/distances.

**NOTE:** the corresponding haplotype data *must* have already been
inserted into the kalis cache by a call to
[`CacheHaplotypes()`](https://kalis.louisaslett.com/reference/CacheHaplotypes.md),
since this function performs checks to confirm the dimensionality
matches.

See page 3 in Supplemental Information for the original ChromoPainter
paper (Lawson et al., 2012) for motivation behind our parameterisation,
which is as follows:

\$\$\rho = 1 - \exp(-s \times cM^\gamma)\$\$

For a complete description, see the main kalis paper, Aslett and Christ
(2024).

## References

Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation
of the Li & Stephens model for local ancestry inference in R", *BMC
Bioinformatics*, **25**(1). Available at:
[doi:10.1186/s12859-024-05688-8](https://doi.org/10.1186/s12859-024-05688-8)
.

Lawson, D.J., Hellenthal, G., Myers, S. and Falush, D. (2012) "Inference
of population structure using dense haplotype data", *PLoS genetics*,
**8**(1). Available at:
[doi:10.1371/journal.pgen.1002453](https://doi.org/10.1371/journal.pgen.1002453)
.

## See also

[`Parameters()`](https://kalis.louisaslett.com/reference/Parameters.md)
to use the resulting recombination probabilities to construct a
`kalisParameters` object.

## Examples

``` r
# Load the mini example data and recombination map from the package built-in #' # dataset
data("SmallHaps")
data("SmallMap")

# Load haplotypes into cache
CacheHaplotypes(SmallHaps)
#> Warning: haplotypes already cached ... overwriting existing cache.

# Compute the recombination probabilities, leaving s and gamma as default
rho <- CalcRho(diff(SmallMap))
rho
#>   [1] 3.072298e-02 1.883221e-03 6.798092e-06 1.197364e-05 6.537176e-05
#>   [6] 1.977829e-06 1.764480e-06 1.928517e-04 3.333315e-05 4.866593e-04
#>  [11] 4.649365e-06 8.483164e-06 6.170919e-04 3.330155e-04 1.082906e-03
#>  [16] 8.186886e-06 7.548185e-07 1.058505e-05 6.291293e-04 6.527013e-04
#>  [21] 1.425614e-04 2.766100e-03 1.809975e-04 1.025515e-03 9.376016e-06
#>  [26] 1.366644e-04 9.586410e-05 1.667846e-06 3.597911e-05 9.346260e-06
#>  [31] 1.386988e-05 9.455863e-06 8.053657e-04 0.000000e+00 1.097848e-04
#>  [36] 1.706280e-04 1.131905e-03 1.698680e-06 5.592445e-05 9.441074e-04
#>  [41] 2.994247e-06 0.000000e+00 1.232414e-05 1.969309e-06 1.635813e-05
#>  [46] 5.673478e-05 2.093181e-05 1.423591e-05 3.106822e-05 6.190222e-06
#>  [51] 8.187811e-06 1.948645e-03 9.612262e-05 0.000000e+00 1.386898e-05
#>  [56] 0.000000e+00 2.374524e-05 1.158938e-02 6.207018e-07 1.245026e-06
#>  [61] 2.714202e-06 5.819492e-06 3.178513e-07 8.887262e-08 5.539737e-07
#>  [66] 2.195520e-05 0.000000e+00 5.382504e-07 2.416272e-05 4.578088e-07
#>  [71] 8.066226e-07 1.780718e-05 0.000000e+00 1.114559e-06 6.387084e-07
#>  [76] 3.355397e-07 3.341899e-05 2.370182e-05 9.361219e-07 3.121628e-06
#>  [81] 0.000000e+00 1.380908e-06 3.981544e-07 1.521412e-05 6.088880e-05
#>  [86] 1.829450e-06 1.745363e-05 2.690787e-04 4.737639e-06 1.236435e-06
#>  [91] 7.144977e-06 4.596637e-05 7.501879e-07 3.207836e-07 8.624930e-05
#>  [96] 2.330897e-06 1.065604e-07 4.962100e-06 1.729355e-06 5.080312e-07
#> [101] 1.034343e-04 4.158098e-07 6.029703e-07 9.577280e-05 2.072616e-04
#> [106] 1.729982e-06 0.000000e+00 1.743994e-04 5.596931e-06 2.077639e-04
#> [111] 0.000000e+00 4.275686e-07 2.199119e-04 4.466601e-06 4.042581e-06
#> [116] 8.883807e-05 2.001860e-06 1.445364e-06 7.091782e-06 1.843242e-06
#> [121] 7.805524e-08 1.987546e-06 1.778252e-05 5.467262e-07 0.000000e+00
#> [126] 2.070603e-04 4.999370e-08 8.507876e-06 1.999806e-05 7.009794e-05
#> [131] 3.312554e-06 4.301399e-07 7.706821e-05 3.026263e-04 6.900658e-05
#> [136] 1.899490e-06 6.199921e-06 2.736548e-06 9.220733e-07 9.314320e-07
#> [141] 1.663724e-06 2.209034e-05 1.352584e-06 4.297543e-07 1.989024e-07
#> [146] 1.033614e-04 2.183316e-04 4.810175e-06 5.215310e-06 1.189485e-05
#> [151] 1.249872e-06 1.312833e-04 1.788202e-06 6.841195e-07 3.491718e-05
#> [156] 5.630717e-07 1.237170e-07 8.678187e-07 0.000000e+00 6.413444e-07
#> [161] 8.816561e-07 9.823999e-07 2.463841e-07 4.682670e-05 1.668122e-06
#> [166] 2.154239e-06 4.765185e-05 1.642860e-04 3.892054e-06 8.037314e-07
#> [171] 4.042199e-06 8.785437e-05 1.235422e-05 5.410768e-07 5.245492e-06
#> [176] 4.105771e-07 0.000000e+00 1.241292e-04 1.060935e-05 5.366691e-06
#> [181] 1.407887e-04 3.162847e-06 1.119661e-05 8.036889e-07 5.227231e-06
#> [186] 3.987610e-07 1.889812e-07 1.043175e-04 2.765477e-05 2.861844e-04
#> [191] 1.287753e-05 1.722533e-06 8.253517e-06 1.182719e-06 7.600410e-07
#> [196] 2.537284e-04 3.357894e-07 7.303595e-07 8.226384e-04 2.071366e-05
#> [201] 1.697994e-04 7.199889e-06 0.000000e+00 3.181682e-05 1.647047e-04
#> [206] 2.340387e-05 8.627927e-07 2.843545e-04 2.263947e-04 1.707134e-04
#> [211] 1.979789e-05 1.045979e-05 4.116453e-04 2.500143e-05 2.502130e-05
#> [216] 7.333162e-04 9.264360e-06 1.104631e-05 1.500975e-05 4.899939e-05
#> [221] 2.003711e-07 7.709842e-05 5.405506e-06 6.700180e-07 4.089772e-05
#> [226] 3.553471e-05 9.620512e-05 1.021574e-04 3.473059e-06 5.032733e-06
#> [231] 3.018856e-05 1.570943e-06 1.538014e-06 3.153386e-05 5.634462e-06
#> [236] 2.271145e-06 2.180938e-05 1.426748e-07 1.200679e-07 1.120902e-05
#> [241] 2.057690e-06 4.751340e-06 9.565998e-07 3.448003e-07 7.817790e-06
#> [246] 6.618692e-05 2.474977e-06 6.714171e-07 1.490395e-04 1.533094e-05
#> [251] 4.255423e-05 2.259110e-06 1.879833e-06 1.136652e-06 6.267087e-05
#> [256] 1.612992e-05 9.213218e-07 1.540307e-04 4.434347e-05 1.713461e-04
#> [261] 5.019734e-05 5.248481e-04 9.425490e-05 3.401882e-06 5.815356e-05
#> [266] 2.681589e-05 3.080475e-05 8.073736e-06 2.251303e-04 2.033318e-04
#> [271] 1.669898e-05 2.291142e-06 2.677244e-06 7.340165e-06 0.000000e+00
#> [276] 7.362316e-06 6.189150e-05 1.219512e-07 1.540678e-06 3.029898e-07
#> [281] 2.478905e-07 1.968149e-05 1.026901e-04 0.000000e+00 4.076519e-06
#> [286] 4.511654e-05 5.770477e-06 3.796285e-05 2.070229e-04 1.900436e-05
#> [291] 1.272606e-04 8.449114e-07 2.451426e-05 2.294817e-04 1.206330e-05
#> [296] 6.293434e-05 5.150321e-04 1.564527e-05 4.724179e-06 1.890515e-06
#> [301] 4.621392e-05 9.506681e-07 4.164448e-06 2.516349e-04 1.268009e-06
#> [306] 5.715582e-07 9.048111e-07 4.275813e-06 8.560334e-08 1.888363e-04
#> [311] 6.132239e-07 2.739148e-07 3.742452e-07 3.828106e-06 1.366745e-05
#> [316] 3.776691e-06 4.477607e-07 2.787425e-06 2.406681e-06 1.070108e-05
#> [321] 2.248604e-05 1.166044e-05 9.383598e-06 5.497293e-05 6.959451e-04
#> [326] 1.024863e-05 3.931442e-05 1.049756e-04 0.000000e+00 7.245571e-07
#> [331] 7.967463e-05 2.471565e-05 9.505406e-05 6.543679e-06 0.000000e+00
#> [336] 1.142243e-04 2.811166e-05 4.018208e-07 3.753972e-05 3.781747e-06
#> [341] 2.033923e-05 2.285410e-07 4.755785e-08 8.195639e-07 2.838863e-04
#> [346] 5.748447e-06 9.565791e-07 4.381465e-05 3.580701e-07 4.789086e-07
#> [351] 3.817816e-06 1.064850e-04 5.310647e-07 0.000000e+00 1.525290e-05
#> [356] 1.726751e-04 6.388522e-05 8.799876e-05 5.533851e-06 4.455855e-05
#> [361] 5.480269e-04 5.857182e-05 1.531392e-04 1.125972e-04 1.014507e-04
#> [366] 2.229192e-06 2.998638e-06 9.328138e-06 5.902549e-06 1.406753e-04
#> [371] 1.843086e-05 3.556065e-05 8.552933e-05 1.286796e-04 3.244456e-06
#> [376] 3.571614e-06 1.137062e-04 5.807801e-06 4.244475e-05 2.822669e-05
#> [381] 2.453923e-05 1.481526e-06 3.326299e-06 1.803912e-04 3.336949e-06
#> [386] 2.781070e-06 5.636435e-05 7.262678e-07 1.604063e-06 4.337679e-05
#> [391] 6.282122e-07 3.600771e-05 2.134311e-06 2.632778e-06 9.781323e-07
#> [396] 1.129465e-06 7.535395e-07 1.568985e-05 1.792468e-05
```
