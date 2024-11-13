#' Small example phased haplotype dataset and recombination map
#'
#' Simulated dataset with \eqn{N=300} haplotypes at \eqn{L=400} variants generated using msprime
#' (Kelleher et al., 2016), together with the recombination map.
#'
#' @docType data
#'
#' @usage data("SmallHaps")
#'
#' @format
#'   For `SmallHaps`, an object of class `matrix` of dimensions \eqn{400 \times 300}{400 x 300} where each column is a simulated haplotype.
#'   This is suitable for passing directly to [CacheHaplotypes()].
#'
#' @keywords datasets
#'
#' @references
#' Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation of the Li & Stephens model for local ancestry inference in R", *BMC Bioinformatics*, **25**(1). Available at: \doi{10.1186/s12859-024-05688-8}.
#'
#' Kelleher, J., Etheridge, A.M. and McVean, G. (2016) "Efficient coalescent simulation and genealogical analysis for large sample sizes", *PLoS Computational Biology*, **12**(5). Available at: \doi{10.1371/journal.pcbi.1004842}.
#'
#' @examples
#' data("SmallHaps")
#'
#' \donttest{
#' # Plot Allele Frequencies
#' hist(rowMeans(SmallHaps),breaks=20)
#' }
#'
#' # Import into kalis cache directly ...
#' CacheHaplotypes(SmallHaps)
#'
#' data("SmallMap")
#'
#' # Find parameters
#' pars <- Parameters(CalcRho(diff(SmallMap)))
#' pars
#'
"SmallHaps"

#' @rdname SmallHaps
#' @usage data("SmallMap")
#' @format
#'   For `SmallMap`, a vector of length 400 representing the recombination map for the `SmallHaps` data.
#'   This can be used with [CalcRho()], by converting to recombination distances using `diff(SmallMap)`.
"SmallMap"
