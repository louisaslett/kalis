#' SmallHaps Example Phased Haplotype Dataset for kalis
#'
#' Simulated dataset with N=100 haplotypes at L=500 loci generated using msprime
#' (Kelleher et al., 2016).
#'
#' @docType data
#'
#' @usage data(SmallHaps)
#'
#' @format
#'   An object of class \code{matrix} where each column is a
#'   simulated haplotype.
#'   This is suitable for passing directly to \code{\link{CacheHaplotypes}}.
#'
#' @keywords datasets
#'
#' @references
#'   Kelleher, J., Etheridge, A. M., & McVean, G. (2016). Efficient coalescent
#'   simulation and genealogical analysis for large sample sizes.
#'   *PLoS computational biology*, **12**(5).
#'
#' @examples
#' \dontrun{
#' data(SmallHaps)
#' # Plot Allele Frequencies
#' hist(rowMeans(SmallHaps),breaks=20)
#'
#' # Store in kalis hdf5 format
#' WriteIndividualHaplotypeH5("SmallHaps.h5", SmallHaps)
#'
#' # Import into kalis cache directly ...
#' CacheHaplotypes(SmallHaps)
#' # ... or via the HDF5 file written out above
#' CacheHaplotypes("SmallHaps.h5")
#' }
#'
"SmallHaps"