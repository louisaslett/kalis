#' SmallHaps Example Phased Haplotype Dataset for kalis
#'
#' Simulated dataset with N=100 haplotypes at L=500 loci generated using msprime (see reference).
#'
#' @docType data
#'
#' @usage data(SmallHaps)
#'
#' @format An object of class \code{"matrix"} where each column is a haplotype.
#'
#' @keywords datasets
#'
#' @references J Kelleher, et al. (2016) PLoS Comput Biol 12(5): e1004842.
#' (\href{https://doi.org/10.1371/journal.pcbi.1004842}{Article})
#'
#'
#' @examples
#' \dontrun{
#' data(SmallHaps)
#' # Plot Allele Frequencies
#' hist(rowMeans(SmallHaps),breaks=20)
#'
#' # Store in kalis hdf5 format
#' WriteIndividualHaplotypeH5("SmallHaps.h5",SmallHaps)
#'
#' # Import into kalis cache
#' CacheAllHaplotypesH5("SmallHaps.h5")
#' }
#'
"SmallHaps"