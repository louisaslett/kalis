% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SmallHaps-data.R
\docType{data}
\name{SmallHaps}
\alias{SmallHaps}
\alias{SmallMap}
\title{Small example phased haplotype dataset and recombination map}
\format{
For \code{SmallHaps}, an object of class \code{matrix} of dimensions \eqn{400 \times 300}{400 x 300} where each column is a simulated haplotype.
This is suitable for passing directly to \code{\link[=CacheHaplotypes]{CacheHaplotypes()}}.

For \code{SmallMap}, a vector of length 400 representing the recombination map for the \code{SmallHaps} data.
This can be used with \code{\link[=CalcRho]{CalcRho()}}, by converting to recombination distances using \code{diff(SmallMap)}.
}
\usage{
data("SmallHaps")

data("SmallMap")
}
\description{
Simulated dataset with \eqn{N=300} haplotypes at \eqn{L=400} variants generated using msprime
(Kelleher et al., 2016), together with the recombination map.
}
\examples{
data("SmallHaps")

\donttest{
# Plot Allele Frequencies
hist(rowMeans(SmallHaps),breaks=20)
}

# Import into kalis cache directly ...
CacheHaplotypes(SmallHaps)

data("SmallMap")

# Find parameters
pars <- Parameters(CalcRho(diff(SmallMap)))
pars

}
\references{
Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation of the Li & Stephens model for local ancestry inference in R", \emph{BMC Bioinformatics}, \strong{25}(1). Available at: \doi{10.1186/s12859-024-05688-8}.

Kelleher, J., Etheridge, A.M. and McVean, G. (2016) "Efficient coalescent simulation and genealogical analysis for large sample sizes", \emph{PLoS Computational Biology}, \strong{12}(5). Available at: \doi{10.1371/journal.pcbi.1004842}.
}
\keyword{datasets}
