% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FB.R
\name{Backward}
\alias{Backward}
\title{Propagate a backward table}
\usage{
Backward(
  bck,
  pars,
  t = bck$l - 1,
  nthreads = min(parallel::detectCores(logical = FALSE), bck$to_recipient -
    bck$from_recipient + 1),
  beta.theta = FALSE
)
}
\arguments{
\item{bck}{a \code{kalisBackwardTable} object, as returned by
\code{\link[=MakeBackwardTable]{MakeBackwardTable()}}.}

\item{pars}{a \code{kalisParameters} object, as returned by
\code{\link[=Parameters]{Parameters()}}.}

\item{t}{a variant to move the backward table to.
Must be less than or equal to current variant of \code{bck}.
By default, it advances to the next variant upstream (or if an uninitialised table, to the last variant).}

\item{nthreads}{if a scalar, the number of CPU cores to use.
If a vector, launch as many threads as the length of the vector and attempt to pin the threads to those CPU cores (requires system to support thread affinity).
By default uses the \code{parallel} package to detect the number of physical cores.}

\item{beta.theta}{logical indicating whether the table should be returned in beta-theta space or in the standard space upon reaching the target variant \code{t}.
See the Details section.}
}
\value{
There is nothing returned.

\strong{NOTE:} for performance reasons, \code{bck} is updated in-place.
}
\description{
Propagates a \code{kalisBackwardTable} to an upstream variant.
The table is updated in-place.
}
\details{
\code{Backward} implements the backward algorithm to propagate the Li and Stephens rescaled hidden Markov model backward probabilities to a new target
variant.
Naturally, this can only propagate a table to variants upstream of its current position.

For mathematical details please see Section 2 of the kalis paper (Aslett and Christ, 2024).
Note that the precise formulation of the backward equation is determined by whether the flag \code{use.spiedel} is set in the parameters provided in \code{pars}.

\strong{Beta-theta space}

The rescaled HMM backward probabilities incorporate all of the haplotype relatedness information downstream from but NOT including the target variant, as is standard in the definition of HMM backward probabilities -- we refer to this as "beta space", or "rescaled probability space."
In contrast, the rescaled forward probabilities, naturally incorporate all of the haplotype relatedness information upstream from AND including the target variant.
Setting \code{beta.theta = TRUE} incorporates the target variant (analogous to the rescaled forward probabilities) -- we refer to this as "beta-theta space."

A backward table in beta-theta space (with \code{beta.theta = TRUE}) can be propagated to an upstream variant without incorporating that variant, thereby moving to beta space (\code{beta.theta = FALSE}), and vice versa.
However, while a backward table in beta space (\code{beta.theta = FALSE}) can be updated to incorporate the current variant, a backward table that is already in beta-theta space can not move to beta space without changing variants -- that would involve "forgetting" the current variant (see Examples).
}
\examples{

# Load the toy haplotype example and set toy parameters
data("SmallHaps")
data("SmallMap")

CacheHaplotypes(SmallHaps)

rho <- CalcRho(diff(SmallMap))
pars <- Parameters(rho)

# Create the backward table we want to propagate
bck <- MakeBackwardTable(pars)

# Calling Backward on this uninitialised table moves it to the last variant
Backward(bck, pars)
bck

# And again moves it to the next variant (etcetera)
Backward(bck, pars)
bck

# Or zoom to a particular variant
Backward(bck, pars, 150)
bck

# Now moving to variant 125 AND specifying "beta space" (though this is the
# default, just being very clear)
Backward(bck, pars, 125, beta.theta = FALSE)
bck

# Now just moving to "beta-theta space" (can be done with or without
# propagating to a new variant)
Backward(bck, pars, 125, beta.theta = TRUE)
bck

# However, attempting to move from "beta-theta space" back to "beta space"
# without propagating is not possible (see Details).
# The following will give an error (hence wrapped in try())
try(Backward(bck, pars, 125, beta.theta = FALSE))
bck

}
\references{
Aslett, L.J.M. and Christ, R.R. (2024) "kalis: a modern implementation of the Li & Stephens model for local ancestry inference in R", \emph{BMC Bioinformatics}, \strong{25}(1). Available at: \doi{10.1186/s12859-024-05688-8}.
}
\seealso{
\code{\link[=MakeBackwardTable]{MakeBackwardTable()}} to generate a backward table;
\code{\link[=Forward]{Forward()}} for the analogous forward propagation function;
\code{\link[=CopyTable]{CopyTable()}} to create a copy of table.
}
