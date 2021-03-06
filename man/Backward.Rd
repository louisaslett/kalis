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
\code{\link{MakeBackwardTable}}.}

\item{pars}{a \code{kalisParameters} object, as returned by
\code{\link{Parameters}}.}

\item{t}{a variant to move the backward table to.
Must be less than or equal to current variant of \code{bck}.
By default, it simply advances to the variant immediately upstream.}

\item{nthreads}{the number of CPU cores to use.
By default uses the \code{parallel} package to detect the number of physical cores.}

\item{beta.theta}{logical indicating whether the table should be returned in
beta-theta space or in the standard space upon reaching the target variant
\code{t}.  See the details section.}
}
\value{
There is nothing returned.
For performance reasons, \code{bck} is updated in-place.
}
\description{
Propagates a \code{kalisBackwardTable} to an upstream variant.
The table is updated in-place.
}
\details{
\code{Backward} implements the backward algorithm to propagate the Li and
Stephens rescaled hidden Markov model backward probabilities to a new target
variant.  Naturally, this can only propagate a table to variants upstream
of its current position.

For mathematical details please see Section 2 of the kalis paper (TODO: ref).
Note that the precise formulation of the backward equation is determined by
whether the flag \code{use.spiedel} is set in the parameters provided in
\code{pars}.

\strong{Beta-theta space}

The rescaled HMM backward probabilities incorporate
all of the haplotype relatedness information downstream from but NOT including the
target variant, as is standard in the definition of HMM backward probabilities -- we refer to this as "beta space."
In contrast, the rescaled forward probabilities, naturally incorporate all of the
haplotype relatedness information upstream from AND including the target variant.
Setting \code{beta.theta = TRUE} incorporates the target variant (analogous to the rescaled forward probabilities) --
we refer to this as "beta-theta space".

A backward table in beta-theta space (with \code{beta.theta = TRUE}) can be propagated
to an upstream variant without incorporating that variant, thereby moving to beta space (\code{beta.theta = FALSE}), and vice versa.
However, while a backward table in beta space (\code{beta.theta = FALSE})
can be updated to incorporate the current variant, a backward table
that is already in beta-theta space can not move to beta space
without changing variants -- that would involve "forgetting" the current variant (see Examples).
}
\examples{
\dontrun{
# Load the toy haplotype example and set toy parameters
CacheHaplotypes(SmallHaplotypes)
pars <- Parameters(rho = CalcRho(cM = SmallMap))

# Create a backward table for the hidden Markov model incorporating all
# recipient and donor haplotypes
bck <- MakeBackwardTable(pars)

# Create a backward table for the hidden Markov model incorporating only
# recipient haplotypes 100 to 200 (inclusive) and all donor haplotypes.
bck <- MakeBackwardTable(pars, 100, 200)

# This table is uninitialised, but ready to pass to the Backward function
# which will trigger initialisation and propagation from the last variant.
# For example, initialise and propagate backward to the 10th variant:
Backward(bck, pars, 10, nthreads = 8)
print(bck)

#### Beta-theta space example ####
# Now moving to beta-theta space
Backward(bck, pars, 8, beta.theta=TRUE)

# Now moving to beta space (can be done with or without propagating to a new variant)
Backward(bck, pars, 8, beta.theta=FALSE)

# Attempting to move from beta space to beta-theta space without propagating
# is not possible (see Details).
Backward(bck, pars, 8, beta.theta=TRUE)
}

}
\seealso{
\code{\link{MakeBackwardTable}} to generate a backward table;
\code{\link{Forward}} for the analogous forward propagation function;
\code{\link{CopyTable}} to create a copy of table.
}
