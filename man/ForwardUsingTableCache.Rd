% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TableCache.R
\name{ForwardUsingTableCache}
\alias{ForwardUsingTableCache}
\title{Use a forward table cache to propagate}
\usage{
ForwardUsingTableCache(
  fwd,
  pars,
  cache,
  t = fwd$l - 1,
  nthreads = min(parallel::detectCores(logical = FALSE), fwd$to_recipient -
    fwd$from_recipient + 1)
)
}
\arguments{
\item{fwd}{a \code{kalisForwardTable} object, as returned by \code{\link[=MakeForwardTable]{MakeForwardTable()}}.}

\item{pars}{a \code{kalisParameters} object, as returned by \code{\link[=Parameters]{Parameters()}}.}

\item{cache}{a cache of forward tables as generated by \code{\link[=CreateForwardTableCache]{CreateForwardTableCache()}} and filled using \code{\link[=FillTableCache]{FillTableCache()}}.}

\item{t}{a locus position to move the forward table to, starting the forward propagation from whatever table in the \code{cache} variable is immediately before locus \code{t}.
By default, it simply advances to the previous locus (which is the natural direction to move when using the cache).}

\item{nthreads}{the number of CPU cores to use.
By default no parallelism is used.}
}
\value{
There is nothing returned.

\strong{NOTE:} for performance reasons, \code{fwd} is updated in-place.
}
\description{
An in-memory cache for forward tables, which has already been filled, can be used to move more quickly to a specified locus.
}
\details{
If the objective is to run the Li and Stephens hidden Markov model both forwards and backwards to the same locus (and to do so for every possible locus), then considerable efficiency can be achieved by first performing a full scan forwards, filling a geometrically spaced cache whilst doing so.
Then, by working backwards, the backward propagation moves one locus at a time and the forward propagation can move backwards by moving forward from a recently cached local table.

Memory for a cache can be allocated using \code{\link[=CreateForwardTableCache]{CreateForwardTableCache()}} and should then be filled with \code{\link[=FillTableCache]{FillTableCache()}}.
To use the cache, then instead of using the \code{\link[=Forward]{Forward()}} function, use this function.

Note that the \code{cache} which is passed to this function will be dynamically updated based on the locus requested: the assumption is that the cache is used to propagate in reverse so any cache entries for a locus position past \code{t} are taken to be no longer needed and that space will redeployed to more densely fill the cache with earlier locus positions.
}
\examples{
\dontrun{
# This code assumes you have already:
#  i) cached the haplotypes using CacheHaplotypes function
#  ii) setup parameters in a variable called pars
#  iii) set the number of loci in a variable called L

# Allocate up to 10GB to a cache, with parameters already setup in pars ...
cache <- CreateForwardTableCache(pars, 10)
# ... and fill it
FillTableCache(cache, pars, nthreads = 8)

# Create forward and backward tables
fwd <- MakeForwardTable(pars)
bck <- MakeBackwardTable(pars)

# Then reach every locus faster by iterating backwards, using the cache to
# move the forward table into position faster
for(l in L:1) {
  Backward(bck, pars, l, nthreads = 8)
  ForwardUsingTableCache(fwd, pars, cache, l, nthreads = 8)
  # Do whatever work is required at
  # every locus here using fwd and bck
}
}

}
\references{
Christ, R.R., Wang, X., Aslett, L.J.M., Steinsaltz, D. and Hall, I. (2024) "Clade Distillation for Genome-wide Association Studies", bioRxiv 2024.09.30.615852. Available at: \doi{10.1101/2024.09.30.615852}.
}
\seealso{
\code{\link[=MakeForwardTable]{MakeForwardTable()}} to make a forward table;
\code{\link[=CreateForwardTableCache]{CreateForwardTableCache()}} to generate a cache;
\code{\link[=FillTableCache]{FillTableCache()}} to fill a cache;
\code{\link[=Forward]{Forward()}} for forward function without using a cache.

Alternatively, see \code{\link[=ForwardIterator]{ForwardIterator()}} to create an iterator which internally creates or uses a table cache.
}
