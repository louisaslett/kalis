# kalis: High Performance Li &amp; Stephens Local Ancestry Inference

kalis provides a high performance implementation of the Li &
Stephens model (<https://www.ncbi.nlm.nih.gov/pubmed/14704198>) for local
ancestry inference (local referring to a region of the genome).  For a set of
N phased haplotypes, kalis computes the posterior marginal probability of each
haplotype copying every other haplotype by running N hidden Markov models in
parallel.  This yields an N x N matrix that summarizes the recent local
ancestry at each locus of interest.  The package provides functionality for
specifying a recombination map, site-specific mutation rates, and differing
prior copying probabilities for each recipient haplotype.  Extensive use is 
made of low level threading and CPU vector instructions.
