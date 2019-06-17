## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
require(kalis)
data("SmallHaps")

## ------------------------------------------------------------------------
WriteIndividualHaplotypeH5("SmallHaps.h5",SmallHaps)
CacheHaplotypes("SmallHaps.h5")

## ---- echo = F, results = 'hide'-----------------------------------------
system2("rm","SmallHaps.h5")

## ------------------------------------------------------------------------
m <- rbeta(500-1,1,10)*1e-6
pars <- Parameters(CalcRho(morgan.dist = m, Ne = 1, gamma = 1), mu = 1e-8)

## ------------------------------------------------------------------------
fwd <- MakeForwardTable(pars)
bck <- MakeBackwardTable(pars)

Forward(fwd, pars, 250)
Backward(bck, pars, 250)

## ------------------------------------------------------------------------
p <- PostProbs(fwd,bck)
d <- DistMat(fwd,bck)

## ---- results='asis'-----------------------------------------------------
plot(d)

