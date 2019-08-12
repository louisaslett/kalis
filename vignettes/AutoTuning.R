## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
)

## ----run.dir-------------------------------------------------------------
run.dir <- "~/Desktop/autotuning_test_files/"

## ----load.data, results = "hide",message=FALSE---------------------------
require(kalis)
require(RSpectra)
require(mlrMBO)
require(future)
require(pryr)
require(QForm)
require(plotly)
require(akima)

source("../../AutoTune2_local.R")

data("SmallHaps")
#haps <- SmallHaps[polymorphic.sites,]
CacheHaplotypes(SmallHaps)
m <- rbeta(nrow(SmallHaps)-1,1,10)*1e-6
pars <- Parameters(CalcRho(morgan.dist = m, Ne = 1e-250, gamma = 1), mu = 1e-250)

polymorphic.sites <- rowSums(SmallHaps)!=0 & rowSums(SmallHaps)!=1
G <- tcrossprod(scale(t(SmallHaps[polymorphic.sites,])))
bg.vecs <- RSpectra::eigs(G,10)$vectors

## ----autotuning_pars-----------------------------------------------------
nthreads <- 2
num.target.loci <- 10
targets <- InvRecombMap(m,num.target.loci = num.target.loci)

