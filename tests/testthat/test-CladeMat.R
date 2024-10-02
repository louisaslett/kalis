test_that("CladeMat Matches GoldMaster CladeMat", {

  CacheHaplotypes(SmallHaps)
  map <- SmallMap

  target.idx <- floor(L()/2)
  pars <- Parameters(CalcRho(diff(map),s = 0.01),mu = 0.0001) # specify HMM parameters
  fwd <- MakeForwardTable(pars)
  bck <- MakeBackwardTable(pars)
  Forward(fwd, pars, target.idx)
  Backward(bck, pars, target.idx)

  unit.dist <- -log(pars$pars$mu)
  thresh <- 0.2

  M <- matrix(0,N()/2,N()/2)
  neigh <- CladeMat(fwd,bck,M,unit.dist = unit.dist, thresh = thresh ,max1var = TRUE)

  M2 <- kalis:::CladeMat.GM(fwd,bck,unit.dist = unit.dist, thresh = thresh)

  expect_equal(M, M2)
})
