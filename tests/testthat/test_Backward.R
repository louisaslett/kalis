make.bck <- expression({
  bck <- MakeBackwardTable(pars, from_recipient, to_recipient)
  bck.scmu <- MakeBackwardTable(pars.scmu, from_recipient, to_recipient)
  bck.scPi <- MakeBackwardTable(pars.scPi, from_recipient, to_recipient)
  bck.scmuPi <- MakeBackwardTable(pars.scmuPi, from_recipient, to_recipient)

  nt.bck <- MakeBackwardTable(pars, from_recipient, to_recipient)
  nt.bck.scmu <- MakeBackwardTable(pars.scmu, from_recipient, to_recipient)
  nt.bck.scPi <- MakeBackwardTable(pars.scPi, from_recipient, to_recipient)
  nt.bck.scmuPi <- MakeBackwardTable(pars.scmuPi, from_recipient, to_recipient)

  Rbck <- MakeBackwardTable(pars, from_recipient, to_recipient)
  Rbck.scmu <- MakeBackwardTable(pars.scmu, from_recipient, to_recipient)
  Rbck.scPi <- MakeBackwardTable(pars.scPi, from_recipient, to_recipient)
  Rbck.scmuPi <- MakeBackwardTable(pars.scmuPi, from_recipient, to_recipient)
})

call.bck <- expression({
  Backward(bck, pars, t = target.l, nthreads = 1, beta.theta = target.beta.theta)
  Backward(bck.scmu, pars.scmu, t = target.l, nthreads = 1, beta.theta = target.beta.theta)
  Backward(bck.scPi, pars.scPi, t = target.l, nthreads = 1, beta.theta = target.beta.theta)
  Backward(bck.scmuPi, pars.scmuPi, t = target.l, nthreads = 1, beta.theta = target.beta.theta)
  expect_identical(bck$l, target.l)
  expect_identical(bck.scmu$l, target.l)
  expect_identical(bck.scPi$l, target.l)
  expect_identical(bck.scmuPi$l, target.l)

  Backward(nt.bck, pars, t = target.l, beta.theta = target.beta.theta)
  Backward(nt.bck.scmu, pars.scmu, t = target.l, beta.theta = target.beta.theta)
  Backward(nt.bck.scPi, pars.scPi, t = target.l, beta.theta = target.beta.theta)
  Backward(nt.bck.scmuPi, pars.scmuPi, t = target.l, beta.theta = target.beta.theta)
  expect_identical(nt.bck$l, target.l)
  expect_identical(nt.bck.scmu$l, target.l)
  expect_identical(nt.bck.scPi$l, target.l)
  expect_identical(nt.bck.scmuPi$l, target.l)

  expect_warning(Backward(Rbck, pars, t = target.l, nthreads = "R", beta.theta = target.beta.theta), "Warning: using gold master R implementation")
  expect_warning(Backward(Rbck.scmu, pars.scmu, t = target.l, nthreads = "R", beta.theta = target.beta.theta), "Warning: using gold master R implementation")
  expect_warning(Backward(Rbck.scPi, pars.scPi, t = target.l, nthreads = "R", beta.theta = target.beta.theta), "Warning: using gold master R implementation")
  expect_warning(Backward(Rbck.scmuPi, pars.scmuPi, t = target.l, nthreads = "R", beta.theta = target.beta.theta), "Warning: using gold master R implementation")
  expect_identical(Rbck$l, target.l)
  expect_identical(Rbck.scmu$l, target.l)
  expect_identical(Rbck.scPi$l, target.l)
  expect_identical(Rbck.scmuPi$l, target.l)
})

test.bck <- expression({
  expect_identical(bck, nt.bck)
  expect_identical(bck.scmu, nt.bck.scmu)
  expect_identical(bck.scPi, nt.bck.scPi)
  expect_identical(bck.scmuPi, nt.bck.scmuPi)

  expect_equal(Rbck$beta, nt.bck$beta)
  e <- err(nt.bck$beta, Rbck$beta)
  expect_lt(e$rel.err, 1e-10)
  expect_equal(e$zero.mismatches, 0)

  expect_equal(Rbck.scmu$beta, nt.bck.scmu$beta)
  e <- err(nt.bck.scmu$beta, Rbck.scmu$beta)
  expect_lt(e$rel.err, 1e-10)
  expect_equal(e$zero.mismatches, 0)

  expect_equal(Rbck.scPi$beta, nt.bck.scPi$beta)
  e <- err(nt.bck.scPi$beta, Rbck.scPi$beta)
  expect_lt(e$rel.err, 1e-10)
  expect_equal(e$zero.mismatches, 0)

  expect_equal(Rbck.scmuPi$beta, nt.bck.scmuPi$beta)
  e <- err(nt.bck.scmuPi$beta, Rbck.scmuPi$beta)
  expect_lt(e$rel.err, 1e-10)
  expect_equal(e$zero.mismatches, 0)


  expect_equal(Rbck$beta.g, nt.bck$beta.g)
  e <- err(nt.bck$beta.g, Rbck$beta.g)
  expect_lt(e$rel.err, 1e-10)
  expect_equal(e$zero.mismatches, 0)

  expect_equal(Rbck.scmu$beta.g, nt.bck.scmu$beta.g)
  e <- err(nt.bck.scmu$beta.g, Rbck.scmu$beta.g)
  expect_lt(e$rel.err, 1e-10)
  expect_equal(e$zero.mismatches, 0)

  expect_equal(Rbck.scPi$beta.g, nt.bck.scPi$beta.g)
  e <- err(nt.bck.scPi$beta.g, Rbck.scPi$beta.g)
  expect_lt(e$rel.err, 1e-10)
  expect_equal(e$zero.mismatches, 0)

  expect_equal(Rbck.scmuPi$beta.g, nt.bck.scmuPi$beta.g)
  e <- err(nt.bck.scmuPi$beta.g, Rbck.scmuPi$beta.g)
  expect_lt(e$rel.err, 1e-10)
  expect_equal(e$zero.mismatches, 0)
})

run.test <- expression({
  for(start.l in c(L(), L()-1L)) {
    for(move.l in as.integer(c(NA,0,1,5,10))) {
      for(initial.beta.theta in c(TRUE, FALSE)) {
        for(final.beta.theta in c(TRUE, FALSE)) {
          # Skip impossible moves
          if(!is.na(move.l) && move.l == 0 && initial.beta.theta && !final.beta.theta) next
          if(is.na(move.l) && !final.beta.theta) next

          eval(make.bck)

          target.l <- start.l
          target.beta.theta <- initial.beta.theta
          eval(call.bck)

          if(!is.na(move.l)) {
            target.l <- start.l - move.l
            target.beta.theta <- final.beta.theta
            eval(call.bck)
          }

          eval(test.bck)
        }
      }
    }
  }
})

vecwidth <- .Call(kalis:::CCall_VectorBitWidth)
Ns <- c(10,(vecwidth-2):(vecwidth+2),(vecwidth+32-2):(vecwidth+32+2),(2*vecwidth-2):(2*vecwidth+2))

for(N in Ns) {
  if(N > Ns[1]) {
    skip_on_cran()
  }

  p <- start.new.example(N, 20)
  for(use.speidel in c(TRUE, FALSE)) {
    eval(make.pars)

    test_that(glue::glue("backward with {N} haps, all recipients, works"), {
      from_recipient <- 1
      to_recipient <- N
      eval(run.test)
    })

    recip <- sample(1:N, 2)
    from_recipient <- min(recip)
    to_recipient <- max(recip)
    test_that(glue::glue("backward with {N} haps, recipients {from_recipient}-{to_recipient}, works"), {
      eval(run.test)
    })
  }
}
