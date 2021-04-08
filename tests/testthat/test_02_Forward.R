make.fwd <- expression({
  fwd <- MakeForwardTable(pars, from_recipient, to_recipient)
  fwd.scmu <- MakeForwardTable(pars.scmu, from_recipient, to_recipient)
  fwd.scPi <- MakeForwardTable(pars.scPi, from_recipient, to_recipient)
  fwd.scmuPi <- MakeForwardTable(pars.scmuPi, from_recipient, to_recipient)

  nt.fwd <- MakeForwardTable(pars, from_recipient, to_recipient)
  nt.fwd.scmu <- MakeForwardTable(pars.scmu, from_recipient, to_recipient)
  nt.fwd.scPi <- MakeForwardTable(pars.scPi, from_recipient, to_recipient)
  nt.fwd.scmuPi <- MakeForwardTable(pars.scmuPi, from_recipient, to_recipient)

  Rfwd <- MakeForwardTable(pars, from_recipient, to_recipient)
  Rfwd.scmu <- MakeForwardTable(pars.scmu, from_recipient, to_recipient)
  Rfwd.scPi <- MakeForwardTable(pars.scPi, from_recipient, to_recipient)
  Rfwd.scmuPi <- MakeForwardTable(pars.scmuPi, from_recipient, to_recipient)
})

call.fwd <- expression({
  Forward(fwd, pars, t = target.l, nthreads = 1)
  Forward(fwd.scmu, pars.scmu, t = target.l, nthreads = 1)
  Forward(fwd.scPi, pars.scPi, t = target.l, nthreads = 1)
  Forward(fwd.scmuPi, pars.scmuPi, t = target.l, nthreads = 1)
  expect_identical(fwd$l, target.l)
  expect_identical(fwd.scmu$l, target.l)
  expect_identical(fwd.scPi$l, target.l)
  expect_identical(fwd.scmuPi$l, target.l)

  Forward(nt.fwd, pars, t = target.l)
  Forward(nt.fwd.scmu, pars.scmu, t = target.l)
  Forward(nt.fwd.scPi, pars.scPi, t = target.l)
  Forward(nt.fwd.scmuPi, pars.scmuPi, t = target.l)
  expect_identical(nt.fwd$l, target.l)
  expect_identical(nt.fwd.scmu$l, target.l)
  expect_identical(nt.fwd.scPi$l, target.l)
  expect_identical(nt.fwd.scmuPi$l, target.l)

  expect_warning(Forward(Rfwd, pars, t = target.l, nthreads = "R"), "Warning: using gold master R implementation")
  expect_warning(Forward(Rfwd.scmu, pars.scmu, t = target.l, nthreads = "R"), "Warning: using gold master R implementation")
  expect_warning(Forward(Rfwd.scPi, pars.scPi, t = target.l, nthreads = "R"), "Warning: using gold master R implementation")
  expect_warning(Forward(Rfwd.scmuPi, pars.scmuPi, t = target.l, nthreads = "R"), "Warning: using gold master R implementation")
  expect_identical(Rfwd$l, target.l)
  expect_identical(Rfwd.scmu$l, target.l)
  expect_identical(Rfwd.scPi$l, target.l)
  expect_identical(Rfwd.scmuPi$l, target.l)
})

test.fwd <- expression({
  expect_identical(fwd, nt.fwd)
  expect_identical(fwd.scmu, nt.fwd.scmu)
  expect_identical(fwd.scPi, nt.fwd.scPi)
  expect_identical(fwd.scmuPi, nt.fwd.scmuPi)

  expect_equal(Rfwd$alpha, nt.fwd$alpha)
  e <- err(nt.fwd$alpha, Rfwd$alpha)
  expect_lt(e$rel.err, 1e-10)
  expect_equal(e$zero.mismatches, 0)

  expect_equal(Rfwd.scmu$alpha, nt.fwd.scmu$alpha)
  e <- err(nt.fwd.scmu$alpha, Rfwd.scmu$alpha)
  expect_lt(e$rel.err, 1e-10)
  expect_equal(e$zero.mismatches, 0)

  expect_equal(Rfwd.scPi$alpha, nt.fwd.scPi$alpha)
  e <- err(nt.fwd.scPi$alpha, Rfwd.scPi$alpha)
  expect_lt(e$rel.err, 1e-10)
  expect_equal(e$zero.mismatches, 0)

  expect_equal(Rfwd.scmuPi$alpha, nt.fwd.scmuPi$alpha)
  e <- err(nt.fwd.scmuPi$alpha, Rfwd.scmuPi$alpha)
  expect_lt(e$rel.err, 1e-10)
  expect_equal(e$zero.mismatches, 0)


  expect_equal(Rfwd$alpha.f, nt.fwd$alpha.f)
  e <- err(nt.fwd$alpha.f, Rfwd$alpha.f)
  expect_lt(e$rel.err, 1e-10)
  expect_equal(e$zero.mismatches, 0)

  expect_equal(Rfwd.scmu$alpha.f, nt.fwd.scmu$alpha.f)
  e <- err(nt.fwd.scmu$alpha.f, Rfwd.scmu$alpha.f)
  expect_lt(e$rel.err, 1e-10)
  expect_equal(e$zero.mismatches, 0)

  expect_equal(Rfwd.scPi$alpha.f, nt.fwd.scPi$alpha.f)
  e <- err(nt.fwd.scPi$alpha.f, Rfwd.scPi$alpha.f)
  expect_lt(e$rel.err, 1e-10)
  expect_equal(e$zero.mismatches, 0)

  expect_equal(Rfwd.scmuPi$alpha.f, nt.fwd.scmuPi$alpha.f)
  e <- err(nt.fwd.scmuPi$alpha.f, Rfwd.scmuPi$alpha.f)
  expect_lt(e$rel.err, 1e-10)
  expect_equal(e$zero.mismatches, 0)
})

run.test <- expression({
  for(start.l in 1:2) {
    for(move.l in as.integer(c(NA,0,1,5,10))) {
      eval(make.fwd)

      target.l <- start.l
      eval(call.fwd)

      if(!is.na(move.l)) {
        target.l <- start.l + move.l
        eval(call.fwd)
      }

      eval(test.fwd)
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

    test_that(glue::glue("forward with {N} haps, all recipients, works"), {
      from_recipient <- 1
      to_recipient <- N
      eval(run.test)
    })

    recip <- sample(1:N, 2)
    from_recipient <- min(recip)
    to_recipient <- max(recip)
    test_that(glue::glue("forward with {N} haps, recipients {from_recipient}-{to_recipient}, works"), {
      eval(run.test)
    })
  }
}
