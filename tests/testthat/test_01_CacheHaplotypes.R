sample.haps <- function(sorted = FALSE) {
  vars <- 1
  haps <- 1
  # Get a sub-matrix with at least 50 mutations/no more than 350 or else we're not really testing anything!
  while(sum(SmallHaps[vars,haps]) < 50 || sum(SmallHaps[vars,haps]) > 350) {
    vars <- sample(1:L(), 20)
    haps <- sample(1:N(), 20)
  }
  if(sorted) {
    vars <- sort(vars)
    haps <- sort(haps)
  }
  list(vars = vars, haps = haps)
}

test_that("caching haplotype R matrix works", {
  x <- CacheHaplotypes(SmallHaps)
  expect_identical(SmallHaps, QueryCache())
  expect_identical(dim(SmallHaps), x)
  expect_output(CacheSummary(), glue::glue("Cache currently loaded with {dim(SmallHaps)[2]} haplotypes, each with {dim(SmallHaps)[1]} variants"))
})

test_that("cache sizing information correct", {
  expect_identical(dim(SmallHaps)[2], N())
  expect_identical(dim(SmallHaps)[1], L())
})

test_that("extracting submatrix from cache works", {
  s <- sample.haps()
  expect_identical(SmallHaps[s$vars,s$haps], QueryCache(s$vars, s$haps))
})

test_that("loading only submatrix to cache works", {
  s <- sample.haps(TRUE)
  expect_warning(x <- CacheHaplotypes(SmallHaps, s$vars, s$haps), "haplotypes already cached ... overwriting existing cache")
  expect_identical(SmallHaps[s$vars,s$haps], QueryCache())
  expect_identical(c(20L,20L), x)
})

test_that("cache can be cleared", {
  ClearHaplotypeCache()
  expect_identical(NULL, N())
  expect_identical(NULL, L())
  expect_output(CacheSummary(), "Cache currently empty")
})

test_that("caching hdf5 haplotype matrix works", {
  x <- CacheHaplotypes(system.file("small_example", "small.h5", package = "kalis"))
  expect_identical(SmallHaps, QueryCache())
  expect_identical(dim(SmallHaps), x)
  expect_output(CacheSummary(), glue::glue("Cache currently loaded with {dim(SmallHaps)[2]} haplotypes, each with {dim(SmallHaps)[1]} variants"))
})

test_that("caching hdf5 haplotype sub-matrix works", {
  s <- sample.haps(TRUE)
  expect_warning(x <- CacheHaplotypes(system.file("small_example", "small.h5", package = "kalis"), s$vars, s$haps), "haplotypes already cached ... overwriting existing cache")
  expect_identical(SmallHaps[s$vars,s$haps], QueryCache())
  expect_identical(c(20L,20L), x)
  expect_output(CacheSummary(), glue::glue("Cache currently loaded with 20 haplotypes, each with 20 variants"))
})

test_that("caching .hap.hz haplotype matrix works", {
  expect_warning(x <- CacheHaplotypes(system.file("small_example", "small.hap.gz", package = "kalis")), "haplotypes already cached ... overwriting existing cache")
  expect_identical(SmallHaps, QueryCache())
  expect_identical(dim(SmallHaps), x)
  expect_output(CacheSummary(), glue::glue("Cache currently loaded with {dim(SmallHaps)[2]} haplotypes, each with {dim(SmallHaps)[1]} variants"))
})

test_that("caching .hap.hz haplotype sub-matrix works", {
  s <- sample.haps(TRUE)
  expect_warning(x <- CacheHaplotypes(system.file("small_example", "small.hap.gz", package = "kalis"), s$vars, s$haps), "haplotypes already cached ... overwriting existing cache")
  expect_identical(SmallHaps[s$vars,s$haps], QueryCache())
  expect_identical(c(20L,20L), x)
  expect_output(CacheSummary(), glue::glue("Cache currently loaded with 20 haplotypes, each with 20 variants"))
})
