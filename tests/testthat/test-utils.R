# test script for utils.R - testcases are NOT comprehensive!

## -----------------------------------------------------------
## 1.  selfoldforCV   ----------------------------------------
## -----------------------------------------------------------
test_that("selfoldforCV chooses folds correctly", {
  # Standard: N divisible by many folds; expect 5
  expect_equal(selfoldforCV(180), 5)
  
  expect_equal(selfoldforCV(96), 4)
  
  # No suitable folds => returns 1
  expect_equal(selfoldforCV(89), 1)
})

## -----------------------------------------------------------
## 2.  soft -----------------------------------------------
## -----------------------------------------------------------
test_that("soft thresholding behaves as expected", {
  expect_equal(soft( 3, 1),  2)
  expect_equal(soft(-3, 1), -2)
  expect_equal(soft( 0.5, 1), 0)
  expect_equal(soft(-0.5, 1), 0)
  # Vectorised safety
  x <- c(-3, -0.5, 0.5, 3)
  expect_equal(soft(x, 1), c(-2, 0, 0, 2))
})

## -----------------------------------------------------------
## 3.  CDfgL  ----------------------------------------------
## -----------------------------------------------------------
skip_if_not(exists("CDfgL"), "CDfgL not present")

test_that("CDfgL solves identity‑matrix sub‑problem exactly", {
  set.seed(1)
  p <- 3
  V <- diag(p)
  u <- c(1, -2, 0.5)
  rho <- 0.3
  sol <- CDfgL(V, beta_i = rep(0, p), u = u, rho = rho,
               print.detail = FALSE, maxIter = 1000)
  
  expect_equal(sol, soft(u, rho), tolerance = 1e-6)
})

## -----------------------------------------------------------
## 4.  update_CorrectCoef  ----------------------------------
## -----------------------------------------------------------
skip_if_not(exists("update_CorrectCoef"), "update_CorrectCoef not present")

test_that("update_CorrectCoef returns sane dimensions and small changes on identical data", {
  set.seed(123)
  G <- 1; n <- 12; p <- 4
  X0 <- list(matrix(rnorm(n * p), n, p))
  X1 <- X0                                   # identical batches
  Theta.list <- list(diag(p))                # identity precision
  res <- update_CorrectCoef(
    X0, X1, Theta.list,
    a.i = rep(1, p), b.i = rep(0, p),
    penal.ksi = 0.1, penal.gamma = 0.1,
    print.detail = FALSE
  )
  
  expect_named(res, c("coef.a", "coef.b"), ignore.order = TRUE)
  expect_length(res$coef.a, p)
  expect_length(res$coef.b, p)
  
  # On identical batches, a ≈ 1 and b ≈ 0
  expect_equal(res$coef.a, rep(0.2, p), tolerance = 0.2)
  expect_equal(res$coef.b, rep(0, p), tolerance = 0.2)
})

## -----------------------------------------------------------
## 5.  findBestPara  ----------------------------------------
## -----------------------------------------------------------
skip_if_not(exists("BEgLasso"), "BEgLasso not available — skipping findBestPara tests")

test_that("findBestPara returns selected penalties", {
  set.seed(11)
  G <- 1; n <- 10; p <- 3
  X0 <- list(matrix(rnorm(n * p), n, p))
  X1 <- list(matrix(rnorm(n * p), n, p))
  
  res <- findBestPara(X0, X1, penal.rho = 0.2, eps = 1e-2, print.detail = FALSE)
  expect_type(res, "list")
  expect_named(res, c("penal.ksi", "penal.gamma", "MinAvedist"), ignore.order = TRUE)
  expect_true(res$penal.ksi %in% c(1, 0.5, 0.3, 0.1))
  expect_true(res$penal.gamma %in% c(1, 0.5, 0.3, 0.1))
  expect_true(is.finite(res$MinAvedist))
})

## -----------------------------------------------------------
## 6.  selrho.useCVBIC  -------------------------------------
## -----------------------------------------------------------
skip_if_not(exists("graphicalLasso"), "graphicalLasso not available")

test_that("selrho.useCVBIC returns rho from candidate grid", {
  set.seed(7)
  X <- matrix(rnorm(120), 20, 6)
  out <- selrho.useCVBIC(X, print.detail = FALSE)
  expect_length(out, 2)
  expect_true(out[1] %in% seq(0.1, 0.9, 0.1))
  expect_true(is.finite(out[2]))
})

## -----------------------------------------------------------
## 7.  DelOutlier  ------------------------------------------
## -----------------------------------------------------------
skip_if_not_installed("mixOmics")

test_that("DelOutlier detects and removes extreme samples", {
  set.seed(21)
  X <- matrix(rnorm(200), 20, 10)
  X[1, ] <- 10                               # big outlier row
  
  res <- DelOutlier(X)
  expect_type(res, "list")
  expect_named(res, c("delsampIdx", "X.out"), ignore.order = TRUE)
  expect_true(1 %in% res$delsampIdx)
  expect_equal(nrow(res$X.out), 19)
})

test_that("DelOutlier keeps data intact when no strong outliers", {
  X <- matrix(rnorm(100), 10, 10)
  res <- DelOutlier(X)
  expect_length(res$delsampIdx, 0)
  expect_equal(nrow(res$X.out), 10)
})

## -----------------------------------------------------------
## 8.  ImputeOutlier  ---------------------------------------
## -----------------------------------------------------------

test_that("ImputeOutlier replaces extreme points and leaves no NAs", {
  set.seed(31)
  X <- matrix(rnorm(150), 15, 10)
  X[5, 3] <- 30                              # extreme high value
  
  Y <- ImputeOutlier(X)
  expect_equal(dim(Y), dim(X))
  expect_false(anyNA(Y))
  # Value should have shrunk toward centre (not remain extreme)
  expect_true(abs(Y[5, 3]) < 10)
})
