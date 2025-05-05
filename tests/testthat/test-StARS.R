# test script for StARS.R - testcases are NOT comprehensive!

# -----------------------------------------------------------
## Skip gracefully if dependencies are missing  --------------
## -----------------------------------------------------------
skip_if_not(
  exists("graphicalLasso"),
  message = "graphicalLasso() not found; skipping StARS tests"
)

## -----------------------------------------------------------
## Helper: simulate quick toy data ---------------------------
## -----------------------------------------------------------
sim_X <- function(n = 24, p = 7) {
  matrix(rnorm(n * p), n, p)
}

beta_target <- 0.05          # hard‑coded inside StARS

## -----------------------------------------------------------
## 1.  Structural integrity & range checks -------------------
## -----------------------------------------------------------
test_that("StARS returns numeric vector of reasonable values", {
  set.seed(11)
  X <- sim_X()
  
  out <- StARS(X, b = 12, M = 25, print.detail = FALSE)
  expect_type(out, "double")
  expect_length(out, 2)
  
  rho <- out[1]; D  <- out[2]
  expect_true(rho >= 0.01 && rho <= 1)
  expect_true(D  >= 0    && D  <= 1)
  # By definition D_var must not exceed beta (≈0.05)
  expect_true(D <= beta_target + 1e-10)
})

## -----------------------------------------------------------
## 2.  Determinism given internal set.seed -------------------
## -----------------------------------------------------------
test_that("StARS result is deterministic for given data", {
  X <- sim_X()
  
  res1 <- StARS(X, b = 10, M = 20, print.detail = FALSE)
  res2 <- StARS(X, b = 10, M = 20, print.detail = FALSE)
  expect_identical(res1, res2)
})

## -----------------------------------------------------------
## 3.  print.detail flag toggles console output --------------
## -----------------------------------------------------------
test_that("print.detail controls console chatter", {
  X <- sim_X()
  # Silent execution
  expect_silent(StARS(X, b = 8, M = 15, print.detail = FALSE))
  # Verbose execution emits at least one message
  expect_message(StARS(X, b = 8, M = 15, print.detail = TRUE),
                 regexp = ".", all = FALSE)
})

## -----------------------------------------------------------
## 4.  Basic input validation / corner cases -----------------
## -----------------------------------------------------------
test_that("StARS errors when b exceeds sample size", {
  X <- sim_X(n = 15, p = 5)
  # b larger than N should trigger an error from sample.int
  expect_error(StARS(X, b = 20, M = 10, print.detail = FALSE),
               regexp = "cannot take a sample|sample.int")
})

test_that("b = 1 still returns valid output", {
  X <- sim_X(n = 10, p = 4)
  out <- StARS(X, b = 1, M = 12, print.detail = FALSE)  # degenerate subsamples
  expect_length(out, 2)
  expect_true(out[2] <= beta_target + 1e-10)
})
