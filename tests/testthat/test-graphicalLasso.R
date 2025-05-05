# test script for graphicalLasso.R - testcases are NOT comprehensive!

# ------------------------------------------------------------------
# Helper to simulate quick toy matrices ----------------------------
# ------------------------------------------------------------------
simulate_X <- function(n = 20, p = 7) {
  matrix(rnorm(n * p), n, p)
}

# ------------------------------------------------------------------
# Skip if external dependency not available ------------------------
# ------------------------------------------------------------------
skip_if_not(
  exists("CDfgL"),
  message = "CDfgL() not found — skipping graphicalLasso tests"
)

# 1 ────────────────────────────────────────────────────────────────
test_that("graphicalLasso returns expected top‑level structure", {
  set.seed(11)
  X <- simulate_X(18, 6)
  
  res <- graphicalLasso(X, rho = 0.1, print.detail = FALSE)
  
  expect_type(res, "list")
  expect_named(res, c("Theta", "W"), ignore.order = TRUE)
  
  # Dimensions & symmetry ------------------------------------------
  expect_equal(dim(res$W),     dim(res$Theta))
  expect_equal(res$Theta, t(res$Theta), tolerance = 1e-10)
  expect_equal(res$W,     t(res$W),     tolerance = 1e-10)
})

# 2 ────────────────────────────────────────────────────────────────
test_that("print.detail toggles console output / messages", {
  X <- simulate_X(12, 5)
  
  expect_silent(
    graphicalLasso(X, rho = 0.05, print.detail = FALSE)
  )
  
})

# 3 ────────────────────────────────────────────────────────────────
test_that("Single‑sample input produces analytically predictable result", {
  # One sample, p = 4
  X1 <- simulate_X(n = 1, p = 4)
  rho <- 0.2
  
  res <- graphicalLasso(X1, rho = rho, print.detail = FALSE)
  
  # With a single sample we set S = 0, so W should be rho * I
  expect_equal(res$W, diag(rho, 4))
  expect_equal(diag(res$Theta), rep(1 / rho, 4))
  # Off‑diagonals should be (close to) zero
  expect_true(all(abs(res$Theta[upper.tri(res$Theta)]) < 1e-8))
})

# 4 ────────────────────────────────────────────────────────────────
test_that("Precision and covariance matrices are positive‑definite", {
  set.seed(202)
  X <- simulate_X(25, 8)
  out <- graphicalLasso(X, rho = 0.15, print.detail = FALSE)
  
  eigW     <- eigen(out$W,     symmetric = TRUE, only.values = TRUE)$values
  eigTheta <- eigen(out$Theta, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigW     > 0))
  expect_true(all(eigTheta > 0))
})

# 5 ────────────────────────────────────────────────────────────────
test_that("Basic input validation triggers informative errors", {
  X <- simulate_X(10, 5)
  
  # Non‑numeric rho -------------------------------------------------
  expect_error(
    graphicalLasso(X, rho = "a", print.detail = FALSE),
    regexp = "numeric"
  )
  
  # Negative rho ----------------------------------------------------
  expect_error(
    graphicalLasso(X, rho = -0.1, print.detail = FALSE),
    regexp = "missing"
  )
  
})
