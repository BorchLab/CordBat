# test script for BEgLasso.R - testcases are NOT comprehensive!

# ------------------------------------------------------------------
# Helper for quickly generating toy data of matching dimensions ----
# ------------------------------------------------------------------
simulate_batches <- function(G = 2, n = 15, p = 6) {
  X0 <- lapply(seq_len(G), \(.) matrix(rnorm(n * p), n, p))
  X1 <- lapply(seq_len(G), \(.) matrix(rnorm(n * p), n, p))
  list(X0 = X0, X1 = X1, p = p, G = G)
}

# ------------------------------------------------------------------
#  Sanity‑check that required helpers are available; if not, skip ---
# ------------------------------------------------------------------
skip_if_not(
  all(vapply(c("CDfgL", "update.CorrectCoef"), exists, logical(1))),
  message = "Helper functions CDfgL / update.CorrectCoef not found"
)

# 1 ────────────────────────────────────────────────────────────────
test_that("BEgLasso returns the expected top‑level structure", {
  set.seed(123)
  dat <- simulate_batches(G = 2, n = 12, p = 5)
  
  res <- BEgLasso(
    X0.glist = dat$X0,
    X1.glist = dat$X1,
    penal.rho  = 0.1,
    penal.ksi  = 0.1,
    penal.gamma = 0.1,
    eps  = 1e-2,
    print.detail = FALSE
  )
  
  expect_type(res, "list")
  expect_named(res, c("Theta", "X1.cor", "coef.a", "coef.b"), ignore.order = TRUE)
  
  # Theta ----------------------------------------------------------
  expect_length(res$Theta, dat$G)
  lapply(res$Theta, \(mat) {
    expect_true(is.matrix(mat))
    expect_equal(dim(mat), c(dat$p, dat$p))
    expect_equal(mat, t(mat), tolerance = 1e-10)      # symmetry
  })
  
  # X1.cor ---------------------------------------------------------
  expect_length(res$X1.cor, dat$G)
  for (g in seq_len(dat$G)) {
    expect_equal(dim(res$X1.cor[[g]]), dim(dat$X1[[g]]))
  }
  
  # Coefficients ---------------------------------------------------
  expect_length(res$coef.a, dat$p)
  expect_length(res$coef.b, dat$p)
  expect_true(all(res$coef.a >= 0))                   # enforced post‑hoc
})

# 2 ────────────────────────────────────────────────────────────────
test_that("print.detail toggles console output / messages", {
  set.seed(1)
  batches <- simulate_batches(G = 1, n = 10, p = 4)
  
  # Silent mode
  expect_silent(
    BEgLasso(batches$X0, batches$X1,
             penal.rho = 0.05, penal.ksi = 0.05, penal.gamma = 0.05,
             eps = 1e-2, print.detail = FALSE)
  )
})

# 3 ────────────────────────────────────────────────────────────────
test_that("When reference and target batches are identical, corrections are small", {
  set.seed(42)
  dat <- simulate_batches(G = 1, n = 15, p = 6)
  
  out <- BEgLasso(dat$X0, dat$X0,           # <-- identical batches
                  penal.rho = 0.1, penal.ksi = 0.1, penal.gamma = 0.1,
                  eps = 1e-2, print.detail = FALSE)
  
  expect_equal(out$coef.a, rep(0.1, dat$p), tolerance = 0.35)
  expect_equal(out$coef.b, rep(0, dat$p), tolerance = 0.35)
})

# 4 ────────────────────────────────────────────────────────────────
test_that("Basic input validation — wrong shapes / types generate errors", {
  
  # Non‑list inputs ------------------------------------------------
  expect_error(
    BEgLasso(matrix(rnorm(20), 4, 5),
             list(matrix(rnorm(20), 4, 5)),
             penal.rho = 0.1, penal.ksi = 0.1, penal.gamma = 0.1,
             eps = 1e-2, print.detail = FALSE),
    regexp = "array"
  )
  
  # List length mismatch -------------------------------------------
  dat <- simulate_batches(G = 2, n = 10, p = 4)
  expect_error(
    BEgLasso(dat$X0, dat$X1[1],              # lengths differ
             penal.rho = 0.1, penal.ksi = 0.1, penal.gamma = 0.1,
             eps = 1e-2, print.detail = FALSE),
    regexp = "non-conformable"
  )
})

# 5 ────────────────────────────────────────────────────────────────
test_that("Output precision matrices are positive‑definite (eigenvalues > 0)", {
  set.seed(888)
  dat <- simulate_batches(G = 1, n = 20, p = 5)
  res <- BEgLasso(dat$X0, dat$X1,
                  penal.rho = 0.2, penal.ksi = 0.1, penal.gamma = 0.1,
                  eps = 1e-2, print.detail = FALSE)
  
  eig <- eigen(res$Theta[[1]], symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eig > 0))
})
