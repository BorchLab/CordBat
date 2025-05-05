# test script for CordBat.R - testcases are NOT comprehensive!

## -----------------------------------------------------------
## 0.  Skip if heavy dependencies are missing ----------------
## -----------------------------------------------------------
needed_fns <- c("CordBat", "BEgLasso", "findBestPara", "getAllCom",
                "StARS", "selrho.useCVBIC", "graphicalLasso")
if (!all(vapply(needed_fns, exists, logical(1)))) {
  skip("Some helper functions for CordBat are not available")
}

## -----------------------------------------------------------
## Helper to quickly simulate tiny expression data -----------
## -----------------------------------------------------------
sim_data <- function(n_batch = 2, per_batch = 6, p = 4, add_qc = FALSE) {
  set.seed(333)
  n  <- n_batch * per_batch + if (add_qc) 2 else 0
  X  <- matrix(rnorm(n * p), n, p)
  batch <- rep(seq_len(n_batch), each = per_batch)
  if (add_qc) {
    group <- c(rep(1, n - 2), "QC", "QC")
    batch[(n - 1):n] <- 1                      # QC in ref batch
  } else {
    group <- rep(1, n)
  }
  list(X = X, batch = batch, group = group)
}

beta_target <- 0.05                    # for StARS stability tests indirectly

## -----------------------------------------------------------
## 1.  Basic structure with skip.impute = TRUE ---------------
## -----------------------------------------------------------
test_that("CordBat returns complete structure and keeps reference batch unchanged", {
  dat <- sim_data()
  res <- CordBat(dat$X, dat$batch,
                 ref.batch = 1,
                 print.detail = FALSE,
                 skip.impute = TRUE)
  
  expect_type(res, "list")
  expect_named(res, c("batch.level", "delsampIdx", "batch.new",
                      "group.new", "X.delout", "X.cor",
                      "X.cor.1", "X.cor.withQC", "Xcor.para"),
               ignore.order = TRUE)
  
  n <- nrow(dat$X); p <- ncol(dat$X)
  expect_equal(dim(res$X.delout), c(n, p))
  expect_equal(dim(res$X.cor),    c(n, p))
  expect_equal(dim(res$X.cor.1),  c(n, p))
  
  # Reference batch samples remain identical in X.cor
  ref_idx <- which(dat$batch == 1)
  expect_equal(res$X.cor[ref_idx, ], dat$X[ref_idx, ], tolerance = 1e-8)
})

## -----------------------------------------------------------
## 2.  Handling of QC samples --------------------------------
## -----------------------------------------------------------
test_that("QC samples bypass correction but are returned in X.cor.withQC", {
  dat <- sim_data(add_qc = TRUE)
  res <- CordBat(dat$X, dat$batch, dat$group,
                 grouping     = TRUE,
                 ref.batch    = 1,
                 print.detail = FALSE,
                 skip.impute  = TRUE)
  
  expect_true(!is.null(res$X.cor.withQC))
  expect_equal(dim(res$X.cor.withQC), dim(dat$X))
  
  qc_idx <- which(dat$group == "QC")
  # QC rows in X.cor.withQC should equal original data
  expect_equal(res$X.cor.withQC[qc_idx, ], dat$X[qc_idx, ], tolerance = 1e-8)
})

## -----------------------------------------------------------
## 3.  print.detail toggles console output -------------------
## -----------------------------------------------------------
test_that("print.detail produces messages", {
  dat <- sim_data()
  expect_silent(
    CordBat(dat$X, dat$batch, ref.batch = 1,
            print.detail = FALSE, skip.impute = TRUE)
  )
  
  expect_message(
    CordBat(dat$X, dat$batch, ref.batch = 1,
            print.detail = TRUE, skip.impute = TRUE),
    regexp = ".", all = FALSE
  )
})

## -----------------------------------------------------------
## 4.  Error on invalid reference batch ----------------------
## -----------------------------------------------------------
test_that("CordBat errors if ref.batch is not present", {
  dat <- sim_data()
  expect_error(
    CordBat(dat$X, dat$batch, ref.batch = "nonexistent",
            print.detail = FALSE, skip.impute = TRUE),
    regexp = "invalid"
  )
})

