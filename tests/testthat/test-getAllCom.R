# test script for getAllCom.R - testcases are NOT comprehensive!

## -----------------------------------------------------------
## Helper: build block‑correlated matrix ---------------------
## -----------------------------------------------------------
sim_block_data <- function(n = 40, blocks = c(3, 3, 3), noise_sd = 0.02) {
  set.seed(999)
  X <- NULL
  grp_labels <- rep(seq_along(blocks), times = blocks)
  for (k in seq_along(blocks)) {
    base <- rnorm(n)                      # same driver for block
    block <- matrix(rep(base, each = blocks[k]), n, blocks[k])
    block <- block + matrix(rnorm(n * blocks[k],
                                  sd = noise_sd), n, blocks[k])
    X <- cbind(X, block)
  }
  list(X = X, labels = grp_labels)
}

## -----------------------------------------------------------
## 1.  Basic structure and coverage --------------------------
## -----------------------------------------------------------
test_that("getAllCom returns complete, non‑overlapping communities", {
  set.seed(1)
  X <- matrix(rnorm(200), 20, 10)
  coms <- getAllCom(X)
  
  # Correct type
  expect_type(coms, "list")
  
  p <- ncol(X)
  all_idx <- sort(unlist(coms))
  # Contains every feature exactly once
  expect_equal(all_idx, seq_len(p))
  
  # No overlaps
  expect_equal(length(all_idx), p)
})

## -----------------------------------------------------------
## 2.  Direct test of ComtyDet on simple adjacency -----------
## -----------------------------------------------------------
test_that("ComtyDet splits a toy graph into two communities", {
  # Build 4‑node graph: nodes 1‑2 tightly connected, 3‑4 tightly connected
  G <- matrix(0.1, 4, 4)
  diag(G) <- 1
  G[1, 2] <- G[2, 1] <- 0.9
  G[3, 4] <- G[4, 3] <- 0.9
  
  inputCOM   <- list(1:4)     # start with one community containing all nodes
  min_sz     <- 2
  out_com    <- ComtyDet(G, inputCOM, min_sz)
  
  # We expect exactly two communities of size 2
  sizes <- vapply(out_com, length, integer(1))
  expect_equal(sort(sizes), c(2, 2))
  
  # Check membership is {1,2} and {3,4} (order‑free)
  expect_true(setequal(out_com[[1]], c(1, 2)) ||
                setequal(out_com[[1]], c(3, 4)))
  expect_true(setequal(out_com[[2]], c(1, 2)) ||
                setequal(out_com[[2]], c(3, 4)))
})
