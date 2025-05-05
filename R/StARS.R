#' Stability Approach to Regularization Selection (StARS)
#'
#' This function selects an optimal regularization parameter (`rho`) for 
#' high-dimensional graphical models using the StARS stability approach.
#'
#' @param X A numeric matrix where rows represent samples and columns represent features.
#' @param b The number of subsamples to draw from the dataset for stability selection.
#' @param M The number of subsampling iterations.
#' @param print.detail Logical. If TRUE, prints detailed progress (default: TRUE).
#'
#' @return A numeric vector containing:
#' \itemize{
#'   \item \code{Sel.rho} - Selected optimal regularization parameter.
#'   \item \code{D_var} - Stability measure at the selected `rho`.
#' }
#'
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(100), 10, 10)
#' StARS(X, b = 7, M = 50)
#'
#' @export
StARS <- function(X, b, M, print.detail = TRUE) {
  set.seed(42)
  N <- nrow(X)
  p <- ncol(X)
  beta <- 0.05
  D_var <- 0
  Sel.rho <- 0.5
  
  # Function to compute stability measure
  compute_stability <- function(X, b, M, rho) {
    sum.fai <- matrix(0, p, p)
    
    for (i in seq_len(M)) {
      subsampIdx <- sample.int(N, b, replace = FALSE)
      Si <- X[subsampIdx, , drop = FALSE]
      Theta <- graphicalLasso(Si, rho, print.detail = print.detail)$Theta
      Theta <- Theta * (abs(Theta) > 1e-6)  # Thresholding
      sum.fai <- sum.fai + (Theta != 0)
    }
    
    thetb <- sum.fai / M
    ksib <- 2 * thetb * (1 - thetb)
    sum(ksib[upper.tri(ksib, diag = FALSE)]) / (p * (p - 1) / 2)
  }
  
  # Initial stability evaluation
  if (print.detail) {
    detailOutput <- capture.output(
      Db.5 <- compute_stability(X, b, M, rho = 0.5)
    , type = "message")
    if(length(detailOutput) > 0) {
      # Print the entire captured output as one message
      message(paste(unique(detailOutput)))
    }
  } else {
    Db.5 <- compute_stability(X, b, M, rho = 0.5)
  }
  
  # Determine search range for `rho`
  if (Db.5 > beta) {
    rhos <- seq(0.9, 0.5, by = -0.1)
  } else {
    Db.1 <- compute_stability(X, b, M, rho = 0.1)
    
    if (Db.1 > beta) {
      rhos <- seq(0.5, 0.1, by = -0.1)
    } else {
      Db.05 <- compute_stability(X, b, M, rho = 0.05)
      
      if (Db.05 > beta) {
        rhos <- seq(0.1, 0.05, by = -0.01)
      } else {
        rhos <- seq(0.05, 0.01, by = -0.01)
      }
    }
  }
  
  # Find the optimal `rho`
  for (rho in rhos) {
    Db <- compute_stability(X, b, M, rho)
    
    if (Db <= beta) {
      D_var <- max(D_var, Db)
    } else {
      Sel.rho <- ifelse(rho > 0.1, rho + 0.1, rho + 0.01)
      break
    }
  }
  
  if (print.detail) {
    message("Selected rho =", Sel.rho, "\n")
  }
  
  return(c(Sel.rho, D_var))
}
