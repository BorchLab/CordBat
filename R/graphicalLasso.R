#' Graphical Lasso Algorithm
#'
#' This function estimates a sparse inverse covariance matrix (precision matrix) 
#' using the graphical lasso algorithm with an L1-penalized covariance estimation approach.
#'
#' @param X A numeric matrix where rows represent samples and columns represent features.
#' @param rho Regularization parameter controlling sparsity.
#' @param print.detail Logical flag to print processing details (default: TRUE).
#' 
#' @return A list containing:
#' \itemize{
#'   \item \code{Theta} - The estimated precision (inverse covariance) matrix.
#'   \item \code{W} - The estimated covariance matrix.
#' }
#' 
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(100), 10, 10)
#' result <- graphicalLasso(X, rho = 0.1)
#' 
#' @importFrom stats cov
#' @importFrom utils capture.output
graphicalLasso <- function(X, rho, print.detail) {
  N <- nrow(X)
  p <- ncol(X)
  
  if (N > 1) { # When there are multiple samples, center and scale normally. 
    # Standardize data: center and scale
    X <- scale(X, center = TRUE, scale = TRUE) 
    # Compute covariance matrix
    S <- cov(X) 
  } else { # With one sample, standardization isn’t possible since the sample variance is undefined.
    # Define S as a zero matrix (since there is no variability)  
    S <- matrix(0, p, p) 
  }
  
  # Initialize variables
  Theta <- matrix(0, p, p)
  W <- S + rho * diag(p)  # Regularized covariance matrix
  B <- matrix(0, p - 1, p)
  
  threshold <- 1e-5  # Convergence threshold
  
  # Iterative update
  repeat {
    W_old <- W  # Store previous iteration
    
    if (print.detail) {
      # Capture all output generated during the loop execution
      detailOutput <- capture.output(
        for (i in seq_len(p)) {
          idx <- setdiff(seq_len(p), i)  # Exclude current index
          
          W_11 <- W[idx, idx]
          s_12 <- S[idx, i]
          
          # Update B using penalized regression
          B[, i] <- CDfgL(W_11, B[, i], s_12, rho, print.detail = print.detail)
          W[idx, i] <- W_11 %*% B[, i]
          W[i, idx] <- W[idx, i]
        }, type = "message")
      if(length(detailOutput) > 0) {
        # Print the entire captured output as one message
        message(paste(unique(detailOutput)))
      }
    } else {
      for (i in seq_len(p)) {
        idx <- setdiff(seq_len(p), i)  # Exclude current index
        
        W_11 <- W[idx, idx]
        s_12 <- S[idx, i]
        
        # Update B using penalized regression
        B[, i] <- CDfgL(W_11, B[, i], s_12, rho, print.detail = print.detail)
        W[idx, i] <- W_11 %*% B[, i]
        W[i, idx] <- W[idx, i]
      }
    }
    
    # Convergence check
    dW <- W - W_old
    S_ndiag <- S[upper.tri(S, diag = FALSE)]
    
    # Avoid division by zero in convergence metric
    if (mean(abs(dW)) / mean(abs(S[upper.tri(S, diag = FALSE)] + 1e-6)) < threshold) break
  }
  
  # Compute precision matrix Theta
  for (i in seq_len(p)) {
    idx <- setdiff(seq_len(p), i)
    Theta[i, i] <- 1 / (W[i, i] - t(W[idx, i]) %*% B[, i])
    Theta[idx, i] <- -B[, i] * Theta[i, i]
    Theta[i, idx] <- Theta[idx, i]
  }
  
  # Return results
  list(Theta = Theta, W = W)
}
