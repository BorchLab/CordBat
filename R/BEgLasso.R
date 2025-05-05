#' Batch Effect Correction using Group Graphical Lasso
#'
#' This function applies a batch effect correction method using the 
#' Group Graphical Lasso (GGM) approach. It estimates correction coefficients
#' and applies them iteratively to adjust for batch effects in multi-group datasets.
#'
#' @param X0.glist A list of matrices where each matrix corresponds to a batch 
#' in group 0 (reference batch).
#' @param X1.glist A list of matrices where each matrix corresponds to a batch 
#' in group 1 (batch to be corrected).
#' @param penal.rho Regularization parameter for the graphical lasso.
#' @param penal.ksi Regularization parameter for coefficient update.
#' @param penal.gamma Additional penalty parameter for coefficient update.
#' @param eps Convergence threshold for stopping criteria.
#' @param print.detail Logical flag to print processing details (default: TRUE).
#' 
#' @return A list containing:
#' \itemize{
#'   \item \code{Theta} - List of precision matrices for each group.
#'   \item \code{X1.cor} - List of corrected matrices for group 1.
#'   \item \code{coef.a} - Corrected scaling coefficients.
#'   \item \code{coef.b} - Corrected offset coefficients.
#' }
#' 
#' @examples
#' set.seed(123)
#' X0.glist <- list(matrix(rnorm(100), 10, 10), matrix(rnorm(100), 10, 10))
#' X1.glist <- list(matrix(rnorm(100), 10, 10), matrix(rnorm(100), 10, 10))
#' res <- BEgLasso(X0.glist, X1.glist, 0.1, 0.1, 0.1, 1e-4, print.detail = FALSE)
#'
#' @importFrom utils capture.output
#' @export
BEgLasso <- function(X0.glist, 
                     X1.glist, 
                     penal.rho, 
                     penal.ksi,
                     penal.gamma, 
                     eps, 
                     print.detail) {
  
  G <- length(X0.glist)
  p <- ncol(X0.glist[[1]])
  
  # Precompute sample sizes
  N0_gvec <- sapply(X0.glist, nrow)
  N1_gvec <- sapply(X1.glist, nrow)
  
  # Compute mean matrices
  X0.m.gmat <- do.call(cbind, lapply(X0.glist, colMeans))
  X1.m.gmat <- do.call(cbind, lapply(X1.glist, colMeans))
  
  # Initialize variables
  Theta.list <- vector("list", G)
  B.list <- vector("list", G)
  coef.as <- X0.m.gmat / X1.m.gmat
  coef.a <- rowMeans(coef.as)
  coef.b <- numeric(p)
  
  coef.A <- diag(coef.a)
  
  # Initialize Theta and B for each group
  for (g in seq_len(G)) {
    Theta.list[[g]] <- matrix(0, p, p)
    B.list[[g]] <- matrix(0, p - 1, p)
  }
  
  # Correct X1 using initial coefficients
  X1.cor.glist <- lapply(seq_len(G), function(g) {
    coef.B.gi <- matrix(rep(coef.b, each = N1_gvec[g]), nrow = N1_gvec[g])
    X1.glist[[g]] %*% coef.A + coef.B.gi
  })
  
  # Compute initial covariance matrices
  W.list <- lapply(seq_len(G), function(g) {
    X.gi <- rbind(X0.glist[[g]], X1.cor.glist[[g]])
    S0_gi <- cov(scale(X.gi, center = TRUE, scale = TRUE))
    S0_gi + penal.rho * diag(1, p)
  })
  
  # Iteration variables
  finished.gmat <- matrix(FALSE, (p - 1) * p / 2, G)
  times0.gvec <- integer(G)
  times1.gvec <- integer(G)
  
  # Iterative correction
  repeat {
    W_old.list <- W.list
    
    S.list <- lapply(seq_len(G), function(g) {
      coef.B.gi <- matrix(rep(coef.b, each = N1_gvec[g]), nrow = N1_gvec[g])
      X1.gi.cor <- X1.glist[[g]] %*% coef.A + coef.B.gi
      X.gi <- rbind(X0.glist[[g]], X1.gi.cor)
      cov(scale(X.gi, center = TRUE, scale = TRUE))
    })
    
    # Update B and W using graphical lasso approach
    if (print.detail) {
      # Capture all output generated during the loop execution
      detailOutput <- capture.output(
        for (j in seq_len(p)) {
          idx <- setdiff(seq_len(p), j)
          
          for (g in seq_len(G)) {
            Wi_11 <- W.list[[g]][idx, idx]
            si_12 <- S.list[[g]][idx, j]
            
            B.list[[g]][, j] <- CDfgL(Wi_11, B.list[[g]][, j], si_12, penal.rho, 
                                      print.detail = print.detail)
            W.list[[g]][idx, j] <- Wi_11 %*% B.list[[g]][, j]
            W.list[[g]][j, idx] <- W.list[[g]][idx, j]
          }
        }
      , type = "message")
      if(length(detailOutput) > 0) {
        # Print the entire captured output as one message
        message(paste(unique(detailOutput)))
      }
    } else {
      for (j in seq_len(p)) {
        idx <- setdiff(seq_len(p), j)
        
        for (g in seq_len(G)) {
          Wi_11 <- W.list[[g]][idx, idx]
          si_12 <- S.list[[g]][idx, j]
          
          B.list[[g]][, j] <- CDfgL(Wi_11, B.list[[g]][, j], si_12, penal.rho, 
                                    print.detail = print.detail)
          W.list[[g]][idx, j] <- Wi_11 %*% B.list[[g]][, j]
          W.list[[g]][j, idx] <- W.list[[g]][idx, j]
        }
      }
    }
    
    
    # Update Theta
    for (g in seq_len(G)) {
      for (j in seq_len(p)) {
        idx <- setdiff(seq_len(p), j)
        Theta.list[[g]][j, j] <- W.list[[g]][j, j] - t(W.list[[g]][idx, j]) %*% B.list[[g]][, j]
        Theta.list[[g]][idx, j] <- -B.list[[g]][, j] * Theta.list[[g]][j, j]
        Theta.list[[g]][j, idx] <- Theta.list[[g]][idx, j]
      }
    }
    
    # Check convergence
    W_new.ndiag <- lapply(seq_len(G), function(g) W.list[[g]][upper.tri(W.list[[g]], diag = FALSE)])
    W_old.ndiag <- lapply(seq_len(G), function(g) W_old.list[[g]][upper.tri(W_old.list[[g]], diag = FALSE)])
    
    zeroIdx.list <- lapply(W_old.ndiag, function(w) w < 1e-5)
    
    for (g in seq_len(G)) {
      zeroIdx_g <- zeroIdx.list[[g]]
      if (any(zeroIdx_g)) {
        if (all(W_new.ndiag[[g]][zeroIdx_g] < 1e-5)) {
          times0.gvec[g] <- times0.gvec[g] + 1
        } else {
          times0.gvec[g] <- 0
        }
        if (times0.gvec[g] >= 5) {
          finished.gmat[zeroIdx_g, g] <- TRUE
        }
      }
      
      if (any(!zeroIdx_g)) {
        dW_g <- max(abs((W_new.ndiag[[g]][!zeroIdx_g] - W_old.ndiag[[g]][!zeroIdx_g]) 
                        / W_old.ndiag[[g]][!zeroIdx_g]))
        
        if (dW_g < eps) {
          times1.gvec[g] <- times1.gvec[g] + 1
        } else {
          times1.gvec[g] <- 0
        }
        
        if (times1.gvec[g] >= 3) {
          finished.gmat[!zeroIdx_g, g] <- TRUE
        }
      }
    }
    
    if (all(as.vector(finished.gmat))) break
    
    # Update coefficients
    coef.update <- update_CorrectCoef(X0.glist, X1.glist, Theta.list, coef.a, 
                                      coef.b, penal.ksi, penal.gamma, print.detail)
    coef.a <- coef.update$coef.a
    coef.b <- coef.update$coef.b
  }
  
  # Ensure non-negative coef.a
  coef.a[coef.a < 0] <- 0.05
  
  list(Theta = Theta.list, X1.cor = X1.cor.glist, coef.a = coef.a, coef.b = coef.b)
}
